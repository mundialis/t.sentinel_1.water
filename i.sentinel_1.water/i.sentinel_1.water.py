#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      i.sentinel_1.water
# AUTHOR(S):   Guido Riembauer
#
# PURPOSE:     Water area mapping based on Sentinel 1 GRD data
#
# COPYRIGHT:    (C) 2020-2022 by mundialis GmbH & Co. KG and the GRASS Development Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#############################################################################
# %Module
# % description: Automatic thresholding and region growing to detect water areas in S1 GRD scenes.
# % keyword: raster
# % keyword: Sentinel-1
# % keyword: Hydrology
# %End

# %option
# % key: input
# % type: string
# % required: yes
# % multiple: no
# % label: S-1 map name
# % description: S-1 data must be calibrated and in logarithmic scale (dB). The second last part of the scene name split by '_' must contain the polarisation (VV or VH)
# %end

# %option
# % key: output
# % type: string
# % required: yes
# % multiple: no
# % label: Output binary water map
# %end

# %option
# % key: hand_rast
# % type: string
# % required: yes
# % multiple: no
# % label: Name of the HAND index raster
# % description: Raster containing Height Above Nearest Drainage (HAND) values
# %end

# %option
# % key: water_core_thresh
# % type: string
# % required: no
# % multiple: no
# % label: Manual threshold for core water areas (in dB)
# % description: All areas below water_core_thresh will be selected as core water areas
# %end

# %option
# % key: water_boundary_thresh
# % type: string
# % required: no
# % multiple: no
# % label: Manual threshold for water boundary areas (in dB)
# % description: All areas below water_boundary_thresh will be selected as water boundary areas
# %end

# %option G_OPT_MEMORYMB
# % description: Available memory for running i.segment
# %end

# %option
# % key: iterations
# % type: integer
# % required: no
# % multiple: no
# % label: Iterations of region growing algorithm (i.segment)
# % answer: 50
# %end

# %option
# % key: segment_threshold
# % type: string
# % required: no
# % multiple: no
# % label: Threshold of region merging for region growing algorithm (i.segment)
# % description: Must be between 0 and 1
# % answer: 0.1
# %end

# %option
# % key: hand_max
# % type: string
# % required: no
# % multiple: no
# % label: Maximum HAND value (meters above drainage network) where water areas are assumed
# % answer: 15.0
# %end

# %option
# % key: polarization
# % type: string
# % required: yes
# % multiple: no
# % label: Polarization of dataset
# % description: Must be VV or VH
# % options: VV,VH
# %end

# %flag
# % key: m
# % description: Don't run automatic threshold identification. Provide thresholds via water_boundary_thresh and water_core_thresh instead.
# %end

import os
import atexit
import grass.script as grass
from grass.pygrass import raster
from grass.pygrass.gis.region import Region
from grass.pygrass.raster.abstract import RasterAbstractBase
from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter1d
from skimage.filters import threshold_minimum
import numpy as np

# define lists to cleanup
rm_regions = []
rm_vectors = []
rm_rasters = []
rm_groups = []


def cleanup():
    """
    Function that loops through global lists with regions, vectors, rasters,
    and groups and removes all items.
    """

    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmr in rm_regions:
        if rmr in [x for x in grass.parse_command("g.list", type="region")]:
            grass.run_command("g.remove", type="region", name=rmr, **kwargs)
    for rmv in rm_vectors:
        if grass.find_file(name=rmv, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmv, **kwargs)
    for rmrast in rm_rasters:
        if grass.find_file(name=rmrast, element="raster")["file"]:
            grass.run_command("g.remove", type="raster", name=rmrast, **kwargs)
    for rmgroup in rm_groups:
        if grass.find_file(name=rmgroup, element="group")["file"]:
            grass.run_command("g.remove", type="group", name=rmgroup, **kwargs)


def raster2numpy(rastname, mapset=""):
    """
    Function that converts a GRASS GIS raster into a numpy array
    Args:
        rastname (string): Name of the input rastermap
        mapset (string): Name of the input mapset
    Returns:
        np.array(rast): rastermap as numpy array
    """
    # we have to use RasterAbstractBase.set_region() to set the region for the
    # GRASS C libs
    curr_reg = Region()
    tmpraster = RasterAbstractBase(rastname)
    tmpraster.set_region(curr_reg)
    with raster.RasterRow(rastname, mapset=mapset, mode="r") as rast:
        return np.array(rast)


def bimodalcheck(array):
    """
    Function that smoothes the histogram of input data and counts the number of
    local maxima
    Args:
        array (numpy array): a 1D numpy array containing data
    Returns:
        maxima (list): a list of indices where local maxima occur in the bins
        firstmax (int): the position of the first maximum in the bins
        hist_bars_smoothed (numpy array): smoothed histogram bars
        bin_centers (numpy array): positions of histogram bars (bins)
    """
    # compute histogram, bins are dependent on input length
    num_bins = int(np.shape(array)[0] / 1000)
    if num_bins < 100:
        num_bins = 100
    tile_hist = np.histogram(array, bins=num_bins)
    hist_bars = tile_hist[0]
    # convert the bins to their center coordinates
    bins = tile_hist[1]
    bin_width = bins[-1] - bins[-2]
    bin_centers = bins + bin_width
    bin_centers = bin_centers[:-1]
    # smooth the histogram, sigma=2 provided good results in tests
    hist_bars_smoothed = gaussian_filter1d(input=hist_bars, sigma=2)
    # compute the local maxima from the smoothed histogram, the prominence is
    # dependent on the bins
    prominence = num_bins
    maxima = find_peaks(hist_bars_smoothed, prominence=prominence)[0]
    firstmax = maxima[0]
    return maxima, firstmax, hist_bars_smoothed, bin_centers


def tiled_autothreshold(
    input, tilesize, nsres, ewres, val_pixels, check_value, polarization
):
    """
    Function that splits an input GRASS raster into regular tiles, checks the
    tiles for bimodality and for each bimodal tile saves the value of the first
    relative maximum and the position of the valley between maxima.
    Args:
        input (string): the input GRASS raster
        tilesize (string): the size of one side of the tile in meters
        nsres (string): north-south resolution of the tiles
        ewres (string): east-west resolution of the tiles
        val_pixels (int): minimum number of no-nan pixels in the tile
        check_value (int): minimum number of occurences of the first maximum in
                           the histogram
        polarizaion (string): VV or VH, used to determine absolute check
                              thresholds
    Returns:
        first_mode_list (list): list of all first relative maxima from bimodal
                                tiles
        threshold_min_list (list): list of all histogram valleys from bimodal
                                   tiles
    """

    # lists for threshold and mode
    threshold_min_list = []
    first_mode_list = []

    tilesize_str = "%s,%s" % (tilesize, tilesize)
    gridname = "grid_%s_%s" % (tilesize, os.getpid())
    # create grid
    grass.use_temp_region()
    grass.run_command("g.region", raster=input, zoom=input)
    grass.run_command("v.mkgrid", box=tilesize_str, map=gridname, quiet=True)
    grass.del_temp_region()
    rm_vectors.append(gridname)
    cats = list(
        grass.parse_command("v.category", input=gridname, option="print").keys()
    )
    cat_ids_bimodal = []
    # setting check values for thresholds depending on polarization
    if polarization == "VV":
        firstmode_max = -17
        minimum_max = -14
    elif polarization == "VH":
        firstmode_max = -24
        minimum_max = -19
    # loop through the tiles
    grass.message(
        _("Looping through tiles of size %s for scene %s..." % (tilesize_str, input))
    )
    cats_done = []
    for cat in cats:
        grass.use_temp_region()
        # set region to tile
        cat_reg = grass.parse_command(
            "v.db.select", flags="r", map=gridname, where="cat=%s" % cat
        )
        grass.run_command(
            "g.region", flags="a", nsres=nsres, ewres=ewres, **cat_reg, quiet=True
        )
        # check if there are at least (val_pixels) in the tile
        tilestats = grass.parse_command("r.univar", map=input, flags="g")
        valid_pixels = int(tilestats["n"])
        if valid_pixels > val_pixels:
            # convert to numpy array
            tiledata = raster2numpy(input)
            tile_data_nonan = tiledata[~np.isnan(tiledata)]
            # check for bimodality and return the smoothed histogram
            maxima, firstmax, hist_bars_smoothed, bin_centers = bimodalcheck(
                tile_data_nonan
            )
            if len(maxima) > 1:
                # do a check that there is at least (check_value) pixels in
                # the first maximum
                if hist_bars_smoothed[firstmax] > check_value:
                    # do a check that the first maximum lies below x dB
                    # (such that its for sure water)
                    if bin_centers[firstmax] < firstmode_max:
                        # try to apply the threshold_minimum function
                        # if that fails, skip
                        try:
                            thresh = threshold_minimum(tile_data_nonan)
                            # check if the minimum is below the pol-dependent
                            # max-minimum and that it is larger than
                            # the first mode
                            if thresh < minimum_max and bin_centers[firstmax] < thresh:
                                threshold_min_list.append(thresh)
                                first_mode_list.append(bin_centers[firstmax])
                                cat_ids_bimodal.append(cat)
                        except Exception:
                            pass
            tiledata = None
            tile_data_nonan = None
        grass.del_temp_region()
        cats_done.append(cat)
        # grass.percent(int(cat)+1, len(cats), 1)
        # percentage_done = str(round(len(cats_done)/len(cats)))
        # print('%s%%' % percentage_done, end="")
    return first_mode_list, threshold_min_list


def adapted_regiongrowing(
    input, output, minsize, threshold, iterations, memory, seeds, bounds
):
    """
    Function that runs region growing within i.segment and only keeps those
    regions that originate from a seed
    Args:
        input (string): the input GRASS raster
        output (string): name of the output GRASS raster
        minsize (string): minimum region size (parameter of i.segment)
        threshold (string): region merging threshold (parameter of i.segment)
        iterations (string): maximum number of iterations
                             (parameter of i.segment)
        memory (string): available memory for region growing
                         (parameter of i.segment)
        seeds (string): name of raster map to be used as seeds
                        (parameter of i.segment)
        bounds (string) name of raster map to be used as boundary
                        (parameter of i.segment)
    """

    # add raster to a group
    groupname = "%s_group_%s" % (input, os.getpid())
    grass.run_command("i.group", input=input, group=groupname, quiet=True)
    rm_groups.append(groupname)
    # run conventional region growing
    reg_grow_outname = "reggrow_tmp_%s" % os.getpid()
    grass.run_command(
        "i.segment",
        group=groupname,
        output=reg_grow_outname,
        minsize=minsize,
        threshold=threshold,
        iterations=iterations,
        method="region_growing",
        memory=memory,
        seeds=seeds,
        bounds=bounds,
        overwrite=True,
        quiet=True,
    )
    # get only the regions that have their origin in the seeds
    core_values_name = "core_values_tmp_%s" % os.getpid()
    grass.run_command(
        "r.mapcalc",
        expression="%s = if(isnull(%s*%s),null(),%s)"
        % (core_values_name, seeds, reg_grow_outname, reg_grow_outname),
        overwrite=True,
        quiet=True,
    )
    rm_rasters.append(core_values_name)
    # get the unique categories from the regions
    region_list = list(
        grass.parse_command("r.stats", input=core_values_name, flags="nl").keys()
    )
    # prepare a string for the reclass file
    reclass_string = ""
    for item in region_list:
        reclass_string += "%s = %s\n" % (item, item)
    reclass_string += "* = NULL"
    # build a file for r.reclass
    rulefile = grass.tempfile()
    with open(rulefile, "w") as file:
        file.write(reclass_string)
    # run r.reclass
    reclass_table = "reclasstable_%s" % os.getpid()
    grass.run_command(
        "r.reclass",
        input=reg_grow_outname,
        output=reclass_table,
        rules=rulefile,
        overwrite=True,
        quiet=True,
    )
    rm_rasters.append(reclass_table)
    rm_rasters.append(reg_grow_outname)
    # remove rulefile
    grass.try_remove(rulefile)
    # convert reclassed to a map
    reclass_map = "reclassmap_%s" % os.getpid()
    grass.run_command(
        "r.mapcalc",
        expression="%s = %s" % (reclass_map, reclass_table),
        overwrite=True,
        quiet=True,
    )
    rm_rasters.append(reclass_map)
    # convert to binary map and save as output
    grass.run_command(
        "r.mapcalc",
        expression="%s = if(isnull(%s),null(),1)" % (output, reclass_map),
        overwrite=True,
        quiet=True,
    )


def main():
    # get the user parameters
    input = options["input"]
    hand = options["hand_rast"]
    memory = int(options["memory"])
    output = options["output"]
    hand_max = options["hand_max"]
    iterations = str(options["iterations"])
    segment_threshold = float(options["segment_threshold"])
    if not 0.0 <= segment_threshold <= 1.0:
        grass.fatal(_("Segment threshold must be between 0.0 and 1.0"))
    segment_threshold = str(segment_threshold)
    manual_thresh = flags["m"]
    if manual_thresh:
        water_boundary_thresh = float(options["water_boundary_thresh"])
        water_core_thresh = float(options["water_core_thresh"])
    polarization = options["polarization"]

    # get current region info
    in_reg = grass.region()
    nsres = in_reg["nsres"]
    ewres = in_reg["ewres"]

    empty_rast = True

    if not manual_thresh:
        #  run image tiling and automatic threshold identification
        grass.message(
            _(
                "Automatic identification of threshold and"
                " water mode of scene %s" % input
            )
        )
        tilesizes = [5000, 4000, 3000]
        for tilesize in tilesizes:
            if empty_rast is True:
                grass.message(_("Setting tilesize to %s m" % tilesize))
                first_mode_list, threshold_min_list = tiled_autothreshold(
                    input, tilesize, nsres, ewres, 50000, 500, polarization
                )
                if len(first_mode_list) > 0 and len(threshold_min_list) > 0:
                    if np.mean(first_mode_list) < np.mean(threshold_min_list):
                        empty_rast = False
        if empty_rast is True:
            grass.warning(
                _(
                    "No water area found for any of the tilesizes."
                    " Consider providing manual thresholds..."
                )
            )
        else:
            first_mode = np.mean(first_mode_list)
            threshold_min = np.mean(threshold_min_list)

    else:
        empty_rast = False
        first_mode = water_core_thresh
        threshold_min = water_boundary_thresh

    grass.use_temp_region()
    grass.run_command("g.region", raster=input)
    if empty_rast is False:
        # apply thresholds to get water core and water boundary area

        print(
            "Setting water core area threshold to %s dB for scene %s "
            % (np.round(first_mode, 2), input)
        )
        print(
            "Setting water boundary threshold to %s dB for scene %s"
            % (np.round(threshold_min, 2), input)
        )

        water_core = "water_core_%s" % os.getpid()
        water_boundary = "water_boundary_%s" % os.getpid()

        grass.run_command(
            "r.mapcalc",
            expression="%s = if(%s < %s,1,null())" % (water_core, input, first_mode),
            quiet=True,
        )
        grass.run_command(
            "r.mapcalc",
            expression="%s = if(%s < %s,1,null())"
            % (water_boundary, input, threshold_min),
            quiet=True,
        )
        rm_rasters.append(water_core)
        rm_rasters.append(water_boundary)

        # run region growing based on the core areas
        grass.message(
            _(
                "Running region growing and correcting water areas "
                "with HAND threshold %s meters for scene %s " % (hand_max, input)
            )
        )
        outname_reg = "water_areas_tmp_%s" % os.getpid()
        adapted_regiongrowing(
            input,
            outname_reg,
            minsize="1",
            threshold=segment_threshold,
            iterations=iterations,
            memory=str(memory),
            seeds=water_core,
            bounds=water_boundary,
        )
        rm_rasters.append(outname_reg)

        # Correcting the water areas based on HAND thresholds
        outname_hand = "water_areas_hand_tmp_%s" % os.getpid()
        grass.run_command(
            "r.mapcalc",
            expression="%s = if(%s == 1 && %s < %s,1, null())"
            % (outname_hand, outname_reg, hand, hand_max),
            quiet=True,
        )
        rm_rasters.append(outname_hand)

        # fill the remaining pixels that are within the scene with 0
        # (not null!)
        grass.run_command(
            "r.mapcalc",
            expression="%s = if(isnull(%s),null(),"
            "if(isnull(%s),0,1))" % (output, input, outname_hand),
            quiet=True,
        )
    else:
        # give an empty map as output
        grass.run_command(
            "r.mapcalc",
            expression="%s = if(isnull(%s),null(),0)" % (output, input),
            quiet=True,
        )

    grass.del_temp_region()


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
