#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      t.sentinel_1.watermasks
# AUTHOR(S):   Guido Riembauer
#
# PURPOSE:     Calculates water masks from Sentinel-1 scenes and registers them
#              in a space time dataset
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
# % description: Calculates water masks from Sentinel-1 scenes and registers them in a space time dataset.
# % keyword: temporal
# % keyword: Sentinel-1
# % keyword: SAR
# % keyword: watermasks
# %End

# %option
# % key: input
# % type: string
# % required: yes
# % multiple: yes
# % label: Sentinel-1 scenes
# % description: List of Sentinel-1 scenes to process
# %end

# %option
# % key: output
# % type: string
# % required: yes
# % multiple: no
# % label: Name of output STRDS
# %end

# %option
# % key: output_title
# % type: string
# % required: no
# % multiple: no
# % label: Title of output STRDS
# %end

# %option
# % key: output_desc
# % type: string
# % required: no
# % multiple: no
# % label: Description of output STRDS
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
# % key: hand_max
# % type: string
# % required: no
# % multiple: no
# % label: Maximum HAND value (meters above drainage network) where water areas are assumed
# % answer: 15.0
# %end

# %option G_OPT_MEMORYMB
# % label: Available memory for running i.segment
# %end

# %option G_OPT_M_NPROCS
# % description: Number of cores for multiprocessing, -2 is n_cores-1
# % answer: 1
# % guisection: Optional
# %end


import grass.script as grass
import os
import psutil
import multiprocessing as mp
from grass.pygrass.modules import Module, ParallelModuleQueue
import atexit

rm_mapsets = []


def cleanup():
    env = grass.gisenv()
    gisdbase = env["GISDBASE"]
    location = env["LOCATION_NAME"]
    for mapset in rm_mapsets:
        grass.utils.try_rmdir(os.path.join(gisdbase, location, mapset))


def test_nprocs_memory():
    # Test nprocs settings
    nprocs = int(options["nprocs"])
    nprocs_real = mp.cpu_count()
    if nprocs > nprocs_real:
        grass.warning(
            "Using %d parallel processes because only %d CPUs available."
            % (nprocs_real, nprocs_real)
        )
        options["nprocs"] = str(nprocs_real)
    # check momory
    memory = int(options["memory"])
    free_ram = freeRAM("MB", 100)
    if free_ram < memory:
        grass.warning("Using %d MB but only %d MB RAM available." % (memory, free_ram))
        options["memory"] = free_ram
        grass.warning("Set used memory to %d MB." % (options["memory"]))


def freeRAM(unit, percent=100):
    """The function gives the amount of the percentages of the installed RAM.
    Args:
        unit(string): 'GB' or 'MB'
        percent(int): number of percent which shoud be used of the free RAM
                      default 100%
    Returns:
        memory_MB_percent/memory_GB_percent(int): percent of the free RAM in
                                                  MB or GB

    """
    # use psutil cause of alpine busybox free version for RAM/SWAP usage
    tot_m = psutil.virtual_memory().total  # in Bytes
    swap_tot_m = psutil.swap_memory().total  # in Bytes
    memory_GB = (tot_m - swap_tot_m) / 1024.0 / 1024.0 / 1024.0
    memory_MB = (tot_m - swap_tot_m) / 1024.0 / 1024.0

    if unit == "MB":
        memory_MB_percent = memory_MB * percent / 100.0
        return int(round(memory_MB_percent))
    elif unit == "GB":
        memory_GB_percent = memory_GB * percent / 100.0
        return int(round(memory_GB_percent))
    else:
        grass.fatal("Memory unit %s not supported" % unit)


def main():
    global rm_mapsets
    # get the user parameters
    input = options["input"].split(",")
    output = options["output"]
    title = options["output_title"]
    desc = options["output_desc"]
    hand_rast = options["hand_rast"]
    hand_max = options["hand_max"]
    test_nprocs_memory()

    if not grass.find_program("i.sentinel_1.water", "--help"):
        grass.fatal(
            _("The 'i.sentinel_1.water' module was not found, install it first:")
            + "\n"
            + "g.extension i.sentinel_1.water url=path/to/addon"
        )

    if not grass.find_program("i.sentinel_1.water.worker", "--help"):
        grass.fatal(
            _("The 'i.sentinel_1.water.worker' module was not found, install it first:")
            + "\n"
            + "g.extension i.sentinel_1.water.worker url=path/to/addon"
        )

    # save current mapset
    env = grass.gisenv()
    start_cur_mapset = env["MAPSET"]

    # create STRDS
    grass.message(_("Creating space time dataset %s..." % (output)))
    grass.run_command(
        "t.create",
        output=output,
        title=title if title else "",
        description=desc if desc else "",
        quiet=True,
        overwrite=True,
    )

    # loop through scenes
    grass.message(_("Creating water masks for input rasters"))

    watermasks = []
    mapsetids = []

    queue_water = ParallelModuleQueue(nprocs=options["nprocs"])
    memory = round(int(options["memory"]) / int(options["nprocs"]))

    for scene in input:
        scene_split = scene.split("_")
        scene_date_raw = scene_split[4].split("T")[0]
        pol = scene_split[10]
        watermask = "S1_watermask_%s_%s_%s" % (scene_date_raw, scene_split[8], pol)
        mapsetid = "tmp_mapset_s1.water_%s_%s_%s_%s" % (
            scene_date_raw,
            scene_split[8],
            pol,
            os.getpid(),
        )
        mapsetids.append(mapsetid)
        rm_mapsets.append(mapsetid)
        grass.message(_("Processing scene %s..." % (scene)))
        params = {
            "input": scene,
            "output": watermask,
            "mapsetid": mapsetid,
            "hand_rast": hand_rast,
            "hand_max": hand_max,
            "memory": memory,
            "pol": pol,
        }

        # copy MASK if there is one
        if grass.find_file(name="MASK", element="raster")["file"]:
            params["flags"] = "m"

        s1_water_worker = Module("i.sentinel_1.water.worker", run_=False, **params)
        watermasks.append(watermask)
        queue_water.put(s1_water_worker)
    queue_water.wait()

    # verify that switching the mapset worked
    env = grass.gisenv()
    cur_mapset = env["MAPSET"]
    if cur_mapset != start_cur_mapset:
        grass.fatal(
            "New mapset is %s, but should be %s" % (cur_mapset, start_cur_mapset)
        )

    # copy maps to current mapset
    for mapsetid in mapsetids:
        for rast in grass.parse_command("g.list", type="raster", mapset=mapsetid):
            grass.run_command(
                "g.copy", raster="%s@%s,%s" % (rast, mapsetid, rast), quiet=True
            )

            scene_date_raw = rast.split("_")[2]
            scene_date = "%s-%s-%s" % (
                scene_date_raw[:4],
                scene_date_raw[4:6],
                scene_date_raw[6:],
            )
            # register in STRDS
            grass.run_command(
                "t.register",
                input=output,
                maps=rast,
                start=scene_date,
                quiet=True,
                overwrite=True,
            )
            grass.message(
                _(
                    "Registered raster map <%s> in space time dataset %s"
                    % (rast, output)
                )
            )


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
