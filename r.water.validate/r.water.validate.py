#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.water.validate
# AUTHOR(S):   Guido Riembauer
#
# PURPOSE:     Compares estimated watermap with raster reference
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
# % description: Compares estimated watermap with raster reference.
# % keyword: hydrology
# % keyword: raster
# %End

# %option
# % key: input_reference
# % type: string
# % required: yes
# % multiple: no
# % label: Input raster dataset containing reference water areas
# % description: Non-water areas must be 0, water areas must be 1
# %end

# %option
# % key: input_watermap
# % type: string
# % required: yes
# % multiple: no
# % label: Estimated watermap to be compared with reference
# % description: Non-water areas must be 0, water areas must be 1
# %end


import os
import atexit
import grass.script as grass
import numpy as np

# define lists to cleanup
rm_rasters = []


def cleanup():
    """
    Function that loops through global lists with regions, vectors, rasters,
    and groups and removes all items.
    """

    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmrast in rm_rasters:
        if grass.find_file(name=rmrast, element="raster")["file"]:
            grass.run_command("g.remove", type="raster", name=rmrast, **kwargs)


def main():
    global rm_rasters
    input_ref = options["input_reference"]
    input_water = options["input_watermap"]

    # get number of total water pixels in reference and watermap
    tot_water_ref = int(
        grass.parse_command("r.univar", map=input_ref, flags="g")["sum"]
    )
    tot_water_map = int(
        grass.parse_command("r.univar", map=input_water, flags="g")["sum"]
    )
    correct_pixels = "correct_water_pixels_%s" % os.getpid()
    rm_rasters.append(correct_pixels)
    grass.run_command(
        "r.mapcalc",
        expression="%s = if(%s+%s==2,1,0)" % (correct_pixels, input_ref, input_water),
        quiet=True,
    )
    tot_correct = int(
        grass.parse_command("r.univar", map=correct_pixels, flags="g")["sum"]
    )

    correctness = np.round(tot_correct / tot_water_map * 100, 2)
    completeness = np.round(tot_correct / tot_water_ref * 100, 2)

    grass.message(_("Correctness of input layer: %s %%" % correctness))
    grass.message(_("Completeness of input layer: %s %%" % completeness))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
