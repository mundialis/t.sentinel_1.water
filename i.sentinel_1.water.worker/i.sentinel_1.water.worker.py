#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      i.sentinel_1.water.worker
# AUTHOR(S):   Guido Riembauer
#
# PURPOSE:     i.sentinel_1.water.worker is a worker addon for the parallel processing
#              of i.sentinel_1.water
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
# % description: i.sentinel_1.water.worker is a worker addon for the parallel processing of i.sentinel_1.water.
# % keyword: temporal
# % keyword: Sentinel-1
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

# %option G_OPT_MEMORYMB
# %end

# %option
# % key: mapsetid
# % type: string
# % required: yes
# % multiple: no
# % key_desc: name
# % description: Name of new mapset to run i.sentinel_1.water
# % guisection: Required
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

# %option
# % key: polarisation
# % type: string
# % required: yes
# % multiple: no
# % label: Polarisation: VV or VH
# % options: VV,VH
# %end

# %flag
# % key: m
# % description: Copy and apply MASK from original mapset
# %end


import os
import shutil
import grass.script as grass


def main():
    # actual mapset, location, ...
    env = grass.gisenv()
    gisdbase = env["GISDBASE"]
    location = env["LOCATION_NAME"]
    old_mapset = env["MAPSET"]

    new_mapset = options["mapsetid"]
    input = options["input"]
    output = options["output"]
    hand_rast = "%s@%s" % (options["hand_rast"], old_mapset)
    hand_max = options["hand_max"]
    memory = options["memory"]
    pol = options["polarisation"]

    # set some common environmental variables, like:
    os.environ.update(
        dict(
            GRASS_COMPRESS_NULLS="1",
            GRASS_COMPRESSOR="LZ4",
            GRASS_MESSAGE_FORMAT="plain",
        )
    )

    if not grass.find_program("i.sentinel_1.water", "--help"):
        grass.fatal(
            _("The 'i.sentinel_1.water' module was not found, install it first:")
            + "\n"
            + "g.extension i.sentinel_1.water url=path/to/addon"
        )

    grass.message("New mapset: %s" % new_mapset)
    grass.utils.try_rmdir(os.path.join(gisdbase, location, new_mapset))

    # create a private GISRC file for each job
    gisrc = os.environ["GISRC"]
    newgisrc = "%s_%s" % (gisrc, str(os.getpid()))
    grass.try_remove(newgisrc)
    shutil.copyfile(gisrc, newgisrc)
    os.environ["GISRC"] = newgisrc

    reg = grass.region()

    # change mapset
    grass.message("GISRC: %s" % os.environ["GISRC"])
    grass.run_command("g.mapset", flags="c", mapset=new_mapset, quiet=True)

    # set region
    del reg["cells"]
    del reg["rows"]
    del reg["cols"]
    del reg["zone"]
    del reg["projection"]
    grass.run_command("g.region", **reg, quiet=True)

    # input raster is copied to current mapset to not mess with group creation
    grass.run_command(
        "g.copy", raster="%s@%s,%s" % (input, old_mapset, input), quiet=True
    )

    # MASK is copied if indicated
    if flags["m"]:
        grass.run_command("g.copy", raster="MASK@%s,MASK" % old_mapset, quiet=True)

    grass.run_command(
        "i.sentinel_1.water",
        input=input,
        output=output,
        hand_rast=hand_rast,
        hand_max=hand_max,
        memory=memory,
        pol=pol,
        overwrite=True,
        quiet=False,
    )

    # cannot remove rasters from different mapset so it is manually removed
    # here
    try:
        grass.run_command("g.remove", type="raster", name=input, flags="f", quiet=True)
        if flags["m"]:
            grass.run_command("r.mask", flags="r", quiet=True)
    except Exception:
        grass.message("Could not remove %s from mapset %s" % (input, new_mapset))

    grass.utils.try_remove(newgisrc)
    os.environ["GISRC"] = gisrc


if __name__ == "__main__":
    options, flags = grass.parser()
    main()
