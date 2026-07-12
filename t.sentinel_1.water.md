[![image-alt](grass_logo.png)](https://grass.osgeo.org/grass-stable/manuals/index.html)

------------------------------------------------------------------------

## NAME

***t.sentinel_1.water*** - Toolset for detecting water areas in
Copernicus Sentinel-1 GRD scenes and registration of the watermasks in a
space time dataset.

## KEYWORDS

[imagery](https://grass.osgeo.org/grass-stable/manuals/imagery.html),
[import](https://grass.osgeo.org/grass-stable/manuals/topic_import.html),
[satellite](keywords.md#satellite), [Sentinel](keywords.md#Sentinel)

## DESCRIPTION

The *t.sentinel_1.water* toolset consists of currently four modules:

- [i.sentinel_1.water](i.sentinel_1.water.md): detects water areas in
  Sentinel-1 GRD scenes using automatic or manual thresholding and
  region growing. Water areas are corrected using a
  Height-Above-Nearest-Drainage (HAND) raster.
- [i.sentinel_1.water.worker](i.sentinel_1.water.worker.md): this is a
  worker Addon to run *i.sentinel_1.water* in parallel in a dedicated
  mapset.
- [r.water.validate](r.water.validate.md): compares the estimated
  watermap with a raster map reference.
- [t.sentinel_1.watermasks](t.sentinel_1.watermasks.md): detects water
  areas in Sentinel-1 GRD scenes using
  [i.sentinel_1.water](i.sentinel_1.water.md) and registers the
  watermasks in a space time dataset.

## AUTHOR

Guido Riembauer, [mundialis](https://www.mundialis.de/), Germany
