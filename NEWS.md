# ignobioR 2.1.0

## New features

- **`sampleboost()`**: Multi-objective optimisation of field sampling design.
  Generates random configurations of non-overlapping circular plots and scores
  each on NDVI variance, mean MRFI, and spatial dispersion (mean
  nearest-neighbour distance). Any objective can be switched off by setting its
  weight to 0.
- Batched point generation and raster extraction for improved performance.
- Automatic working CRS detection from the NDVI raster (no mandatory
  `CRS.new`).
- Optional file output (CSV field sheet, score table, PDF maps) controlled by
  `output_dir`.

## Dependencies

- Added `tidyterra`, `dplyr`, `rlang` to Imports.

# ignobioR 2.0.0

## Breaking changes

- Switched from `raster`/`sp` to `terra`/`sf` throughout.
- `ignorance_map()` now returns a `terra` `SpatRaster` instead of a `raster`
  `RasterLayer`.
- `virtual_list()` uses `sf` geometries internally.

## New features

- Enhanced PDF reports with quantile and continuous colour scales.
- Coverage-weighted rasterization for more accurate MRFI calculation.
- Added `floratus`, `park`, and `unsuitablezone` example datasets.

## Dependencies

- Removed `raster`, `sp`, `rgdal`, `rgeos`.
- Added `sf`, `terra`, `ggplot2`, `gridExtra`, `scales`, `RColorBrewer`.

# ignobioR 1.0.0

- Initial CRAN-style release.
- `ignorance_map()`: Map of Relative Floristic Ignorance (MRFI).
- `virtual_list()`: Virtual Floristic List (VFL).
