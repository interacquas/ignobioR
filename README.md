# ignobioR

**Version 2.1.0 — Next Generation Floristics toolkit for R**

[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

The `ignobioR` package implements Next Generation Floristics (NGF) methodology to explicitly account for spatial and temporal uncertainties in botanical occurrence records. It provides three core tools:

1. **Map of Relative Floristic Ignorance (MRFI)** — Identifies knowledge gaps across your study area
2. **Virtual Floristic List (VFL)** — Estimates occurrence probabilities for all potentially present taxa
3. **SampleBoost** — Multi-objective optimisation of field sampling design

## What's New in 2.1.0

- **`sampleboost()`**: New function for optimised placement of sampling plots, maximising environmental heterogeneity (NDVI variance), spatial dispersion (mean nearest-neighbour distance), and optionally prioritising understudied areas (MRFI). Objectives can be switched on/off via weights.
- **Batched sampling engine**: Configurations are generated and scored in blocks for improved performance.

## Installation

```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("interacquas/ignobioR")
```

## Quick Start

### MRFI and VFL

```r
library(ignobioR)

data(floratus)
data(park)
data(unsuitablezone)

# Map of Relative Floristic Ignorance
mrfi <- ignorance_map(
  data_flor = floratus,
  site = park,
  excl_areas = unsuitablezone,
  tau = 20,
  cellsize = 2000
)

# Virtual Floristic List
vfl <- virtual_list(
  data_flor = floratus,
  site = park,
  excl_areas = unsuitablezone,
  tau = 20
)
```

### SampleBoost

```r
# Optimise sampling design using NDVI, ignorance and spatial dispersion
res <- sampleboost(
  ndvi = my_ndvi,
  ignorance = mrfi,
  site = park,
  nplot = 50,
  plot_radius = 5.64,
  perm = 1000,
  seed = 1
)

# Inspect results
res$best_scores
res$plots$ndvi

# Use NDVI + distance only (switch off ignorance)
res2 <- sampleboost(
  ndvi = my_ndvi,
  ignorance = mrfi,
  site = park,
  nplot = 50,
  plot_radius = 5.64,
  perm = 1000,
  seed = 1,
  igno.weight = 0
)

# Write field-ready outputs to disk
res3 <- sampleboost(
  ndvi = my_ndvi,
  ignorance = mrfi,
  site = park,
  nplot = 50,
  plot_radius = 5.64,
  perm = 1000,
  seed = 1,
  output_dir = "output"
)
```

## Documentation

Full documentation and vignettes: https://interacquas.github.io/ignobioR/

## Citation

D'Antraccoli, M., Bedini, G., & Peruzzi, L. (2022). Maps of relative floristic ignorance and virtual floristic lists: An R package to incorporate uncertainty in mapping and analysing biodiversity data. *Ecological Informatics*, 67, 101512. https://doi.org/10.1016/j.ecoinf.2021.101512

## Authors

- Marco D'Antraccoli ([University of Pisa](https://people.unipi.it/marco_dantraccoli/))
- Giuseppe Antonelli ([ResearchGate](https://www.researchgate.net/profile/Giuseppe-Antonelli))
- Gianni Bedini ([University of Pisa](https://people.unipi.it/gianni_bedini/))
- Lorenzo Peruzzi ([University of Pisa](https://people.unipi.it/lorenzo_peruzzi/))

## License

GPL-3