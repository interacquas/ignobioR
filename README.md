# ignobioR

**Version 2.0.0 - Next Generation Floristics toolkit for R**

[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

The `ignobioR` package implements Next Generation Floristics (NGF) methodology to explicitly account for spatial and temporal uncertainties in botanical occurrence records. It provides two core tools:

1. **Map of Relative Floristic Ignorance (MRFI)** - Identifies knowledge gaps across your study area
2. **Virtual Floristic List (VFL)** - Estimates occurrence probabilities for all potentially present taxa

## What's New in 2.0.0

- **Modern spatial packages**: Now uses `sf` and `terra` (faster, more reliable)
- **Enhanced visualizations**: Automatic PDF reports with quantile and continuous scales
- **Improved accuracy**: Coverage-weighted rasterization for precise MRFI calculation

## Installation
```r
# Install from GitHub
if(!require(devtools)) install.packages("devtools")
devtools::install_github("interacquas/ignobioR")
```

## Quick Start
```r
library(ignobioR)

# Load example data
data(floratus)
data(park)
data(unsuitablezone)

# Create Map of Relative Floristic Ignorance
mrfi <- ignorance_map(
  data_flor = floratus,
  site = park,
  excl_areas = unsuitablezone,
  tau = 20,
  cellsize = 2000
)

# Generate Virtual Floristic List
vfl <- virtual_list(
  data_flor = floratus,
  site = park,
  excl_areas = unsuitablezone,
  tau = 20
)
```

## Documentation

Full documentation and vignettes: https://interacquas.github.io/ignobioR/

## Citation

D'Antraccoli, M., Bedini, G., & Peruzzi, L. (2022). Maps of relative floristic ignorance and virtual floristic lists: An R package to incorporate uncertainty in mapping and analysing biodiversity data. *Ecological Informatics*, 67, 101512. https://doi.org/10.1016/j.ecoinf.2021.101512

## Authors

- Marco D'Antraccoli ([University of Pisa](https://people.unipi.it/marco_dantraccoli/))
- Giuseppe Antonelli ([ResearchGate](https://www.researchgate.net/profile/Giuseppe-Antonelli))

## License

GPL-3
