# ignobioR

The huge amount of biodiversity records currently available offers increasing opportunities of data analyses, allowing to improve our knowledge of natural systems and their dynamics. This amount of data poses new challenges for reliable analyses and correct interpretation of results. Indeed, to safely deal with occurrence records we must consider their uncertainty, which can introduce biases within analyses. 

This R package provides an objective framework to explicitly include spatial and temporal uncertainties during the mapping and listing of plant occurrence records (both current and historical) for a given study area; ignobioR package returns the here defined (1) ‘Map of Relative Floristic Ignorance’ (MRFI), which represents the spatial distribution of the lack of floristic knowledge, and a (2) ‘Virtual Floristic List’ (VFL), i.e. a list of taxa potentially occurring in the area, showing a probability of occurrence for each taxon.

# How to download and install the R package
if(!require(devtools)){install.packages("devtools"); library(devtools)} 

install_github("interacquas/ignobioR")

# How to use the R package (in brief):

library(ignobioR)

data(floratus)

data(park)

### MAP OF RELATIVE FLORISTIC IGNORANCE (MRFI)
#### Short example
set.seed(123)

mrfi <- ignorance_map(data_flor= floratus[sample(nrow(floratus), 2000), ],  site=park, tau= 20, cellsize= 2000)

#### Extended example
mrfi <- ignorance_map(data_flor = floratus, excl_areas = unsuitablezone, site = park, tau = 20, cellsize = 2000)

### VIRTUAL FLORISTIC LIST (VFL)

#### Short example
set.seed(123)

vfl <- virtual_list(data_flor= floratus[sample(nrow(floratus), 2000), ], site = park, excl_areas = unsuitablezone, tau = 20, upperlimit = 25)

#### Extended example 
vfl <- virtual_list(data_flor = floratus, site = park, excl_areas = unsuitablezone, tau = 20, upperlimit = 25)


# The Authors

Marco D'Antraccoli, Botanic Garden and Museum, University of Pisa (marco.dantraccoli@unipi.it; https://people.unipi.it/marco_dantraccoli/)

Giuseppe Antonelli

Gianni Bedini, Department of Biology, University of Pisa

Lorenzo Peruzzi, Department of Biology, University of Pisa


# Please check our Vignette
https://interacquas.github.io/ignobioR/
