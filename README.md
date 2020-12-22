# ignobioR

The huge amount of biodiversity records currently available offers increasing opportunities of data analyses, allowing to improve our knowledge of natural systems and their dynamics. This amount of data poses new challenges for reliable analyses and correct interpretation of results. Indeed, to safely deal with occurrence records we must consider their uncertainty, which can introduce biases within analyses. We developed an objective framework coded in R programming language, to explicitly include spatial and temporal uncertainties during the mapping and listing of plant occurrence records (both current and historical) for a given study area. Our workflow returns the here defined (i) ‘Map of Floristic Ignorance’ (MFI), which represents the spatial distribution of the lack of floristic knowledge, and a (ii) ‘Virtual Floristic List’ (VFL), i.e. a list of taxa potentially occurring in the area, showing a probability of occurrence for each taxon. The method here presented can manage a huge amount of occurrence data and to represent floristic ignorance across a study area with a sustainable computational effort. Several parameters can be set up by the user, conferring high flexibility to the method. Uncertainty cannot be avoided, but it may be incorporated into biodiversity analyses through appropriate methodological approaches and innovative spatial representations. This contribution introduces a workflow which pushes forward the analytical capacities to deal with uncertainty in biological occurrence records, allowing to produce more reliable outputs.

# How to download and install the R package
if(!require(devtools)){install.packages("devtools"); library(devtools)} 

install_github("interacquas/ignobioR")

# How to use the R package (in brief):

library(ignobioR)

data(floratus)

data(park)


mfi <- ignorance_map(data_flor=floratus, site=park, tau=20, cellsize=2000) # draft the Map of Floristic Ignorance

vfl <- virtual_list(data_flor=floratus, site=park, tau=20) # draft the Virtual floristic List 


# The Authors

Marco D'Antraccoli, University of Pisa,  Botanic Garden and Museum (marco.dantraccoli@unipi.it)

Gianni Bedini, University of Pisa, Department of Biology

Lorenzo Peruzzi, University of Pisa, Department of Biology


# Please check our Vignette
https://interacquas.github.io/ignobioR/
