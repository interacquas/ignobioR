# ignobioR

Spiegare in poche righe cosa fa questo pacchetto

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

Marco D'Antraccoli, Pisa Botanic garden and Museum (marco.dantraccoli@unipi.it)

Gianni Bedini, Department of Biology

Lorenzo Peruzzi, Department of Biology


# Please check our Vignette
https://interacquas.github.io/ignobioR/
