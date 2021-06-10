#### This is the submitted version of the code; please check on GiThub repository for the updated version ###

if(!require(sf)){install.packages("sf"); library(sf)}
if(!require(raster)){install.packages("raster"); library(raster)}
library(ignobioR)

data(park)
data(floratus)

# Create a grid of a study area
grid <- sf::st_make_grid(st_as_sf(park), cellsize=0.05)
grid <- as_Spatial(grid)
plot(grid)
plot(park, add=TRUE, border="red")

dataset <- list() # create an empty list

for (i in 1:length(grid@plotOrder)) {
  message()
  message(paste0("Cell number  ", i, "of", length(grid@plotOrder))
  
  list <- virtual_list(data_flor=floratus[1:5000,], excl_areas = unsuitablezone,  site=grid[i], tau=20, verbose=FALSE) # draft the Virtual floristic List

  dataset[[i]] <- cbind(list$VFL[,1:2], i)
  
}

df <- do.call(rbind.data.frame, dataset)
colnames(df) <- c("Taxon", "Probability", "Cell number")





