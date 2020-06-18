#' @title Draft the virtual floristic list of a given territory
#'
#' @description Bla bla bla
#' @param site,
#'        tau
#' @return A dataframe
#' @export
#' @examples
#' \dontrun{virtual_list()}


virtual_list <- function(data_flor, site, year_study, excl_areas=NULL, CRS.new, tau, upperlimit) {

start_time <- Sys.time() ## starting time
raster::crs(site) <- sp::CRS("+init=epsg:4326")
CRS.new <- paste0("+init=epsg:",CRS.new)
print(paste0("Chosen Coordinate Reference System:", " ", CRS.new))

probsptemp_species <- function(specie) {
  only_spatiotemporal_nintersect <- c(NA)
  site_species <- overlayXYT[which(overlayXYT$Taxon == specie),]
  v <- site_species$p_occurrence_spatiotemporal
  l <- length(v)

  if (l == 1) {
    data <- c(v[1], l, max(v), min(v))

    return(data)
  }

  if (l == 0) {
    data <- c(0, l, 0, 0)
    return(data)}

  if(max(v) >= 0.99) {
    data <- c(1, l)
    data <- c(1, l, max(v), min(v))
    return(data)}

  if (l > 1 & l<= upperlimit) {

    only_spatiotemporal_sum <- sum(v)
    only_spatiotemporal_nintersect <- c()

    for (i in 2:l) {

      only_spatiotemporal_nintersect[[i-1]] <- sum(combn(v, m = i, FUN = prod))
    }

    factor_inclexcl <- c() #### apply the 'inclusione-exclusion’ principle

    for(i in 1:length(only_spatiotemporal_nintersect)) {
      if(floor(i/2)==i/2) {factor_inclexcl[i] = 1} else {factor_inclexcl[i]=-1}
    }

    only_spatiotemporal_nintersect_INCLEXCL <- factor_inclexcl * only_spatiotemporal_nintersect
    elements_intersect <- sum(only_spatiotemporal_nintersect_INCLEXCL)

    only_spatiotemporal_INCLEXCL <- sum(only_spatiotemporal_sum, elements_intersect)

    data <- c(only_spatiotemporal_INCLEXCL, length(v), max(v), min(v))
    return(data)
  }

  v2 <- sort(v, decreasing = TRUE)
  v2 <- v[1:upperlimit]

  l <- length(v2)
  only_spatiotemporal_sum <- sum(v2)

  for (i in 2:l) {
    only_spatiotemporal_nintersect[[i-1]] <- sum(combn(v2, m = i, FUN = prod))
  }

  factor_inclexcl <- c() #### apply the 'inclusion-exclusion’ principle

  for(i in 1:length(only_spatiotemporal_nintersect)) {
    if(floor(i/2)==i/2) {factor_inclexcl[i] = 1} else {factor_inclexcl[i]=-1}
  }

  only_spatiotemporal_nintersect_INCLEXCL <- factor_inclexcl * only_spatiotemporal_nintersect
  elements_intersect <- sum(only_spatiotemporal_nintersect_INCLEXCL)
  only_spatiotemporal_INCLEXCL <- sum(only_spatiotemporal_sum, elements_intersect)
  data <- c(only_spatiotemporal_INCLEXCL, length(v), max(v2), min(v2))
  return(data)
} # Define the core function

if(is.null(excl_areas)==TRUE) {print("No unsuitable areas provided")
                               cont <- 0} else
                                 {print("Unsiuitable areas provided")
                                  cont <- 1}

print("Preparing spatial objects!")
# Create a ‘SpatialPointsdataframe’
data_flor_planar <- data_flor
coordinates(data_flor_planar) <- ~ Long+Lat
proj4string(data_flor_planar) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

if(cont==1)
{
  crs(excl_areas) <- CRS("+init=epsg:4326")
  # Crop the shapefile of the seas with study area
  excl_areas <- crop(excl_areas, extent(data_flor_planar))
  excl_areas_3035 <- spTransform(excl_areas, CRS.new) # CRS conversion to new CRS
  plot(excl_areas_3035, main="Area of exclusion uploaded!")
  }


data_flor_planar <- spTransform(data_flor_planar, CRS.new)
points_3035 <- data_flor_planar
site_3035 <- spTransform(site, CRS.new)
data_flor_planar <- spTransform(data_flor_planar, CRS.new)
data_flor_planar$lat <- data_flor_planar@coords[,2]
data_flor_planar$long <- data_flor_planar@coords[,1]

# Apply for cycle to taxa having buffer intersecting with the polygon of the study area
data_flor_buffer <- gBuffer(data_flor_planar, width=(data_flor_planar$uncertainty), byid=TRUE)
result <- raster::intersect(data_flor_buffer, site_3035)
DF <- as.data.frame(result)

# Intersection between the ‘SpatialPointsDataframe’ and the site layer

points_INS <- raster::intersect(points_3035, site_3035)


# Create a vector with taxa present in the input dataframe
list <- unique(DF$Taxon)
list<- as.vector(list)

included_species <- GISTools::poly.counts(data_flor_planar, site_3035)
number_included_species <- max(included_species)
sapply(over(site_3035, geometry(data_flor_planar), returnList = FALSE), length)


# Drafting the VFL
print("Drafting the Virtual Floristic List")
listing_time_START <- Sys.time() # record the starting time of the analysis

# Subsetting the ‘SpatialPolygonDataframe’ with buffers using dataframe 'result' (i.e. select occurrence records which intersect the study area)
subset_vector <- result$Taxon
df_spt <- subset(data_flor_buffer, data_flor_buffer$Taxon %in% subset_vector)

# Excluding the areas covered by sea (only work if the layer is provided)
if (cont==1) {
  df_spt <- df_spt - excl_areas_3035
}


#Plot buffers and study area
plot(site_3035)
plot(df_spt, add=TRUE)

# Measure the area of the buffers
x <- as.vector(unique(df_spt$Taxon))
df_spt$area_buffer <- gArea(df_spt, byid=TRUE)
overlayXYT <- raster::intersect(site_3035, df_spt)
overlayXYT$area_intersection = sapply(slot(overlayXYT, "polygons"), slot, "area")

#Plot buffers and study area
plot(site_3035)
plot(overlayXYT, add=TRUE)

# Spatial probability
overlayXYT$p_occurrence_spatial <- overlayXYT$area_intersection/overlayXYT$area_buffer

# Temporal probability
overlayXYT$p_occurrence_temporal <- (1-(tau/100))^((year_study- overlayXYT$year)/100)

# Spatiotemporal probability of occurrence
overlayXYT$p_occurrence_spatiotemporal <- (overlayXYT$p_occurrence_spatial * overlayXYT$p_occurrence_temporal)

cl2 <- makeCluster(detectCores()/2) #### to perform a parallel computing
registerDoSNOW(cl2)

pb2 <- txtProgressBar(min = 0, max = length(x), style = 3)
progress <- function(n) setTxtProgressBar(pb2, n)
opts <- base::list(progress = progress)

output <- foreach(i = 1:length(x), .combine = rbind,
                  .options.snow = opts, .packages= "sp") %dopar% {
                    probsptemp_species(x[i])
                  }

close(pb2)
stopCluster(cl2)

output <- as.data.frame(output)
output$taxon <- x
colnames(output)[1:4] <- c("Estimated_Spatiotemporal_probability", "Number_of_records", "Max_probability", "Min_probability")
output <- output[,c(5,1,2,3,4)]
rownames(output) <- NULL
output[c("Estimated_Spatiotemporal_probability", "Max_probability", "Min_probability")] <- lapply(output[c("Estimated_Spatiotemporal_probability", "Max_probability", "Min_probability")], function(x) 100 * x)

listing_time_FIN <- Sys.time()

print(paste0("Virtual floristic List drafting time:", round(as.numeric(difftime(time1 = listing_time_FIN, time2 = listing_time_START, units = "mins")), 5), " ????"))

# FINAL STEPS

#1 Remove rows with 'Estimated spatiotemporal probability' equal to 0
output2<-output[!(output$Estimated_Spatiotemporal_probability==0),]

#2 Order by decreasing 'Estimated spatiotemporal probability'
output3 <- output2[order(-output2$Estimated_Spatiotemporal_probability, -output2$Max_probability, output2$taxon),]
is.num <- sapply(output3, is.numeric)
output3[is.num] <- lapply(output3[is.num], round, 1) # round the values

#3 Store into environment and save the .csv file

to <- list(VFL = output3, Statistics= c(Study_year= year_study, CRS=CRS.new))
write.csv(output3, row.names=FALSE, "Virtual floristic list.csv")
return(to)

}

