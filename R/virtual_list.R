#' @title Virtual floristic list (VFL)
#'
#' @description A list of species potentially occurring within a study site, in which a probability of occurrence is computed for every taxon
#'
#' @param data_flor dataframe having 5 columns, namely ‘Taxon’ (species identity), ‘Long’ (longitude coordinates), ‘Lat’ (latitude coordinates), ‘uncertainty’ (radius of uncertainty, in metres), and ‘year’ (year of the record)
#' @param site a layer object of class ‘SpatialPolygonsDataFrame’ representing the study area, having CRS: +init=epsg:4326
#' @param year_study the present-year in which you perform the analysis
#' @param excl_areas a layer object of class ‘SpatialPolygonsDataFrame’ to delimit certainly unsuitable areas adjacent or within the study area, having CRS: +init=epsg:4326
#' @param CRS.new the new Coordinate Reference System. Note: must be in XXXXXXXXXXXXlll
#' @param tau percentual value of taxa loss in 100 years time-span (see below for further details)
#' @param upperlimit maximum number of floristic records for a given taxon managed by the algorithm to compute the VFL. You can increase this number; a number between 20 and 30 is suggested. Then, the algorithm considers only the n most relevant, in terms of probability, occurrence records
#' 
#' @return an object of class list with X slots:
#' \itemize{
##'  \item{"The virtual floristic list"}{A dataframe}
##'  \item{"Statistics"}{bla bla bla}
##' }
#' @export
#' @examples \dontrun{
#'
#' data(floratus)
#' data(park)
#' data(unsuitablezone)
#' 
#' virtual_list(data_flor=datashort, site=site, excl_areas=unsuitablezone, tau = 30, upperlimit = 25)
#' 
#' }

virtual_list<-function(data_flor,site,year_study=NULL,excl_areas=NULL, CRS.new = 3035, tau, upperlimit=20){


  ################## Check for settings ##################
  if (length(year_study) == 0) 
  {
    year_study <- Sys.Date()
    year_study <- as.numeric(substr(year_study, start = 1, stop = 4))
  }
  
  if (max(data_flor$year) > year_study) 
  {
    stop("Some occurrence dates are more recent than the year of the study")
  }

  if (tau < 0 | tau >= 100) 
  {
    stop(" 0 <= tau < 100 is FALSE. Please set up another tau value")
  }
  
# Preliminary steps
  
start_time <- Sys.time() ## starting time
raster::crs(site) <- sp::CRS("+init=epsg:4326")
CRS.new <- paste0("+init=epsg:",CRS.new)
message(paste0("Chosen Coordinate Reference System:", " ", CRS.new))

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
      only_spatiotemporal_nintersect <- unlist(only_spatiotemporal_nintersect)
      
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
    only_spatiotemporal_nintersect <- unlist(only_spatiotemporal_nintersect)
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

if(is.null(excl_areas)==TRUE) {message("No unsuitable areas provided")
                               cont <- 0} else
                                 {message("Unsiuitable areas provided")
                                  cont <- 1}

message("Preparing spatial objects!")

# Create a ‘SpatialPointsdataframe’
data_flor_planar <- data_flor

xy <- data_flor_planar[,c(2,3)]
ttt <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
data_flor_planar <- sp::SpatialPointsDataFrame(coords = xy, data = data_flor_planar,
                               proj4string = ttt)


site_3035 <- sp::spTransform(site, CRS.new)


if(cont==1)
{
  raster::crs(excl_areas) <- sp::CRS("+init=epsg:4326")
  # Crop the shapefile of the seas with study area
  excl_areas <- raster::crop(excl_areas, raster::extent(data_flor_planar))
  excl_areas_3035 <- sp::spTransform(excl_areas, CRS.new) # CRS conversion to new CRS
  sp::plot(site_3035, lwd = 2, border="black", main="Buffers (spatial uncertainty)")
  sp::plot(excl_areas_3035, add=TRUE, col=rgb(1,0,0, 0.2), main="Area of exclusion uploaded!")
   }

data_flor_planar <- sp::spTransform(data_flor_planar, CRS.new)
points_3035 <- data_flor_planar
data_flor_planar <- sp::spTransform(data_flor_planar, CRS.new)
data_flor_planar$lat <- data_flor_planar@coords[,2]
data_flor_planar$long <- data_flor_planar@coords[,1]

# Apply for cycle to taxa having buffer intersecting with the polygon of the study area
data_flor_buffer <- rgeos::gBuffer(data_flor_planar, width=(data_flor_planar$uncertainty), byid=TRUE)
result <- raster::intersect(data_flor_buffer, site_3035)
DF <- as.data.frame(result)

# Intersection between the ‘SpatialPointsDataframe’ and the site layer

points_INS <- raster::intersect(points_3035, site_3035)

# Create a vector with taxa present in the input dataframe
list <- unique(DF$Taxon)
list<- as.vector(list)

included_species <- GISTools::poly.counts(data_flor_planar, site_3035)
number_included_species <- max(included_species)
TA<- sp::geometry(data_flor_planar)
sapply(sp::over(site_3035, TA, returnList = FALSE), length)

# Drafting the VFL
message("Drafting the Virtual Floristic List")
listing_time_START <- Sys.time() # record the starting time of the analysis

# Subsetting the ‘SpatialPolygonDataframe’ with buffers using dataframe 'result' (i.e. select occurrence records which intersect the study area)
subset_vector <- result$Taxon
df_spt <- subset(data_flor_buffer, data_flor_buffer$Taxon %in% subset_vector)

# Excluding the areas covered by sea (only work if the layer is provided)
if (cont==1) {
  df_spt <- df_spt - excl_areas_3035
}

# Measure the area of the buffers
x <- as.vector(unique(df_spt$Taxon))
df_spt$area_buffer <- rgeos::gArea(df_spt, byid=TRUE)
overlayXYT <- raster::intersect(site_3035, df_spt)
overlayXYT$area_intersection = sapply(slot(overlayXYT, "polygons"), slot, "area")

#Plots

#Plot buffers and study area

if(cont==0) {    
  sp::plot(site_3035, lwd = 2, border="red", main="Buffers (spatial uncertainty)")
  sp::plot(df_spt, add=TRUE, col=rgb(0,0,1, 0.2), lty=2) } 
  else {
    par(mfrow=c(1,2))
    sp::plot(site_3035, main="Buffers (spatial uncertainty)")
    sp::plot(data_flor_buffer, add=TRUE, col=rgb(0,0,1, 0.2), lty=2) 
    
    sp::plot(site_3035, main="Intersection (+ area of axclusion)")
    sp::plot(overlayXYT, add=TRUE, col=rgb(0,0,1, 0.2), lty=2)
    sp::plot(excl_areas_3035, add=TRUE, col=rgb(1,0,0, 0.2), lty=2)
    
    }
  par(mfrow=c(1,1))  

  # Spatial probability
  overlayXYT$p_occurrence_spatial <- overlayXYT$area_intersection/overlayXYT$area_buffer

  # Temporal probability
  overlayXYT$p_occurrence_temporal <- (1-(tau/100))^((year_study- overlayXYT$year)/100)

  # Spatiotemporal probability of occurrence
  overlayXYT$p_occurrence_spatiotemporal <- (overlayXYT$p_occurrence_spatial * overlayXYT$p_occurrence_temporal)

cl2 <- parallel::makeCluster(parallel::detectCores()/2) #### to perform a parallel computing
doSNOW::registerDoSNOW(cl2)

pb2 <- txtProgressBar(min = 0, max = length(x), style = 3)
progress <- function(n) setTxtProgressBar(pb2, n)
opts <- base::list(progress = progress)
`%dopar%` <- foreach::`%dopar%`

output <- foreach::foreach(i = 1:length(x), .combine = rbind,
                  .options.snow = opts, .packages= "sp") %dopar% {
                    probsptemp_species(x[i])
                  }

close(pb2)
parallel::stopCluster(cl2)

output <- as.data.frame(output)
output$taxon <- x
colnames(output)[1:4] <- c("Estimated_Spatiotemporal_probability", "Number_of_records", "Max_probability", "Min_probability")
output <- output[,c(5,1,2,3,4)]
rownames(output) <- NULL
output[c("Estimated_Spatiotemporal_probability", "Max_probability", "Min_probability")] <- lapply(output[c("Estimated_Spatiotemporal_probability", "Max_probability", "Min_probability")], function(x) 100 * x)

listing_time_FIN <- Sys.time()

message(paste0("Virtual floristic List drafting time:", round(as.numeric(difftime(time1 = listing_time_FIN, time2 = listing_time_START, units = "mins")), 5), " minutes"))

# FINAL STEPS

#1 Remove rows with 'Estimated spatio-temporal probability' equal to 0
is.num <- sapply(output, is.numeric)
output[is.num] <- lapply(output[is.num], round, 1) # round the values
output2 <- output[(output$Estimated_Spatiotemporal_probability>0),]

#2 Order taxa by alphabetical order
output3 <- output2[order(output2$taxon, -output2$Estimated_Spatiotemporal_probability, -output2$Max_probability),]

#3 Store into environment and save the .csv file

to <- list(VFL = output3, Statistics= data.frame(Study_year= year_study, CRS=CRS.new))
write.csv(output3, row.names=FALSE, "Virtual floristic list.csv")
return(to)

}


