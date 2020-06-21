#' @title Ignorance map
#'
#' @description A list of species potentially occurring within a study site, in which a probability of occurrence is computed for every taxon
#'
#' @param data_flor dataframe having 5 columns, namely ‘Taxon’ (species identity), ‘Long’ (longitude coordinates), ‘Lat’ (latitude coordinates), ‘uncertainty’ (radius of uncertainty, in metres), and ‘year’ (year of the record)
#' @param site a layer object of class ‘SpatialPolygonsDataFrame’ representing the study area, having CRS: +init=epsg:4326
#' @param year_study the present-year in which you perform the analysis
#' @param excl_areas a layer object of class ‘SpatialPolygonsDataFrame’ to delimit certainly unsuitable areas adjacent or within the study area, having CRS: +init=epsg:4326
#' @param CRS.new he new Coordinate Reference System. Note: must be in XXXXXXXXXXXXlll
#' @param tau percentual value of taxa loss in 100 years time-span (see below for further details)
#' @param cellsize the resolution of the ignorance map in meters
#' 
#' @return A list with 4 objects:
#' \itemize{
##'  \item{"MFI"}{ The Map of Floristic Ignorance}
##'  \item{"RICH"}{ The corresponding map computed without taking into account spatial and temporal uncertainties}
##'  \item{"uncertainties"}{ The corresponding map computed without taking into account spatial and temporal uncertainties}
##'  \item{"Statistics"}{ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx}
##' }
#' @export
#' @examples \dontrun{
#' data(datashort)
#' data(site)
#' data(exclareas)
#' 
#' ignorance_map(data_flor=datashort, site=site, tau=20, cellsize=10000)
#' }


ignorance_map <- function(data_flor, site, year_study=NULL, excl_areas=NULL, CRS.new=3035, tau, cellsize) {

  ################## Check for settings #############
  if (length(year_study) ==0) 
  {
    year_study <- Sys.Date()
    year_study <- as.numeric(substr(year_study, start = 1, stop = 4))
  }
  
  if (max(data_flor$year) > year_study) 
  {
    stop("Some occurrence dates are more recent than the year of the study")
  }
  
  if (class(site) != "SpatialPolygonsDataFrame" | class(excl_areas) != "SpatialPolygonsDataFrame") 
  {
    print("Layers must be of class SpatialPolygonsDataFrame")
  }
  
  if (tau < 0 | tau >= 100) 
  {
    stop(" 0 <= tau < 100 is FALSE. Please set up another tau value")
  }
  
  if (length(which(2*data_flor$uncertainty < (cellsize/20))) > 1) 
  {print(data_flor[(data_flor$uncertainty * 2) < (cellsize/20),])
    stop("There are ", paste(length(which(2*data_flor$uncertainty < (cellsize/20)))), " occurrence records having an uncertainty value too small in respect to the cell size. They could be lost during the rasterisation process. Digit ‘help(rasterize)' for more details")
  }
  
  # Preliminary steps
  start_time <- Sys.time() ## starting time
  raster::crs(site) <- sp::CRS("+init=epsg:4326")
  CRS.new <- paste0("+init=epsg:", CRS.new)
  print(paste0("Chosen Coordinate Reference System:", " ", CRS.new))
  print(paste0("Chosen tau:", " ", tau))
  
ignorance_species <- function(dfOBJ) {
  ###### Create buffers having radius = 'uncertainty'
  DDF_buffer <- rgeos::gBuffer(dfOBJ, width=(dfOBJ $uncertainty), byid=TRUE)

  if (cont==1) {
    DDF_buffer <- DDF_buffer - excl_areas_3035
    }

  ##### Create the Dataframe

  #### Spatial probability
  xy <- raster::extract(r, DDF_buffer)
  spt_scor <- c()
  for (i in 1:length(xy)) {
    spt_scor[[i]] <- length(xy[[i]])
  }
  DDF_buffer@data$spatial_score <- 1/spt_scor

  #### Temporal probability
  DDF_buffer@data$time_score <- (1-(tau/100))^((year_study - DDF_buffer@data$year)/100)

  ### Spatio-temporal probability
  DDF_buffer@data$st_ignorance <- DDF_buffer@data$spatial_score * DDF_buffer@data$time_score

  ##### Modify the extent
  DDF_buffer@bbox <- as.matrix(raster::extent(r))

  ####### Rasterise
  x <- raster::stack()
  for(i in 1:nrow(DDF_buffer)){
    xxx <- raster::rasterize(DDF_buffer[i,], r2, 'st_ignorance', update=TRUE, getCover=TRUE)
    xxx[xxx > 0] <- DDF_buffer@data[i,"st_ignorance"]
    x <- raster::stack(x, xxx)

  }

  r.polys <- raster::stackApply(x, indices = rep(1, raster::nlayers(x)), fun = max)
  r.polys[is.na(r.polys[])] <- 0
  return(r.polys)
} # define the core function

if(is.null(excl_areas)==TRUE)
  {print("No unsuitable areas provided")
  cont <- 0} else {print("Unsuitable areas provided. Plotting")
                   sp::plot(site)
                   cont <- 1}

print("Creating spatial objects")

# Create a ‘SpatialPointsdataframe’
data_flor_planar <- data_flor

xy <- data_flor_planar[,c(2,3)]
ttt <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
data_flor_planar <- sp::SpatialPointsDataFrame(coords = xy, data = data_flor_planar, proj4string = ttt)

if(cont==1)
{
  raster::crs(excl_areas) <- sp::CRS("+init=epsg:4326")
  # Crop the shapefile of the seas with study area
  excl_areas <- raster::crop(excl_areas, raster::extent(data_flor_planar))
  excl_areas_3035 <- sp::spTransform(excl_areas, CRS.new) # CRS conversion to new CRS
  }

data_flor_planar <- sp::spTransform(data_flor_planar, CRS.new)
points_3035 <- data_flor_planar
site_3035 <- sp::spTransform(site, CRS.new)
print("ok1")
data_flor_planar$lat <- data_flor_planar@coords[,2]
data_flor_planar$long <- data_flor_planar@coords[,1]

# Apply for cycle to taxa having buffer intersecting with the polygon of the study area
data_flor_buffer <- rgeos::gBuffer(data_flor_planar, width=(data_flor_planar$uncertainty), byid=TRUE)
print("ok2")

##### Plot intermediate steps

if(cont==1)
{
  sp::plot(data_flor_buffer)
  sp::plot(site_3035, add=TRUE)
  sp::plot(excl_areas_3035, add =TRUE, border="red",lty =2)
  sp::plot(points_3035, add=TRUE, col="blue")


} else {
  sp::plot(data_flor_buffer)
  sp::plot(site_3035, add=TRUE)
  sp::plot(points_3035, add=TRUE, col="blue")
}

print("Filtering occurrence records with buffers intersecting the study area")
result <- raster::intersect(data_flor_buffer, site_3035)
DF <- as.data.frame(result)

# Intersection between the ‘SpatialPointsDataframe’ and the site layer

points_INS <- raster::intersect(points_3035, site_3035)

# Create a vector with taxa present in the input dataframe
list <- unique(DF$Taxon)
list <- as.vector(list)

# Create an empty raster
r <- raster::raster()
raster::xmin(r) <- min(result@bbox[1,1]) - cellsize
raster::xmax(r) <- max(result@bbox[1,2]) + cellsize
raster::ymin(r) <- min(result@bbox[2,1]) - cellsize
raster::ymax(r) <- max(result@bbox[2,2]) + cellsize
raster::res(r) <- cellsize
raster::crs(r) <- CRS.new
raster::values(r) <- 1:raster::ncell(r)
r2 <- r
r2[]<-NA
print("ok4")
rich <- raster::rasterize(data_flor_planar, r, 'Taxon', function(x, ...) length(unique(na.omit(x))))
rich[is.na(rich)] <- 0
print("ok5")

included_species <- GISTools::poly.counts(data_flor_planar, site_3035)
number_included_species <- max(included_species)
TA2<- sp::geometry(data_flor_planar)
sapply(sp::over(site_3035, TA2, returnList = FALSE), length)


cl <- parallel::makeCluster(parallel::detectCores()-1) #### to perform a parallel computing
doSNOW::registerDoSNOW(cl)

print("Starting to draft the Map of Floristic Ignorance!")
pb <- utils::txtProgressBar(min = 0, max = length(list), style = 3) # progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- base::list(progress = progress)

`%dopar%` <- foreach::`%dopar%`

raster_stack <- foreach::foreach(i= 1:length(list), .options.snow = opts, .packages="raster", .combine=raster::stack) %dopar% {
  yo <- DF[which(DF$Taxon == list[i]),]
  xy <- yo[,c(7,6)]
  ppp <- raster::crs(site_3035)
  df_species <- sp::SpatialPointsDataFrame(coords = xy, data = yo, proj4string = ppp)
  ignorance_species(df_species)
}

close(pb)
parallel::stopCluster(cl)

#####
print("Almost done. Preparing outputs")

#### Sum the single rasters
base::names(raster_stack) <- list
raster_sum <- sum(raster_stack, na.rm = TRUE)

############ Rescale the raster to show IFI
r.max = raster::cellStats(raster_sum, "max")
raster_sum <- r.max - raster_sum

############ Plot the map after a 'mask' operation

raster_sum2 <- raster::crop(raster_sum, r)

x_crop <- raster::crop(raster_sum2, r)

e <- raster::extent(site_3035)
x_crop <- raster::extend(x_crop, raster::extent(e[1]-cellsize, e[2] + cellsize, e[3] - cellsize, e[4]+cellsize), value=raster::maxValue(raster_sum2))

rgdal::writeOGR(site_3035, tempdir(), f <- basename(tempfile()), 'ESRI Shapefile')
gdalUtils::gdal_rasterize(sprintf('%s/%s.shp', tempdir(), f),
               f2 <- tempfile(fileext='.tif'), at=T,
               tr=raster::res(x_crop), te=c(sp::bbox(x_crop)), burn=1,
               init=0, a_nodata=0, ot='Byte')

raster_new <- x_crop*raster::raster(f2)

### Record the time elapsed
end_time <- Sys.time()

#### Create the dataframe storing the descriptive statistics
names <-c("Started", "Finished", "Elapsed time", "CRS", "Cell size (Km)", "100 years % loss ratio (tau)",
          "Occurrence whithin", "Total occurrence computed", "Occurrence uncertanties (metres; median value)", "Occurrence dates (year; median value)")
values <- c(as.character(start_time), as.character(end_time), end_time-start_time, CRS.new, cellsize/1000, tau,
            nrow(points_INS), nrow(DF), round(median(DF$uncertainty)), median(DF$year))
statistics <- as.data.frame(cbind(names, values))
statistics

#### Producing the images to export #####

### Plot n° 1
sp::plot(raster_new, main="raster_new")

test_spdf <- as(raster_new, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

p1 <- ggplot2::ggplot(test_df)+ ggplot2::coord_equal()+ ggplot2::theme_classic()+
  ggplot2::labs(fill="IFI")+
  ggplot2::theme(legend.position="right",legend.direction='vertical')+ ggplot2::theme(legend.key.width=grid::unit(0.6, "cm"))+
  ggplot2::xlab("Latitude") + ggplot2::ylab("Longitude")+
  ggplot2::scale_fill_distiller(palette = "Spectral", direction = -1, guide = ggplot2::guide_legend(),breaks=rev(seq(0, raster::maxValue(raster_new), raster::maxValue(raster_new)/10)),
                                labels=round(rev(seq(0, raster::maxValue(raster_new), raster::maxValue(raster_new)/10))), limits = c(0, raster::maxValue(raster_new)))+
  ggplot2::geom_tile(test_df, mapping = ggplot2::aes(x=.data$x, y=.data$y, fill=.data$value), alpha=0.8)+
  ggplot2::geom_polygon(site_3035, mapping= ggplot2::aes(x=long, y=lat, group=group), fill=NA, color="black", size=1)+
  ggplot2::geom_polygon(data=raster::rasterToPolygons(raster_new), mapping=ggplot2::aes(x=long, y=lat, group=group), color="black", alpha=0)+
  ggplot2::ggtitle("Map of Floristic Ignorance (MFI)")

# Plot n° 2
x_crop_rich <- raster::extend(rich, raster_new, value=0)
sp::plot(x_crop_rich, main="rich")
sp::plot(site_3035, add=TRUE)

rgdal::writeOGR(site_3035, tempdir(), f <- basename(tempfile()), 'ESRI Shapefile')
gdalUtils::gdal_rasterize(sprintf('%s/%s.shp', tempdir(), f),
               f3 <- tempfile(fileext='.tif'), at=T,
               tr=raster::res(x_crop_rich), te=c(sp::bbox(x_crop_rich)), burn=1,
               init=0, a_nodata=0, ot='Byte')

raster_new_rich <- x_crop_rich*raster::raster(f3) # multiply the raster by 1 or NA

test_spdf2 <- as(raster_new_rich, "SpatialPixelsDataFrame")
test_df2 <- as.data.frame(test_spdf2)
colnames(test_df2) <- c("value", "x", "y")

p2 <- ggplot2::ggplot(test_df2) + ggplot2::coord_equal() + ggplot2::theme_classic() +
  ggplot2::theme(legend.position="right", legend.direction='vertical', legend.key.width=grid::unit(0.6, "cm"))+
  ggplot2::geom_tile(data=test_df2, mapping=ggplot2::aes(x=.data$x, y=.data$y, fill=.data$value), alpha=0.8) +
  ggplot2::geom_polygon(data=site_3035, mapping=ggplot2::aes(x=long, y=lat, group=group),fill=NA, color="black", size=1) +
  ggplot2::geom_polygon(data=raster::rasterToPolygons(raster_new_rich), mapping = ggplot2::aes(x=long, y=lat, group=group), color="black", alpha=0)+
  ggplot2::scale_fill_distiller(palette = "Spectral", direction = +1, guide = ggplot2::guide_legend(),breaks=rev(seq(0, raster::maxValue(raster_new_rich), raster::maxValue(raster_new_rich)/10)),
                                labels=round(rev(seq(0, raster::maxValue(raster_new_rich), raster::maxValue(raster_new_rich)/10))), limits = c(0, raster::maxValue(raster_new_rich)))+
  ggplot2::ggtitle("Traditional richness map (without uncertainties)")+
  ggplot2::xlab("Latitude") + ggplot2::ylab("Longitude")


# Plot n° 3
p3 <-  ggplot2::ggplot(DF)+
    ggplot2::aes(x = year, y = ..count../sum(..count..))+
    ggplot2::geom_histogram(alpha=.6, fill="#FF6666", binwidth = diff(range(DF$year))/30)+
    ggplot2::coord_cartesian(xlim = c(min(DF$year), year_study))+
    ggplot2::scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
    ggplot2::ggtitle("Occurrence dates distribution")+
    ggplot2::xlab("Year") + ggplot2::ylab("Density")+
    ggplot2::labs(fill="Number of taxa")+
    ggplot2::theme_classic()


# Plot n° 4

p4 <- ggplot2::ggplot(DF) +
  ggplot2::aes(x = uncertainty, y = ..count../sum(..count..)) +
  ggplot2::geom_histogram(alpha=.6, fill="#FF6666", binwidth = diff(range(DF$uncertainty))/30)+
  ggplot2::coord_cartesian(xlim = c(min(DF$uncertainty), max(DF$uncertainty)))+
  ggplot2::scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  ggplot2::ggtitle("Occurrence uncertainties distribution")+
  ggplot2::xlab("Uncertainty") + ggplot2::ylab("Density")+
  ggplot2::theme_classic()


# Plot n° 5

ss <-  gridExtra::tableGrob(statistics, theme=gridExtra::ttheme_minimal())

# Creating the .pdf file
grDevices::pdf("IgnoranceMap.pdf", onefile = TRUE)
p1
p2
p3
p4
plot(ss)
dev.off()

# Write to file the raster of the ‘Map of Floristic Ignorance’ and a .csv file listing the taxa considered to draft the map
raster::writeRaster(raster_new, filename = "MAPignorance", format="GTiff", overwrite=TRUE)
write.csv(list, row.names=FALSE, "Taxa considered to compute the Floristic Ignorance Map.csv")
print(paste0("Done! The files have been saved here", getwd()))

### Print images
print(p1)
print(p2)
print(p3)
print(p4)
plot(ss)

# Save into a list

to2 <- list(MFI = raster_new, RICH = raster_new_rich, Uncertainties = data.frame(uncertainty= DF$uncertainty, year= DF$year) , Statistics= statistics)

}
