% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/virtual_list.R
\name{virtual_list}
\alias{virtual_list}
\title{Virtual floristic list (VFL)}
\usage{
virtual_list(
  data_flor,
  site,
  year_study = NULL,
  excl_areas = NULL,
  CRS.new = 3035,
  tau,
  upperlimit = 20
)
}
\arguments{
\item{data_flor}{dataframe having 5 columns, namely ‘Taxon’ (species identity), ‘Long’ (longitude coordinates), ‘Lat’ (latitude coordinates), ‘uncertainty’ (radius of uncertainty, in metres), and ‘year’ (year of the record)}

\item{site}{a layer object of class ‘SpatialPolygonsDataFrame’ representing the study area, having CRS: +init=epsg:4326}

\item{year_study}{the present-year in which you perform the analysis}

\item{excl_areas}{a layer object of class ‘SpatialPolygonsDataFrame’ to delimit certainly unsuitable areas adjacent or within the study area, having CRS: +init=epsg:4326}

\item{CRS.new}{the new Coordinate Reference System. Note: must be in a projected CRS. Default = 3035}

\item{tau}{percentual value of taxa loss in 100 years time-span (see below for further details)}

\item{upperlimit}{maximum number of floristic records for a given taxon managed by the algorithm to compute the VFL. You can increase this number; a number between 20 and 30 is suggested. Then, the algorithm considers only the n most relevant, in terms of probability, occurrence records}
}
\value{
an object of class list with 2 slots:
\itemize{
 \item{"The virtual floristic list"}{A dataframe with five columns}
 \item{"Statistics"}{A set of summarising statistics of the VFL drafting process }
}
}
\description{
A list of species potentially occurring within a study site, in which a probability of occurrence is computed for every taxon
}
\examples{
\dontrun{

data(floratus)
data(park)
data(unsuitablezone)

virtual_list(data_flor = floratus, site = park, excl_areas = unsuitablezone, tau = 30, upperlimit = 25)

}
}
