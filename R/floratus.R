#' Floristic occurrence records for Tuscany region
#'
#' A dataset containing 80,000 florisitic occurrence records for Tuscany
#' 
#'
#' @format A data frame with 80,000 rows and 5 variables:
#' \describe{
#'   \item{Taxon}{Taxonomic Identity}
#'   \item{Long}{Longitude (in WGS84)}
#'   \item{Lat}{Latitude (in WGS84)}
#'   \item{uncertainty}{radius of the uncertainty buffer rising from the occurrence geographic position, in metres}
#'   \item{year}{the year of the occurrence}
#' }
#' @source \url{http://bot.biologia.unipi.it/wpb/italia/index.html; https://www.gbif.org/}
"floratus"
