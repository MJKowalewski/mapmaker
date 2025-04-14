#' mapratio
#'
#' Computes map dimensions for drawing maps based on GSHHS-NOAA Data
#'
#' Computes map aspect ratio based on the x-axis (longitudinal) and y-axis
#' (latitudinal) ranges of the map, while also correcting for
#' longitude/latitude aspect ratio using mid-point latitude of the map.
#' Designed for GSHHS-NOAA .bna files (Atlas Boundary Files) with
#' shoreline coordinates defined as hierarchically arranged closed polygons.
#' Files accessible at \url{https://gnome.orr.noaa.gov/goods/tools/GSHHS/coast_subset}
#'
#' Data need to be uploaded into an R data.frame to create an object
#' acceptable by the "mapratio" and "mapmaker" function.
#' Use this: object <- read.delim('filename.bna', sep=',', header=TRUE)
#'
#' The resolution of a map drawn by mapmaker will depend on
#' the resolution chosen when downloading the .bna file from NOAA website.
#' For maps crossing the dateline longitude west values are > 180.
#' The map is always produced in the northward orientation.
#'
#' @param coords a data.frame with 3 columns based on a
#' GSHHS-NOAA .bna file
#'
#' @param map.width a numerical value defining
#' map width in inches (default = 3)
#'
#' @param map.height a numerical value defining map height.
#' If map height is provided, map width value is ignored and
#' the map dimensions are calculated based on height value.
#'
#' @returns a vector with two values representing map dimensions (width, height)
#'
#'
#' @examples
#' mapratio(bahamas)
#' mapratio(bahamas, map.height = 4)
#'
#' @export

mapratio <- function(coords, map.width = 3, map.height = NULL) {
  if (!is.data.frame(coords)) stop('"coords" is not a data.frame')
  if (ncol(coords) != 3) stop('"coords" must have three columns')
  if (!is.numeric(coords[,1]))
    stop('longitude coordinates ("coords[,1]") must be numeric')
  if (!is.numeric(coords[,2]))
    stop('latitude coordinates ("coords[,2]") must be numeric')
  if (sum(is.na(coords[,1:2]) > 0))
    stop('no missing values are alllowed for "coords" coordinates (coords[,1:2])')
  mapratio <- cos((pi/180) * mean(coords[1:4, 2])) # map width coefficient
  long.rng <- diff(range(coords[1:4, 1])) # longitude range
  lat.rng <- diff(range(coords[1:4, 2])) # latitude range
  mapdim <- long.rng / lat.rng # map dimension ratio
  map.dims <- c(map.width, map.width/(mapdim*mapratio)) # map dimensions
  if (length(map.height>0)) {
    map.dims <- c(map.height*(mapdim*mapratio), map.height)
    message(paste0('map.height=', map.height, ' was provided (map.width value will be ignored)'))
  }
  return(map.dims)
}
