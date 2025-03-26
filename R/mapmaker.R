#' mapmaker
#'
#' Returns a map with coordinates sourced from NOAA ".bna" data file
#'
#' Generates a map based on .bna file downloaded from:
#' https://gnome.orr.noaa.gov/goods/tools/GSHHS/coast_subset.
#' A NOAA file can be uploaded using: read.delim('filename.bna', sep=',', header=TRUE)
#' If argument pdf = TRUE, a PDF file is generated with aspect ratios
#' based on mid-point latitude of the map.
#' A map will have a resolution given by the .bna file downloaded
#' by the user. For maps crossing the dateline longitude west
#' values are > 180. The map is always produced in the
#' northward orientation. Map aspect ratio is scaled only
#' for pdf outputs (pdf = TRUE) and local R plots (pdf = FALSE)
#' should be used for preview only. Map width < map height and
#' the difference in dimensions (aspect ratio)
#' increases with latitude.
#'
#' @param coords a data.frame or matrix (ncol = 3) with geographic
#'  coordinates based on a NOAA .bna file. The file
#' @param pdf logical (default = FALSE), determines if map is
#' plotted locally or saved as pdf file
#' @param filename  character string (default = 'mymap.pdf')
#' defining name of the pdf file (ignored when pdf = FALSE)
#' @param map.height numerical value (default = 5) defining
#' the height of the map in inches (ignored if pdf = FALSE)
#' @param scale.bar logical (default = TRUE), determines if
#' the scale bar is plotted
#' @param scale.width numerical value (default = NULL),
#' the length of the scale bar (km)
#' @param scale.lwd numerical (default = 1.5), thickness of the
#'  scale bar line
#' @param scale.lat numerical (default = NULL), latitude of the
#' scale bar (ignored if scale.long is not provided)
#' @param scale.long numerical (default = NULL), west longitude
#'  of the scale bar (ignored if scale.long is not provided)
#' @param arr.north logical (default = TRUE), indicates
#' if N arrow should be plotted
#' @param arr.lwd numerical (default = 1.5), thickness of N arrow
#' @param arr.cex numerical (default = 0.8), the font size
#' of 'N' symbol of N arrow
#' @param arr.coord numerical vector (lenght=3) (default = NULL)
#' defining N arrow location (long, start lat, end lat)
#' @param sea.col color (default = 'lightskyblue1') for sea
#' @param coast.col color (default = 'khaki3') for coastline
#' @param land.col color name (default = 'khaki1') for land
#' @param grayscale logical (default = FALSE) overrides
#' all color definitions to plot gray-scale map
#' @param rscript R script (default = NULL) with additional
#' plot statements (e.g., mtext(), points(), etc.)
#'
#' @return a plot (map)
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
#' @examples
#' mapmaker(bahamas)
#'
#' @export

mapmaker <- function(coords, pdf = FALSE, filename='mymap.pdf',
                     map.height = 5, scale.bar = TRUE,
                     scale.width = NULL, scale.lwd = 1.5, scale.lat = NULL,
                     scale.long = NULL, arr.north=TRUE, arr.lwd = 1.5,
                     arr.cex=0.8, arr.coord = NULL, sea.col = 'lightskyblue1',
                     coast.col = 'khaki3', land.col = 'khaki1',
                     grayscale = FALSE, rscript = NULL) {

  if (!(is.data.frame(coords) | is.matrix(coords)))
    stop('error: coords should be data.frame or matrix')
  if (ncol(coords) != 3)
    stop(paste('error: incorrect dimensions, coords have',
               ncol(coords), 'columns, while 3 columns are required'))
  if (sum(is.na(coords[,1:2])) > 0)
    stop('missing values not allowed for long-lat coordinates
         (the first two columns of coord)')

  plgs <- which(!is.na(coords[,3])) # find polygons
  coords <- coords[,-3] # keep coordinates only
  mapratio <- cos((pi/180) * mean(coords[1:4, 2])) # map width coefficient
  long.rng <- diff(range(coords[1:4, 1])) # longitude range
  lat.rng <- diff(range(coords[1:4, 2])) # latitude range
  mapdim <- long.rng / lat.rng # map dimension ratio

  if (grayscale) {
    sea.col <- 'white'; coast.col <- 'gray75'; land.col <- 'gray85'
    col1 <- 'black'; col2 <- 'black'; bg1 <-  'white'; bg2 = 'white'
    collab1 = 'gray30'; collab2 = 'black'
  }

  if (scale.bar) {
    if (length(scale.width) == 0) {
      if (long.rng <= 0.1) scale.width <- 1
      if (0.1 < long.rng & long.rng <= 5) scale.width <- 20
      if (5 < long.rng & long.rng <= 10) scale.width <- 100
      if (10 < long.rng & long.rng <= 20) scale.width <- 200
      if (20 < long.rng & long.rng <= 50) scale.width <- 500
      if (50 < long.rng) scale.width <- 1000
    }
    if (length(scale.long) + length(scale.lat) == 2) {
      bar.l <- scale.width / (111.3 * cos((pi/180) * scale.lat)) # bar length in degrees longitude
      bar.coords <- cbind(c(scale.long, scale.long + bar.l),
                          rep(scale.lat, 2)) # scale bar coordinates
    }
    if (length(scale.long) + length(scale.lat) < 2) {
      scale.lat <- max(coords[1:4, 2]) - 0.06 * lat.rng
      scale.long <- max(coords[1:4, 1]) - 0.15 * lat.rng
      bar.l <- scale.width / (111.3 * cos((pi/180) * scale.lat)) # bar length in degrees longitude
      bar.coords <- cbind(c(scale.long, scale.long + bar.l),
                          rep(scale.lat, 2)) # scale bar coordinates
    }
  }

  if (arr.north) {
    if (length(arr.coord) == 0) {
      arr.x <- max(coords[1:4, 1]) - 0.03*long.rng
      arr.y1 <- max(coords[1:4, 2]) - 0.2*lat.rng
      arr.y2 <- max(coords[1:4, 2]) - 0.08*lat.rng
      arr.coord <- c(arr.x, arr.y1, arr.y2)
    }
  }

  if (pdf) pdf(filename, width=map.height*mapdim*mapratio, height=map.height) # scale plot height using map ratio

  graphics::par(oma = c(0, 0, 0, 0), mar = c(3, 3, 1, 1)) # define figure margins
  graphics::plot(coords[1:4,], type = 'n', xlab = '', ylab = '', axes = F)
  graphics::polygon(coords[1:4,], col = sea.col, border = NA)
  for (i in 1:length(plgs)) {
    if (i < length(plgs)) {
      a1 <- plgs[i] + 1
      a2 <- plgs[i+1] - 1
      graphics::polygon(coords[a1:a2, ], col = land.col, border = coast.col, lwd=0.3)
    }
  }
  graphics::rect(coords[1,1], coords[1,2], coords[3,1], coords[2,2]) # add border
  graphics::axis(2, pos = coords[1,1], las = 1, cex.axis = 0.7,
       tck = -0.02, hadj = 0.6) # y axis
  graphics::axis(1, pos = coords[1,2], cex.axis = 0.7,
       tck = -0.02, padj = -1.5) # x axis
  graphics::mtext(side = 1, line = 1, expression(Longitude~degree))
  graphics::mtext(side = 2, line = 1.5, expression(Latitude~degree))

  if (scale.bar) {
    graphics::points(bar.coords, type='l', lwd=scale.lwd, lend=3, col='black')
    graphics::text(mean(bar.coords[,1]), bar.coords[1,2], paste(scale.width, 'km'),
         cex=0.6, pos=3, offset=0.25)
  }

  if (arr.north) {
    graphics::arrows(arr.coord[1], arr.coord[2] + 0.666*diff(arr.coord[2:3]),
           arr.coord[1], arr.coord[3],
           lend = 3, length = 0.1, lwd = arr.lwd)
    graphics::points(rep(arr.coord[1],2),
           c(arr.coord[2], arr.coord[2] + 0.333*diff(arr.coord[2:3])),
           type = 'l', lwd = arr.lwd, lend = 3)
    graphics::text(arr.coord[1], arr.coord[2] + 0.5*diff(arr.coord[2:3]),
         'N', cex = arr.cex)
  }
  if (length(rscript) > 0) rscript
  if (pdf) grDevices::dev.off()
}
