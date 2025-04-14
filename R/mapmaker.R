#' mapmaker
#'
#' Draw Maps using GSHHS-NOAA Data
#'
#' Generates a map using GSHHS-NOAA data that can be downloaded at
#' \url{https://gnome.orr.noaa.gov/goods/tools/GSHHS/coast_subset}.
#' The data download as .bna files (Atlas Boundary Files) with
#' shoreline coordinates defined as hierarchically arranged closed polygons.
#'
#'
#' Data need to be uploaded into an R data.frame to create an object
#' acceptable by the mapmaker function.
#'
#' One way to do it is to use read.delim:
#' object <- read.delim('filename.bna', sep=',', header=TRUE)
#'
#' The resolution of a map drawn by mapmaker will depend on
#' the resolution chosen when downloading the .bna file from NOAA website.
#' For maps crossing the dateline longitude west values are > 180.
#' The map is always produced in the northward orientation.
#'
#' @param coords a data.frame or matrix (ncol = 3) based on a
#' GSHHS-NOAA .bna file
#' @param scale.bar logical, determines if the scale bar is plotted
#' (default = TRUE). If scale.bar = FALSE, other scale parameters are ignored
#' @param scale.width numerical value, the length of the scale bar
#' in kilometers (default = NULL)
#' @param scale.lwd numerical value, the scale bar line
#' width (default = 1.5)
#' @param scale.cex numerical value, the font size of
#' scale label (default = 0.7)
#' @param scale.lat numerical value, latitude of the scale bar
#' (default = NULL). Ignored if scale.long not provided
#' @param scale.long numerical value, west longitude endpoint
#' of the scale bar (default = NULL). Ignored if scale.lat not provided.
#' @param arr.north logical, indicates if N arrow is plotted (default = TRUE).
#' When arr.north = FALSE all arrow parameters are ignored.
#' @param arr.lwd numerical value, thickness of N arrow (default = 1.5)
#' @param arr.cex numerical value, the font size of 'N' symbol
#' (default = 0.7)
#' @param arr.coord numerical vector with three values defining
#' N arrow location (long, start lat, end lat) (default = NULL)
#' If not provided, N arrow is plotted in the upper right corner of the map.
#' @param cex.lab numerical value for scaling size of axis labels (default = 0.8)
#' @param axes.lab logical, whether the axes labels should be plotted (default = TRUE)
#' @param xcex numerical value for x-axis labels (default = 0.7)
#' @param xtck numerical value for x-axis ticks (default = -0.02)
#' @param xadj numerical value for x-axis adjustment (default = 0.6)
#' @param ycex numerical value for x-axis labels (default = 0.7)
#' @param ytck numerical value for x-axis ticks (default = -0.02)
#' @param yadj numerical value for x-axis adjustment (default = -1.5)
#' @param sea.col the color used for sea areas (default = 'lightskyblue1')
#' @param coast.col the color used for coastline (default = 'khaki3')
#' @param land.col the color used for land areas (default = 'khaki1')
#' @param grayscale logical (default = FALSE), if the grayscale map
#'  should be drawn (overrides other color definitions)
#' @param transp numerical value for land color transparency (default = 0.5)
#' @param rscript R script (default = NULL) with additional
#' plot statements (e.g., mtext(), points(), etc.)
#'
#' @return a single plot
#'
#' @examples
#' mapmaker(bahamas) # draws a map with arbitrary dimensions
#'
#' map.dim <- mapratio(coords = bahamas, map.width = 4)
#' map.dim # find correct map dimensions
#' par(pin=map.dim)
#' mapmaker(coords = bahamas, grayscale = TRUE, rscript = {
#' points(bahamas.towns[, 2:3], pch = 24, cex = 0.6, col = 'black', bg = 'blue')
#' text(bahamas.towns[, 2:3], cex = 0.5, col = 'blue', bahamas.towns[,1], pos = 2, offset=0.15)
#' points(bahamas.sites[, 2:3], pch = 21, cex = 1.3, col = 'black', bg = 'white')
#' text(bahamas.sites[, 2:3], cex = 0.4, col = 'black', bahamas.sites[,1])
#' })
#'
#' @export

mapmaker <- function(coords, scale.bar = TRUE, scale.width = NULL,
                     scale.lwd = 1.5, scale.cex = 0.7,
                     scale.lat = NULL, scale.long = NULL,
                     arr.north=TRUE, arr.lwd = 1.5,
                     arr.cex=0.7, arr.coord = NULL, cex.lab=0.8,
                     axes.lab = TRUE,
                     xcex = 0.7,  xtck = -0.02, xadj = 0.6,
                     ycex = 0.7, ytck = -0.02, yadj = -1.5,
                     sea.col = 'lightskyblue1',
                     coast.col = 'khaki3', land.col = 'khaki1',
                     grayscale = FALSE, transp = 0.5,
                     rscript = NULL) {

  if (!(is.data.frame(coords) | is.matrix(coords)))
    stop('error: coords should be data.frame or matrix')
  if (ncol(coords) != 3)
    stop(paste('error: incorrect dimensions, coords have',
               ncol(coords), 'columns, while 3 columns are required'))
  if (sum(is.na(coords[ , 1:2])) > 0)
    stop('missing values not allowed for long-lat coordinates
         (the first two columns of coord)')

  plgs <- which(!is.na(coords[ , 3])) # find polygons
  coords <- coords[,-3] # keep coordinates only

  long.rng <- diff(range(coords[1:4, 1])) # longitude range
  lat.rng <- diff(range(coords[1:4, 2])) # latitude range

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
      scale.lat <- max(coords[1:4, 2]) - 0.1 * lat.rng
      scale.long <- max(coords[1:4, 1]) - 0.5 * long.rng
      bar.l <- scale.width / (111.3 * cos((pi/180) * scale.lat)) # bar length in degrees longitude
      bar.coords <- cbind(c(scale.long, scale.long + bar.l),
                          rep(scale.lat, 2)) # scale bar coordinates
    }
  }

  if (arr.north) {
    if (length(arr.coord) == 0) {
      arr.x <- max(coords[1:4, 1]) - 0.03*long.rng
      arr.y1 <- max(coords[1:4, 2]) - 0.15*lat.rng
      arr.y2 <- max(coords[1:4, 2]) - 0.02*lat.rng
      arr.coord <- c(arr.x, arr.y1, arr.y2)
    }
  }

 # graphics::par(pin=map.dims)
  graphics::plot(coords[1:4,], type = 'n', xlab = '', ylab = '', axes = F)
  graphics::polygon(coords[1:4,], col = sea.col, border = NA)
  for (i in 1:length(plgs)) {
    if (i < length(plgs)) {
      a1 <- plgs[i] + 1
      a2 <- plgs[i+1] - 1
      graphics::polygon(coords[a1:a2, ],
                        col = grDevices::adjustcolor(land.col, transp),
                        border = coast.col, lwd=0.3)
    }
  }
  graphics::rect(coords[1,1], coords[1,2], coords[3,1], coords[2,2]) # add border
  graphics::axis(2, pos = coords[1,1], las = 1, cex.axis = xcex,
       tck = xtck, hadj = xadj) # y axis
  graphics::axis(1, pos = coords[1,2], cex.axis = ycex,
       tck = ytck, padj = yadj) # x axis

  if (axes.lab) {
  graphics::mtext(side = 1, line = 1, cex = cex.lab,
                  expression(Longitude~degree))
  graphics::mtext(side = 2, line = 1.5, cex = cex.lab,
                  expression(Latitude~degree))
  }

  if (scale.bar) {
    graphics::points(bar.coords, type='l', lwd=scale.lwd,
                     lend=3, col='black')
    graphics::text(mean(bar.coords[,1]), bar.coords[1,2],
                   paste(scale.width, 'km'),
                   cex=scale.cex, pos=3, offset=0.3)
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
}
