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
#'
#' object <- read.delim('filename.bna', sep=',', header=TRUE)
#'
#'
#' The resolution of a map drawn by mapmaker will depend on
#' the resolution chosen when downloading the .bna file from NOAA website.
#' For maps crossing the dateline longitude west values are > 180.
#' The map is always produced in the northward orientation.
#' The function computes map aspect ratio based on the mid-point
#' latitude of a given map. The aspect ratio only applies to pdf-formatted
#' figures produced by mapmaker function when argument pdf = TRUE.
#' The scaling is not applied When the map is previewed in R plot window
#' (pdf = FALSE).
#'
#' @param coords a data.frame or matrix (ncol = 3) based on a
#' GSHHS-NOAA .bna file
#' @param pdf logical (default = FALSE), if map is saved to an
#' external pdf file or plotted locally
#' @param filename character string (default = 'mymap.pdf')
#' defining the name of the pdf output file (ignored when pdf = FALSE)
#' @param map.height numerical value (default = 5) defining
#' the height of the map in inches (ignored if pdf = FALSE)
#' @param lakes logical (default = FALSE), plot inland water
#' bodies (lakes, etc.) with a distinct color (WARNING: will
#' slow down the function for large files)
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
#' @param sea.col the color used for sea areas (default = 'lightskyblue1')
#' @param coast.col the color used for coastline (default = 'khaki3')
#' @param land.col the color used for land areas (default = 'khaki1')
#' @param lake.col the color used for lakes (default = 'seagreen1')
#' (ignored if lakes = FALSE)
#' @param grayscale logical (default = FALSE), if the grayscale map
#'  should be drawn (overrides other color definitions)
#' @param rscript R script (default = NULL) with additional
#' plot statements (e.g., mtext(), points(), etc.)
#'
#' @return either a plot (pdf = FALSE) or a pdf file (pdf = TRUE)
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
#' @examples
#' mapmaker(bahamas)
#'
#' @export

mapmaker <- function(coords, pdf = FALSE, filename='mymap.pdf',
                     map.height = 5, lakes = FALSE,
                     scale.bar = TRUE, scale.width = NULL,
                     scale.lwd = 1.5, scale.lat = NULL,
                     scale.long = NULL, arr.north=TRUE,
                     arr.lwd = 1.5, arr.cex=0.8, arr.coord = NULL,
                     sea.col = 'lightskyblue1', coast.col = 'khaki3',
                     land.col = 'khaki1', lake.col = 'seagreen1',
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

  if (lakes) {
    pol.rng <- NULL
    for (i in 1:(length(plgs)-1)) {
      ipol <- as.vector(apply(coords[(plgs[i] + 1):(plgs[i+1] - 1),], 2, range))
      pol.rng <- rbind(pol.rng, ipol)
    }
    out2 <- NULL
    for (i in 1:(length(plgs)-1)) {
      out1 <- 0
      apol <- colMeans(coords[(plgs[i] + 1):(plgs[i+1] - 1),])
      for (j in 1:(length(plgs)-1)) {
        if (j != i) {
          longT <- apol[1] < pol.rng[j,2] & apol[1] > pol.rng[j,1]
          latT <- apol[2] < pol.rng[j,4] & apol[2] > pol.rng[j,3]
          if (longT & latT) out1 <- out1 + 1
        }
      }
      out2 <- rbind(out2, cbind(i, out1))
    }
  }

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
      scale.long <- max(coords[1:4, 1]) - 0.5 * lat.rng
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

  if (pdf) pdf(filename, width=map.height*mapdim*(1/mapratio),
               height=map.height) # scale plot height using map ratio
  graphics::par(oma = c(0, 0, 0, 0),
                mar = c(3, 3, 1, 1)) # define figure margins
  graphics::plot(coords[1:4,], type = 'n',
                 xlab = '', ylab = '', axes = F,)
  graphics::polygon(coords[1:4,], col = sea.col, border = NA)
  for (i in 1:length(plgs)) {
    if (i < length(plgs)) {
      a1 <- plgs[i] + 1
      a2 <- plgs[i+1] - 1
      if (!lakes) graphics::polygon(coords[a1:a2, ],
                                   col = land.col,
                                   border = coast.col,
                                   lwd=0.3)
      if (lakes) {
      if (out2[i,2] == 0) graphics::polygon(coords[a1:a2, ],
                                    col = land.col,
                                    border = coast.col,
                                    lwd=0.3)
      if (out2[i,2] > 0) graphics::polygon(coords[a1:a2, ],
                                                    col = lake.col,
                                                    border = coast.col,
                                                    lwd=0.3)
      }
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
