\name{bearing.straits}
\alias{bearing.straits}
\docType{data}
\title{bearing.straits map coordinates}
\description{a dataframe based on a .bna file downloaded from NOAA}
\format{A dataframe.
\describe{a dataframe (ncol=3) with longitude coordinates [,1],
latitude coordinates [,2], and polygon attributes for i polygons [ith,1:3] }
}
\details{A dataframe with longitude and latitude coordinates of polygons
for the Bearing Straits based on a .bna file downloaded from NOAA.
The first column is longitude value, the second column is latitude,
the third column is NA except for rows describing polygon attributes.
Because the map crosses the dateline, longitudes West are positive values
> 180}
\keyword{datasets}
\examples{str(bearing.straits)}
