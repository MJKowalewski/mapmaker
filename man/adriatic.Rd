\name{adriatic}
\alias{adriatic}
\docType{data}
\title{adriatic map coordinates}
\description{a dataframe based on a .bna file downloaded from NOAA}
\format{A dataframe.
\describe{a dataframe (ncol=3) with longitude coordinates [,1],
latitude coordinates [,2], and polygon attributes for i polygons [ith,1:3] }
}
\details{A dataframe with longitude and latitude coordinates of polygons
for the Adriatic Sea based on a .bna file downloaded from NOAA.
The first column is longitude value, the second column is latitude,
the third column is NA except for rows describing polygon attributes.}
\keyword{datasets}
\examples{str(adriatic)}
