% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualizations.R
\name{diversityStack}
\alias{diversityStack}
\title{Diversity stack}
\usage{
diversityStack(rasterList, template)
}
\arguments{
\item{rasterList}{A \code{list} of \code{SpatRaster} objects, which
are interpreted as species distributions (1 = presence,
0 = absence).}

\item{template}{A \code{SpatRaster} with the desired extent}
}
\value{
A \code{SpatRaster}
}
\description{
Takes list of rasters of species distributions
(interpreted as 1 = presence, 0 = absence), which do not
have to have the same extents, and stack them to create an
estimate of species richness that matches the extent and
resolution of a template.
}
\examples{
library(terra)
rast1 <- rast(ncol=10, nrow=10)
values(rast1) <- rep(0:1, 50)

rast2 <- rast(ncol=10, nrow=10)
values(rast2) <- c(rep(0, 50), rep(1,50))

rastList <- list(rast1, rast2)
result <- diversityStack(rasterList = rastList,
                         template = rast2)
result
plot(result)

}
\keyword{plotting}
