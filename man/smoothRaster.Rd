% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolateRaster.R
\name{smoothRaster}
\alias{smoothRaster}
\title{Smooth rasters}
\usage{
smoothRaster(inputRaster, fast = FALSE, ...)
}
\arguments{
\item{inputRaster}{An object of class \code{SpatRaster}}

\item{fast}{A logical operator. Setting to \code{TRUE} triggers use
of \code{fastTps} instead of \code{Tps}.}

\item{...}{For any additional arguments passed to \code{Tps} or \code{fastTps}}
}
\value{
An object of class \code{SpatRaster}
}
\description{
Uses thin plate spline regression from
\code{fields} package to smooth raster values.
}
\details{
Original raster is smoothed using a thin
plate spline. This may be desirable in cases where
the user has a reasonable expectation of spatial autocorrelation,
but observes putative measurement errors in a raster. The user has
the option of choosing \code{fastTps} to speed calculation,
but be advised that this is only an approximation
of a true thin plate spline.
}
\examples{
library(terra)
library(fields)
# Create sample raster
r <- rast(ncol=100, nrow=100)
values(r) <- 1:10000

# Introduce a "bubble"
values(r)[720:725] <- 9999
plot(r)

# Smooth bubble with smoothRaster
fastSmooth <- smoothRaster(r, fast = TRUE, aRange = 10.0)
plot(fastSmooth)

}
\seealso{
\code{\link[fields]{Tps}}, \code{\link[fields]{fastTps}}
}
\keyword{dataPrep}
