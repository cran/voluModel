% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualizations.R
\name{areColors}
\alias{areColors}
\title{Are Colors}
\usage{
areColors(col)
}
\arguments{
\item{col}{A vector of anything to be interpreted by \code{rgb}
as a color.}
}
\value{
A logical vector stating whether inputs
can be interpreted as colors.
}
\description{
Checks to see if a given vector can be
interpreted by R as a color or colors
}
\examples{

areColors(col = c("red", "prairie_chicken", 2))

}
\keyword{internal}
\keyword{plotting}
