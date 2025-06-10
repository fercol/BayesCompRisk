\name{BayesCRout}
\alias{BayesCRout}
\docType{data}
\title{
Output from a BayesCR analysis on a simulated dataset.
}
\description{
This dataset is the output of a BayesCR analysis on the simulated data \code{\link{BayesCRdata}}.
}
\usage{data("BayesCRout")}
\format{
  An object of class BayesCR (see values in \code{\link{BayesCR}}).
}
\examples{
## Load output example:
data("BayesCRout", package = "BayesCR")

## Print results:
summary(BayesCRout, digits = 3)

## Plot traces for survival parameters:
plot(BayesCRout)
}
\keyword{datasets}
