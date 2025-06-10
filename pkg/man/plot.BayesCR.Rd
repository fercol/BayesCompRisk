\name{plot.BayesCR}
\alias{plot.BayesCR}
\alias{plot.BayesCR}
\alias{print.BayesCR}
\alias{summary.BayesCR}

\title{
 Summarizing and plotting Bayesian competing risk model outputs.
}

\description{
 These functions are all generic methods for class \code{BayesCR}.
}

\usage{
 \method{plot}{BayesCR} (x, type = "traces", noCIs = FALSE, logMort = FALSE,  \dots)

\method{print}{BayesCR} (x, \dots)

\method{summary}{BayesCR} (object, \dots)

}

 \arguments{
	\item{x }{An object of class \code{BayesCR}.}

	\item{object }{An object of class \code{BayesCR}.}

	\item{type }{Character string specifying the type of plot, with options \dQuote{\code{traces}}, \dQuote{\code{densities}}, \dQuote{\code{densComp}}, 
	\dQuote{\code{demorates}}, and \dQuote{\code{gof}} (see \code{details}).}

	\item{noCIs }{Logical indicating whether the 95\% credible intervals of survival and mortality should be drawn for plot of \code{type = } \dQuote{\code{demorates}}.}

	\item{logMort }{Logical indicating whether mortalities should be drawn in log scale for plot of \code{type = } \dQuote{\code{demorates}}.}

	\item{\dots }{Additional arguments  passed to functions \code{print}, \code{summary} and \code{plot} (see details).}

}

\details{
FILL UP 
}

\value{
	\item{fill up }{FILL UP}
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{FILL UP}

\examples{
## ----------------------- #
## Check BayesCR outputs:
## ----------------------- #
## Load output example:
data("BayesCRout", package = "BayesCompRisk")

## Print results:
summary(BayesCRout, digits = 3)

## Plot traces for survival parameters:
plot(BayesCRout)

## Plot posterior densities of cause of death parameters:
plot(BayesCRout, type = "densities")

## Plot comparision of parameter posterior between causes of death:
plot(BayesCRout, type = "densComp")

## Plot survival and mortality curves:
plot(BayesCRout, type = "demorates")

## Plot survival and mortality without credible intervals:
plot(BayesCRout, type = 'demorates', noCIs = FALSE)

## Plot survival and mortality without credible intervals and log-mortality:
plot(BayesCRout, type = 'demorates', noCIs = TRUE, logMort = TRUE)

## Goodness of fit plot:
plot(BayesCRout, type = "gof")

## Plot cumulative incidences:
plot(BayesCRout, type = "cumInc")
}

\keyword{methods}

