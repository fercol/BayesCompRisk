\name{BayesCRdata}
\alias{BayesCRdata}
\docType{data}
\title{
Example of cause of death data for BayesCR analysis.
}
\description{
Simulated data on causes of death.
}
\usage{data("BayesCRdat")}
\format{
  A data frame with 400 observations on the following 2 variables:
  \describe{
    \item{\code{age}}{a numeric vector of ages at death.}
    \item{\code{cause}}{a string vector of causes of death.}
  }
}
\details{
The data was simulated assuming three different causes of death with bathtub-logistic cause specific hazards. For each simulated individual, three competing ages at death were randomly drawn for each of the three cause-specific hazards. The invidual age at death and corresponding cause was that of the earliest cause-specific age at death.
}
\examples{
\donttest{
## Load data:
data("BayesCRdata", package = "BayesCR")

## Run analysis:
out <- BayesCR(BayesCRdata, nsim = 4, ncpus = 4, 
               UPDJUMP = TRUE, jumpSD = NULL)
}
}
\keyword{datasets}
