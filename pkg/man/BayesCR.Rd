\name{BayesCR}
\alias{BayesCR}
\alias{BayesCR.default}

\title{
 Parametric Bayesian estimation of age-specific competing risks.
}

\description{
 This function performs multiple Markov Chain Monte Carlo (MCMC) simulations for the Bayesian estimation of age-specific competing risks for cause of death analyses.
}

\usage{
 BayesCR(object, \dots)

\method{BayesCR}{default} (object, niter = 11000, burnin = 1001, thinning = 20,  
         nsim = 4, ncpus = 4, UPDJUMP = TRUE, jumpSD = NULL, \dots)

}

 \arguments{
	\item{object }{A \code{data.frame} to be used as input for the analysis (see \code{details})}

	\item{niter }{The total number of MCMC steps.}

	\item{burnin }{The number of iterations for the burn in (see \code{details}).}

	\item{thinning }{The number of skipped MCMC steps to minimize serial autocorrelation (see \code{details}).}

	\item{nsim }{A numerical value for the number of simulations to be run.}

	\item{ncpus }{a numerical value that indicates the number of cpus to be used for analyses in parallel.}

	\item{UPDJUMP }{A logical argument indicating wheter to update jump standard deviations (adaptive independent Metropolis) until an update rate of 0.25 is achieved (see \code{details}).}

	\item{jumpSD }{Numerical vector of jump standard deviations in case the analysis was carried out previously.}

	\item{\dots }{Additional arguments to be passed to function \code{basta} (see \code{details})}
	
}

\details{
FILL UP 
}

\value{
	  \item{coefficients }{a matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, update rates per parameter, serial autocorrelation on the thinned sequences and the potential scale reduction factor for convergence (see \code{Convergence} value below).}
	  
    \item{x }{a vector of the ages used to calculate \code{mort} and \code{surv}.}
    
    \item{mort }{median and 95\% predictive intervals for the estimated mortality per cause of death and combined (i.e., \dQuote{\code{All}}), called by functions \code{plot} or \code{\link{plot.BayesCR}}.}

    \item{surv }{median and 95\% predictive intervals for the estimated survival per cause of death and combined (i.e., \dQuote{\code{All}}), called by functions \code{plot} or \code{\link{plot.BayesCR}}.}

    \item{cuts }{numerical vector of the index of ages in \code{x} for which the \code{All} survival is larger than or equal to 0.01.}

    \item{theta }{a matrix of the sampled mortality parameters after convergence and thinning.}

    \item{likePost }{a matrix with the log likelihood and log posterior after convergence and thinning.}

    \item{DIC }{basic deviance information criterion (DIC) calculations to be used for model choice (Spiegelhalter \emph{et al.} 2002, Celeux \emph{et al.} 2006). If any of the parameter chains did not converge, then the returned value is \code{NA}.}

    \item{PLE }{list object with the product limit estimator per cause of death and for all combined, called by functions \code{plot} or \code{\link{plot.BayesCR}}.}

    \item{runs }{list object with the outputs of each individual MCMC run, called by functions \code{plot} or \code{\link{plot.BayesCR}}.}

    \item{data }{list object with details of the data provided by the user.}

    \item{settings }{a vector indicating the number of iterations for each MCMC, the burn in sequence, the thinning interval, and the number of simulations that were run, called by \code{summary}.}

    \item{keep }{numerical vector of indices indicating the converged MCMC traces in \code{runs}.}

    \item{params }{nist object with details of the parameters, including initial parameters, prior means and standard deviations, jump standard deviations, lower bounds of the support for each parameter's prior distributions, etc.}

    \item{convergence }{logical indicating whether all parameters converged}
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{summary.BayesCR}}, \code{\link{plot.BayesCR}}.}

\examples{
\donttest{
## Load data:
data("BayesCRdata", package = "BayesCompRisk")

## Run analysis:
out <- BayesCR(BayesCRdata, nsim = 4, ncpus = 4, 
               UPDJUMP = TRUE, jumpSD = NULL)
}
\dontshow{
data("BayesCRout", package = "BayesCompRisk")
out <- BayesCRout
}
## ----------------------- #
## Check BayesCR outputs:
## ----------------------- #
## Print results:
summary(out, digits = 3)

## Plot traces for survival parameters:
plot(out)

## Plot posterior densities of cause of death parameters:
plot(out, type = "densities")

## Plot survival and mortality curves:
plot(out, type = "demorates")
}

\keyword{htest}
\concept{Survival analysis}
\concept{Bayesian inference}

