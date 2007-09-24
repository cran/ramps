\name{summary.ramps}
\alias{summary.ramps}

\title{Posterior Summaries of georamps Model Fits}

\description{
   Posterior summaries of the sampled model parameters returned by \code{\link{georamps}}.
}

\usage{
   \method{summary}{ramps}(object, iter, \dots)
}

\arguments{
   \item{object}{fitted object returned by \code{georamps}.}
   \item{iter}{optional numerical vector specifying the MCMC iterations to be included in the summaries.  All are included by default.}
   \item{\dots}{additional arguments to be passed to \code{summary.mcmc}.}
}

\value{
   Two sets of summary statistics for each model parameter.  Sample mean, standard deviation, naive standard error of the mean, and time-series-based standard error are included in the first set.  Quantiles are included in the second.
}

\author{
   Brian Smith \email{brian-j-smith@uiowa.edu}
}

\seealso{
   \code{\link{georamps}}
   \code{\link[coda]{summary.mcmc}}
}

\examples{
## Posterior summaries for georamps example results

\dontrun{
summary(NURE.fit)
}
}

\keyword{models}
