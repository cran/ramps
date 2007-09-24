\name{DIC}
\alias{DIC}
\alias{DIC.ramps}

\title{Deviance Information Criterion (DIC)}

\description{
   Function to calculate DIC for fitted model objects returned by \code{\link{georamps}}.
}

\usage{
   \method{DIC}{ramps}(object, iter, \dots)
}

\arguments{
   \item{object}{fitted model object from georamps.}
   \item{iter}{optional numerical vector specifying the MCMC iterations to be included in the calculation.  All are included by default.}
   \item{\dots}{some methods for this generic require additional arguments.  None are used in this method.}
}

\details{
}

\value{The numerical DIC value.}

\references{
   Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der Linde, A. (2002) \dQuote{Bayesian Measures of Model Complexity and Fit}, Journal of the Royal Statistical Society - Series B, 64, 583-639.
}

\author{
   Brian Smith \email{brian-j-smith@uiowa.edu}
}

\seealso{
   \code{\link{georamps}}
}

\examples{
## DIC calculation for georamps example results

\dontrun{
DIC(NURE.fit)
}
}

\keyword{models}
