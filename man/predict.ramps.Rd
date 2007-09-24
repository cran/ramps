\name{predict.ramps}
\alias{predict.ramps}

\title{Prediction Method for georamps Model Fits}

\description{
   Prediction of latent spatial parameters at unmeasured geographic sites based on a fitted \code{\link{georamps}} model object.
}

\usage{
   \method{predict}{ramps}(object, newdata, \dots)
}

\arguments{
   \item{object}{fitted object returned by \code{georamps}.}
   \item{newdata}{data frame containing unmeasured coordinates at which to predict latent spatial parameters.}
   \item{\dots}{some methods for this generic require additional arguments.  None are used in this method.}
}

\details{
   Prediction will be performed only at the coordinates in \code{newdata} that differ from those used in the initial \code{georamps} model fitting.  In particular, overlapping coordinates will be excluded automatically in the prediction.  Latent spatial parameters are sampled from the posterior predictive distribution with spatial variance set equal to 1.  The scale may be changed by multiplying resulting parameters by the square-root of the desired spatial variance.
}


\value{
   \code{"predict.ramps"} object, inheriting from class \code{"matrix"}, of sampled latent parameters at the new set of coordinates.  Parameter labels are supplied in the returned column names and MCMC iteration numbers in the row names.  A matrix containing the new coordinates is supplied in the \code{coords} attribute of the object.
}

\author{
   Brian Smith \email{brian-j-smith@uiowa.edu}
}

\seealso{
   \code{\link{georamps}}
}

\examples{
## Prediction for georamps example results

\dontrun{
newsites <- data.frame(lon = c(-75.1, -74.7, -74.3),
                       lat = c(39.6, 40.2, 40.8))
predict(NURE.fit, newsites)
}
}

\keyword{models}

