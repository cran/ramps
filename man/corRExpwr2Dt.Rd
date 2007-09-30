\name{corRExpwr2Dt}
\alias{corRExpwr2Dt}

\title{Non-Separable Temporally Integrated Powered Exponential Spatial Correlation Structure}

\description{
   This function is a constructor for the \code{'corRExpwr2Dt'} class, representing a non-separable spatial correlation structure for temporally integrated measurements. Letting \eqn{rs} denote the spatial range, \eqn{ps} the spatial shape, \eqn{rt} the temporal range, \eqn{lambda} the space-time interaction, and \eqn{n} the nugget effect, the correlation between two observations a distance \eqn{d} apart in space and \eqn{t} in time is \eqn{\exp(-(d/rs)^ps - t/rt - \lambda (d/rs)^ps (t/rt))}{exp(-(d/rs)^ps - t/rt - lambda * (d/rs)^ps * (t/rt))} when no nugget effect is present and \eqn{(1-n) \exp(-(d/rs)^ps - t/rt - \lambda (d/rs)^ps (t/rt))}{(1-n)*exp(-(d/rs)^ps - t/rt - lambda * (d/rs)^ps * (t/rt))} when a nugget effect is assumed.
}

\usage{
   corRExpwr2Dt(value = numeric(0), form = ~ 1, nugget = FALSE,
                metric = c("euclidean", "maximum", "manhattan", "haversine"),
                radius = 3956, fixed = FALSE)
}

\arguments{
   \item{value}{optional vector with the parameter values in constrained form. If \code{nugget} is \code{FALSE}, \code{value} can have only four elements, corresponding to the \dQuote{spatial range}, \dQuote{spatial shape}, \dQuote{temporal range}, and \dQuote{space-time interaction} of the powered exponential correlation structure, all of which must be greater than zero. If \code{nugget} is \code{TRUE}, meaning that a nugget effect is present, \code{value} can contain four or five elements, the first four as described previously and the fifth the \dQuote{nugget effect} (one minus the correlation between two observations taken arbitrarily close together); the first four must be greater than zero and the fifth must be between zero and one. Defaults to \code{numeric(0)}, which results in a range of 90\% of the minimum distance and a nugget effect of 0.1 being assigned to the parameters when \code{object} is initialized.}

   \item{form}{one sided formula of the form \code{~ S1+...+Sp+T1+T2}, specifying spatial covariates \code{S1} through \code{Sp} and the times \code{(T1, T2)} at which measurement periods begin and end, respectively.}

   \item{nugget}{optional logical value indicating whether a nugget effect is present. Defaults to \code{FALSE}.  This argument exists for consistency with the \code{nlme} library and should be left set at its default value when used in \code{georamps} since the associated model includes a separate measurement error variance parameter.}

   \item{metric}{optional character string specifying the distance metric to be used. The currently available options are \code{"euclidean"} for the root sum-of-squares of distances;   \code{"maximum"} for the maximum difference; \code{"manhattan"} for the sum of the absolute differences; and \code{"haversine"} for the great circle distance (miles) between latitude/longitude coordinates. Partial matching of arguments is used, so only the first three characters need to be provided. Defaults to \code{"euclidean"}.}

   \item{radius}{radius to be used in the haversine formula for great circle distance.  Defaults to the Earth's radius of 3,956 miles.}

   \item{fixed}{optional logical value indicating whether the coefficients should be allowed to vary or be kept fixed at their initial value. This argument exists for consistency with the \code{nlme} library and is ignored in the \code{ramps} algorithm.}
}

\value{
   Object of class \code{'corRExpwr2Dt'}, also inheriting from class \code{'corSpatial'}, representing a non-separable spatial correlation structure.
}

\references{
   Cressie, N. and Huang, H.-C. (1993) \dQuote{Classes of Nonseperable, Spatio-Temporal Stationary Covariance Functions}, Journal of the American Statistical Association, 94, 1330-1340.

   Smith, B.J. and Oleson, J.J. (2007) \dQuote{Geostatistical Hierarchical Model for Temporally Integrated Radon Measurements}, Jounal of Agricultural, Biological, and Environmental Statistics, in press.
}

\author{
   Brian Smith \email{brian-j-smith@uiowa.edu}
}

\seealso{
   \code{\link{corClasses}},
   \code{\link[nlme]{Initialize.corStruct}},
   \code{\link[nlme]{summary.corStruct}}
}

\examples{
sp1 <- corRExpwr2Dt(form = ~ x + y + t1 + t2)

spatDat <- data.frame(x = (0:4)/4, y = (0:4)/4, t1=(0:4)/4, t2=(1:5)/4)

cs1ExpwrDt <- corRExpwr2Dt(c(1, 1, 1, 1), form = ~ x + y + t1 + t2)
cs1ExpwrDt <- Initialize(cs1ExpwrDt, spatDat)
corMatrix(cs1ExpwrDt)

cs2ExpwrDt <- corRExpwr2Dt(c(1, 1, 1, 1), form = ~ x + y + t1 + t2, metric = "man")
cs2ExpwrDt <- Initialize(cs2ExpwrDt, spatDat)
corMatrix(cs2ExpwrDt)

cs3ExpwrDt <- corRExpwr2Dt(c(1, 1, 1, 1, 0.2), form = ~ x + y + t1 + t2, nugget = TRUE)
cs3ExpwrDt <- Initialize(cs3ExpwrDt, spatDat)
corMatrix(cs3ExpwrDt)
}

\keyword{models}
