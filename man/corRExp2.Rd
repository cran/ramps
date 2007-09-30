\name{corRExp2}
\alias{corRExp2}

\title{Non-Separable Exponential Spatio-Temporal Correlation Structure}

\description{
   This function is a constructor for the \code{'corRExp2'} class, representing a non-separable spatial correlation structure. Letting \eqn{rs} denote the spatial range, \eqn{rt} the temporal range, \eqn{lambda} the space-time interaction, and \eqn{n} the nugget effect, the correlation between two observations a distance \eqn{d} apart in space and \eqn{t} in time is \eqn{\exp(-d/rs - t/rt - \lambda (d/rs) (t/rt))}{exp(-d/rs - t/rt - lambda * (d/rs) * (t/rt))} when no nugget effect is present and \eqn{(1-n) \exp(-d/rs - t/rt - \lambda (d/rs) (t/rt))}{(1-n)*exp(-d/rs - t/rt - lambda * (d/rs) * (t/rt))} when a nugget effect is assumed.
}

\usage{
   corRExp2(value = numeric(0), form = ~ 1, nugget = FALSE,
            metric = c("euclidean", "maximum", "manhattan", "haversine"),
            radius = 3956, fixed = FALSE)
}

\arguments{
   \item{value}{optional vector with the parameter values in constrained form. If \code{nugget} is \code{FALSE}, \code{value} can have only three elements, corresponding to the \dQuote{spatial range}, \dQuote{temporal range}, and \dQuote{space-time interaction} of the exponential correlation structure, all of which must be greater than zero. If \code{nugget} is \code{TRUE}, meaning that a nugget effect is present, \code{value} can contain three or four elements, the first three as described previously and the fourth the \dQuote{nugget effect} (one minus the correlation between two observations taken arbitrarily close together); the first three must be greater than zero and the fourth must be between zero and one. Defaults to \code{numeric(0)}, which results in a range of 90\% of the minimum distance and a nugget effect of 0.1 being assigned to the parameters when \code{object} is initialized.}

   \item{form}{one sided formula of the form \code{~ S1+...+Sp+T}, specifying spatial covariates \code{S1} through \code{Sp} and the times \code{T} at which measurement were taken.}

   \item{nugget}{optional logical value indicating whether a nugget effect is present. Defaults to \code{FALSE}.  This argument exists for consistency with the \code{nlme} library and should be left set at its default value when used in \code{georamps} since the associated model includes a separate measurement error variance parameter.}

   \item{metric}{optional character string specifying the distance metric to be used. The currently available options are \code{"euclidean"} for the root sum-of-squares of distances;   \code{"maximum"} for the maximum difference; \code{"manhattan"} for the sum of the absolute differences; and \code{"haversine"} for the great circle distance (miles) between latitude/longitude coordinates. Partial matching of arguments is used, so only the first three characters need to be provided. Defaults to \code{"euclidean"}.}

   \item{radius}{radius to be used in the haversine formula for great circle distance.  Defaults to the Earth's radius of 3,956 miles.}

   \item{fixed}{optional logical value indicating whether the coefficients should be allowed to vary or be kept fixed at their initial value. This argument exists for consistency with the \code{nlme} library and is ignored in the \code{ramps} algorithm.}
}

\value{
   Object of class \code{'corRExp2'}, inheriting from class \code{'corSpatioTemporal'}, representing a non-separable spatial correlation structure.
}

\references{
   Cressie, N. and Huang, H.-C. (1993) \dQuote{Classes of Nonseperable, Spatio-Temporal Stationary Covariance Functions}, Journal of the American Statistical Association, 94, 1330-1340.
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
sp1 <- corRExp2(form = ~ x + y + t)

spatDat <- data.frame(x = (0:4)/4, y = (0:4)/4, t=(0:4)/4)

cs1Exp <- corRExp2(c(1, 1, 1), form = ~ x + y + t)
cs1Exp <- Initialize(cs1Exp, spatDat)
corMatrix(cs1Exp)

cs2Exp <- corRExp2(c(1, 1, 1), form = ~ x + y + t, metric = "man")
cs2Exp <- Initialize(cs2Exp, spatDat)
corMatrix(cs2Exp)

cs3Exp <- corRExp2(c(1, 1, 1, 0.2), form = ~ x + y + t, nugget = TRUE)
cs3Exp <- Initialize(cs3Exp, spatDat)
corMatrix(cs3Exp)
}

\keyword{models}
