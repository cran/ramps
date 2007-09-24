\name{corRExpwr}
\alias{corRExpwr}

\title{Powered Exponential Spatial Correlation Structure}

\description{
   This function is a constructor for the \code{corRExpwr} class, representing a powered exponential spatial correlation structure. Letting \eqn{r} denote the range, \eqn{p} the shape, and \eqn{n} the nugget effect, the correlation between two observations a distance \eqn{d} apart is \eqn{\exp(-(d/r)^p)}{exp(-(d/r)^p)} when no nugget effect is present and \eqn{(1-n) \exp(-(d/r)^p)}{(1-n)*exp(-(d/r)^p)} when a nugget effect is assumed.
}

\usage{
   corRExpwr(value = numeric(0), form = ~ 1, nugget = FALSE,
             metric = c("euclidean", "maximum", "manhattan", "haversine"),
             radius = 3956, fixed = FALSE)
}

\arguments{
   \item{value}{optional vector with the parameter values in constrained form. If \code{nugget} is \code{FALSE}, \code{value} can have only two elements, corresponding to the \dQuote{range} and \dQuote{shape} of the powered exponential correlation structure, both of which must be greater than zero. If \code{nugget} is \code{TRUE}, meaning that a nugget effect is present, \code{value} can contain two or three elements, the first two as described previously and the third the \dQuote{nugget effect} (one minus the correlation between two observations taken arbitrarily close together); the first two must be greater than zero and the third must be between zero and one. Defaults to \code{numeric(0)}, which results in a range of 90\% of the minimum distance and a nugget effect of 0.1 being assigned to the parameters when \code{object} is initialized.}

   \item{form}{one sided formula of the form \code{~ S1+...+Sp}, specifying spatial covariates \code{S1} through \code{Sp}. Defaults to \code{~ 1}, which corresponds to using the order of the observations in the data as a covariate, and no groups.}

   \item{nugget}{optional logical value indicating whether a nugget effect is present. Defaults to \code{FALSE}.  This argument exists for consistency with the \code{nlme} library and should be left set at its default value when used in \code{georamps} since the associated model includes a separate measurement error variance parameter.}

   \item{metric}{optional character string specifying the distance metric to be used. The currently available options are \code{"euclidean"} for the root sum-of-squares of distances;   \code{"maximum"} for the maximum difference; \code{"manhattan"} for the sum of the absolute differences; and \code{"haversine"} for the great circle distance (miles) between latitude/longitude coordinates. Partial matching of arguments is used, so only the first three characters need to be provided. Defaults to \code{"euclidean"}.}

   \item{radius}{radius to be used in the haversine formula for great circle distance.  Defaults to the Earth's radius of 3,956 miles.}

   \item{fixed}{optional logical value indicating whether the coefficients should be allowed to vary or be kept fixed at their initial value. This argument exists for consistency with the \code{nlme} library and is ignored in the \code{ramps} algorithm.}
}

\value{
   Object of class \code{corRExpwr}, also inheriting from class \code{corSpatial}, representing a powered exponential spatial correlation structure.
}

\references{
   Cressie, N.A.C. (1993), \dQuote{Statistics for Spatial Data}, J. Wiley & Sons.

   Venables, W.N. and Ripley, B.D. (1997) \dQuote{Modern Applied Statistics with S-plus}, 2nd Edition, Springer-Verlag.
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
sp1 <- corRExpwr(form = ~ x + y + z)

spatDat <- data.frame(x = (0:4)/4, y = (0:4)/4)

cs1Expwr <- corRExpwr(c(1, 1), form = ~ x + y)
cs1Expwr <- Initialize(cs1Expwr, spatDat)
corMatrix(cs1Expwr)

cs2Expwr <- corRExpwr(c(1, 1), form = ~ x + y, metric = "man")
cs2Expwr <- Initialize(cs2Expwr, spatDat)
corMatrix(cs2Expwr)

cs3Expwr <- corRExpwr(c(1, 1, 0.2), form = ~ x + y, nugget = TRUE)
cs3Expwr <- Initialize(cs3Expwr, spatDat)
corMatrix(cs3Expwr)
}

\keyword{models}
