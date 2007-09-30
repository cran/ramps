\name{corRMatern}
\alias{corRMatern}

\title{Matern Spatial Correlation Structure}

\description{
   This function is a constructor for the \code{'corRMatern'} class, representing a Matern spatial correlation structure. Letting \eqn{r} denote the range, \eqn{s} the shape, and \eqn{n} the nugget effect, the correlation between two observations a distance \eqn{d} apart is \eqn{1/(2^{s-1} \Gamma(s)) (d/r)^s K_s(d/r)}{1/(2^(s-1)*gamma(s))*(d/r)^s*besselK(d/r, s)} when no nugget effect is present and \eqn{(1-n) 1/(2^{s-1} \Gamma(s)) (d/r)^s K_s(d/r)}{(1-n)*1/(2^(s-1)*gamma(s))*(d/r)^s*besselK(d/r, s)} when a nugget effect is assumed.
}

\usage{
   corRMatern(value = numeric(0), form = ~ 1, nugget = FALSE,
              metric = c("euclidean", "maximum", "manhattan", "haversine"),
              radius = 3956, fixed = FALSE)
}

\arguments{
   \item{value}{optional vector with the parameter values in constrained form. If \code{nugget} is \code{FALSE}, \code{value} can have only two elements, corresponding to the \dQuote{range} and \dQuote{scale} of the Matern correlation structure, both of which must be greater than zero. If \code{nugget} is \code{TRUE}, meaning that a nugget effect is present, \code{value} can contain two or three elements, the first two described previously and the third the \dQuote{nugget effect} (one minus the correlation between two observations taken arbitrarily close together); the first two must be greater than zero and the third must be between zero and one. Defaults to \code{numeric(0)}, which results in a range of 90\% of the minimum distance and a nugget effect of 0.1 being assigned to the parameters when \code{object} is initialized.}

   \item{form}{one sided formula of the form \code{~ S1+...+Sp}, specifying spatial covariates \code{S1} through \code{Sp}. Defaults to \code{~ 1}, which corresponds to using the order of the observations in the data as a covariate, and no groups.}

   \item{nugget}{optional logical value indicating whether a nugget effect is present. Defaults to \code{FALSE}.  This argument exists for consistency with the \code{nlme} library and should be left set at its default value when used in \code{georamps} since the associated model includes a separate measurement error variance parameter.}

   \item{metric}{optional character string specifying the distance metric to be used. The currently available options are \code{"euclidean"} for the root sum-of-squares of distances;   \code{"maximum"} for the maximum difference; \code{"manhattan"} for the sum of the absolute differences; and \code{"haversine"} for the great circle distance (miles) between latitude/longitude coordinates. Partial matching of arguments is used, so only the first three characters need to be provided. Defaults to \code{"euclidean"}.}

   \item{radius}{radius to be used in the haversine formula for great circle distance.  Defaults to the Earth's radius of 3,956 miles.}

   \item{fixed}{optional logical value indicating whether the coefficients should be allowed to vary or be kept fixed at their initial value. This argument exists for consistency with the \code{nlme} library and is ignored in the \code{ramps} algorithm.}
}

\value{
   Object of class \code{'corRMatern'}, also inheriting from class \code{'corSpatial'}, representing a Matern spatial correlation structure.
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
sp1 <- corRMatern(form = ~ x + y + z)

spatDat <- data.frame(x = (0:4)/4, y = (0:4)/4)

cs1Matern <- corRMatern(c(1, 1), form = ~ x + y)
cs1Matern <- Initialize(cs1Matern, spatDat)
corMatrix(cs1Matern)

cs2Matern <- corRMatern(c(1, 1), form = ~ x + y, metric = "man")
cs2Matern <- Initialize(cs2Matern, spatDat)
corMatrix(cs2Matern)

cs3Matern <- corRMatern(c(1, 1, 0.2), form = ~ x + y, nugget = TRUE)
cs3Matern <- Initialize(cs3Matern, spatDat)
corMatrix(cs3Matern)
}

\keyword{models}