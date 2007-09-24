\name{corRExp}
\alias{corRExp}

\title{Exponential Spatial Correlation Structure}

\description{
   This function is a constructor for the \code{corRExp} class, representing an exponential spatial correlation structure. Letting \eqn{r} denote the range and \eqn{n} the nugget effect, the correlation between two observations a distance \eqn{d} apart is \eqn{\exp(-d/r)}{exp(-d/r)} when no nugget effect is present and \eqn{(1-n) \exp(-d/r)}{(1-n)*exp(-d/r)} when a nugget effect is assumed.
}

\usage{
   corRExp(value = numeric(0), form = ~ 1, nugget = FALSE,
           metric = c("euclidean", "maximum", "manhattan", "haversine"),
           radius = 3956, fixed = FALSE)
}

\arguments{
   \item{value}{optional vector with the parameter values in constrained form. If \code{nugget} is \code{FALSE}, \code{value} can have only one element, corresponding to the \dQuote{range} of the exponential correlation structure, which must be greater than zero. If \code{nugget} is \code{TRUE}, meaning that a nugget effect is present, \code{value} can contain one or two elements, the first being the \dQuote{range} and the second the \dQuote{nugget effect} (one minus the correlation between two observations taken arbitrarily close together); the first must be greater than zero and the second must be between zero and one. Defaults to \code{numeric(0)}, which results in a range of 90\% of the minimum distance and a nugget effect of 0.1 being assigned to the parameters when \code{object} is initialized.}

   \item{form}{one sided formula of the form \code{~ S1+...+Sp}, specifying spatial covariates \code{S1} through \code{Sp}. Defaults to \code{~ 1}, which corresponds to using the order of the observations in the data as a covariate, and no groups.}

   \item{nugget}{optional logical value indicating whether a nugget effect is present. Defaults to \code{FALSE}.  This argument exists for consistency with the \code{nlme} library and should be left set at its default value when used in \code{georamps} since the associated model includes a separate measurement error variance parameter.}

   \item{metric}{optional character string specifying the distance metric to be used. The currently available options are \code{"euclidean"} for the root sum-of-squares of distances;   \code{"maximum"} for the maximum difference; \code{"manhattan"} for the sum of the absolute differences; and \code{"haversine"} for the great circle distance between latitude/longitude coordinates. Partial matching of arguments is used, so only the first three characters need to be provided. Defaults to \code{"euclidean"}.}

   \item{radius}{radius to be used in the haversine formula for great circle distance.  Defaults to the Earth's radius of 3,956 miles.}

   \item{fixed}{optional logical value indicating whether the coefficients should be allowed to vary or be kept fixed at their initial value. This argument exists for consistency with the \code{nlme} library and is ignored in the \code{ramps} algorithm.}
}

\value{
   Object of class \code{corRExp}, also inheriting from class \code{corSpatial}, representing an exponential spatial correlation structure.
}

\references{
   Cressie, N.A.C. (1993), \dQuote{Statistics for Spatial Data}, J. Wiley & Sons.

   Venables, W.N. and Ripley, B.D. (1997) \dQuote{Modern Applied Statistics with S-plus}, 2nd Edition, Springer-Verlag.
}

\author{
   Brian Smith \email{brian-j-smith@uiowa.edu} and Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com}, and Douglas Bates \email{bates@stat.wisc.edu} 
}

\seealso{
   \code{\link{corClasses}},
   \code{\link[nlme]{Initialize.corStruct}},
   \code{\link[nlme]{summary.corStruct}}
}

\examples{
sp1 <- corRExp(form = ~ x + y + z)

spatDat <- data.frame(x = (0:4)/4, y = (0:4)/4)

cs1Exp <- corRExp(1, form = ~ x + y)
cs1Exp <- Initialize(cs1Exp, spatDat)
corMatrix(cs1Exp)

cs2Exp <- corRExp(1, form = ~ x + y, metric = "man")
cs2Exp <- Initialize(cs2Exp, spatDat)
corMatrix(cs2Exp)

cs3Exp <- corRExp(c(1, 0.2), form = ~ x + y, nugget = TRUE)
cs3Exp <- Initialize(cs3Exp, spatDat)
corMatrix(cs3Exp)
}

\keyword{models}
