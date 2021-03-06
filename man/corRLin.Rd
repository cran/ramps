\name{corRLin}
\alias{corRLin}

\title{Linear Spatial Correlation Structure}

\description{
   This function is a constructor for the \code{'corRLin'} class, representing a linear spatial correlation structure. Letting \eqn{r} denote the range, the correlation between two observations a distance \eqn{d < r} apart is \eqn{1-(d/r)}. If \eqn{d \geq r}{d >= r} the correlation is zero.
}

\usage{
   corRLin(value = numeric(0), form = ~ 1,
           metric = c("euclidean", "maximum", "manhattan", "haversine"),
           radius = 3956)
}

\arguments{
   \item{value}{optional numeric \dQuote{range} parameter value for the linear correlation structure, which must be greater than zero.  Defaults to \code{numeric(0)}, which results in a range of 90\% of the minimum distance being assigned to the parameter when \code{object} is initialized.}

   \item{form}{one-sided formula of the form \code{~ S1+...+Sp}, specifying spatial covariates \code{S1} through \code{Sp}. Defaults to \code{~ 1}, which corresponds to using the order of the observations in the data as a covariate.}

   \item{metric}{optional character string specifying the distance metric to be used. The currently available options are \code{"euclidean"} for the root sum-of-squares of distances;   \code{"maximum"} for the maximum difference; \code{"manhattan"} for the sum of the absolute differences; and \code{"haversine"} for the great-circle distance (miles) between longitude/latitude coordinates. Partial matching of arguments is used, so only the first three characters need to be provided. Defaults to \code{"euclidean"}.}

   \item{radius}{radius to be used in the haversine formula for great-circle distance.  Defaults to the Earth's radius of 3,956 miles.}
}

\note{
   When \code{"haversine"} is used as the distance metric, longitude and latitude coordinates must be given as the first and second covariates, respectively, in the formula specification for the \code{form} argument.
}

\value{
   Object of class \code{'corRLin'}, also inheriting from class \code{'corRSpatial'}, representing a linear spatial correlation structure.
}

\references{
   Cressie, N.A.C. (1993), \dQuote{Statistics for Spatial Data}, J. Wiley & Sons.

   Venables, W.N. and Ripley, B.D. (1997) \dQuote{Modern Applied Statistics with S-plus}, 2nd Edition, Springer-Verlag.
}

\author{
   Brian Smith \email{brian-j-smith@uiowa.edu} and Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com}, and Douglas Bates \email{bates@stat.wisc.edu} 
}

\seealso{
   \code{\link{corRClasses}}
}

\examples{
sp1 <- corRLin(form = ~ x + y + z)

spatDat <- data.frame(x = (0:4)/4, y = (0:4)/4)

cs1Lin <- corRLin(1, form = ~ x + y)
cs1Lin <- Initialize(cs1Lin, spatDat)
corMatrix(cs1Lin)

cs2Lin <- corRLin(1, form = ~ x + y, metric = "man")
cs2Lin <- Initialize(cs2Lin, spatDat)
corMatrix(cs2Lin)
}

\keyword{models}
