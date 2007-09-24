DIC <- function(object, ...)
{
   UseMethod("DIC")
}

DIC.ramps <- function(object, iter, ...)
{
   ## Vector of iterations for the calculation
   iter <- if (missing(iter)) 1:nrow(object$params)
           else match.iter(iter, object$params)

	## Posterior parameter means
   w <- rep(1 / length(iter), length(iter))
   params <- as.vector(crossprod(w, object$params[iter,]))
   beta <- params2beta(params, object$control)
   phi <- params2phi(params, object$control)
   sigma2 <- params2kappa(params, object$control)
   sigma2.e <- kappa2kappa.e(sigma2, object$control)
   sigma2.re <- kappa2kappa.re(sigma2, object$control)
   sigma2.z <- kappa2kappa.z(sigma2, object$control)
   loglik <- as.numeric(crossprod(w, object$loglik[iter]))

   ## Mean structure and measurement variance
   MU <- object$y - object$xmat %*% beta
   SIGMA <- Diagonal(x = sigma2.e[object$etype] / object$weights)

   ## Spatial variance
   KMAT <- object$kmat %*%
                 as(Diagonal(x = sqrt(sigma2.z)[object$ztype]), "sparseMatrix")
   object$correlation[] <- unconstrained(object$correlation, phi) 
   SIGMA <- SIGMA + tcrossprod(KMAT %*% corMatrix(object$correlation), KMAT)

   ## Random effects variance
   if(ncol(object$wmat) > 0) {
      SIGMA <- SIGMA + tcrossprod(object$wmat %*%
               as(Diagonal(x = sqrt(sigma2.re)[object$retype]), "sparseMatrix"))
   }

   ## Log-likelihood evaluated at posterior parameter means
   uiSIGMA <- solve(chol(SIGMA))
   loglik.mu <- -length(object$y) / 2.0 * log(2.0 * pi) +
                   sum(log(diag(uiSIGMA))) -
                   as.numeric(crossprod(crossprod(uiSIGMA, MU))) / 2.0

   ## Deviance information criterion
   -2.0 * (2.0 * loglik - loglik.mu)
}
