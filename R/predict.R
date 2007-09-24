predict.ramps <- function(object, newdata, ...)
{
   ## Extract coordinates from the supplied data frame
   spf <- terms(getCovariateFormula(object$correlation))
   attr(spf, "intercept") <- 0
   newcoords <- model.matrix(spf, newdata)

   ## Construct matrix of measured and unmeasured sites
   ## Initialize the correlation structure using all sites
   n <- ncol(newcoords)
   val <- merge(cbind(newcoords, 1:nrow(newcoords)),
                cbind(object$coords, 1:nrow(object$coords)),
                by=1:n, all.x=TRUE, all.y=FALSE)
   val <- val[order(val[, n+1]),]
   monitor <- is.na(val[, n+2])

   if (!any(monitor))
      stop("No new sites supplied in 'newdata'.  Spatial prediction is not",
           " supported\n\tat point-source measurement sites or at areal grid",
           " sites used in 'georamps'.")

   sites <- unique.sites(newcoords[monitor,,drop=FALSE])
   correlation <- Initialize(
              eval(object$call[match("correlation", names(object$call))][[1]]),
              data = as.data.frame(rbind(object$coords, sites$coords)))

   ## Extract sampled correlation and variance parameters 
   idx <- 1:ncol(object$params)
   phi <- object$params[, params2phi(idx, object$control), drop=FALSE]
   idx <- params2kappa(idx, object$control)
   sigma.z <- sqrt(object$params[, kappa2kappa.z(idx, object$control), drop=FALSE])

   ## K matrix components needed for mpdensity
   val <- unique.sites(as.matrix(object$kmat))
   xk1mat <- cBind(object$xmat,
                   as(model.matrix(~ factor(val$idx) - 1), "dgCMatrix"))
   k2mat <- as(val$coords, "dgCMatrix")

   ## Constants for the sampling algorithm
   allz <- all(object$control$z$monitor)
   nz <- nrow(k2mat)
   n1 <- nrow(object$coords)
   n2 <- nrow(sites$coords)
   p <- ncol(object$xmat)

   znew <- matrix(0, nrow(object$params), n2)

   cat("MCMC Sampler Progress (N = ", nrow(znew), "):\n", sep="")
   for(i in 1:nrow(znew)) {
      if (allz) {
         z <- k2mat %*% object$z[i,]
      } else {
         mpd <- mpdensity(object$params[i,], object$y, xk1mat, k2mat, object$wmat,
                          object$correlation, object$etype, object$ztype,
                          object$retype, object$weights, object$control)
         BETA <- mpd$betahat + solve(mpd$uXtSiginvX, rnorm(p + nz)) 
         z <- BETA[seq(p + 1, length.out = nz)]
      }

      correlation[] <- unconstrained(correlation, phi[i,])
      R <- corMatrix(correlation)
      KMAT <- k2mat %*% as(Diagonal(x = sigma.z[i, object$ztype]), "sparseMatrix")
      R11 <- as.matrix(tcrossprod(KMAT %*% R[1:n1, 1:n1], KMAT))
      R21 <- R[(n1+1):(n1+n2), 1:n1] %*% t(KMAT)
      R22 <- R[(n1+1):(n1+n2), (n1+1):(n1+n2)]

      uiR11 <- as(backsolve(chol(R11), diag(nrow(R11))), "dtCMatrix")
      val <- R21 %*% uiR11
      znew[i,] <- rmvnorm2(1, tcrossprod(val, uiR11) %*% z,
                              R22 - tcrossprod(val))[1,]

      print.iter(i)
   }
   cat("\n")

   val <- structure(znew, dimnames = list(rownames(object$params),
                                          paste("zp_", 1:n2, sep="")))
   attr(val, "coords") <- sites$coords
   class(val) <- c("predict.ramps", "matrix")

   val
}


print.predict.ramps <- function(x, ...)
{
   attr(x, "coords") <- NULL
   print.default(x[], ...)
}
