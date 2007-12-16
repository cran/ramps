predict.ramps <- function(object, newdata, ...)
{
   ## Create data frame containing all relevant variables
   mt <- delete.response(object$terms)
   spf <- getCovariateFormula(object$correlation)
   variance <- eval(object$call[["variance"]])
   val <- reformulate(c(labels(mt), all.vars(spf), all.vars(variance$spatial)))
   mfdata <- model.frame(val, newdata, xlev = object$xlevels)

   ## Extract spatial coordinates
   spt <- terms(spf)
   attr(spt, "intercept") <- 0
   newcoords <- model.matrix(spt, mfdata)

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

   # Extract main effects design matrix and spatial variances indices
   idx <- as.numeric(rownames(sites$coords))
   xmat <- as(model.matrix(mt, mfdata[idx,,drop=FALSE],
                           attr(object$xmat, "contrasts")), "dgCMatrix")
   ztype <- if (is.null(variance$spatial)) factor(rep(1, length(idx)))
            else factor(getCovariate(mfdata[idx,,drop=FALSE], variance$spatial),
                        levels = levels(object$ztype))

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

   y <- structure(matrix(0, nrow(object$params), n2),
           dimnames = list(object$control$iter, paste("yp_", 1:n2, sep="")),
           coords = sites$coords,
           class = c("predict.ramps", "matrix"))

   cat("MCMC Sampler Progress (N = ", max(object$control$iter), "):\n", sep="")
   for(i in 1:nrow(y)) {
      if (allz) {
         beta <- params2beta(object$params[i,], object$control)
         z <- k2mat %*% object$z[i,]
      } else {
         mpd <- mpdbetaz(object$params[i,], object$y, xk1mat, k2mat,
                         object$wmat, object$correlation, object$etype,
                         object$ztype, object$retype, object$weights,
                         object$control)
         BETA <- mpd$betahat + solve(mpd$uXtSiginvX, rnorm(p + nz))
         beta <- BETA[seq(length.out = p)] 
         z <- BETA[seq(p + 1, length.out = nz)]
      }

      coef(correlation) <- phi[i,]
      R <- corMatrix(correlation)
      KMAT <- k2mat %*% as(Diagonal(x = sigma.z[i, object$ztype]), "sparseMatrix")
      R11 <- as.matrix(tcrossprod(KMAT %*% R[1:n1, 1:n1], KMAT))
      R21 <- R[(n1+1):(n1+n2), 1:n1] %*% t(KMAT)
      R22 <- R[(n1+1):(n1+n2), (n1+1):(n1+n2)]

      uiR11 <- as(backsolve(chol(R11), diag(nrow(R11))), "dtCMatrix")
      val <- R21 %*% uiR11
      y[i,] <- as.vector(xmat %*% beta) + sqrt(sigma.z[i, ztype]) *
                  rmvnorm2(1, tcrossprod(val, uiR11) %*% z,
                              R22 - tcrossprod(val))[1,]

      print.iter(object$control$iter[i])
   }
   cat("\n")

   y
}


print.predict.ramps <- function(x, ...)
{
   attr(x, "coords") <- NULL
   print.default(x[], ...)
}
