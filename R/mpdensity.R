################################################################################
## mpdensity/mpdensity2 - Marginalized posterior density, up to a constant scale
##
## Arguments:
##    params  - numerical vector of phi and kappa parameter values
##    y       - numerical response vector (n x 1)
##    xk1mat  - design Matrix cbind(xmat, k1mat) where 'xmat' is the covariate
##              Matrix and 'k1mat' is such that kmat = k1mat %*% k2mat
##    k2mat   - design Matrix giving the unique rows of kmat
##    wmat    - design Matrix of non-spatial random effect (n by q)
##    spcor   - initialized (nlme) spatial correlation structure
##    etype   - factor indexing the measurement variances (n x 1)
##    ztype   - factor indexing the spatial variances (nz x 1)
##    retype  - factor indexing the random effect variances (q x 1)
##    weights - numerical vector by which to weight the measurement error
##              variance
##    control - ramps.control object
################################################################################

mpdensity <- function(params, y, xk1mat, k2mat, wmat, spcor, etype, ztype, retype,
                      weights, control)
{
   spcor[] <- unconstrained(spcor, params2phi(params, control))
   kappa <- params2kappa(params, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   p <- length(control$beta)
   nz <- nrow(k2mat)

   if (ncol(wmat) == 0) {
      uiSIGMA.11 <- as(as(Diagonal(x = sqrt(weights / kappa.e[etype])),
                          "sparseMatrix"), "dtCMatrix")
   } else {
      SIGMA.11 <- Diagonal(x = kappa.e[etype] / weights) + tcrossprod(wmat %*%
                     as(Diagonal(x = sqrt(kappa.re)[retype]), "sparseMatrix"))
      uiSIGMA.11 <- solve(chol(SIGMA.11))
   }

   KMAT <- k2mat %*% as(Diagonal(x = sqrt(kappa.z)[ztype]), "sparseMatrix")
   uSIGMA.22 <- chol(as.matrix(tcrossprod(KMAT %*% corMatrix(spcor), KMAT)))

   logsqrtdet <- sum(log(diag(uSIGMA.22))) - sum(log(diag(uiSIGMA.11)))

   linvX.r1 <- crossprod(uiSIGMA.11, xk1mat)
   linvX.22 <- as(backsolve(uSIGMA.22, diag(-1, nrow(uSIGMA.22)),
                            transpose = TRUE), "dtCMatrix")
   linvX <- rBind(linvX.r1, cBind(Matrix(0, nz, p), linvX.22))
   linvY <- rBind(as(crossprod(uiSIGMA.11, y), "dgCMatrix"), Matrix(0, nz, 1))

   XtSiginvX <- crossprod(linvX)
   uXtSiginvX <- chol(XtSiginvX)
   logsqrtdetXtSiginvX <- sum(log(diag(uXtSiginvX)))
   XtSiginvY <- crossprod(linvX, linvY)

   betahat <- solve(XtSiginvX, XtSiginvY)

   resids <- linvY - linvX %*% betahat
   quadform <- as.numeric(crossprod(resids))

   shape <- sigma2shape(control)
   loglik <- -logsqrtdet - logsqrtdetXtSiginvX -
               (sum(shape) + (length(y) - p) / 2.0) *
               log(quadform / 2.0 + sum(sigma2scale(control) / kappa)) -
               sum((shape + 1.0) * log(kappa))

   list(value = loglik, betahat = betahat, quadform = quadform,
        uXtSiginvX = uXtSiginvX, uSig11inv = uiSIGMA.11)
}


mpdensity2 <- function(params, y, xk1mat, k2mat, wmat, spcor, etype, ztype, retype,
                       weights, control)
{
   spcor[] <- unconstrained(spcor, params2phi(params, control))
   kappa <- params2kappa(params, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   p <- length(control$beta)
   nz <- nrow(k2mat)

   if (ncol(wmat) == 0) {
      uSIGMA.11 <- diag(sqrt(kappa.e[etype] / weights))
   } else {
      SIGMA.11 <- Diagonal(x = kappa.e[etype] / weights) + tcrossprod(wmat %*%
                     as(Diagonal(x = sqrt(kappa.re)[retype]), "sparseMatrix"))
      uSIGMA.11 <- as.matrix(chol(SIGMA.11))
   }

   tKMAT <- t(k2mat) * sqrt(kappa.z)[ztype]
   uSIGMA.22 <- chol(crossprod(tKMAT, corMatrix(spcor) %*% tKMAT))

   logsqrtdet <- sum(log(c(diag(uSIGMA.11), diag(uSIGMA.22)))) 

   linvX.r1 <- backsolve(uSIGMA.11, xk1mat, transpose = TRUE)
   linvX.22 <- backsolve(uSIGMA.22, diag(-1, nrow(uSIGMA.22)), transpose = TRUE)
   linvX <- rbind(linvX.r1, cbind(matrix(0, nz, p), linvX.22))
   linvY <- c(backsolve(uSIGMA.11, y, transpose = TRUE), rep(0, nz))

   XtSiginvX <- crossprod(linvX)
   uXtSiginvX <- chol(XtSiginvX)
   logsqrtdetXtSiginvX <- sum(log(diag(uXtSiginvX)))

   XtSiginvY <- crossprod(linvX, linvY)

   betahat <- chol2inv(uXtSiginvX) %*% XtSiginvY

   resids <- linvY - linvX %*% betahat
   quadform <- as.numeric(crossprod(resids))

   shape <- sigma2shape(control)
   loglik <- -logsqrtdet - logsqrtdetXtSiginvX -
               (sum(shape) + (length(y) - p) / 2.0) *
               log(quadform / 2.0 + sum(sigma2scale(control) / kappa)) -
               sum((shape + 1.0) * log(kappa))

   list(value = loglik, betahat = betahat, quadform = quadform,
        uXtSiginvX = uXtSiginvX, uSig11inv = uiSIGMA.11)
}
