################################################################################
## mpdpred - Marginalized posterior density with flexible prediction sites
##
## Arguments:
##    params  - numerical vector of phi and kappa parameter values
##    Y       - concatenated response Matrix Y = c(y, rep(0, nzp))
##    X       - concatenated model Matrix ((n + nzp) x (p + nzp))
##    k22mat  - design Matrix of nonprediction spatial random effect (ny2 by nzu)
##    wmat    - design Matrix of non-spatial random effect (n by q)
##    spcor   - initialized (nlme) spatial correlation structure
##    etype   - factor indexing the measurement variances (n x 1)
##    ztype   - factor indexing the spatial variances (nz x 1)
##    retype  - factor indexing the random effect variances (q x 1)
##    weights - numerical vector by which to weight the measurement error
##              variance
##    control - ramps.control object
##
## Notes:
##    need dimension ny1, ny2, nzp to split the rows in Y = c(y1, y2, 0)
##    note that if ny1 = 0 then nzp = 0 and vice versa
##           this is because ny1 is caused by nzp;
##           if nzp = 0, then everything is in ny2.
##    note that nzp = length(0) = length(zp)
##    note that ny2 and nzp cannot be 0 at the same time, i.e., ny2 + nzp > 0
################################################################################

mpdpred <- function(params, Y, X, k22mat, wmat, spcor, etype, ztype, retype,
                    weights, control)
{
   ## Get dimensions
   n <- length(weights)
   nz <- length(ztype)
   p <- length(control$beta)
   ny2 <- nrow(k22mat)
   nzu <- ncol(k22mat)
   ny1 <- n - ny2
   nzp <- nz - nzu

   ## Extract current parameter values
   phi <- params2phi(params, control)
   kappa <- params2kappa(params, control)
   kappa.e <- kappa2kappa.e(kappa, control)
   kappa.z <- kappa2kappa.z(kappa, control)
   kappa.re <- kappa2kappa.re(kappa, control)

   SIGMA.e <- as(as(Diagonal(x = kappa.e[etype] / weights), "sparseMatrix"),
                 "dtCMatrix")

   ## BEGIN: uiSIGMA.11
   if (ny1 == 0) {
      uiSIGMA.11 <- Matrix(0, 0, 0)
   } else {
      y1Idx <- 1:ny1
      SIGMA.11 <- SIGMA.e[y1Idx, y1Idx]
      if (ncol(wmat) > 0) {
        SIGMA.11 <- SIGMA.11 + tcrossprod(wmat[y1Idx, , drop = FALSE] %*%
                       as(Diagonal(x = sqrt(kappa.re)[retype]), "sparseMatrix"))
      }
      uiSIGMA.11 <- solve(chol(SIGMA.11))
   }
   ## END: uiSIGMA.11

   ## BEGIN: uiSIGMA.22
   spcor[] <- unconstrained(spcor, phi)
   vroot <- sqrt(kappa.z)[ztype]
   SigmaZ <- as(vroot * t(corMatrix(spcor) * vroot), "dgCMatrix")

   y2Idx <- ZuIdx <- ZpIdx <- numeric(0)
   S11 <- S22 <- Matrix(0, 0, 0)
   if (ny2 > 0) {
      y2Idx <- (ny1 + 1):n
      ZuIdx <- (nzp + 1):nz
      S11 <- k22mat %*% tcrossprod(SigmaZ[ZuIdx, ZuIdx], k22mat) +
                SIGMA.e[y2Idx, y2Idx]
      if (ncol(wmat) > 0) {
         S11 <- S11 + tcrossprod(wmat[y2Idx, , drop = FALSE] %*%
                       as(Diagonal(x = sqrt(kappa.re)[retype]), "sparseMatrix"))
      }
   }
   if (nzp > 0) {
      ZpIdx <- 1:nzp
      S22 <- SigmaZ[ZpIdx, ZpIdx]
   }
   S12 <- if (ny2 > 0 && nzp > 0) k22mat %*% SigmaZ[ZuIdx, ZpIdx]
          else Matrix(0, ny2, nzp)
   uiSIGMA.22 <- solve(chol(rBind(cBind(S11, S12), cBind(t(S12), S22))))
   ## END: uiSIGMA.22

   logsqrtdet <- -1 * sum(log(c(diag(uiSIGMA.11), diag(uiSIGMA.22))))

   linvX.r1 <- if (ny1 > 0) crossprod(uiSIGMA.11, X[1:ny1,])
               else Matrix(0, 0, ncol(X))
   linvX.r2 <- crossprod(uiSIGMA.22, X[(ny1 + 1):(ny1 + ny2 + nzp),])
   linvX <- rBind(as(linvX.r1, "dgCMatrix"), as(linvX.r2, "dgCMatrix"))

   linvY.r1 <- if (ny1 > 0) crossprod(uiSIGMA.11, Y[1:ny1])
               else Matrix(0, 0, 1)
   linvY.r2 <- crossprod(uiSIGMA.22, Y[(ny1 + 1):(ny1 + ny2 + nzp)])
   linvY <- rBind(as(linvY.r1, "dgCMatrix"), as(linvY.r2, "dgCMatrix"))

   XtSiginvX <- crossprod(linvX)
   XtSiginvY <- crossprod(linvX, linvY)

   list(betahat = solve(XtSiginvX, XtSiginvY), uXtSiginvX = chol(XtSiginvX))

   #############################################################################
   ## Legacy code used previously for calling mpdpred in the slice sampler.
   ## mpdensity is now used exclusively for slice sampling, so this block is no
   ## longer needed. 
   #############################################################################
   ## betahat <- solve(XtSiginvX, XtSiginvY)
   ## uXtSiginvX <- chol(XtSiginvX)
   ##
   ## logsqrtdetXtSiginvX <- sum(log(diag(uXtSiginvX)))
   ## resids <- linvY - linvX %*% betahat  
   ## quadform <- as.numeric(crossprod(resids))
   ##
   ## shape <- sigma2shape(control)
   ## loglik <- -logsqrtdet - logsqrtdetXtSiginvX -
   ##              (sum(shape) + (length(Y) - p) / 2) *
   ##              log(quadform / 2 + sum(sigma2scale(control) / kappa)) -
   ##              sum((shape + 1) * log(kappa))
   ##
   ## list(value = loglik, betahat = betahat, quadform = quadform,
   ##      uXtSiginvX = uXtSiginvX)
   #############################################################################

}
