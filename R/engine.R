ramps.engine <- function(y, xmat, kmat, wmat, spcor, etype, ztype, retype,
                         weights, control)
################################################################################
# y:        the response, n by 1
# xmat:     the covariate matrix, n by p
# kmat:     the design matrix of spatial random effect, n by NS;
# wmat:     the design matrix of non-spatial random effect, n by q (number of
#              random effects)
# spcor:    initialized (nlme) spatial correlation structure
# etype:    integer vector giving the type of measurement variances, n by 1
# ztype:    integer vector giving the type of spatial variances, nz by 1
# retype:   integer vector giving the type of random effect variances, q by 1
# weights:  number of observations each observed value is based upon
# control:  ramps.control object
################################################################################
{
  ## Global constants
  n <- length(y)
  p <- length(control$beta)
  nzp <- sum(control$z$monitor)
  iter <- control$iter

  ## Return objects
  val <- c(name.ext("phi", 1:length(control$phi)),
           name.ext("sigma2.e", levels(etype)),
           name.ext("sigma2.z", levels(ztype)),
           name.ext("sigma2.re", levels(retype)), colnames(xmat))
  params <- matrix(NA, length(iter), length(val), dimnames = list(iter, val))

  params.z <- matrix(NA, length(iter), nzp, dimnames = list(iter, NULL))
  if (nzp > 0) colnames(params.z) <- paste("z_", 1:nzp, sep="") 

  loglik <- structure(rep(NA, length(iter)), names = iter)
  evals <- 0

  ## Initialize external output files
  write.header(colnames(params), control$file$params)
  write.header(colnames(params.z), control$file$z)

  ## Construct matrix components for mpdensity
  val <- unique.sites(as.matrix(kmat))
  k2mat <- as(val$coords, "dgCMatrix")
  xk1mat <- cBind(xmat, as(model.matrix(~ factor(val$idx) - 1), "dgCMatrix"))

  zidx <- seq(length.out = nzp)
  idx <- match(zidx, as.vector((val$coords == 1) %*% 1:ncol(kmat)))
  if (any(is.na(idx))) {
    pred <- TRUE
    zidx <- p + zidx

    ## Extract matrix components for mpdpred
    val <- as.vector(kmat %*% as.numeric(!control$z$monitor)) == 0
    idx <- order(val, decreasing = TRUE)
    ny1 <- sum(val)
    ny2 <- n - ny1
    r1 <- seq(length.out = ny1)
    r2 <- seq(length(r1) + 1, length.out = ny2) 
    c1 <- seq(length.out = sum(control$z$monitor))
    c2 <- seq(length(c1) + 1, length.out = sum(!control$z$monitor))
    k11mat <- kmat[idx[r1], c1]
    k21mat <- kmat[idx[r2], c1]
    k22mat <- kmat[idx[r2], c2]

    ## Reorder remaining data structures by idx
    y <- y[idx]
    if (ncol(xmat) > 0) xmat <- xmat[idx, , drop = FALSE]
    xk1mat <- xk1mat[idx, , drop = FALSE]
    if (ncol(wmat) > 0) wmat <- wmat[idx, , drop = FALSE]
    weights <- weights[idx]
    etype <- etype[idx]

    ## Construct X and Y structures for mpdpred
    Y <- c(y, rep(0, nzp))
    X <- Matrix(0, n + nzp, p + nzp)
    X[1:n, seq(length.out=p)] <- xmat
    if (nzp > 0) {
      X[1:ny1, (p + 1):(p + nzp)] <- k11mat
      X[(n + 1):(n + nzp), (p + 1):(p + nzp)] <-
                                 as(Diagonal(x = rep(-1, nzp)), "sparseMatrix")
    }
    if (nzp > 0 && ny2 > 0) X[(ny1 + 1):n, (p + 1):(p + nzp)] <- k21mat
  } else {
    pred <- FALSE
    zidx <- p + idx
  }

  ## First mpdensity evaluation
  theta <- sigma2init(control)
  theta <- c(control$phi$init, theta / sum(theta))
  curreval <- mpdensity(theta, y, xk1mat, k2mat, wmat, spcor, etype, ztype,
                        retype, weights, control)

  idx <- 1
  cat("MCMC Sampler Progress (N = ", max(iter), "):\n", sep="")
  for (i in 1:max(iter)) {
    ## Draw phi and kappa
    curreval <- sliceSimplex(x=theta, mpdfun=mpdensity, log=TRUE,
                             fx=curreval$value, y=y, xk1mat=xk1mat, k2mat=k2mat,
                             wmat=wmat, spcor=spcor, etype=etype, ztype=ztype,
                             retype=retype, weights=weights, control=control)

    theta <- curreval$params
    evals <- evals + curreval$newevals

    if (i == iter[idx]) {
      ## Draw variance parameters sigma2 
      kappa <- params2kappa(theta, control)
      as2 <- sum(sigma2shape(control)) + (n - p) / 2.0
      bs2 <- sum(sigma2scale(control) / kappa) + curreval$quadform / 2.0
      sigma2.tot <- 1.0 / rgamma(1, as2, bs2)
      sigma2 <- sigma2.tot * kappa

      ## Likelihood evaluation
      sigma2.e <- kappa2kappa.e(sigma2, control)
      sigma2.re <- kappa2kappa.re(sigma2, control)

      BETA <- curreval$betahat + sqrt(sigma2.tot) *
                solve(curreval$uXtSiginvX, rnorm(length(curreval$betahat)))
      MU <- y - xk1mat %*% BETA
      SIGMA <- Diagonal(x = sigma2.e[etype] / weights)
      if(ncol(wmat) > 0) {
        SIGMA <- SIGMA + tcrossprod(wmat %*%
                   as(Diagonal(x = sqrt(sigma2.re)[retype]), "sparseMatrix"))
      }
      uiSIGMA <- solve(chol(SIGMA))
      loglik[idx] <- sum(log(diag(uiSIGMA))) -
                       as.numeric(crossprod(crossprod(uiSIGMA, MU))) / 2.0

      ## Draw beta and z parameters
      if (pred) {
        val <- mpdpred(theta, Y, X, k22mat, wmat, spcor, etype, ztype, retype,
                       weights, control)
        BETA <- val$betahat + sqrt(sigma2.tot) *
                  solve(val$uXtSiginvX, rnorm(p + nzp))
      }
      beta <- BETA[seq(length.out = p)]
      z <- BETA[zidx]

      ## Save sampled parameters
      params[idx, ] <- c(params2phi(theta, control), sigma2, beta)
      write.params(i, params[idx,], control$file$params)
      params.z[idx, ] <- z
      write.params(i, z, control$file$z)

      idx <- idx + 1
    }

    print.iter(i)
  }
  cat("\n")

  list(params = params, z = params.z,
       loglik = loglik - (n / 2.0 * log(2.0 * pi)), evals = evals)
}
