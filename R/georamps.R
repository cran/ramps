georamps <- function(fixed, random, correlation, data = sys.frame(sys.parent()),
   subset, weights, variance = list(fixed = ~ 1, random = ~ 1, spatial = ~ 1),
   aggregate = list(grid = NULL, blockid = ""), control = ramps.control(...),
   contrast = NULL, ...)
{
   ## Create data frame containing all relevant variables, except random which
   ##    is added later
   Call <- match.call()
   val <- list(fixed, variance$fixed, variance$random, variance$spatial)
   if (!is.null(aggregate$grid)) {
      val <- c(val, asOneSidedFormula(aggregate$blockid))
   }
   spvars <- all.vars(getCovariateFormula(correlation))
   val <- c(unlist(lapply(val, "all.vars")), spvars)
   mfargs <- list(formula = paste("~", paste(val, collapse = "+")),
                  data = data, weights = Call[["weights"]],
                  subset = Call[["subset"]], na.action = na.pass)
   mfdata <- do.call("model.frame", mfargs)

   ## Zero out point-source coordinates for aggregate measurements
   val <- mfdata[,aggregate$blockid]
   if (!is.null(val)) {
      idx <- c(aggregate$blockid, spvars)
      if (!all(idx %in% colnames(aggregate$grid)))
         stop("coordinates and 'blockid' must be given in 'grid'")
      aggregate$grid <- na.omit(aggregate$grid[,idx])
      mfdata[is.element(val, aggregate$grid[,aggregate$blockid]), spvars] <- 0
   }

   ## Remove incomplete records from the data frame
   mfdata <- na.omit(mfdata)

   ## Extract weights
   weights <- model.weights(mfdata)

   ## Initial model frame attributes are no longer needed
   attr(mfdata, "terms") <- NULL

   # Extract response variable and main effects design matrix
   mf <- model.frame(fixed, data = mfdata, drop.unused.levels = TRUE)
   y <- model.response(mf, "numeric")
   xmat <- Matrix(model.matrix(fixed, mf, contrast), sparse = TRUE)

   ## Indices to map measurement error variances
   variance$fixed <- if (is.null(variance$fixed)) factor(rep(1, nrow(mfdata)))
                     else factor(getCovariate(mfdata, variance$fixed))

   ## Structures for latent spatial parameters
   if (missing(correlation)) {
      stop("Unspecified correlation structure")
   } else {
      ## Create matrix of unique coordinates for latent parameters
      spf <- terms(getCovariateFormula(correlation))
      attr(spf, "intercept") <- 0

      if (is.null(aggregate$grid)) {
         idx1 <- rep(TRUE, nrow(mfdata))
         val <- model.matrix(spf, mfdata)
         idx2 <- NULL
      } else {
         idx1 <- !is.element(mfdata[,aggregate$blockid],
                             aggregate$grid[,aggregate$blockid])
         val <- model.matrix(spf, mfdata[idx1,,drop=FALSE])
         idx2 <- is.element(aggregate$grid[,aggregate$blockid],
                            mfdata[,aggregate$blockid])
         val <- rbind(val, model.matrix(spf, aggregate$grid[idx2,,drop=FALSE]))
      }
      sites <- unique.sites(val)

      ## Logical vector indicating z values to be monitored
      if (is.logical(control$z$monitor)) {
         control$z$monitor <- rep(control$z$monitor,
                                  length.out = nrow(sites$coords))
      } else {
         idx <- colnames(sites$coords)
         if (!all(idx %in% colnames(control$z$monitor)))
            stop("'z' monitor must be a logical value or matrix of coordinates")
         val <- unique.sites(control$z$monitor[,idx,drop=FALSE])
         val <- merge(cbind(sites$coords, 1:nrow(sites$coords)),
                      cbind(val$coords, 1), by=idx, all.x = TRUE)
         n <- length(idx)
         control$z$monitor <- !is.na(val[order(val[, n+1]), n+2])
      }
      ## Order latent parameters as (z$monitor == T, z$monitor == F)
      idx <- order(control$z$monitor, decreasing = TRUE)
      control$z$monitor <- control$z$monitor[idx]
      sites$coords <- sites$coords[idx, ]
      sites$idx <- order(idx)[sites$idx]

      ## Initialize correlation structure
      correlation <- Initialize(correlation, data = as.data.frame(sites$coords))

      ## Matrix to map latent parameters to observed data
      k <- model.matrix(~ factor(sites$idx) - 1)
      kmat <- matrix(0, nrow(mfdata), nrow(sites$coords))
      kmat[idx1,] <- k[seq(length.out = sum(idx1)),]
      if (length(idx2) > 0) {
         idx <- aggregate$grid[idx2, aggregate$blockid]
         kmat[match(sort(unique(idx)), mfdata[, aggregate$blockid]),] <-
            t(model.matrix(~ factor(idx) - 1)) %*%
            k[seq(sum(idx1) + 1, length.out = sum(idx2)),] / as.vector(table(idx))
      }

      ## Indices to map spatial variances
      val <- if (is.null(variance$spatial)) factor(rep(1, nrow(mfdata)))
             else factor(getCovariate(mfdata, variance$spatial))
      idx <- unlist(apply(as.numeric(val) * (kmat > 0), 2, unique))
      idx <- idx[idx > 0]
      if (length(idx) != ncol(kmat))
         stop("Differing spatial variances specified for measurements from the",
              " same site")
      else variance$spatial <- as.factor(levels(val)[idx])

      kmat <- as(kmat, "dgCMatrix")
   }

   ## Structures for random effects parameters
   if (missing(random)) {
      wmat <- matrix(numeric(0), 0, 0)
   } else {
      ## Matrix to map random effects to observed data
      g <- factor(getGroups(data, random)[as.numeric(rownames(mfdata))])
      w <- model.matrix(~ g - 1)
      wmat <- matrix(0, nrow(mfdata), ncol(w))
      wmat[as.numeric(rownames(w)), ] <- w

      ## Indices to map random effects variances
      val <- if (is.null(variance$random)) factor(rep(1, nrow(mfdata)))
             else factor(getCovariate(mfdata, variance$random))
      idx <- unlist(apply(as.numeric(val) * (wmat > 0), 2, unique))
      idx <- idx[idx > 0]
      if (length(idx) != ncol(wmat))
         stop("Differing random effects variances specified for measurements",
              " within the same group")
      else variance$random <- as.factor(levels(val)[idx])

      wmat <- as(wmat, "dgCMatrix")
   }

   ## Default values for weights if not supplied
   if (is.null(weights)) weights <- rowSums(as(kmat, "lgCMatrix"))

   ## Check parameter specifications against supplied data
   if (length(control$beta) != (n <- ncol(xmat)))
      stop("'beta' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$sigma2.e) != (n <- nlevels(variance$fixed)))
      stop("'sigma2.e' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$phi) != (n <- length(correlation)))
      stop("'phi' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$sigma2.z) != (n <- nlevels(variance$spatial)))
      stop("'sigma2.z' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$sigma2.re) != (n <- nlevels(variance$random)))
      stop("'sigma2.re' parameter specification in 'ramps.control' must be of",
           " length ", n)

   ## Set a single tuning parameter for the sigma2 parameters
   val <- sigma2tuning(control)
   idx <- val != 1
   val <- ifelse(any(idx), min(val[idx]), 1)
   if (length(control$sigma2.e)) control$sigma2.e$tuning[] <- val
   if (length(control$sigma2.z)) control$sigma2.z$tuning[] <- val
   if (length(control$sigma2.re)) control$sigma2.re$tuning[] <- val

   ## Obtain MCMC samples from ramps engine
   val <- ramps.engine(y, xmat, kmat, wmat, correlation, variance$fixed,
                       variance$spatial, variance$random, weights, control)

   structure(
      list(params = as.mcmc(val$params),
           z = if (any(control$z$monitor)) as.mcmc(val$z) else val$z, 
           coords = sites$coords, loglik = val$loglik,
           y = y, xmat = xmat, etype = variance$fixed, weights = weights,
           correlation = correlation, kmat = kmat, ztype = variance$spatial,
           wmat = wmat, retype = variance$random,
           control = control, evals = val$evals, call = Call),
      class = "ramps")
}


print.ramps <- function(x, ...)
{
   cat("\nCall: ", paste(deparse(x$call), collapse = "\n"), "\n")

   params <- colnames(x$params)
   cat("\nCoefficients:\n")
   if (length(tmp <- params2beta(params, x$control)) > 0)
      print.default(tmp, print.gap = 2, quote = FALSE)

   sigma2 <- params2kappa(params, x$control)

   n <- sum(rowSums(as(x$kmat, "lgCMatrix")) > 1)
   cat("\nMeasurements\n",
       " N = ", length(x$y), "\n",
       " Point Source = ", length(x$y) - n, "\n",
       " Areal = ", n, "\n",
       " Error Variance: ",
       paste(kappa2kappa.e(sigma2, x$control), collapse = " "), "\n", sep = "")

   cat("\nLatent Spatial Process\n",
       " Sites = ", ncol(x$kmat), "\n",
       " Correlation: ", class(x$correlation)[1], "(",
          paste(params2phi(params, x$control), collapse = ", "), ")\n",
       " Variance: ",
       paste(kappa2kappa.z(sigma2, x$control), collapse = " "), "\n", sep = "")

   if (ncol(x$wmat) > 0) {
      cat("\nRandom Effects\n",
          " N = ", ncol(x$wmat), "\n",
          " Variance: ",
          paste(kappa2kappa.re(sigma2, x$control), collapse = " "), "\n",
          sep = "")
   }

   rn <- rownames(x$params)
   cat("\nMCMC Output\n",
       " Saved Samples = ", length(rn), " (", paste(rn[1:3], collapse = ", "),
          ", ..., ", rn[length(rn)], ")\n",
       " Slice Evaluations = ", x$evals, "\n", sep = "")

   invisible(x)
}


summary.ramps <- function(object, iter, ...)
{
   ## Vector of iterations for the calculation
   iter <- if (missing(iter)) 1:nrow(object$params)
           else match.iter(iter, object$params)

   summary(as.mcmc(object$params[iter, , drop=F]), ...)
}
