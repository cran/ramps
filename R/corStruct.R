################################################################################
# Generic methods to be used with correlation structures
################################################################################

unconstrained <- function(object, ...)
{
   UseMethod("unconstrained")
}


################################################################################
# corStruct method functions
################################################################################

Initialize.corStruct <- function(object, data, ...)
{
   form <- formula(object)
   ## obtaining the groups information, if any
   if (!is.null(getGroupsFormula(form))) {
      attr(object, "groups") <- getGroups(object, form, data = data)
      attr(object, "Dim") <- Dim(object, attr(object, "groups"))
   } else {                              # no groups
      attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
   }
   ## obtaining the covariate(s)
   attr(object, "covariate") <- getCovariate(object, data = data)

   object
}


Dim.corStruct <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   ugrp <- unique(groups)
   groups <- factor(groups, levels = ugrp)
   len <- table(groups)
   list(N = length(groups),
        M = length(len),
        maxLen = max(len),
        sumLenSq = sum(len^2),
        len = len,
        start = match(ugrp, groups) - 1)
}


################################################################################
# corSpatial method functions
################################################################################

Initialize.corSpatial <- function(object, data, ...)
{
   if (!is.null(attr(object, "minD"))) { #already initialized
      return(object)
   }
   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) > 0) {		# initialized
      if (val[1] <= 0) {
         stop("Range must be > 0 in \"corSpatial\" initial value")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 1) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else {
            if (length(val) != 2) {
               stop("Initial value for \"corSpatial\" parameters of wrong dimension")
            }
         }
         if ((val[2] <= 0) || (val[2] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else {				# only range parameter
         if (length(val) != 1) {
            stop("Initial value for \"corSpatial\" parameters of wrong dimension")
         }
      }
   } else {
      val <- min(unlist(attr(object, "covariate"))) * 0.9
      if (nug) val <- c(val, 0.1)
   }
   val[1] <- log(val[1])
   if (nug) val[2] <- log(val[2] / (1 - val[2]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "minD") <- min(unlist(attr(object, "covariate")))
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


getCovariate.corSpatial <- function(object, form = formula(object), data)
{
   if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
      if (missing(data)) {
         stop("Need data to calculate covariate")
      }
      covForm <- getCovariateFormula(form)
      if (length(all.vars(covForm)) > 0) { # covariate present
         if (attr(terms(covForm), "intercept") == 1) {
            covForm <-
               eval(parse(text = paste("~", deparse(covForm[[2]]),"-1",sep="")))
         }
         covar <-
            as.data.frame(unclass(model.matrix(covForm,
                          model.frame(covForm, data, drop.unused.levels = TRUE))))
      } else {
         covar <- NULL
      }

      if (!is.null(getGroupsFormula(form))) { # by groups
         grps <- getGroups(object, data = data)
         if (is.null(covar)) {
            covar <- lapply(split(grps, grps),
                            function(x) as.vector(dist2(1:length(x))))
         } else {
            covar <- lapply(split(covar, grps),
                            function(el, metric, radius) {
                               el <- as.matrix(el)
                               if (nrow(el) > 1) {
                                  as.vector(dist2(el, metric, r = radius))
                               } else {
                                  numeric(0)
                               }
                            }, metric = attr(object, "metric"),
                               radius = attr(object, "radius"))
         }
         covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups
      } else {				# no groups
         if (is.null(covar)) {
            covar <- as.vector(dist2(1:nrow(data)))
         } else {
            covar <- as.vector(dist2(as.matrix(covar),
                                     method = attr(object, "metric"),
                                     r = attr(object, "radius")))
         }
      }
      if (any(unlist(covar) == 0)) {
         stop("Cannot have zero distances in \"corSpatial\"")
      }
   }

   covar
}


corFactor.corSpatial <- function(object, ...)
{
   if (!is.null(aux <- attr(object, "factor"))) {
      return(aux)
   }
   val <- corMatrix(object, corr = FALSE, ...)
   lD <- attr(val, "logDet")
   if (is.list(val)) val <- unlist(val)
   else val <- as.vector(val)
   names(val) <- NULL
   attr(val, "logDet") <- lD

   val
}


corMatrix.corSpatial <- function(object, covariate = getCovariate(object),
   corr = TRUE, ...)
{
   par <- coef(object, unconstrained = FALSE)
   cor0 <- 1
   switch(class(object)[1],
      corRExp = {
         val <- cor.exp(unlist(covariate), par[1])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[2])
            val <- cor0 * val
         }
      },
      corRExpwr = {
         val <- cor.exp(unlist(covariate), par[1], par[2])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[3])
            val <- cor0 * val
         }
      },
      corRGaus = {
         val <- cor.exp(unlist(covariate), par[1], 2)
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[2])
            val <- cor0 * val
         }
      },
      corRLin = {
         val <- cor.lin(unlist(covariate), par[1])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[2])
            val <- cor0 * val
         }
      },
      corRMatern = {
         val <- cor.matern(unlist(covariate), par[1], par[2])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[3])
            val <- cor0 * val
         }
      },
      corRRatio = {
         val <- cor.ratio(unlist(covariate), par[1])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[2])
            val <- cor0 * val
         }
      },
      corRSpher = {
         val <- cor.spher(unlist(covariate), par[1])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[2])
            val <- cor0 * val
         }
      },
      corRWave = {
         val <- cor.wave(unlist(covariate), par[1])
         if (attr(object, "nugget")) {
            cor0 <- (1 - par[2])
            val <- cor0 * val
         }
      }
   )

   if (data.class(covariate) == "list") {
      if (is.null(names(covariate))) {
         names(covariate) <- 1:length(covariate)
      }
      corD <- Dim(object, rep(names(covariate),
                  unlist(lapply(covariate,
                  function(el) round((1 + sqrt(1 + 8 * length(el)))/2)))))
   } else {
      corD <- Dim(object, rep(1, round((1 + sqrt(1 + 8* length(covariate)))/2)))
   }
   len <- corD[["len"]]
   val <- split(val, rep(1:corD[["M"]], (len * (len - 1)) / 2))
   lD <- NULL
   for(i in seq(val)) {
      x <- matrix(0, len[i], len[i])
      x[lower.tri(x)] <- val[[i]]
      if (corr) {
         val[[i]] <- x + t(x)
         diag(val[[i]]) <- cor0
      } else {
         diag(x) <- cor0
         l <- chol(t(x))
         val[[i]] <- t(backsolve(l, diag(len[i])))
         lD <- c(lD, diag(l))
      }
   }
   if (i < 2) val <- val[[1]]
   else names(val) <- names(len)

   if (!is.null(lD)) lD <- -1 * sum(log(lD))
   attr(val, "logDet") <- lD

   val
}


unconstrained.corSpatial <- function(object, value, inv = FALSE, ...)
{
   n <- length(object)
   if(!inv) {
      val <- log(value)
      if (attr(object, "nugget")) val[n] <- log(value[n] / (1 - value[n]))
   } else {
      val <- exp(value)
      if (attr(object, "nugget")) val[n] <- val[n] / (1 + val[n])
   }
   val
}


coef.corSpatial <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)
   nm <- "range"
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm
   val
}


"coef<-.corSpatial" <- function(object, ..., value)
{
   if (length(value) != length(object)) {
      stop("Cannot change the length of the parameter after initialization")
   }

   object[] <- value

   ## updating the factor list and logDet
   attr(object, "factor") <- NULL
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


Dim.corSpatial <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   val <- Dim.corStruct(object, groups)
   val[["start"]] <-
      c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
   ## will use third component of Dim list for spClass
   names(val)[3] <- "spClass"
   val[[3]] <- match(class(object)[1], c("corRExp", "corRExpwr", "corRGaus",
                     "corRLin", "corRMatern", "corRRatio", "corRSpher"), 0)

   val
}


################################################################################
# corRExp - exponential spatial correlation structure
################################################################################

corRExp <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRExp", "corSpatial", "corStruct")

   value
}


################################################################################
# corRExpwr - Powered exponential spatial correlation structure
################################################################################

corRExpwr <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRExpwr", "corSpatial", "corStruct")
   value
}


Initialize.corRExpwr <- function(object, data, ...)
{
   if (!is.null(attr(object, "minD"))) { #already initialized
      return(object)
   }
   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) >= 2) {		# initialized
      if (any(val[1:2] <= 0)) {
         stop("Initial values for \"corSpatial\" parameters must be > 0")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 2) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else {
            if (length(val) != 3) {
               stop("Initial values for \"corSpatial\" parameters of wrong dimension")
            }
         }
         if ((val[3] <= 0) || (val[3] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else {				# only correlation parameters
         if (length(val) != 2) {
            stop("Initial values for \"corSpatial\" parameters of wrong dimension")
         }
      }
   } else {
      val <- c(min(unlist(attr(object, "covariate"))) * 0.9, 1)
      if (nug) val <- c(val, 0.1)
   }
   val[1:2] <- log(val[1:2])
   if (nug) val[3] <- log(val[3] / (1 - val[3]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "minD") <- min(unlist(attr(object, "covariate")))
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


coef.corRExpwr <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)
   nm <- c("range", "shape")
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm

   val
}


################################################################################
#corRGaus - Gaussian spatial correlation structure
################################################################################

corRGaus <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRGaus", "corSpatial", "corStruct")
   value
}


################################################################################
# corRLin - Linear spatial correlation structure
################################################################################

corRLin <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRLin", "corSpatial", "corStruct")
   value
}


################################################################################
# corRMatern - Matern spatial correlation structure
################################################################################

corRMatern <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRMatern", "corSpatial", "corStruct")
   value
}


Initialize.corRMatern <- function(object, data, ...)
{
   if (!is.null(attr(object, "minD"))) { #already initialized
      return(object)
   }
   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) >= 2) {		# initialized
      if (any(val[1:2] <= 0)) {
         stop("Initial values for \"corSpatial\" parameters must be > 0")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 2) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else {
            if (length(val) != 3) {
               stop("Initial values for \"corSpatial\" parameters of wrong dimension")
            }
         }
         if ((val[3] <= 0) || (val[3] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else {				# only correlation parameters
         if (length(val) != 2) {
            stop("Initial values for \"corSpatial\" parameters of wrong dimension")
         }
      }
   } else {
      val <- c(min(unlist(attr(object, "covariate"))) * 0.9, 0.5)
      if (nug) val <- c(val, 0.1)
   }
   val[1:2] <- log(val[1:2])
   if (nug) val[3] <- log(val[3] / (1 - val[3]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "minD") <- min(unlist(attr(object, "covariate")))
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


coef.corRMatern <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)
   nm <- c("range", "shape")
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm

   val
}


################################################################################
# corRRatio - rational quadratic spatial correlation structure
################################################################################

corRRatio <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRRatio", "corSpatial", "corStruct")
   value
}


################################################################################
# corRSpher - spherical spatial correlation structure
################################################################################

corRSpher <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRSpher", "corSpatial", "corStruct")
   value
}


################################################################################
# corRWave - sine wave spatial correlation structure
################################################################################

corRWave <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
   attr(value, "formula") <- form
   attr(value, "nugget") <- nugget
   attr(value, "metric") <- match.arg(metric)
   attr(value, "radius") <- radius
   attr(value, "fixed") <- fixed
   class(value) <- c("corRWave", "corSpatial", "corStruct")
   value
}


################################################################################
# corSpatioTemporal method functions
################################################################################

Initialize.corSpatioTemporal <- function(object, data, ...)
{
   if (!is.null(attr(object, "logDet"))) { #already initialized
      return(object)
   }
   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) >= 2) {		# initialized
      if (any(val[1:2] <= 0)) {
         stop("Initial values for \"corSpatioTemporal\" parameters must be > 0")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 2) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else {
            if (length(val) != 3) {
               stop("Initial values for \"corSpatioTemporal\" parameters of wrong dimension")
            }
         }
         if ((val[3] <= 0) || (val[3] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else {				# only correlation parameters
         if (length(val) != 2) {
            stop("Initial values for \"corSpatioTemporal\" parameters of wrong dimension")
         }
      }
   } else {
      covar <- attr(object, "covariate")
      x <- if (is.list(covar)) unlist(sapply(covar, attr, which = "dist"))
           else attr(covar, "dist")
      y <- if (is.list(covar)) unlist(sapply(covar, attr, which = "period"))
           else attr(covar, "period")
      val <- c(ifelse(any(x > 0), min(x[x > 0]), 1) * 0.9,
               ifelse(any(y > 0), min(y[y > 0]), 1) * 0.9)
      if (nug) val <- c(val, 0.1)
   }
   val[1:2] <- log(val[1:2])
   if (nug) val[3] <- log(val[3] / (1 - val[3]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


getCovariate.corSpatioTemporal <- function(object, form = formula(object), data)
{
   if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
      if (missing(data)) {
         stop("Need data to calculate covariate")
      }
      covForm <- getCovariateFormula(form)
      tcovar <- length(all.vars(covForm))
      if (tcovar >= 3) { # covariates present
         if (attr(terms(covForm), "intercept") == 1) {
            covForm <- eval(parse(text = paste("~", deparse(covForm[[2]]),
                          "-1", sep = "")))
         }
         covar <- as.data.frame(unclass(model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))))
         if (nrow(covar) > nrow(unique(covar))) {
            stop("Cannot have duplicate sites in \"corSpatioTemporal\"")
         }
      } else {
         covar <- NULL
      }

      if (!is.null(getGroupsFormula(form))) { # by groups
         if (!is.null(covar)) {
            grps <- getGroups(object, data = data)
            covar <- lapply(split(covar, grps),
                        function(el, metric, radius) {
                            el <- as.matrix(el)
                            r <- nrow(el)
                            if (r > 1) {
                               d <- as.vector(dist2(el[, -tcovar],
                                                    metric, r = radius))
                               x <- matrix(0, r, r)
                               idx <- lower.tri(x)
                               period <- abs(el[col(x)[idx], tcovar] -
                                             el[row(x)[idx], tcovar])
                            } else {
                               d <- numeric(0)
                               period <- numeric(0)
                            }
                            attr(el, "dist") <- d
                            attr(el, "period") <- period
                            el
                         }, metric = attr(object, "metric"),
                            radius = attr(object, "radius"))
         }
         covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups
      } else {				# no groups
         if (!is.null(covar)) {
            covar <- as.matrix(covar)
            attr(covar, "dist") <- as.vector(
                                      dist2(covar[, -tcovar],
                                            method = attr(object, "metric"),
                                            r = attr(object, "radius")))
            x <- matrix(0, nrow(covar), nrow(covar))
            idx <- lower.tri(x)
            attr(covar, "period") <- abs(covar[col(x)[idx], tcovar] -
                                         covar[row(x)[idx], tcovar])
         }
      }
   }

   covar
}


coef.corSpatioTemporal <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)
   nm <- c("spatial range", "temporal range")
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm

   val
}


Dim.corSpatioTemporal <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   val <- Dim.corStruct(object, groups)
   val[["start"]] <-
      c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
   ## will use third component of Dim list for spClass
   names(val)[3] <- "spClass"
   val[[3]] <- match(class(object)[1], c("corRExp2", "corRExpwr2"), 0)

   val
}


corMatrix.corSpatioTemporal <- function(object, covariate = getCovariate(object),
   corr = TRUE, ...)
{
   lD <- NULL
   par <- coef(object, unconstrained = FALSE)

   switch(class(object)[1],
     corRExp2 = corf <- function(dist, period, par, nugget = FALSE) {
        val <- cor.exp2(dist, period, par[1], 1, par[2], 1, par[3])
        if (nugget) val <- (1 - par[4]) * val
        val
     },
     corRExpwr2 = corf <- function(dist, period, par, nugget = FALSE) {
        val <- cor.exp2(dist, period, par[1], par[2], par[3], par[4], par[5])
        if (nugget) val <- (1 - par[6]) * val
        val
     }
   )
   nug <- attr(object, "nugget")

   if (data.class(covariate) == "list") {
      if (is.null(names(covariate))) {
         nm <- seq(covariate)
         names(covariate) <- nm
      } else {
         nm <- names(covariate)
      }
      len <- unlist(lapply(covariate, nrow))
      val <- list()
      for(i in nm) {
         el <- covariate[[i]]
         x <- matrix(0, len[i], len[i])
         x[lower.tri(x)] <- corf(attr(el, "dist"), attr(el, "period"), par, nug)
         if (corr) {
            val[[i]] <- x + t(x)
            diag(val[[i]]) <- corf(0, 0, par, nug)
         } else {
            diag(x) <- corf(0, 0, par, nug)
            l <- chol(t(x))
            val[[i]] <- t(backsolve(l, diag(len[i])))
            lD <- c(lD, diag(l))
         }
      }
   } else {
      len <- nrow(covariate)
      x <- matrix(0, len, len)
      x[lower.tri(x)] <- corf(attr(covariate, "dist"), attr(covariate, "period"),
                              par, attr(object, "nugget"))
      if (corr) {
         val <- x + t(x)
         diag(val) <- corf(0, 0, par, nug)
      } else {
         diag(x) <- corf(0, 0, par, nug)
         l <- chol(t(x))
         val <- t(backsolve(l, diag(len)))
         lD <- c(lD, diag(l))
      }
   }

   if (!is.null(lD)) lD <- -1 * sum(log(lD))
   attr(val, "logDet") <- lD

   val
}


corFactor.corSpatioTemporal <- function(object, ...)
{
   if (!is.null(aux <- attr(object, "factor"))) {
      return(aux)
   }
   val <- corMatrix(object, corr = FALSE, ...)
   lD <- attr(val, "logDet")
   if (is.list(val)) val <- unlist(val)
   else val <- as.vector(val)
   names(val) <- NULL
   attr(val, "logDet") <- lD

   val
}


################################################################################
# corRExp2 - Exponential spatio-temporal correlation structure
################################################################################

corRExp2 <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
  attr(value, "formula") <- form
  attr(value, "nugget") <- nugget
  attr(value, "metric") <- match.arg(metric)
  attr(value, "radius") <- radius
  attr(value, "fixed") <- fixed
  class(value) <- c("corRExp2", "corSpatioTemporal", "corSpatial", "corStruct")

  value
}


Initialize.corRExp2 <- function(object, data, ...)
{
   if (!is.null(attr(object, "logDet"))) { #already initialized
      return(object)
   }
   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) >= 3) {		# initialized
      if (any(val[1:3] <= 0)) {
         stop("Initial values for \"corSpatioTemporal\" parameters must be > 0")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 3) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else {
            if (length(val) != 4) {
               stop("Initial values for \"corSpatioTemporal\" parameters of wrong dimension")
            }
         }
         if ((val[4] <= 0) || (val[4] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else {				# only correlation parameters
         if (length(val) != 3) {
            stop("Initial values for \"corSpatioTemporal\" parameters of wrong dimension")
         }
      }
   } else {
      covar <- attr(object, "covariate")
      x <- if (is.list(covar)) unlist(sapply(covar, attr, which = "dist"))
           else attr(covar, "dist")
      y <- if (is.list(covar)) unlist(sapply(covar, attr, which = "period"))
           else attr(covar, "period")
      val <- c(ifelse(any(x > 0), min(x[x > 0]), 1) * 0.9,
               ifelse(any(y > 0), min(y[y > 0]), 1) * 0.9, 1)
      if (nug) val <- c(val, 0.1)
   }
   val[1:3] <- log(val[1:3])
   if (nug) val[4] <- log(val[4] / (1 - val[4]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


coef.corRExp2 <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)
   nm <- c("spatial range", "temporal range", "interaction")
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm

   val
}


################################################################################
# corRExpwr2 - Powered exponential spatio-temporal correlation structure
################################################################################

corRExpwr2 <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
  attr(value, "formula") <- form
  attr(value, "nugget") <- nugget
  attr(value, "metric") <- match.arg(metric)
  attr(value, "radius") <- radius
  attr(value, "fixed") <- fixed
  class(value) <- c("corRExpwr2", "corSpatioTemporal", "corSpatial", "corStruct")

  value
}


Initialize.corRExpwr2 <- function(object, data, ...)
{
   if (!is.null(attr(object, "logDet"))) { #already initialized
      return(object)
   }
   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) >= 5) {		# initialized
      if (any(val[1:5] <= 0)) {
         stop("Initial values for \"corSpatioTemporal\" parameters must be > 0")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 5) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else {
            if (length(val) != 6) {
               stop("Initial values for \"corSpatioTemporal\" parameters of wrong dimension")
            }
         }
         if ((val[6] <= 0) || (val[6] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else {				# only correlation parameters
         if (length(val) != 5) {
            stop("Initial values for \"corSpatioTemporal\" parameters of wrong dimension")
         }
      }
   } else {
      covar <- attr(object, "covariate")
      x <- if (is.list(covar)) unlist(sapply(covar, attr, which = "dist"))
           else attr(covar, "dist")
      y <- if (is.list(covar)) unlist(sapply(covar, attr, which = "period"))
           else attr(covar, "period")
      val <- c(ifelse(any(x > 0), min(x[x > 0]), 1) * 0.9, 1,
               ifelse(any(y > 0), min(y[y > 0]), 1) * 0.9, 1, 1)
      if (nug) val <- c(val, 0.1)
   }
   val[1:5] <- log(val[1:5])
   if (nug) val[6] <- log(val[6] / (1 - val[6]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


coef.corRExpwr2 <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)
   nm <- c("spatial range", "spatial shape", "temporal range",
           "temporal shape", "interaction")
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm

   val
}


################################################################################
# corRExpwrDt - Joint powered exponential spatial and integrated exponential
#               temporal correlation structure
################################################################################

corRExpwr2Dt <- function(value = numeric(0), form = ~ 1, nugget = FALSE,
   metric = c("euclidean", "maximum", "manhattan", "haversine"), radius = 3956,
   fixed = FALSE)
{
  attr(value, "formula") <- form
  attr(value, "nugget") <- nugget
  attr(value, "metric") <- match.arg(metric)
  attr(value, "radius") <- radius
  attr(value, "fixed") <- fixed
  class(value) <- c("corRExpwr2Dt", "corSpatial", "corStruct")

  value
}


Initialize.corRExpwr2Dt <- function(object, data, ...)
{
   if (!is.null(attr(object, "logDet"))) { #already initialized
      return(object)
   }

   object <- Initialize.corStruct(object, data)
   nug <- attr(object, "nugget")

   val <- as.vector(object)
   if (length(val) >= 4) {		# initialized
      if (any(val[1:4] <= 0)) {
         stop("Initial values for \"corRExpwr2Dt\" parameters must be > 0")
      }
      if (nug) {				# with nugget effect
         if (length(val) == 4) {		# assuming nugget effect not given
            val <- c(val, 0.1)		# setting it to 0.1
         } else if (length(val) != 5) {
             stop("Initial values for \"corRExpwr2Dt\" parameters of wrong dimension")
         }
         if ((val[5] <= 0) || (val[5] >= 1)) {
            stop("Initial value of nugget ratio must be between 0 and 1")
         }
      } else if (length(val) != 4) {
         stop("Initial values for \"corRExpwr2Dt\" parameters of wrong dimension")
      }
   } else {
      covar <- attr(object, "covariate")
      x <- if (is.list(covar)) unlist(sapply(covar, attr, which = "dist"))
           else attr(covar, "dist")
      t1 <- if (is.list(covar)) unlist(sapply(covar, attr, which = "t1"))
            else attr(covar, "t1")
      t2 <- if (is.list(covar)) unlist(sapply(covar, attr, which = "t2"))
            else attr(covar, "t2")
      y <- abs(t2 - t1)
      val <- c(ifelse(any(x > 0), min(x[x > 0]), 1) * 0.9, 1,
               ifelse(any(y > 0), min(y[y > 0]), 1) * 0.9, 1)
      if (nug) val <- c(val, 0.1)
   }
   val[1:4] <- log(val[1:4])
   if (nug) val[5] <- log(val[5] / (1 - val[5]))
   oldAttr <- attributes(object)
   object <- val
   attributes(object) <- oldAttr
   attr(object, "factor") <- corFactor(object)
   attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

   object
}


getCovariate.corRExpwr2Dt <- function(object, form = formula(object), data)
{
   if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
      if (missing(data)) {
         stop("Need data to calculate covariate")
      }
      covForm <- getCovariateFormula(form)
      tcovar <- length(all.vars(covForm)) + c(-1, 0)
      if (tcovar[1] >= 3) { # covariates present
         if (attr(terms(covForm), "intercept") == 1) {
            covForm <- eval(parse(text = paste("~", deparse(covForm[[2]]),
                          "-1", sep = "")))
         }
         covar <- as.data.frame(unclass(model.matrix(covForm,
                     model.frame(covForm, data, drop.unused.levels = TRUE))))
         if (nrow(covar) > nrow(unique(covar))) {
            stop("Cannot have duplicate sites in \"corRExpwr2Dt\"")
         } else if (any(covar[,3] > covar[,4])) {
            stop("Temporal limits must be ascending in \"corRExpwr2Dt\"")
         }
      } else {
         covar <- NULL
      }

      if (!is.null(getGroupsFormula(form))) { # by groups
         if (!is.null(covar)) {
            grps <- getGroups(object, data = data)
            covar <- lapply(split(covar, grps),
                        function(el, metric, radius) {
                            el <- as.matrix(el)
                            r <- nrow(el)
                            if (r > 1) {
                               d <- as.vector(dist2(el[, -tcovar],
                                                    metric, r = radius))
                               x <- matrix(0, r, r)
                               idx <- lower.tri(x)
                               t1 <- el[col(x)[idx], tcovar]
                               t2 <- el[row(x)[idx], tcovar]
                            } else {
                               d <- numeric(0)
                               t1 <- numeric(0)
                               t2 <- numeric(0)
                            }
                            attr(el, "dist") <- d
                            attr(el, "t1") <- t1
                            attr(el, "t2") <- t2
                            el
                         }, metric = attr(object, "metric"),
                            radius = attr(object, "radius"))
         }
         covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups
      } else {				# no groups
         if (!is.null(covar)) {
            covar <- as.matrix(covar)
            attr(covar, "dist") <- as.vector(
                                      dist2(covar[, -tcovar],
                                            method = attr(object, "metric"),
                                            r = attr(object, "radius")))
            x <- matrix(0, nrow(covar), nrow(covar))
            idx <- lower.tri(x)
            attr(covar, "t1") <- covar[col(x)[idx], tcovar]
            attr(covar, "t2") <- covar[row(x)[idx], tcovar]
         }
      }
   }

   covar
}


coef.corRExpwr2Dt <- function(object, unconstrained = TRUE, ...)
{
   if (attr(object, "fixed") && unconstrained) {
      return(numeric(0))
   }
   val <- as.vector(object)
   if (length(val) == 0) {               # uninitialized
      return(val)
   }
   if (!unconstrained) val <- unconstrained(object, val, inv = TRUE)

   nm <- c("spatial range", "spatial shape", "temporal range", "interaction")
   if (attr(object, "nugget")) nm <- c(nm, "nugget")
   names(val) <- nm

   val
}


Dim.corRExpwr2Dt <- function(object, groups, ...)
{
   if (missing(groups)) return(attr(object, "Dim"))
   val <- Dim.corStruct(object, groups)
   val[["start"]] <-
      c(0, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
   ## will use third component of Dim list for spClass
   names(val)[3] <- "spClass"
   val[[3]] <- match(class(object)[1], c("corRExpwr2Dt"), 0)

   val
}


corMatrix.corRExpwr2Dt <- function(object, covariate = getCovariate(object),
   corr = TRUE, ...)
{
   lD <- NULL
   par <- coef(object, unconstrained = FALSE)

   corf <- function(dist, t1, t2, par, nugget = FALSE) {
      val <- cor.exp2dt(dist, t1, t2, par[1], par[2], par[3], par[4])
      if (nugget) val <- (1 - par[5]) * val
      val
   }
   nug <- attr(object, "nugget")

   if (data.class(covariate) == "list") {
      if (is.null(names(covariate))) {
         nm <- seq(covariate)
         names(covariate) <- nm
      } else {
         nm <- names(covariate)
      }
      len <- unlist(lapply(covariate, nrow))
      val <- list()
      for(i in nm) {
         el <- covariate[[i]]
         x <- matrix(0, len[i], len[i])
         x[lower.tri(x)] <- corf(attr(el, "dist"), attr(el, "t1"), attr(el, "t2"),
                                 par, nug)
         if (corr) {
            val[[i]] <- x + t(x)
            diag(val[[i]]) <- corf(0, el[,3:4], el[,3:4], par, nug)
         } else {
            diag(x) <- corf(0, el[,3:4], el[,3:4], par, nug)
            l <- chol(t(x))
            val[[i]] <- t(backsolve(l, diag(len[i])))
            lD <- c(lD, diag(l))
         }
      }
   } else {
      len <- nrow(covariate)
      x <- matrix(0, len, len)
      x[lower.tri(x)] <- corf(attr(covariate, "dist"), attr(covariate, "t1"),
                              attr(covariate, "t2"), par, attr(object, "nugget"))
      if (corr) {
         val <- x + t(x)
         diag(val) <- corf(0, covariate[,3:4], covariate[,3:4], par, nug)
      } else {
         diag(x) <- corf(0, covariate[,3:4], covariate[,3:4], par, nug)
         l <- chol(t(x))
         val <- t(backsolve(l, diag(len)))
         lD <- c(lD, diag(l))
      }
   }

   if (!is.null(lD)) lD <- -1 * sum(log(lD))
   attr(val, "logDet") <- lD

   val
}


corFactor.corRExpwr2Dt <- function(object, ...)
{
   if (!is.null(aux <- attr(object, "factor"))) {
      return(aux)
   }
   val <- corMatrix(object, corr = FALSE, ...)
   lD <- attr(val, "logDet")
   if (is.list(val)) val <- unlist(val)
   else val <- as.vector(val)
   names(val) <- NULL
   attr(val, "logDet") <- lD

   val
}


################################################################################
# Distance and correlation functions
################################################################################

dist2 <- function(x, method = c("euclidean", "maximum", "manhattan", "canberra",
   "binary", "minkowski", "haversine"), diag = FALSE, upper = FALSE,
   p = 2, r = 3956)
{
   METHOD <- match.arg(method)
   switch(METHOD,
      haversine = {
              m <- matrix(NA, nrow(x), nrow(x))
              idx <- lower.tri(m)
              m[idx] <- haversine(x[col(m)[idx],1:2], x[row(m)[idx],1:2], r)
              d <- as.dist(m, diag = diag, upper = upper)
            },
            {
              f <- get("dist", envir = as.environment("package:stats"))
              d <- f(x, method = METHOD, diag = diag, upper = upper, p = p)
            }
   )

   d
}

# Great circle distance
haversine <- function(x, y, r = 3956)
{
   if(is.vector(x)) x <- matrix(x, 1, 2)
   if(is.vector(y)) y <- matrix(y, 1, 2)

   rad <- pi / 180
   z <- sin((y - x) * (rad / 2))^2
   a <- z[,1] + cos(rad * x[,1]) * cos(rad * y[,1]) * z[,2]
   2 * r * atan2(sqrt(a), sqrt(1 - a))
}

# Powered exponential correlation function
cor.exp <- function(x, range, p = 1)
{
   if (range <= 0 || p <= 0)
      stop("Exponential correlation parameter must be > 0")

   if (p == 1) exp(x / (-1 * range))
   else exp(-1 * (x / range)^p)
}

# Linear correlation function
cor.lin <- function(x, range)
{
   if (range <= 0)
      stop("Linear correlation parameter must be > 0")

   r <- x < range
   r[r] <- 1 - x[r] / range
   r
}

# Matern correlation function
cor.matern <- function(x, range, scale)
{
   if(range <= 0 || scale <= 0)
      stop("Matern correlation parameters must be > 0")

   a <- 1 / (2^(scale - 1) * gamma(scale))
   r <- x == 0
   x0 <- x[idx <- !r] / range

   r[idx] <- a * x0^scale * besselK(x0, scale)
   r
}

# Rational quadratic correlation function
cor.ratio <- function(x, range)
{
   if (range <= 0)
      stop("Rational quadratic correlation parameter must be > 0")

   1 / ((x / range)^2 + 1)
}

# Sperical correlation function
cor.spher <- function(x, range)
{
   if (range <= 0)
      stop("Spherical correlation parameter must be > 0")

   r <- x < range
   x0 <- x[r] / range
   r[r] <- 1 - 1.5 * x0 + 0.5 * x0^3
   r
}

# Sine wave correlation function
cor.wave <- function(x, range)
{
   if (range <= 0)
      stop("Sine wave correlation parameter must be > 0")

   x0 <- x / range
   sin(x0) / x0
}

# Non-separable exponential spatio-temporal correlation function
cor.exp2 <- function(x, t, x.range, x.p = 1, t.range, t.p = 1, lambda = 0)
{
   if (t.range <= 0 || x.range <= 0 || x.p <= 0 || lambda < 0)
      stop("Exponential correlation parameters must be > 0")

   x0 <- if (x.p == 1) x / (-1 * x.range)
         else -1 * (x / x.range)^x.p
   t0 <- if (t.p == 1) t / (-1 * t.range)
         else -1 * (t / t.range)^t.p

   exp(x0 - lambda * x0 * t0 + t0) 
}

# Non-separable temporally integrated exponential spatial correlation function
cor.exp2dt <- function(x, t1, t2, x.range, x.p = 1, t.range, lambda = 0)
{
   if (t.range <= 0 || x.range <= 0 || x.p <= 0 || lambda < 0)
      stop("Exponential correlation parameters must be > 0")

   if (is.vector(t1)) t1 <- matrix(t1, 1, 2)
   if (is.vector(t2)) t2 <- matrix(t2, 1, 2)

   x0 <- if (x.p == 1) x / (-1 * x.range)
         else -1 * (x / x.range)^x.p

   overlap <- pmin(t1[,2], t2[,2]) - pmax(t1[,1], t2[,1])
   overlap[overlap < 0] <- 0
   norm <- (t1 %*% c(-1, 1)) * (t2 %*% c(-1, 1))

   if (lambda == 0) theta <- t.range
   else theta <- t.range / (1 - lambda * x0)

   val <- (theta^2 * exp(abs(t1[, c(1,1,2,2)] - t2[, c(2,1,2,1)]) /
             (-1 * theta)) %*% c(1, -1, -1, 1) + 2 * theta * overlap) / norm
   exp(x0) * as.vector(val)
}
