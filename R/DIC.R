DIC <- function(object, ...)
{
   UseMethod("DIC")
}

DIC.ramps <- function(object, ...)
{
   ## Posterior parameter means
   w <- rep(1 / nrow(object$params), nrow(object$params))
   params <- as.vector(crossprod(w, object$params))
   loglik <- as.numeric(crossprod(w, object$loglik))

   switch(object$control$mpdfun,
      mpdbeta = {
         val <- mpdbeta(params, object$y, object$xmat, object$kmat, object$wmat,
                        object$correlation, object$etype, object$ztype,
                        object$retype, object$weights, object$control)
      },
      mpdbetaz = {
         sites <- unique.sites(as.matrix(object$kmat))
         xk1mat <- cBind(object$xmat,
                      as(model.matrix(~ factor(sites$idx) - 1), "dgCMatrix"))
         k2mat <- as(sites$coords, "dgCMatrix")

         val <- mpdbetaz(params, object$y, xk1mat, k2mat, object$wmat,
                         object$correlation, object$etype, object$ztype,
                         object$retype, object$weights, object$control)
      }
   )
   loglik.mu <-  -0.5 * length(object$y) * log(2.0 * pi) - val$logsqrtdet -
                    val$quadform[1] / 2.0

   ## Deviance information criterion
   pD <- -2.0 * (loglik - loglik.mu)
   c(DIC = pD - 2.0 * loglik, pD = pD)
}
