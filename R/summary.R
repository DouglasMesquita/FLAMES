#' Summary for crobit class (coda is necessary)
#'
#' @param object crobit object to summarise
#' @param HPD HPD intervals or quantiles (TRUE or FALSE)
#' @param ... Other parameters for summary
#'
#' @method summary crobit
#'
#' @export

summary.crobit <- function(object, HPD = T, ...){

  ans <- list()
  ans$'Call' <- object$call

  object$beta <- as.matrix(object$beta)

  if(HPD){
    hpd_interval <- apply(X = object$beta, MARGIN = 2, function(x) coda::HPDinterval(coda::as.mcmc(x)))

    coef_info <- data.frame('mean' = colMeans(object$beta),
                            'std_error' = apply(X = object$beta, MARGIN = 2, stats::sd),
                            'lower_95' = hpd_interval[1,],
                            'upper_95' = hpd_interval[2,])
  } else{
    coef_info <- data.frame('mean' = colMeans(object$beta),
                            'std_error' = apply(X = object$beta, MARGIN = 2, stats::sd),
                            'lower_95' = apply(X = object$beta, MARGIN = 2, function(x) stats::quantile(x = x, probs = 0.025)),
                            'upper_95' = apply(X = object$beta, MARGIN = 2, function(x) stats::quantile(x = x, probs = 0.975)))
  }

  ans$'Coeficients' <- coef_info

  if(HPD){
    hpd_c <- coda::HPDinterval(coda::as.mcmc(object$c))
    hpd_df <- coda::HPDinterval(coda::as.mcmc(object$df))

    other_info <- data.frame('mean' = c(mean(object$c), mean(object$df)),
                             'std_error' = c(stats::sd(object$c), stats::sd(object$df)),
                             'lower_95' = c(hpd_c[1],
                                            hpd_df[1]),
                             'upper_95' = c(hpd_c[2],
                                            hpd_df[2]))
  } else{
    other_info <- data.frame('mean' = c(mean(object$c), mean(object$df)),
                             'std_error' = c(stats::sd(object$c), stats::sd(object$nu)),
                             'lower_95' = c(stats::quantile(x = object$c, probs = 0.025),
                                            stats::quantile(x = object$df, probs = 0.025)),
                             'upper_95' = c(stats::quantile(x = object$c, probs = 0.975),
                                            stats::quantile(x = object$df, probs = 0.975)))
  }

  rownames(other_info) <- c("c parameter", "Degrees of freedom")
  ans$'Other parameters' <- other_info

  return(ans)
}
