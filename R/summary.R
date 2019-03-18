#' Summary for robit class (coda is necessary)
#'
#' @param object robit object to summarise
#' @param HPD HPD intervals or quantiles (TRUE or FALSE)
#' @param ... Other parameters for summary
#'
#' @method summary robit
#'
#' @export

summary.robit <- function(object, HPD = T, ...){

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
    if(!is.null(object[["c"]])){
      hpd_c <- coda::HPDinterval(coda::as.mcmc(object$c))
      hpd_df <- coda::HPDinterval(coda::as.mcmc(object$df))

      lower <- c(hpd_c[1], hpd_df[1])
      upper <- c(hpd_c[2], hpd_df[2])

      mean_vec <- c(mean(object$c), mean(object$df))
      sd_vec <- c(stats::sd(object$c), stats::sd(object$df))
    } else{
      hpd_df <- coda::HPDinterval(coda::as.mcmc(object$df))

      lower <- hpd_df[1]
      upper <- hpd_df[1]

      mean_vec <- mean(object$df)
      sd_vec <- stats::sd(object$df)
    }

    other_info <- data.frame('mean' = mean_vec,
                             'std_error' = sd_vec,
                             'lower_95' = lower,
                             'upper_95' = upper)
  } else{
    if(!is.null(object[["c"]])){
      quantile_c <- stats::quantile(x = object$c, probs = c(0.025, 0.975))
      quantile_df <- stats::quantile(x = object$df, probs = c(0.025, 0.975))

      lower <- c(quantile_c[1], quantile_df[1])
      upper <- c(quantile_c[2], quantile_df[2])

      mean_vec <- c(mean(object$c), mean(object$df))
      sd_vec <- c(stats::sd(object$c), stats::sd(object$df))
    } else{
      quantile_df <- stats::quantile(x = object$df, probs = c(0.025, 0.975))

      lower <- quantile_df[1]
      upper <- quantile_df[2]

      mean_vec <- mean(object$df)
      sd_vec <- stats::sd(object$df)
    }

    other_info <- data.frame('mean' = mean_vec,
                             'std_error' = sd_vec,
                             'lower_95' = lower,
                             'upper_95' = upper)
  }

  if(!is.null(object[["c"]])){
    rownames(other_info) <- c("c parameter", "Degrees of freedom")
  } else{
    rownames(other_info) <- "Degrees of freedom"
  }

  ans$'Other parameters' <- other_info
  ans$'Fit measures' <- object$fit_measures

  return(ans)
}
