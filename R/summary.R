#' Summary for FLAMES class (coda is necessary)
#'
#' @param object FLAMES object to summarise
#' @param HPD HPD intervals or quantiles (TRUE or FALSE)
#' @param ... Other parameters for summary
#'
#' @method summary FLAMES
#'
#' @export

summary.FLAMES <- function(object, HPD = T, ...){

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

  if(!is.null(object[["c"]]) | !is.null(object[["d"]]) | !is.null(object[["df"]])){
    other_info <- data.frame()

    if(!is.null(object[["c"]])){
      mean_c <- mean(object$c)
      sd_c <- stats::sd(object$c)

      if(HPD){
        interval_c <- coda::HPDinterval(coda::as.mcmc(object$c))
      } else{
        interval_c <- stats::quantile(x = object$c, probs = c(0.025, 0.975))
      }

      data_c <- data.frame(mean = mean_c, sd = sd_c, lower_95 = interval_c[1], upper_95 = interval_c[2], row.names = "c parameter")
      other_info <- rbind.data.frame(other_info, data_c)
    }

    if(!is.null(object[["d"]])){
      mean_d <- mean(object$d)
      sd_d <- stats::sd(object$d)

      if(HPD){
        interval_d <- coda::HPDinterval(coda::as.mcmc(object$d))
      } else{
        interval_d <- stats::quantile(x = object$d, probs = c(0.025, 0.975))
      }

      data_d <- data.frame(mean = mean_d, sd = sd_d, lower_95 = interval_d[1], upper_95 = interval_d[2], row.names = "d parameter")
      other_info <- rbind.data.frame(other_info, data_d)
    }

    if(!is.null(object[["df"]])){
      mean_df <- mean(object$df)
      sd_df <- stats::sd(object$df)

      if(HPD){
        interval_df <- coda::HPDinterval(coda::as.mcmc(object$df))
      } else{
        interval_df <- stats::quantile(x = object$df, probs = c(0.025, 0.975))
      }

      data_df <- data.frame(mean = mean_df, sd = sd_df, lower_95 = interval_df[1], upper_95 = interval_df[2], row.names = "Degrees of freedom")
      other_info <- rbind.data.frame(other_info, data_df)
    }

    ans$'Other parameters' <- other_info
  }

  if(!is.null( object$fit_measures)) ans$'Fit measures' <- object$fit_measures

  return(ans)
}
