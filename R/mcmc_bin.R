#' @title MCMC for binary regression
#'
#' @description MCMC for binary regression
#'
#' @param data Dataset to be used
#' @param formula Simple formula
#' @param nsim Sample size required for MCMC
#' @param burnin Burn in for MCMC
#' @param lag Lag for MCMC
#' @param type "logit", "probit", "cauchit", "tobit" or "cloglog"
#' @param sample_c Should c be sampled?
#' @param sigma_beta Variance of beta prior
#' @param a_c Shape1 for c prior (beta)
#' @param b_c Shape2 for c prior (beta)
#' @param a_lambda Inferior limit for lambda
#' @param b_lambda Superior limit for lambda
#' @param var_df Variance to sample log(df)
#' @param var_c Variance to sample c
#' @param var_lambda Variance to sample lambda
#' @param bound_beta Limits to sample beta for ARMS
#' @param method "metropolis" or "ARMS"
#'
#' @examples if(interactive()){
#'  set.seed(1)
#'
#'  ##-- Data ----
#'  n <- 1000
#'  n_cov <- 2
#'
#'  ##-- Covariates
#'  X <- matrix(rnorm(n*n_cov), ncol = n_cov)
#'
#'  ##-- Coefficients
#'  betas <- c(-1, 0.5)
#'  XBeta <- X%*%betas
#'
#'  ##-- c parameter
#'  c1 <- 0.2
#'  c2 <- 1 - c1
#'
#'  type_data = "tobit"
#'  df <- 10
#'
#'  ##-- p and y
#'  p <- tbreg:::inv_link(x = XBeta, type = type_data, df = df)*c2 + c1
#'  y <- rbinom(n = n, size = 1, prob = p)
#'
#'  bd <- data.frame(y = y, X)
#'
#'  ##-- Hyperparameters (prioris)
#'  sigma_beta <- 100
#'
#'  a_c <- 0.1  ## non informative beta
#'  b_c <- 0.1  ## non informative beta
#'
#'  a_lambda <- 0.01
#'  b_lambda <- 1.00
#'
#'  ##-- MCMC
#'  nsim <- 1000
#'  burnin <- 5000
#'  lag <- 10
#'
#'  type <- "tobit"
#'  bound_beta <- c(-10, 10)
#'  var_df <- 0.5
#'  n_points <- 6
#'
#'  f <- y ~ X1 + X2
#'
#'  ##-- ARMS ~ 3 minutes
#'  out_arms <- mcmc_bin(data = bd, formula = f,
#'                      nsim = nsim, burnin = burnin, lag = lag,
#'                      type = type, sample_c = TRUE,
#'                      sigma_beta = sigma_beta, a_c = a_c, b_c = b_c,
#'                      a_lambda = a_lambda, b_lambda = b_lambda,
#'                      var_df = var_df, bound_beta = bound_beta,
#'                      method = "ARMS", force_intercept = TRUE)
#'
#'  ##-- ARMS ~ 45 seconds
#'  out_met <- mcmc_bin(data = bd, formula = f,
#'                     nsim = nsim, burnin = burnin, lag = lag,
#'                     type = type, sample_c = TRUE,
#'                     sigma_beta = sigma_beta, a_c = a_c, b_c = b_c,
#'                     a_lambda = a_lambda, b_lambda = b_lambda,
#'                     var_df = var_df, bound_beta = bound_beta,
#'                     method = "metropolis")
#'
#'  ##-- GLM
#'  out_glm <- glm(formula = f, data = bd, family = "binomial")
#'
#'  summary(out_glm)
#'  summary(out_arms)
#'  summary(out_met)
#'
#'  coef(out_glm)
#'  coef(out_arms)
#'  coef(out_met)
#'  betas
#'
#'  plot(out_arms, ask = T)
#'  plot(out_met, ask = T)
#' }
#'
#' @return Chains of all parameters
#' @export

mcmc_bin <- function(data, formula,
                     nsim = 1000, burnin = round(0.1*nsim), lag = 10,
                     type = "logit", sample_c = TRUE,
                     sigma_beta = 100, a_c = 0.01, b_c = 0.01,
                     a_lambda = 0.01, b_lambda = 0.99,
                     var_df = 0.04, var_c = 0.02, var_lambda = 0.2,
                     bound_beta,
                     method = "ARMS", force_intercept = FALSE){

  if(!is.data.frame(data))
    stop("data must be a data.frame")

  if(burnin + nsim*lag < 1000)
    stop("Please consider to increase the nsim, burnin and/or lag.")

  ##-- Getting call, y and X
  call_tbreg <- match.call()

  model_fr <- match.call(expand.dots = FALSE)
  match_strings <- match(c("formula", "data"), names(model_fr), 0L)
  model_fr <- model_fr[c(1L, match_strings)]
  model_fr[[1L]] <- quote(stats::model.frame)
  model_fr <- eval(model_fr, parent.frame())
  model_types <- attr(model_fr, "terms")
  y <- stats::model.response(model_fr, "numeric")
  X <- stats::model.matrix(model_types, model_fr)

  if(sample_c & "(Intercept)" %in% colnames(X)){
    if(!force_intercept){
      X <- X[, -1]
      warning("Intercept was removed since you want to sample c")
    }

    if(ncol(X) == 0) stop("You need at least one covariate. We are working on this limitation.")
  }

  ##-- Link function
  inv_link_f <- function(x, df) inv_link(x = x, type = type, df = df)
  link_f <- function(x, df) link(x = x, type = type, df = df)

  sample_size <- burnin + lag*nsim
  n_cov <- ncol(X)
  n <- nrow(X)

  p_c <- rep(0.5, sample_size)
  if(!sample_c) p_c <- rep(0, sample_size)

  p_prop <- matrix(0.5, nrow = sample_size, ncol = n)
  p_beta <- matrix(0, nrow = sample_size, ncol = n_cov)
  p_df <- rep(5, sample_size)
  p_lambda <- rep(0.2, sample_size)

  time_1 <- Sys.time()
  if(method == "metropolis"){
    samp <- mcmc_bin_metropolis(y = y, X = X,
                                sample_size = sample_size, burnin = burnin,
                                inv_link_f = inv_link_f,
                                type = type, sample_c = sample_c,
                                sigma_beta = sigma_beta,
                                a_c = a_c, b_c = b_c,
                                a_lambda = a_lambda, b_lambda = b_lambda,
                                var_df = var_df, var_c = var_c, var_lambda = var_lambda,
                                p_c = p_c, p_prop = p_prop, p_beta = p_beta, p_df = p_df, p_lambda = p_lambda)

    p_c <- samp$p_c
    p_prop <- samp$p_prop
    p_beta <- samp$p_beta
    p_df <- samp$p_df
    p_lambda <- samp$p_lambda

  } else{
    if(method == "ARMS"){
      samp <- mcmc_bin_arms(y = y, X = X,
                            sample_size = sample_size,
                            inv_link_f = inv_link_f,
                            type = type, sample_c = sample_c,
                            sigma_beta = sigma_beta,
                            a_c = a_c, b_c = b_c,
                            a_lambda = a_lambda, b_lambda = b_lambda,
                            bound_beta = bound_beta,
                            p_c = p_c, p_prop = p_prop, p_beta = p_beta, p_df = p_df, p_lambda = p_lambda)

      p_c <- samp$p_c
      p_prop <- samp$p_prop
      p_beta <- samp$p_beta
      p_df <- samp$p_df
      p_lambda <- samp$p_lambda

    } else{
      stop(paste("The", method, "algorithm is not implemented."))
    }
  }
  time_2 <- Sys.time()
  time_elapsed <- time_2 - time_1

  ##-- Outputs
  pos <- seq((burnin + lag), sample_size, by = lag)

  p_beta <- p_beta[pos, ]
  p_prop <- p_prop[pos, ]
  p_c <- p_c[pos]
  p_df <- p_df[pos]
  p_lambda <- p_lambda[pos]

  colnames(p_beta) <- colnames(X)

  out <- list(p = p_prop, beta = p_beta,
              c = p_c, df = p_df, lambda = p_lambda,
              time = time_elapsed,
              call = call_tbreg)

  class(out) <- "tbreg"

  return(out)
}
