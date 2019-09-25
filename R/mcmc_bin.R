#' @title MCMC for binary regression
#'
#' @description MCMC for binary regression
#'
#' @param data Dataset to be used
#' @param formula Simple formula
#' @param nsim Sample size required for MCMC
#' @param burnin Burn in for MCMC
#' @param lag Lag for MCMC
#' @param type "logit", "probit", "cauchit", "robit", "cloglog" or "loglog"
#' @param sample_c Should c be sampled?
#' @param sample_d Should d be sampled?
#' @param sigma_beta Variance of beta prior
#' @param mean_c Mean for c a priori
#' @param sd_c Standard deviation for c a priori
#' @param mean_d Mean for d a priori
#' @param sd_d Standard deviation for d a priori
#' @param a_lambda Inferior limit for lambda
#' @param b_lambda Superior limit for lambda
#' @param var_df Variance to sample 1-exp(-df/const)
#' @param var_c Variance to sample from c (if sample_d = TRUE) otherwise log(c/(1-c))
#' @param var_d Variance to sample d (if sample_c = TRUE) otherwise log(c/(1-c))
#' @param var_lambda Variance to sample lambda
#' @param method "metropolis" or "ARMS"
#' @param const A constant to help on sampling degrees of freedom \eqn{\tilde{df} = df/c}
#' @param const_beta A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_c A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_d A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_df A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_lambda A constant to tunning the acceptance rate (default = 2.38^2)
#' @param fitm Should return fit measures? "full" for measures based on full dataset or "bootstrap" to use a bootstrap technique
#'
#' @examples \dontrun{
#' ##-- Seed ----
#' set.seed(123456)
#'
#' ##-- Data ----
#' n <- 1000
#' n_cov <- 2
#'
#' ##-- Covariates
#' X <- matrix(rt(n*n_cov, df = 2.5), ncol = n_cov)
#' X <- scale(X, scale = FALSE)
#'
#' ##-- Coefficients
#' betas <- c(0, -1, 0.5)
#' XBeta <- cbind(1, X)%*%betas
#'
#' ##-- c parameter
#' c1 <- 0.25
#' d1 <- 0.95
#'
#' type_data = "cloglog"
#'
#' ##-- p and y
#' p <- robit:::inv_link(x = XBeta, type = type_data)*(d1-c1) + c1
#' y <- rbinom(n = n, size = 1, prob = p)
#'
#' bd <- data.frame(y = y, X)
#'
#' ##-- Hyperparameters (prioris)
#' sigma_beta <- 3
#'
#' mean_c <- robit:::mean_sd_beta(a = 1, b = 2)$mean
#' sd_c <- robit:::mean_sd_beta(a = 1, b = 2)$sd
#'
#' mean_d <- robit:::mean_sd_beta(a = 2, b = 1)$mean
#' sd_d <- robit:::mean_sd_beta(a = 2, b = 1)$sd
#'
#' a_lambda <- 0.01
#' b_lambda <- 0.99
#'
#' ##-- MCMC
#' nsim <- 1000
#' burnin <- 5000
#' lag <- 10
#'
#' f <- y ~ X1 + X2
#' type <- "cloglog"
#'
#' ##-- ARMS
#' out_arms <- mcmc_bin(data = bd, formula = f,
#'                      nsim = nsim, burnin = burnin, lag = lag,
#'                      type = type, sample_c = TRUE, sample_d = TRUE,
#'                      sigma_beta = sigma_beta,
#'                      mean_c = mean_c, sd_c = sd_c, mean_d = mean_d, sd_d = sd_d,
#'                      a_lambda = a_lambda, b_lambda = b_lambda,
#'                      method = "ARMS", const = 50, fitm = "full")
#'
#' ##-- ARMS
#' out_met <- mcmc_bin(data = bd, formula = f,
#'                     nsim = nsim, burnin = burnin, lag = lag,
#'                     type = type, sample_c = TRUE, sample_d = TRUE,
#'                     sigma_beta = sigma_beta,
#'                     mean_c = mean_c, sd_c = sd_c, mean_d = mean_d, sd_d = sd_d,
#'                     a_lambda = a_lambda, b_lambda = b_lambda,
#'                     method = "metropolis", const = 50, fitm = "full")
#'
#' ##-- GLM
#' out_glm <- glm(formula = f, data = bd, family = "binomial")
#'
#' summary(out_glm)
#' summary(out_arms)
#' summary(out_met)
#'
#' coef(out_glm)
#' coef(out_arms)
#' coef(out_met)
#' betas
#'
#' plot(out_arms, ask = T)
#' plot(out_met, ask = T)
#' }
#'
#' @return Chains of all parameters
#' @export

mcmc_bin <- function(data, formula,
                     nsim = 1000, burnin = round(0.1*nsim), lag = 10,
                     type = "logit", sample_c = TRUE, sample_d = TRUE,
                     sigma_beta = 3,
                     mean_c = FLAMES:::mean_sd_beta(a = 1, b = 2)$mean, sd_c = FLAMES:::mean_sd_beta(a = 1, b = 2)$sd,
                     mean_d = FLAMES:::mean_sd_beta(a = 2, b = 1)$mean, sd_d = FLAMES:::mean_sd_beta(a = 2, b = 1)$sd,
                     a_lambda = 0.01, b_lambda = 0.99,
                     var_df = 0.02, var_c = ifelse(sample_d, 0.005, 0.02), var_d = ifelse(sample_c, 0.005, 0.02), var_lambda = 0.05,
                     method = "ARMS",
                     const = 1, const_beta = 2.38^2, const_c = 2.38^2, const_d = 2.38^2, const_df = 2.38^2, const_lambda = 2.38^2,
                     fitm = FALSE){

  if(!is.data.frame(data)) stop("data must be a data.frame")
  #if(burnin + nsim*lag < 1000) stop("Please consider to increase the nsim, burnin and/or lag")
  if(a_lambda < 0) stop("a_lambda must be > 0")
  if(b_lambda < a_lambda) stop("a_lambda must be < b_lambda")

  ##-- c and d
  par_c <- mean_sd_beta(mean = mean_c, sd = sd_c, show_warnings = TRUE)
  a_c <- par_c$a
  b_c <- par_c$b

  par_d <- mean_sd_beta(mean = mean_d, sd = sd_d, show_warnings = TRUE)
  a_d <- par_d$a
  b_d <- par_d$b

  ##-- Getting call, y and X
  call_FLAMES <- match.call()

  model_fr <- match.call(expand.dots = FALSE)
  match_strings <- match(c("formula", "data"), names(model_fr), 0L)
  model_fr <- model_fr[c(1L, match_strings)]
  model_fr[[1L]] <- quote(stats::model.frame)
  model_fr <- eval(model_fr, parent.frame())
  model_types <- attr(model_fr, "terms")
  y <- stats::model.response(model_fr, "numeric")
  X <- stats::model.matrix(model_types, model_fr)
  rm(data)

  ##-- Link function
  inv_link_f <- function(x, df) inv_link(x = x, type = type, df = df)
  link_f <- function(x, df) link(x = x, type = type, df = df)

  n_cov <- ncol(X)
  n <- nrow(X)

  p_prop <- p_beta <- p_lambda <- p_df <- vector(mode = 'list', length = nsim)
  if(sample_c) p_c <- vector(mode = 'list', length = nsim) else p_c <- NULL
  if(sample_d) p_d <- vector(mode = 'list', length = nsim) else p_d <- NULL
  if(type == "robit") p_df <- vector(mode = 'list', length = nsim) else p_df <- NULL
  if(type == "robit") p_lambda <- vector(mode = 'list', length = nsim) else p_lambda <- NULL

  ##-- Start values (as input soon)
  p_prop[[1]] <- rep(0.5, n)
  p_beta[[1]] <- rep(0, n_cov)

  if(sample_c) p_c[[1]] <- ifelse(!sample_c, 0, 0.25)
  if(sample_d) p_d[[1]] <- ifelse(!sample_d, 1, 0.95)
  if(type == "robit") p_lambda[[1]] <- mean(c(a_lambda, b_lambda))
  if(type == "robit") p_df[[1]] <- 10

  ##-- MCMC
  time_1 <- Sys.time()
  if(method == "metropolis"){
    samp <- mcmc_bin_metropolis(y = y, X = X,
                                nsim = nsim, burnin = burnin, lag = lag,
                                inv_link_f = inv_link_f,
                                type = type, sample_c = sample_c, sample_d = sample_d,
                                sigma_beta = sigma_beta,
                                a_c = a_c, b_c = b_c, a_d = a_d, b_d = b_d,
                                a_lambda = a_lambda, b_lambda = b_lambda,
                                var_df = var_df, var_c = var_c, var_d = var_d, var_lambda = var_lambda,
                                p_c = p_c, p_d = p_d, p_prop = p_prop, p_beta = p_beta, p_df = p_df, p_lambda = p_lambda,
                                const = const,
                                const_beta = const_beta, const_c = const_c, const_d = const_d, const_df = const_df, const_lambda = const_lambda)

  } else{
    if(method == "ARMS"){
      samp <- mcmc_bin_arms(y = y, X = X,
                            nsim = nsim, burnin = burnin, lag = lag,
                            inv_link_f = inv_link_f,
                            type = type, sample_c = sample_c, sample_d = sample_d,
                            sigma_beta = sigma_beta,
                            a_c = a_c, b_c = b_c, a_d = a_d, b_d = b_d,
                            a_lambda = a_lambda, b_lambda = b_lambda,
                            p_c = p_c, p_d = p_d, p_prop = p_prop, p_beta = p_beta, p_df = p_df, p_lambda = p_lambda,
                            const = const)

    } else{
      stop(paste("The", method, "algorithm is not implemented."))
    }
  }

  samp$p_prop <- do.call(args = samp$p_prop, what = "rbind")
  samp$p_beta <- do.call(args = samp$p_beta, what = "rbind")
  if(sample_c) samp$p_c <- do.call(args = samp$p_c, what = "c")
  if(sample_d) samp$p_d <- do.call(args = samp$p_d, what = "c")
  if(type == "robit") samp$p_df <- do.call(args = samp$p_df, what = "c")
  if(type == "robit") samp$p_lambda <- do.call(args = samp$p_lambda, what = "c")

  time_2 <- Sys.time()
  time_elapsed <- time_2 - time_1

  ##-- Outputs
  colnames(samp$p_beta) <- colnames(X)

  if(!is.null(fitm)){
    if(fitm == "full") nrep <- NULL else nrep <- 100
    fit <- fit_measures(y = y, p = as.matrix(samp$p_prop), nrep = nrep, nsamp = 100)
  } else{
    fit <- NULL
  }

  invisible(gc(reset = T, verbose = FALSE, full = TRUE))

  out <- list(p = samp$p_prop, beta = samp$p_beta,
              c = samp$p_c, d = samp$p_d, df = samp$p_df, lambda = samp$p_lambda,
              time = time_elapsed,
              call = call_FLAMES,
              fit_measures = fit)

  class(out) <- "FLAMES"

  return(out)
}
