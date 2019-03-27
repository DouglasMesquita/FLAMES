#' @title MCMC for binary regression - ARMS version
#'
#' @description Perform MCMC with ARMS algorithm
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param nsim Sample size required for MCMC
#' @param burnin Burn in for MCMC
#' @param lag Lag for MCMC
#' @param inv_link_f Inverse link function
#' @param type "logit", "probit", "cauchit", "robit" or "cloglog"
#' @param sample_c Should c be sampled?
#' @param sigma_beta Variance of beta prior
#' @param a_c Shape1 for c prior
#' @param b_c Shape2 for c prior
#' @param a_lambda Inferior limit for lambda
#' @param b_lambda Superior limit for lambda
#' @param bound_beta Limits to sample beta
#' @param p_c To restore c
#' @param p_prop To restore p
#' @param p_beta To restore beta
#' @param p_df To restore df
#' @param p_lambda To restore lambda
#' @param const A constant to help on sampling degrees of freedom \eqn{\tilde{df} = df/c}
#'
#' @return Chains of all parameters
#'
#' @import HI
#' @import progress
#' @import stats
#'
#' @export

mcmc_bin_arms <- function(y, X,
                          nsim, burnin, lag,
                          inv_link_f,
                          type, sample_c,
                          sigma_beta, a_c, b_c, a_lambda, b_lambda,
                          bound_beta,
                          p_c, p_prop, p_beta, p_df, p_lambda,
                          const){

  sample_size <- burnin + lag*nsim

  ##-- Indicator functions for HI package
  lim_beta <- exp(bound_beta)/(1+exp(bound_beta))

  ind_fun_beta <- function(x) (x > lim_beta[1])*(x < lim_beta[2])
  ind_fun_c <- function(x) (x > 0)*(x < 1)
  ind_fun_df <- function(x) (x > 0)*(x < 1)
  ind_fun_lambda <- function(x) (x > a_lambda)*(x < b_lambda)

  n_cov <- ncol(X)

  beta_aux <- p_beta[[1]]
  c_aux <- ifelse(!is.null(p_c), p_c[[1]], 0)
  if(!is.null(p_df)) df_aux <- p_df[[1]]
  if(!is.null(p_lambda)) lambda_aux <- p_lambda[[1]]

  ##-- Progress bar
  pb <- progress::progress_bar$new(total = sample_size-1)

  ##-- Loop ----
  for(i in 2:sample_size){
    pb$tick()

    ##-- Regression coefficients
    for(j in 1:n_cov){
      y_start <- exp(beta_aux[j])/(1+exp(beta_aux[j]))
      beta_aux_j <- HI::arms(y.start = y_start,
                             myldens = function(x) beta_posterior(y = y, X = X,
                                                                  p_beta = beta_aux, p_beta_element = x, element = j,
                                                                  p_c = c_aux, p_df = df_aux,
                                                                  inv_link_f = inv_link_f,
                                                                  sigma_beta = sigma_beta,
                                                                  log = TRUE, method = "ARMS"),
                             indFunc = ind_fun_beta,
                             n.sample = 1)

      beta_aux[j] <- log(beta_aux_j/(1-beta_aux_j))
    }

    ##-- c parameter
    if(sample_c){
      y_start <- c_aux
      c_aux <- HI::arms(y.start = y_start,
                        myldens = function(x) c_posterior(y = y, X = X,
                                                          p_beta = beta_aux,
                                                          p_c = x, p_df = df_aux,
                                                          inv_link_f = inv_link_f,
                                                          a_c = a_c, b_c = b_c,
                                                          log = TRUE),
                        indFunc = ind_fun_c,
                        n.sample = 1)
    }

    if(type == "robit"){
      ##-- lambda parameter
      y_start <- lambda_aux
      lambda_aux <- HI::arms(y.start = y_start,
                             myldens = function(x) lambda_posterior(p_df = df_aux, p_lambda = x,
                                                                    inv_link_f = inv_link_f,
                                                                    log = TRUE),
                             indFunc = ind_fun_lambda,
                             n.sample = 1)

      ##-- df parameter
      y_start <- 1-exp(-df_aux/const)
      df_aux <- HI::arms(y.start = y_start ,
                         myldens = function(x) df_posterior(y = y, X = X,
                                                            p_beta = beta_aux, p_c = c_aux, p_df = x, p_lambda = lambda_aux,
                                                            inv_link_f = inv_link_f,
                                                            log = TRUE, method = "ARMS", const = const),
                         indFunc = ind_fun_df,
                         n.sample = 1)

      df_aux <- -const*log(1-df_aux)
    }

    ##-- Saving the iteration i
    if((i-burnin)%%lag == 0 & i > burnin){
      pos <- (i-burnin)/lag

      p_beta[[pos]] <- beta_aux
      if(!is.null(p_c)) p_c[[pos]] <- c_aux
      if(!is.null(p_df)) p_df[[pos]] <- df_aux
      if(!is.null(p_lambda)) p_lambda[[pos]] <- lambda_aux

      p_prop[[pos]] <- make_p(p_c = c_aux, X = X, p_beta = beta_aux, p_df = df_aux, inv_link_f = inv_link_f)

      invisible(gc(reset = T, verbose = FALSE, full = TRUE))
    }
  }

  ##-- Outputs
  return(list(p_c = p_c, p_beta = p_beta, p_prop = p_prop, p_df = p_df, p_lambda = p_lambda))
}
