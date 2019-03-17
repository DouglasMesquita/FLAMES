#' @title MCMC for binary regression - ARMS version
#'
#' @description Perform MCMC with ARMS algorithm
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param sample_size Sample size required for MCMC
#' @param inv_link_f Inverse link function
#' @param type "logit", "probit", "cauchit", "tobit" or "cloglog"
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
#'
#' @return Chains of all parameters
#'
#' @import HI
#' @import progress
#'
#' @export

mcmc_bin_arms <- function(y, X,
                          sample_size,
                          inv_link_f,
                          type, sample_c,
                          sigma_beta, a_c, b_c, a_lambda, b_lambda,
                          bound_beta,
                          p_c, p_prop, p_beta, p_df, p_lambda){

  ##-- Indicator functions for HI package
  lim_beta <- exp(bound_beta)/(1+exp(bound_beta))

  ind_fun_beta <- function(x) (x > lim_beta[1])*(x < lim_beta[2])
  ind_fun_c <- function(x) (x > 0)*(x < 1)
  ind_fun_df <- function(x) (x > 0)*(x < 1)
  ind_fun_lambda <- function(x) (x > a_lambda)*(x < b_lambda)

  n_cov <- ncol(X)

  pb <- progress::progress_bar$new(total = sample_size-1)
  ##-- Loop ----
  for(i in 2:sample_size){
    pb$tick()

    ##-- Regression coefficients
    aux_beta <- p_beta[(i-1), ]

    for(j in 1:n_cov){
      y_start <- exp(p_beta[(i-1), j])/(1+exp(p_beta[(i-1), j]))

      aux_beta_j <- HI::arms(y.start = y_start,
                             myldens = function(x) beta_posterior(y = y, X = X,
                                                                  p_beta = aux_beta, p_beta_element = x, element = j,
                                                                  p_c = p_c[i-1], p_df = p_df[i-1],
                                                                  inv_link_f = inv_link_f,
                                                                  sigma_beta = sigma_beta,
                                                                  log = TRUE, method = "ARMS"),
                             indFunc = ind_fun_beta,
                             n.sample = 1)

      aux_beta[j] <- log(aux_beta_j/(1-aux_beta_j))

      p_beta[i, j] <- aux_beta[j]
    }

    ##-- c parameter
    if(sample_c){
      y_start <- p_c[(i-1)]
      p_c[i] <- HI::arms(y.start = y_start,
                         myldens = function(x) c_posterior(y = y, X = X,
                                                           p_beta = p_beta[i, ],
                                                           p_c = x, p_df = p_df[i-1],
                                                           inv_link_f = inv_link_f,
                                                           a_c = a_c, b_c = b_c,
                                                           log = TRUE),
                         indFunc = ind_fun_c,
                         n.sample = 1)
    }

    if(type == "tobit"){
      ##-- lambda parameter
      y_start <- p_lambda[i-1]
      aux_lambda <- HI::arms(y.start = y_start,
                             myldens = function(x) lambda_posterior(p_df = p_df[i-1], p_lambda = x,
                                                                    inv_link_f = inv_link_f,
                                                                    log = TRUE),
                             indFunc = ind_fun_lambda,
                             n.sample = 1)

      p_lambda[i] <- aux_lambda

      ##-- df parameter
      const <- 50

      y_start <- 1-exp(-p_df[i-1]/const)
      aux_df <- HI::arms(y.start = y_start ,
                         myldens = function(x) df_posterior(y = y, X = X,
                                                            p_beta = p_beta[i, ], p_c = p_c[i], p_df = x, p_lambda = p_lambda[i],
                                                            inv_link_f = inv_link_f,
                                                            log = TRUE, method = "ARMS", const = const),
                         indFunc = ind_fun_df,
                         n.sample = 1)

      p_df[i] <- -const*log(1-aux_df)
    }

    p_prop[i, ] <- make_p(p_c = p_c[i], X = X, p_beta = p_beta[i, ], p_df = p_df[i], inv_link_f = inv_link_f)
  }

  ##-- Outputs
  return(list(p_c = p_c, p_beta = p_beta, p_prop = p_prop, p_df = p_df, p_lambda = p_lambda))
}
