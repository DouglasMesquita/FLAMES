#' @title MCMC for binary regression - Metropolis version
#'
#' @description Perform MCMC with Metropolis Hastings algorithm
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
#' @param var_df Variance to sample log(df)
#' @param var_c Variance to sample c
#' @param var_lambda Variance to sample lambda
#' @param p_c To restore c
#' @param p_prop To restore p
#' @param p_beta To restore beta
#' @param p_df To restore df
#' @param p_lambda To restore lambda
#'
#' @return Chains of all parameters
#'
#' @import progress
#' @import stats
#'
#' @export

mcmc_bin_metropolis <- function(y, X,
                                nsim, burnin, lag,
                                inv_link_f,
                                type, sample_c,
                                sigma_beta, a_c, b_c, a_lambda, b_lambda,
                                var_df, var_c, var_lambda,
                                p_c, p_prop, p_beta, p_df, p_lambda){

  sample_size <- burnin + lag*nsim

  n_cov <- ncol(X)

  beta_aux <- p_beta[1, ]
  c_aux <- p_c[1]
  df_aux <- p_df[1]
  lambda_aux <- p_lambda[1]

  ##-- For adapting MCMC
  it <- 1
  beta_adapt <- matrix(beta_aux, nrow = 100, ncol = n_cov, byrow = TRUE)
  c_adapt <- rep(c_aux, 100)
  df_adapt <- rep(df_aux, 100)
  lambda_adapt <- rep(lambda_aux, 100)

  sigma_beta_met <- sqrt(diag(solve(t(X)%*%X)))
  sigma_df_met <- sqrt(var_df)
  sigma_c_met <- sqrt(var_c)
  sigma_lambda_met <- sqrt(var_lambda)

  ##-- Stop points
  stop_pts <- seq(100, 1000, 100)
  pb <- progress::progress_bar$new(total = sample_size-1)

  ##-- Loop ----
  for(i in 2:sample_size){
    pb$tick()

    ##-- Regression coefficients
    if(i %in% stop_pts){
      sigma_beta_met <- sqrt(diag(stats::var(beta_adapt)))
      sigma_df_met <- stats::sd(log(df_adapt))
      sigma_c_met <- stats::sd(c_adapt)
      sigma_lambda_met <- stats::sd(lambda_adapt)
    }

    beta_at <- beta_aux
    for(j in 1:n_cov){
      post_beta_current <- beta_posterior(y = y, X = X,
                                          p_beta = beta_at, p_beta_element = beta_at[j], element = j,
                                          p_c = c_aux, p_df = df_aux,
                                          inv_link_f = inv_link_f,
                                          sigma_beta = sigma_beta,
                                          log = TRUE, method = "metropolis")

      ##-- + Proposal
      beta_prop <- stats::rnorm(n = 1, mean = beta_at[j], sd = sigma_beta_met[j])
      beta_at[j] <- beta_prop

      post_beta_sampled <- beta_posterior(y = y, X = X,
                                          p_beta = beta_at, p_beta_element = beta_at[j], element = j,
                                          p_c = c_aux, p_df = df_aux,
                                          inv_link_f = inv_link_f,
                                          sigma_beta = sigma_beta,
                                          log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_beta_sampled - post_beta_current)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        beta_aux[j] <- beta_prop
      } else{
        beta_at[j] <- beta_aux[j]
      }
    }

    ##-- c parameter
    if(sample_c){
      post_c_current <- c_posterior(y = y, X = X,
                                    p_beta = beta_aux, p_c = c_aux, p_df = df_aux,
                                    inv_link_f = inv_link_f,
                                    a_c = a_c, b_c = b_c,
                                    log = TRUE)

      ##-- + Proposal
      c_prop <- rtnorm(n = 1, mean = c_aux, sd = sigma_c_met, truncA = 0, truncB = 1)

      prop_1 <- dtnorm(x = c_aux, mean = c_aux, sd = sigma_c_met, truncA = 0, truncB = 1)
      prop_2 <- dtnorm(x = c_prop, mean = c_aux, sd = sigma_c_met, truncA = 0, truncB = 1)

      post_c_sampled <- c_posterior(y = y, X = X,
                                    p_beta = beta_aux, p_c = c_prop, p_df = df_aux,
                                    inv_link_f = inv_link_f,
                                    a_c = a_c, b_c = b_c,
                                    log = TRUE)

      ##-- + Metropolis step
      metropolis <- exp(post_c_sampled - post_c_current + prop_1 - prop_2)

      if(is.nan(metropolis)) metropolis <- 0
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        c_aux <- c_prop
      }
    }

    if(type == "robit"){
      ##-- lambda parameter
      post_lambda_current <- lambda_posterior(p_df = df_aux, p_lambda = lambda_aux,
                                              inv_link_f = inv_link_f,
                                              log = TRUE)

      ##-- + Proposal
      lambda_prop <- rtnorm(n = 1, mean = lambda_aux, sd = sigma_lambda_met, truncA = a_lambda, truncB = b_lambda)

      prop_1 <- dtnorm(x = lambda_aux, mean = lambda_aux, sd = sigma_c_met, truncA = a_lambda, truncB = b_lambda)
      prop_2 <- dtnorm(x = lambda_prop, mean = lambda_aux, sd = sigma_c_met, truncA = a_lambda, truncB = b_lambda)

      post_lambda_sampled <- lambda_posterior(p_df = df_aux, p_lambda = lambda_prop,
                                              inv_link_f = inv_link_f,
                                              log = TRUE)

      ##-- + Metropolis step
      metropolis <- exp(post_lambda_sampled - post_lambda_current + prop_1 - prop_2)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        lambda_aux <- lambda_prop
      }

      ##-- df parameter
      y_star <- log(df_aux)
      post_df_current <- df_posterior(y = y, X = X,
                                      p_beta = beta_aux, p_c = c_aux, p_df = y_star, p_lambda = lambda_aux,
                                      inv_link_f = inv_link_f,
                                      log = TRUE, method = "metropolis")

      ##-- + Proposal
      df_prop <- stats::rnorm(n = 1, mean = y_star, sd = sigma_df_met)

      post_df_sampled <- df_posterior(y = y, X = X,
                                      p_beta = beta_aux, p_c = c_aux, p_df = df_prop, p_lambda = lambda_aux,
                                      inv_link_f = inv_link_f,
                                      log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_df_sampled - post_df_current)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        df_aux <- exp(df_prop)
      }
    }

    if(i <= 1000){
      it <- ifelse(i %in% stop_pts, 1, it + 1)

      beta_adapt[it, ] <- beta_aux
      c_adapt[it] <- c_aux
      df_adapt[it] <- df_aux
      lambda_adapt[it] <- lambda_aux
    }

    ##-- Saving the iteration i
    if((i-burnin)%%lag == 0 & i > burnin){
      pos <- (i-burnin)/lag

      p_beta[pos, ] <- beta_aux
      p_c[pos] <- c_aux
      p_df[pos] <- df_aux
      p_lambda[pos] <- lambda_aux

      p_prop[pos, ] <- make_p(p_c = p_c[pos], X = X, p_beta = p_beta[pos, ], p_df = p_df[pos], inv_link_f = inv_link_f)
    }
  }

  ##-- Outputs ----
  return(list(p_c = p_c, p_beta = p_beta, p_prop = p_prop, p_df = p_df, p_lambda = p_lambda))
}
