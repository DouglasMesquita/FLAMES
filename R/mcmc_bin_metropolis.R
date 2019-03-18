#' @title MCMC for binary regression - Metropolis version
#'
#' @description Perform MCMC with Metropolis Hastings algorithm
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param sample_size Sample size required for MCMC
#' @param burnin Burn in period (to adapt metropolis)
#' @param inv_link_f Inverse link function
#' @param type "logit", "probit", "cauchit", "tobit" or "cloglog"
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
                                sample_size, burnin,
                                inv_link_f,
                                type, sample_c,
                                sigma_beta, a_c, b_c, a_lambda, b_lambda,
                                var_df, var_c, var_lambda,
                                p_c, p_prop, p_beta, p_df, p_lambda){

  sigma_beta_met <- sqrt(diag(solve(t(X)%*%X)))
  sigma_df_met <- sqrt(var_df)
  sigma_c_met <- sqrt(var_c)
  sigma_lambda_met <- sqrt(var_lambda)

  n_cov <- ncol(X)

  ##-- Stop points
  stop_pts <- seq(100, 1000, 100)

  pb <- progress::progress_bar$new(total = sample_size-1)
  ##-- Loop ----
  for(i in 2:sample_size){
    pb$tick()

    ##-- Regression coefficients
    aux_beta <- p_beta[(i-1),]

    if(i %in% stop_pts){
      sigma_beta_met <- sqrt(diag(stats::var(p_beta[(i-100):i, ])))
      sigma_df_met <- stats::sd(log(p_df[(i-100):i]))
      sigma_c_met <- stats::sd(p_c[(i-100):i])
      sigma_lambda_met <- stats::sd(p_lambda[(i-100):i])
    }

    for(j in 1:n_cov){
      post_beta_current <- beta_posterior(y = y, X = X,
                                          p_beta = aux_beta, p_beta_element = aux_beta[j], element = j,
                                          p_c = p_c[(i-1)], p_df = p_df[i-1],
                                          inv_link_f = inv_link_f,
                                          sigma_beta = sigma_beta,
                                          log = TRUE, method = "metropolis")

      ##-- + Proposal
      beta_prop <- stats::rnorm(n = 1, mean = p_beta[(i-1), j], sd = sigma_beta_met[j])
      aux_beta[j] <- beta_prop

      post_beta_sampled <- beta_posterior(y = y, X = X,
                                          p_beta = aux_beta, p_beta_element = aux_beta[j], element = j,
                                          p_c = p_c[(i-1)], p_df = p_df[i-1],
                                          inv_link_f = inv_link_f,
                                          sigma_beta = sigma_beta,
                                          log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_beta_sampled - post_beta_current)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        p_beta[i, j] <- beta_prop
      } else{
        p_beta[i, j] <- p_beta[(i-1), j]
        aux_beta[j] <- p_beta[(i-1), j]
      }
    }

    ##-- c parameter
    if(sample_c){
      post_c_current <- c_posterior(y = y, X = X,
                                    p_beta = p_beta[i,], p_c = p_c[(i-1)], p_df = p_df[(i-1)],
                                    inv_link_f = inv_link_f,
                                    a_c = a_c, b_c = b_c,
                                    log = TRUE)

      ##-- + Proposal
      c_prop <- rtnorm(n = 1, mean = p_c[(i-1)], sd = sigma_c_met, truncA = 0, truncB = 1)

      prop_1 <- dtnorm(x = p_c[(i-1)], mean = p_c[(i-1)], sd = sigma_c_met, truncA = 0, truncB = 1)
      prop_2 <- dtnorm(x = c_prop, mean = p_c[(i-1)], sd = sigma_c_met, truncA = 0, truncB = 1)

      post_c_sampled <- c_posterior(y = y, X = X,
                                    p_beta = p_beta[i,], p_c = c_prop, p_df = p_df[(i-1)],
                                    inv_link_f = inv_link_f,
                                    a_c = a_c, b_c = b_c,
                                    log = TRUE)

      ##-- + Metropolis step
      metropolis <- exp(post_c_sampled - post_c_current + prop_1 - prop_2)

      if(is.nan(metropolis)) metropolis <- 0
      exp_val <- stats::rexp(n = 1, rate = 1)
      if(metropolis > exp_val){
        p_c[i] <- c_prop
      } else{
        p_c[i] <- p_c[(i-1)]
      }
    }

    if(type == "tobit"){
      ##-- lambda parameter
      post_lambda_current <- lambda_posterior(p_df = p_df[i-1], p_lambda = p_lambda[i-1],
                                              inv_link_f = inv_link_f,
                                              log = TRUE)

      ##-- + Proposal
      lambda_prop <- rtnorm(n = 1, mean = p_lambda[(i-1)], sd = sigma_lambda_met, truncA = 0, truncB = 1)

      prop_1 <- dtnorm(x = p_lambda[(i-1)], mean = p_lambda[(i-1)], sd = sigma_c_met, truncA = 0, truncB = 1)
      prop_2 <- dtnorm(x = lambda_prop, mean = p_lambda[(i-1)], sd = sigma_c_met, truncA = 0, truncB = 1)

      post_lambda_sampled <- lambda_posterior(p_df = p_df[i-1], p_lambda = lambda_prop,
                                              inv_link_f = inv_link_f,
                                              log = TRUE)

      ##-- + Metropolis step
      metropolis <- exp(post_lambda_sampled - post_lambda_current + prop_1 - prop_2)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        p_lambda[i] <- lambda_prop
      } else{
        p_lambda[i] <- p_lambda[(i-1)]
      }

      ##-- df parameter
      y_star <- log(p_df[i-1])
      post_df_current <- df_posterior(y = y, X = X,
                                      p_beta = p_beta[i, ], p_c = p_c[i], p_df = y_star, p_lambda = p_lambda[i],
                                      inv_link_f = inv_link_f,
                                      log = TRUE, method = "metropolis")

      ##-- + Proposal
      df_prop <- stats::rnorm(n = 1, mean = y_star, sd = sigma_df_met)

      post_df_sampled <- df_posterior(y = y, X = X,
                                      p_beta = p_beta[i, ], p_c = p_c[i], p_df = df_prop, p_lambda = p_lambda[i],
                                      inv_link_f = inv_link_f,
                                      log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_df_sampled - post_df_current)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        p_df[i] <- exp(df_prop)
      } else{
        p_df[i] <- p_df[(i-1)]
      }
    }

    p_prop[i,] <- make_p(p_c = p_c[i], X = X, p_beta = p_beta[i,], p_df = p_df[i], inv_link_f = inv_link_f)
  }

  ##-- Outputs ----
  return(list(p_c = p_c, p_beta = p_beta, p_prop = p_prop, p_df = p_df, p_lambda = p_lambda))
}
