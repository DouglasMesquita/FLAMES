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
#' @param type "logit", "probit", "cauchit", "robit", "cloglog" or "loglog"
#' @param sample_c Should c be sampled?
#' @param sample_d Should c be sampled?
#' @param sigma_beta Variance of beta prior
#' @param a_c Shape1 for c prior
#' @param b_c Shape2 for c prior
#' @param a_d Shape1 for d prior
#' @param b_d Shape2 for d prior
#' @param a_lambda Inferior limit for lambda
#' @param b_lambda Superior limit for lambda
#' @param var_df Variance to sample 1-exp(-df/const)
#' @param var_c Variance to sample from c (if sample_d = TRUE) otherwise log(c/(1-c))
#' @param var_d Variance to sample d (if sample_c = TRUE) otherwise log(c/(1-c))
#' @param var_lambda Variance to sample lambda
#' @param p_c To restore c
#' @param p_d To restore d
#' @param p_prop To restore p
#' @param p_beta To restore beta
#' @param p_df To restore df
#' @param p_lambda To restore lambda
#' @param const A constant to help on sampling degrees of freedom \eqn{\tilde{df} = df/c}
#' @param const_beta A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_c A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_d A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_df A constant to tunning the acceptance rate (default = 2.38^2)
#' @param const_lambda A constant to tunning the acceptance rate (default = 2.38^2)
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
                                type, sample_c, sample_d,
                                sigma_beta,
                                a_c, b_c, a_d, b_d, a_lambda, b_lambda,
                                var_df, var_c, var_d, var_lambda,
                                p_c, p_d, p_prop, p_beta, p_df, p_lambda,
                                const, const_beta, const_c, const_d, const_df, const_lambda){

  sample_size <- burnin + lag*nsim

  n_cov <- ncol(X)

  beta_aux <- p_beta[[1]]
  c_aux <- ifelse(!is.null(p_c), p_c[[1]], 0)
  d_aux <- ifelse(!is.null(p_d), p_d[[1]], 1)
  if(!is.null(p_df)) df_aux <- p_df[[1]]
  if(!is.null(p_lambda)) lambda_aux <- p_lambda[[1]]

  ##-- For adapting MCMC
  it <- 1
  size_adapt <- 500
  beta_adapt <- matrix(beta_aux, nrow = size_adapt, ncol = n_cov, byrow = TRUE)
  if(!is.null(p_c)) c_adapt <- rep(c_aux, size_adapt)
  if(!is.null(p_d)) d_adapt <- rep(d_aux, size_adapt)
  if(!is.null(p_df)) df_adapt <- rep(df_aux, size_adapt)
  if(!is.null(p_lambda)) lambda_adapt <- rep(lambda_aux, size_adapt)

  sigma_beta_met <- sqrt(diag(solve(t(X)%*%X)))
  if(!is.null(p_c)) sigma_c_met <- sqrt(var_c)
  if(!is.null(p_d)) sigma_d_met <- sqrt(var_d)
  if(!is.null(p_df)) sigma_df_met <- sqrt(var_df)
  if(!is.null(p_lambda)) sigma_lambda_met <- sqrt(var_lambda)

  ##-- Stop points
  stop_pts <- seq(100, sample_size, size_adapt)
  pb <- progress::progress_bar$new(total = sample_size-1)

  ##-- Loop ----
  for(i in 2:sample_size){
    pb$tick()

    ##-- Regression coefficients
    if(i %in% stop_pts){
      sigma_beta_met <- sqrt(diag(stats::var(beta_adapt)))*const_beta + 1e-05
      if(!is.null(p_c)) sigma_c_met <- stats::sd(c_adapt)*const_c + 1e-05
      if(!is.null(p_d)) sigma_d_met <- stats::sd(d_adapt)*const_d + 1e-05
      if(!is.null(p_df)) sigma_df_met <- stats::sd(df_adapt)*const_df + 1e-05
      if(!is.null(p_lambda)) sigma_lambda_met <- stats::sd(lambda_adapt)*const_lambda + 1e-05
    }

    beta_at <- beta_aux
    for(j in 1:n_cov){
      post_beta_current <- beta_fullcond(y = y, X = X,
                                         p_beta = beta_at, p_beta_element = beta_at[j], element = j,
                                         p_c = c_aux, p_d = d_aux, p_df = df_aux,
                                         inv_link_f = inv_link_f,
                                         sigma_beta = sigma_beta,
                                         log = TRUE, method = "metropolis")

      ##-- + Proposal
      beta_prop <- stats::rnorm(n = 1, mean = beta_at[j], sd = sigma_beta_met[j])
      beta_at[j] <- beta_prop

      post_beta_sampled <- beta_fullcond(y = y, X = X,
                                         p_beta = beta_at, p_beta_element = beta_at[j], element = j,
                                         p_c = c_aux, p_d = d_aux, p_df = df_aux,
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
    if(sample_c & !sample_d){
      c_star <- log(c_aux/(1-c_aux))

      post_c_current <- c_fullcond(y = y, X = X,
                                   p_beta = beta_aux, p_c = c_star, p_d = d_aux, p_df = df_aux,
                                   inv_link_f = inv_link_f,
                                   a_c = a_c, b_c = b_c,
                                   log = TRUE, method = "metropolis")

      ##-- + Proposal
      c_prop <- rnorm(n = 1, mean = c_star, sd = sigma_c_met)

      post_c_sampled <- c_fullcond(y = y, X = X,
                                   p_beta = beta_aux, p_c = c_prop, p_d = d_aux, p_df = df_aux,
                                   inv_link_f = inv_link_f,
                                   a_c = a_c, b_c = b_c,
                                   log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_c_sampled - post_c_current)

      if(is.nan(metropolis)) metropolis <- 0
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        c_aux <- exp(c_prop)/(1+exp(c_prop))
      }
    }

    ##-- d parameter
    if(sample_d & !sample_c){
      d_star <- log(d_aux/(1-d_aux))

      post_d_current <- d_fullcond(y = y, X = X,
                                   p_beta = beta_aux, p_c = c_aux, p_d = d_star, p_df = df_aux,
                                   inv_link_f = inv_link_f,
                                   a_d = a_d, b_d = b_d,
                                   log = TRUE, method = "metropolis")

      ##-- + Proposal
      d_prop <- rnorm(n = 1, mean = d_star, sd = sigma_d_met)

      post_d_sampled <- d_fullcond(y = y, X = X,
                                   p_beta = beta_aux, p_c = c_aux, p_d = d_prop, p_df = df_aux,
                                   inv_link_f = inv_link_f,
                                   a_d = a_d, b_d = b_d,
                                   log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_d_sampled - post_d_current)

      if(is.nan(metropolis)) metropolis <- 0
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        d_aux <- exp(d_prop)/(1+exp(d_prop))
      }
    }

    ##-- c and d parameters
    if(sample_d & sample_c){

      post_cd_current <- cd_fullcond(y = y, X = X,
                                     p_beta = beta_aux, p_c = c_aux, p_d = d_aux, p_df = df_aux,
                                     inv_link_f = inv_link_f,
                                     a_c = a_c, b_c = b_c,
                                     a_d = a_d, b_d = b_d,
                                     log = TRUE, method = "metropolis")

      ##-- + Proposal
      lim_cd <- c(0.0001, 0.9999)
      c_prop <- rtnorm(n = 1, mean = c_aux, sd = sigma_c_met, truncA = lim_cd[1], truncB = lim_cd[2])
      d_prop <- rtnorm(n = 1, mean = d_aux, sd = sigma_d_met, truncA = lim_cd[1], truncB = lim_cd[2])

      dens1c <- dtnorm(x = c_aux, mean = c_prop, sd = sigma_c_met, truncA = 0, truncB = 1)
      dens2c <- dtnorm(x = c_prop, mean = c_aux, sd = sigma_c_met, truncA = 0, truncB = 1)
      dens1d <- dtnorm(x = d_aux, mean = d_prop, sd = sigma_d_met, truncA = c_prop, truncB = 1)
      dens2d <- dtnorm(x = d_prop, mean = d_aux, sd = sigma_d_met, truncA = c_prop, truncB = 1)

      post_cd_sampled <- cd_fullcond(y = y, X = X,
                                     p_beta = beta_aux, p_c = c_prop, p_d = d_prop, p_df = df_aux,
                                     inv_link_f = inv_link_f,
                                     a_c = a_c, b_c = b_c,
                                     a_d = a_d, b_d = b_d,
                                     log = TRUE, method = "metropolis")

      ##-- + Metropolis step
      metropolis <- exp(post_cd_sampled - post_cd_current + dens1c + dens1d - dens2c - dens2d)

      if(is.nan(metropolis)) metropolis <- 0
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        c_aux <- c_prop
        d_aux <- d_prop
      }
    }

    if(type == "robit"){
      ##-- lambda parameter
      lambda_01 <- (lambda_aux-a_lambda)/(b_lambda-a_lambda)
      lambda_star <- log(lambda_01/(1-lambda_01))
      post_lambda_current <- lambda_fullcond(p_df = df_aux, p_lambda = lambda_star,
                                             inv_link_f = inv_link_f,
                                             log = TRUE, method = "metropolis",
                                             a_lambda = a_lambda, b_lambda = b_lambda)

      ##-- + Proposal
      lambda_prop <- rnorm(n = 1, mean = lambda_star, sd = sigma_lambda_met)

      post_lambda_sampled <- lambda_fullcond(p_df = df_aux, p_lambda = lambda_prop,
                                             inv_link_f = inv_link_f,
                                             log = TRUE, method = "metropolis",
                                             a_lambda = a_lambda, b_lambda = b_lambda)

      ##-- + Metropolis step
      metropolis <- exp(post_lambda_sampled - post_lambda_current)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        lambda_aux <- exp(lambda_prop)*(b_lambda-a_lambda)/(1+exp(lambda_prop)) + a_lambda
      }

      ##-- df parameter
      df_star <- 1-exp(-df_aux/const)
      post_df_current <- df_fullcond(y = y, X = X,
                                     p_beta = beta_aux, p_c = c_aux, p_d = d_aux, p_df = df_star, p_lambda = lambda_aux,
                                     inv_link_f = inv_link_f,
                                     log = TRUE, method = "metropolis", const = const)

      ##-- + Proposal
      df_prop <- rtnorm(n = 1, mean = df_star, sd = sigma_df_met, truncA = 0.0001, truncB = 0.9999)

      dens1 <- dtnorm(x = df_star, mean = df_prop, sd = sigma_df_met, truncA = 0, truncB = 1)
      dens2 <- dtnorm(x = df_prop, mean = df_star, sd = sigma_df_met, truncA = 0, truncB = 1)

      post_df_sampled <- df_fullcond(y = y, X = X,
                                     p_beta = beta_aux, p_c = c_aux, p_d = d_aux, p_df = df_prop, p_lambda = lambda_aux,
                                     inv_link_f = inv_link_f,
                                     log = TRUE, method = "metropolis", const = const)

      ##-- + Metropolis step
      metropolis <- exp(post_df_sampled - post_df_current + dens1 - dens2)
      exp_val <- stats::rexp(n = 1, rate = 1)

      if(metropolis > exp_val){
        df_aux <- -const*log(1-df_prop)
      }
    }

    if(i <= max(stop_pts)){
      it <- ifelse(i %in% stop_pts, 1, it + 1)

      beta_adapt[it, ] <- beta_aux
      if(!is.null(p_c) & !is.null(p_d)){
        c_adapt[it] <- c_aux
        d_adapt[it] <- d_aux
      }
      if(!is.null(p_c) & is.null(p_d)) c_adapt[it] <- log(c_aux/(1-c_aux))
      if(!is.null(p_d) & is.null(p_c)) d_adapt[it] <- log(d_aux/(1-d_aux))
      if(!is.null(p_df)) df_adapt[it] <- 1-exp(-df_aux/const)
      if(!is.null(p_lambda)){
        lambda_01 <- (lambda_aux-a_lambda)/(b_lambda-a_lambda)
        lambda_adapt[it] <- log(lambda_01/(1-lambda_01))
      }
    }

    ##-- Saving the iteration i
    if((i-burnin)%%lag == 0 & i > burnin){
      pos <- (i-burnin)/lag

      p_beta[[pos]] <- beta_aux
      if(!is.null(p_c)) p_c[[pos]] <- c_aux
      if(!is.null(p_d)) p_d[[pos]] <- d_aux
      if(!is.null(p_df)) p_df[[pos]] <- df_aux
      if(!is.null(p_lambda)) p_lambda[[pos]] <- lambda_aux

      p_prop[[pos]] <- make_p(p_c = c_aux, p_d = d_aux, X = X, p_beta = beta_aux, p_df = df_aux, inv_link_f = inv_link_f)
    }
  }

  ##-- Outputs ----
  return(list(p_c = p_c, p_d = p_d, p_beta = p_beta, p_prop = p_prop, p_df = p_df, p_lambda = p_lambda))
}
