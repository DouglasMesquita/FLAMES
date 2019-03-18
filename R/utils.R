#' @title Inverse link function for several models (cumulative distributions)
#'
#' @param x Value to evaluate the cumulative distribution
#' @param type "logit", "probit", "cauchit", "tobit" or "cloglog"
#' @param df Degrees of freedom (if type = "tobit")
#'
#' @return Cumulative probability in x

inv_link <- function(x, type = "logit", df = 1){

  if(type == "logit"){
    cum_prob <- stats::plogis(q = x)
  }
  if(type == "probit"){
    cum_prob <- stats::pnorm(q = x)
  }
  if(type == "cauchit"){
    cum_prob <- stats::pcauchy(q = x)
  }
  if(type == "tobit"){
    cum_prob <- stats::pt(q = x, df = df)
  }
  if(type == "cloglog"){
    cum_prob <- 1-exp(-exp(x))
  }

  return(cum_prob)
}

#' @title Link function for several models
#'
#' @param x Probability to evaluate the quantile
#' @param type "logit", "probit", "cauchit", "tobit" or "cloglog"
#' @param df Degrees of freedom (if type = "tobit")
#'
#' @return Quantile of x

link <- function(x, type = "logistic", df = 1){

  if(type == "logistic"){
    quantile_x <- stats::qlogis(p = x)
  }
  if(type == "probit"){
    quantile_x <- stats::qnorm(p = x)
  }
  if(type == "cauchit"){
    quantile_x <- stats::qcauchy(p = x)
  }
  if(type == "tobit"){
    quantile_x <- stats::qt(p = x, df = df)
  }
  if(type == "cloglog"){
    quantile_x <- log(-log(1-x))
  }

  return(quantile_x)
}

#' @title Bernoulli likelihood
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param p_beta Regression coefficients
#' @param p_c c parameter
#' @param p_df Degrees of freedom
#' @param inv_link_f Inverse link function
#' @param log For log-scale values
#'
#' @return Likelihood value

bernoulli_like <- function(y, X,
                           p_beta, p_c, p_df,
                           inv_link_f,
                           log = TRUE){

  p <- make_p(p_c = p_c, X = X, p_beta = p_beta, p_df = p_df, inv_link_f = inv_link_f)

  if(sum(p == 1) != 0){
    p <- ifelse(p == 1, 0.999999, p)
  }
  if(sum(p == 0) != 0){
    p <- ifelse(p == 0, .Machine$double.eps, p)
  }

  like <- sum(y*log(p) + (1-y)*log(1-p))

  if(!log){
    like <- exp(like)
  }

  return(like)
}

#' @title Beta posterior distribution
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param p_beta Regression coefficients
#' @param p_beta_element Specific coefficient element
#' @param element Element you are sampling
#' @param p_c c parameter
#' @param p_df Degrees of freedom
#' @param inv_link_f Inverse link function
#' @param sigma_beta Variance of beta[element] prior
#' @param log For log-scale values
#' @param method "ARMS" or "metropolis"
#'
#' @return Likelihood value

beta_posterior <- function(y, X,
                           p_beta, p_beta_element, element, p_c, p_df,
                           inv_link_f,
                           sigma_beta = 100,
                           log = TRUE,
                           method = "ARMS"){

  if(method == "ARMS"){
    p_beta[element] <- log(p_beta_element/(1-p_beta_element))
  } else{
    p_beta[element] <- p_beta_element
  }

  like <- bernoulli_like(y = y, X = X, p_beta = p_beta, p_c = p_c, p_df = p_df, inv_link_f = inv_link_f, log = log)

  post_val <- like - (1/(2*sigma_beta^2))*p_beta[element]^2                                # prior
  if(method == "ARMS")  post_val <- post_val - log(p_beta_element) - log(1-p_beta_element) # jacobian

  if(!log){
    post_val <- exp(post_val)
  }

  return(post_val)
}

#' @title c posterior distribution
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param p_beta Regression coefficients
#' @param p_c c parameter
#' @param p_df Degrees of freedom
#' @param inv_link_f Inverse link function
#' @param a_c Shape1 for c prior
#' @param b_c Shape2 for c prior
#' @param log For log-scale values
#'
#' @return Likelihood value

c_posterior <- function(y, X,
                        p_beta, p_c, p_df,
                        inv_link_f,
                        a_c = 1, b_c = 1,
                        log = TRUE){

  like <- bernoulli_like(y = y, X = X, p_beta = p_beta, p_c = p_c, p_df = p_df, inv_link_f = inv_link_f, log = log)

  post_val <- like + (a_c-1)*log(p_c) + (b_c-1)*log(1-p_c) # prior

  if(!log){
    post_val <- exp(post_val)
  }

  return(post_val)
}

#' @title Degree of freedom posterior distribution
#'
#' @param y Bernoulli observed values
#' @param X Covariate matrix
#' @param p_beta Regression coefficients
#' @param p_c c parameter
#' @param p_df Degrees of freedom
#' @param p_lambda Lambda hyperparameter for p_df
#' @param inv_link_f Inverse link function
#' @param log For log-scale values
#' @param method "ARMS" or "metropolis"
#' @param const For numerical problems
#'
#' @return Likelihood value

df_posterior <- function(y, X,
                         p_beta,  p_c, p_df, p_lambda,
                         inv_link_f,
                         log = TRUE, method = "ARMS", const = 50){

  p_df_star <- p_df

  if(method == "ARMS"){
    p_df <- -const*log(1-p_df)
  } else{
    p_df <- exp(p_df)
  }

  like <- bernoulli_like(y = y, X = X, p_beta = p_beta, p_c = p_c, p_df = p_df, inv_link_f = inv_link_f, log = log)

  post_val <- like - (p_lambda*p_df)                                        # prior
  if(method == "ARMS") post_val <- post_val + log(const) - log(1-p_df_star) # jacobian
  if(method == "metropolis") post_val <- post_val + p_df_star               # jacobian

  if(!log){
    post_val <- exp(post_val)
  }

  return(post_val)
}

#' @title Lambda for degree of freedom posterior distribution
#'
#' @param p_df Degrees of freedom
#' @param p_lambda Lambda hyperparameter for p_df
#' @param inv_link_f Inverse link function
#' @param log For log-scale values
#'
#' @return Likelihood value

lambda_posterior <- function(p_df, p_lambda,
                             inv_link_f,
                             log = TRUE){

  post_val <- log(p_lambda) - p_lambda*p_df

  if(!log){
    post_val <- exp(post_val)
  }

  return(post_val)
}

#' @title Calculate p
#'
#' @param X Covariate matrix
#' @param p_beta Coefficients
#' @param p_c c parameter
#' @param p_df Degrees of freedom
#' @param inv_link_f Inverse link function
#'
#' @return Likelihood value

make_p <- function(X, p_beta, p_c, p_df, inv_link_f){

  out <- inv_link_f(x = X%*%p_beta, df = p_df)*(1-p_c) + p_c

  return(out)
}

rtnorm <- function(n, mean, sd, truncA = -Inf, truncB = Inf){

  u <- stats::runif(n = n, min = 0, max = 1)
  prob = u*(stats::pnorm(q = truncB, mean = mean, sd = sd)-
              stats::pnorm(q = truncA, mean = mean, sd = sd)) + stats::pnorm(q = truncA, mean = mean, sd = sd)
  quant <- stats::qnorm(p = prob, mean = mean, sd = sd)
  return(quant)
}
dtnorm <- function(x, mean, sd, truncA = -Inf, truncB = Inf){
  prob <- ifelse(x < truncA | x > truncB,
                 0,
                 stats::dnorm(x = x, mean = mean, sd = sd)/
                   (stats::pnorm(q = truncB, mean = mean, sd = sd) - stats::pnorm(q = truncA, mean = mean, sd = sd)))

  return(prob)
}
mean_sd_beta <- function(mean = NULL, sd = NULL, a = NULL, b = NULL){

  if(all(is.null(c(mean, var, a, b)))) stop("You must to define mean and sd or a and b.")

  if(all(is.null(c(mean, sd)))){
    var <- (a*b)/((a+b+1)*(a+b)^2)
    sd <- sqrt(var)
    mean <- a/(a+b)
  }

  if(all(is.null(c(a, b)))){
    var_lim <- mean*(1-mean)

    if(sd >= sqrt(var_lim)) stop(sprintf("sd must be smaller than %s.", round(sqrt(var_lim), 2)))

    nu <- var_lim/sd^2 - 1

    a <- mean*nu
    b <- (1-mean)*nu
  }

  return(list(a = a,
              b = b,
              mean = mean,
              sd = sd))
}
