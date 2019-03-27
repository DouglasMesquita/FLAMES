#' @title Inverse link function for several models (cumulative distributions)
#'
#' @param x Value to evaluate the cumulative distribution
#' @param type "logit", "probit", "cauchit", "robit" or "cloglog"
#' @param df Degrees of freedom (if type = "robit")
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
  if(type == "robit"){
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
#' @param type "logit", "probit", "cauchit", "robit" or "cloglog"
#' @param df Degrees of freedom (if type = "robit")
#'
#' @return Quantile of x

link <- function(x, type = "logit", df = 1){

  if(type == "logit"){
    quantile_x <- stats::qlogis(p = x)
  }
  if(type == "probit"){
    quantile_x <- stats::qnorm(p = x)
  }
  if(type == "cauchit"){
    quantile_x <- stats::qcauchy(p = x)
  }
  if(type == "robit"){
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
#' @param method "ARMS" or "metropolis"
#'
#' @return Likelihood value

c_posterior <- function(y, X,
                        p_beta, p_c, p_df,
                        inv_link_f,
                        a_c = 1, b_c = 1,
                        log = TRUE, method = "ARMS"){

  p_c_star <- p_c

  if(method == "metropolis") p_c <- exp(p_c)/(1+exp(p_c))

  like <- bernoulli_like(y = y, X = X, p_beta = p_beta, p_c = p_c, p_df = p_df, inv_link_f = inv_link_f, log = log)

  post_val <- like + (a_c-1)*log(p_c) + (b_c-1)*log(1-p_c)                            # prior
  if(method == "metropolis") post_val <- post_val + p_c_star - 2*log(1+exp(p_c_star)) # Jacobian

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
                         log = TRUE, method = "ARMS", const){

  p_df_star <- p_df

  if(method == "ARMS"){
    p_df <- -const*log(1-p_df)
  } else{
    p_df <- exp(p_df)/const
  }

  like <- bernoulli_like(y = y, X = X, p_beta = p_beta, p_c = p_c, p_df = p_df, inv_link_f = inv_link_f, log = log)

  post_val <- like - (p_lambda*p_df)                                              # prior
  if(method == "ARMS") post_val <- post_val + log(const) - log(1-p_df_star)       # jacobian
  if(method == "metropolis") post_val <- post_val + p_df_star - log(const)        # jacobian

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
#' @param method "ARMS" or "metropolis"
#' @param a_lambda Inferior limit for lambda
#' @param b_lambda Superior limit for lambda
#'
#' @return Likelihood value

lambda_posterior <- function(p_df, p_lambda,
                             inv_link_f,
                             log = TRUE, method = "ARMS",
                             a_lambda = NULL, b_lambda = NULL){
  p_lambda_star <- p_lambda

  if(method == "metropolis") p_lambda <- exp(p_lambda)*(b_lambda - a_lambda)/(1+exp(p_lambda)) + a_lambda

  post_val <- log(p_lambda) - p_lambda*p_df
  if(method == "metropolis") post_val <- post_val + p_lambda_star - 2*log(1+exp(p_lambda_star)) + log(b_lambda - a_lambda) # Jacobian

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
  out <- as.numeric(out)

  return(out)
}

#' @title Conditional Predictive Ordinate for Poisson data
#'
#' @param y Vector of response variable
#' @param p Matrix of probabilities for each individual y in each MCMC iteration
#'
#' @return cpo Conditional Predictive Ordinate
#'

CPO_bernoulli <- function(y, p){
  niter <- length(p)
  px <- dbinom(x = y, size = 1, prob = p)
  px <- ifelse(px < .Machine$double.eps, .Machine$double.eps, px)
  cpo <- niter/sum(1/px, na.rm = T)
  return(cpo)
}

#' @title Log Pseudo Marginal Likelihood
#'
#' @param y Vector of response variable
#' @param p Matrix of probabilities for each individual y in each MCMC iteration
#'
#' @return LPML -2 * Log Pseudo Marginal Likelihood measure

LPML <- function(y, p){
  aux <- data.frame(y, t(p))
  CPO <- apply(X = aux, MARGIN = 1,
               FUN = function(x) CPO_bernoulli(y = x[1], p = x[-1]))
  lpml <- sum(log(CPO), na.rm = T)
  return(-2*lpml)
}

#' @title Deviance Information Criterion
#'
#' @param y Vector of response variable
#' @param p Matrix of probabilities for each individual y in each MCMC iteration
#'
#' @return DIC Deviance Information Criterion measure

DIC <- function(y, p){
  niter <- nrow(p)
  thetaBayes <- colMeans(p)
  log_px <- dbinom(x = y, size = 1, prob = thetaBayes, log = TRUE)
  log_px <- sum(log_px, na.rm = T)
  p <- cbind(y, t(p))
  esp_log_px <- apply(p, MARGIN = 1,
                      FUN = function(x) dbinom(x = x[1], size = 1, prob = x[-1], log = TRUE))
  esp_log_px <- sum(esp_log_px, na.rm = T)/niter

  pDIC <- 2*(log_px - esp_log_px)
  dic <- -2*log_px + 2*pDIC
  return(dic)
}

#' @title Widely Applicable Information Criterion
#'
#' @param y Vector of response variable
#' @param p Matrix of probabilities for each individual y in each MCMC iteration
#'
#' @return WAIC_1 Widely Applicable Information Criterion measure

WAIC <- function(y, p){
  aux <- data.frame(y, t(p))
  px <- apply(aux, MARGIN = 1,
              FUN = function(x) dbinom(x = x[1], size = 1, prob = x[-1]))
  px <- ifelse(px < .Machine$double.eps, .Machine$double.eps, px)
  log_px <- log(px)
  lppd <- sum(log(colMeans(px, na.rm = T)))
  mean_log_px <- colMeans(log_px, na.rm = T)
  log_mean_px <- log(colMeans(px, na.rm = T))
  pWAIC <- 2*sum(log_mean_px - mean_log_px)
  waic <- -2*(lppd - pWAIC)
  return(waic)
}

#' @title LPML, DIC and WAIC measures
#'
#' @param y Vector of response variable
#' @param p Matrix of probabilities for each individual y in each MCMC iteration
#'
#' @return measures A data.frame with LMPL, DIC and WAIC measures

fit_measures <- function(y, p, nrep = NULL, nsamp = 100){

  if(!is.null(nrep)){
    n <- length(y)
    lpml <- dic <- waic <- vector(mode = "list", length = nrep)

    for(i in 1:nrep){
      samp_pos <- sample(x = 1:n, size = nsamp, replace = TRUE)

      lpml[[i]] <- LPML(y[samp_pos], p[, samp_pos])
      dic[[i]] <- DIC(y[samp_pos], p[, samp_pos])
      waic[[i]] <- WAIC(y[samp_pos], p[, samp_pos])
    }

    lpml <- mean(unlist(lpml))*n/nsamp
    dic <- mean(unlist(dic))*n/nsamp
    waic <- mean(unlist(waic))*n/nsamp
  } else{
    lpml <- LPML(y, p)
    dic <- DIC(y, p)
    waic <- WAIC(y, p)
  }

  invisible(gc(reset = TRUE, verbose = FALSE, full = TRUE))

  measures <- data.frame('DIC' = dic, '-2*LPML' = lpml, 'WAIC' = waic, check.names = F)

  return(measures = measures)
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
mean_sd_beta <- function(mean = NULL, sd = NULL, a = NULL, b = NULL, show_warnings = FALSE){

  if(all(is.null(c(mean, var, a, b)))) stop("You must to define mean and sd or a and b.")

  if(all(is.null(c(mean, sd)))){
    var <- (a*b)/((a+b+1)*(a+b)^2)
    sd <- sqrt(var)
    mean <- a/(a+b)
  }

  if(all(is.null(c(a, b)))){
    var_lim <- mean*(1-mean)

    if(sd >= sqrt(var_lim)){
      sd <- sqrt(var_lim)*0.99
      if(show_warnings) warning(sprintf("sd must be smaller than %s. We set it as %s", round(sqrt(var_lim), 2), round(sqrt(var_lim), 2)))
    }

    nu <- var_lim/sd^2 - 1

    a <- mean*nu
    b <- (1-mean)*nu
  }

  return(list(a = a,
              b = b,
              mean = mean,
              sd = sd))
}
