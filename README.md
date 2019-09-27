
# FLAMES

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

[![Travis build
status](https://travis-ci.org/DouglasMesquita/FLAMES.svg?branch=master)](https://travis-ci.org/DouglasMesquita/FLAMES)
<!-- [![Codecov test coverage](https://codecov.io/gh/DouglasMesquita/FLAMES/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/FLAMES?branch=master) -->
<!-- badges: end -->

## Overview

FLAMES: Flexible Link function with AssyMptotES is a package able to fit
binary regression under several link functions and two possible
assymptotes parameters, c and d, in the form

<center>

<img src="https://raw.githubusercontent.com/DouglasMesquita/FLAMES/master/man/figures/gif.gif" alt="formula">

</center>

## Installation

``` r
# Install from CRAN (when available)
install.packages("FLAMES")
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("DouglasMesquita/FLAMES")
```

## Usage

``` r
library(FLAMES)
```

`library(FLAMES)` will load **mcmc\_bin**, a function to fit binary
regression models. The user can choose for MCMC based on ARMS or
Metropolis Hastings algorithms. Also it is possible to fit **logit**,
**probit**, **cauchit**, **cloglog**, **loglog** and **robit**
regressions.

``` r
##-- Seed
set.seed(123456)

##-- Data
n <- 1000
n_cov <- 1

##-- Covariates
X <- matrix(rnorm(n*n_cov), ncol = n_cov)

##-- Coefficients
betas <- c(0, -1)
XBeta <- cbind(1, X)%*%betas

##-- c parameter
c1 <- 0.20
d1 <- 0.95

type_data = "cloglog"

##-- p and y
p <- FLAMES:::inv_link(x = XBeta, type = type_data, df = df)*(d1-c1) + c1
y <- rbinom(n = n, size = 1, prob = p)

bd <- data.frame(y = y, X)

##-- MCMC
nsim <- 2000
burnin <- 10000
lag <- 10

f <- y ~ X
type <- type_data

##-- ARMS ~ 4 minutes (soon in c++)
out_arms <- mcmc_bin(data = bd, formula = f,
                     nsim = nsim, burnin = burnin, lag = lag,
                     type = type, sample_c = TRUE, sample_d = TRUE,
                     method = "ARMS")

##-- ARMS ~ 55 seconds (soon in c++)
out_met <- mcmc_bin(data = bd, formula = f,
                    nsim = nsim, burnin = burnin, lag = lag,
                    type = type, sample_c = TRUE, sample_d = TRUE,
                    method = "metropolis")

##-- GLM
out_glm <- glm(formula = f, data = bd, family = "binomial")
```

``` r
cbind(glm = coef(out_glm),
      arms = coef(out_arms),
      metropolis = coef(out_met),
      real = betas)
#>                   glm        arms  metropolis real
#> (Intercept)  1.005911  0.03550848  0.03918668    0
#> X           -1.071551 -1.58034100 -1.54308670   -1
```

``` r
par(mfrow = c(2, 4), mar = c(3, 4, 1, 1))
plot(out_arms, ask = F)
plot(out_met, ask = F)
```

<center>

<img src="https://raw.githubusercontent.com/DouglasMesquita/FLAMES/master/man/figures/chains.png" alt="chains" style="width:80%">

</center>

``` r
summary(out_arms)
#> $Call
#> mcmc_bin(data = bd, formula = f, nsim = nsim, burnin = burnin, 
#>     lag = lag, type = type, sample_c = TRUE, sample_d = TRUE, 
#>     method = "ARMS")
#> 
#> $Coeficients
#>                    mean std_error   lower_95   upper_95
#> (Intercept)  0.03550848 0.1914609 -0.3369406  0.3833066
#> X           -1.58034100 0.5113920 -2.5337717 -0.7536014
#> 
#> $`Other parameters`
#>                  mean  std_error  lower_95  upper_95
#> c parameter 0.3120890 0.08972683 0.1215077 0.4628683
#> d parameter 0.9322188 0.02085929 0.8912744 0.9722881
#> 
#> $`Fit measures`
#>        DIC  -2*LPML     WAIC
#> 1 1041.917 1042.101 1041.917
```

``` r
summary(out_met)
#> $Call
#> mcmc_bin(data = bd, formula = f, nsim = nsim, burnin = burnin, 
#>     lag = lag, type = type, sample_c = TRUE, sample_d = TRUE, 
#>     method = "metropolis")
#> 
#> $Coeficients
#>                    mean std_error   lower_95   upper_95
#> (Intercept)  0.03918668 0.1991352 -0.3443703  0.4388762
#> X           -1.54308670 0.5257698 -2.5815248 -0.7613898
#> 
#> $`Other parameters`
#>                  mean  std_error  lower_95  upper_95
#> c parameter 0.3024514 0.09427923 0.1067364 0.4676210
#> d parameter 0.9328742 0.02038635 0.8929056 0.9712623
#> 
#> $`Fit measures`
#>        DIC -2*LPML     WAIC
#> 1 1047.701 1047.88 1047.701
```
