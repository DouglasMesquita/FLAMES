
robit
=====

<!-- badges: start -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/robit)](https://cran.r-project.org/package=robit) -->
[![Travis build status](https://travis-ci.org/DouglasMesquita/robit.svg?branch=master)](https://travis-ci.org/DouglasMesquita/robit) <!-- [![Codecov test coverage](https://codecov.io/gh/DouglasMesquita/robit/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/robit?branch=master) --> <!-- badges: end -->

Overview
--------

A package to fit binary regression under several link functions and two possible fraction parameters, c and d, in the form

*p*<sub>*i*</sub> = *c* + (*d* − *c*)*F*<sub>*v*</sub>(X*β*)

Installation
------------

``` r
# Install from CRAN (when available)
install.packages("robit")
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("DouglasMesquita/robit")
```

Usage
-----

`library(robit)` will load **mcmc\_bin**, a function to fit binary regression models. The user can choose for MCMC based on ARMS or Metropolis Hastings algorithms. Also it is possible to fit **logit**, **probit**, **cauchit**, **cloglog** and **robit** regressions.

See `?mcmc_bin` for examples.
