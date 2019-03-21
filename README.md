
robit
=====

<!-- badges: start -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/robit)](https://cran.r-project.org/package=robit) -->
[![Travis build status](https://travis-ci.org/DouglasMesquita/robit.svg?branch=master)](https://travis-ci.org/DouglasMesquita/robit) [![Codecov test coverage](https://codecov.io/gh/DouglasMesquita/robit/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/robit?branch=master) <!-- badges: end -->

Overview
--------

A package to fit binary regression under several link functions and a possible cure fraction term in the form *p*<sub>*i*</sub> = *c* + (1 − *c*)*F*<sub>*v*</sub>(X*β*)

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

`library(robit)` will load the following function:

-   **mcmc\_bin**, for binary regression. The user can choose for MCMC based on ARMS or Metropolis Hastings.

See `?mcmc_bin` for a complete example of how to use this package and for examples.
