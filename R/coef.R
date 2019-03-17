#' Coefficients for tbreg class
#'
#' @param object tbreg object to restore the coefficients
#' @param ... Other parameters for coef
#'
#' @method coef tbreg
#'
#' @export

coef.tbreg <- function(object, ...){
  coefficients <- colMeans(object$beta)
  return(coefficients)
}
