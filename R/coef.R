#' Coefficients for FLAMES class
#'
#' @param object FLAMES object to restore the coefficients
#' @param ... Other parameters for coef
#'
#' @method coef FLAMES
#'
#' @export

coef.FLAMES <- function(object, ...){
  coefficients <- colMeans(object$beta)
  return(coefficients)
}
