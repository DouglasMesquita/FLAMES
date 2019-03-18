#' Coefficients for crobit class
#'
#' @param object crobit object to restore the coefficients
#' @param ... Other parameters for coef
#'
#' @method coef crobit
#'
#' @export

coef.crobit <- function(object, ...){
  coefficients <- colMeans(object$beta)
  return(coefficients)
}
