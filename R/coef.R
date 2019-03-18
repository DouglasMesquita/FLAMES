#' Coefficients for robit class
#'
#' @param object robit object to restore the coefficients
#' @param ... Other parameters for coef
#'
#' @method coef robit
#'
#' @export

coef.robit <- function(object, ...){
  coefficients <- colMeans(object$beta)
  return(coefficients)
}
