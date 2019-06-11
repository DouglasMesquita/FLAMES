#' Print for FLAMES class
#'
#' @param x FLAMES object to print
#' @param ... Other parameters for print
#'
#' @method print FLAMES
#'
#' @export

print.FLAMES <- function(x, ...){
  cat('Title: \n')
  cat('\n')
  cat('Coefficients:')
  cat('\n')
  print(coef.FLAMES(x))
}
