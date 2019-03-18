#' Print for crobit class
#'
#' @param x crobit object to print
#' @param ... Other parameters for print
#'
#' @method print crobit
#'
#' @export

print.crobit <- function(x, ...){
  cat('Title: \n')
  cat('\n')
  cat('Coefficients:')
  cat('\n')
  print(coef.crobit(x))
}
