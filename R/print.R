#' Print for robit class
#'
#' @param x robit object to print
#' @param ... Other parameters for print
#'
#' @method print robit
#'
#' @export

print.robit <- function(x, ...){
  cat('Title: \n')
  cat('\n')
  cat('Coefficients:')
  cat('\n')
  print(coef.robit(x))
}
