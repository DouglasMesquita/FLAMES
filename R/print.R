#' Print for tbreg class
#'
#' @param x tbreg object to print
#' @param ... Other parameters for print
#'
#' @method print tbreg
#'
#' @export

print.tbreg <- function(x, ...){
  cat('Title: \n')
  cat('\n')
  cat('Coefficients:')
  cat('\n')
  print(coef.tbreg(x))
}
