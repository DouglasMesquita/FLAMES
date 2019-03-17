#' Plot for tbreg class
#'
#' @param x tbreg object to plot
#' @param bty Box style
#' @param ask To ask for next plot
#' @param ... Other parameters for plot
#'
#' @method plot tbreg
#'
#' @import coda
#'
#' @export

plot.tbreg <- function(x, bty = "l", ask = T, ...){

  data <- data.frame(x$beta, c = x$c, df = x$df)

  p <- ncol(x$beta)
  n <- nrow(x$beta)

  ##-- Coefficients ----
  for(i in 1:p){
    xlab <- 'Iteration'
    ylab <- colnames(x$beta)[i]
    graphics::plot(x = 1:n, y = x$beta[, i],
                   col = 'grey20', type = 'l', cex = 1.5,
                   xlab = xlab, ylab = ylab,
                   bty = bty, ...)
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
  }

  ##-- Other parameters ----
  ##-- + c parameter ----
  graphics::plot(x = 1:n, y = data$c,
                 col = 'grey20', type = 'l', cex = 1.5,
                 xlab = 'Iteration', ylab = 'c parameter',
                 bty = bty, ...)
  if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))

  ##-- + Temporal parameter ----
  graphics::plot(x = 1:n, y = data$df,
                 col = 'grey20', type = 'l', cex = 1.5,
                 xlab = 'Iteration', ylab = 'Degrees of freedom',
                 bty = bty, ...)
}
