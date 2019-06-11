#' Plot for FLAMES class
#'
#' @param x FLAMES object to plot
#' @param bty Box style
#' @param ask To ask for next plot
#' @param ... Other parameters for plot
#'
#' @method plot FLAMES
#'
#' @import coda
#'
#' @export

plot.FLAMES <- function(x, bty = "l", ask = T, ...){

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
    if(ask & (!is.null(x$c | x$df))) invisible(readline(prompt = "Press [enter] to see the next plot..."))
  }

  ##-- Other parameters ----
  ##-- + c parameter ----
  if(!is.null(x$c)){
    graphics::plot(x = 1:n, y = x$c,
                   col = 'grey20', type = 'l', cex = 1.5,
                   xlab = 'Iteration', ylab = 'c parameter',
                   bty = bty, ...)
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
  }

  ##-- + d parameter ----
  if(!is.null(x$d)){
    graphics::plot(x = 1:n, y = x$d,
                   col = 'grey20', type = 'l', cex = 1.5,
                   xlab = 'Iteration', ylab = 'd parameter',
                   bty = bty, ...)
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
  }

  ##-- + Temporal parameter ----
  if(!is.null(x$df)){
  graphics::plot(x = 1:n, y = x$df,
                 col = 'grey20', type = 'l', cex = 1.5,
                 xlab = 'Iteration', ylab = 'Degrees of freedom',
                 bty = bty, ...)
  }
}
