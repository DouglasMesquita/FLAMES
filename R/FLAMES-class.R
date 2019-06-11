#' An object of 'FLAMES' class
#'
#' @slot p A matrix with estimated probability
#' @slot beta A matrix with coefficients samples
#' @slot c A vector with c parameter sample
#' @slot df A vector with df parameter sample
#' @slot lambda A vector with lambda parameter sample
#' @slot time Time elapsed
#' @slot call Function call
#' @slot fit_measures DIC, -2*LPML and WAIC measures
#'
#' @importFrom methods new
#'
#' @export

FLAMES_class <- setClass(Class = 'FLAMES',
                         representation(p = 'matrix',
                                        beta = 'matrix',
                                        c = 'numeric',
                                        df = 'numeric',
                                        lambda = 'numeric',
                                        time = 'numeric',
                                        call = 'call',
                                        fit_measures = 'data.frame'))
