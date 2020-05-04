#' log sum of exponentials
#'
#' this function exaluates a log sum of exponential
#' @param x is the input
#' @export
logsumexp <- function (x) {
    ## from http://tr.im/hH5A
    y <- max(x)
  y + log(sum(exp(x - y)))
}
