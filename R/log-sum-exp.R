#' log sum of exponentials
#'
#' this function evaluates a log sum of exponential
#' @param x is the input
#' @export
#' @keywords internal
  
logsumexp <- function (x) {
    if (!is_numeric_vector(x, length(x)))
        stop("x must be a numeric vector")
    y <- max(x)
    y + log(sum(exp(x - y)))
}
