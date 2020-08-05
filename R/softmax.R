#' Softmax transformation 
#'
#' this function exaluates a softmax transformation of a vector on the real number line to the simplex using the log sum of exponentials function
#' @param x is the input
#' @export

softmax <- function (x) {
    if (!is_numeric_vector(x, length(x)))
        stop("x must be a numeric vector")
    exp(x - logsumexp(x))
}