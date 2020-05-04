#' Convert compsitional count data to proportions
#'
#'  this function takes a matrix of compositional count data and converts the data to a proportion
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data where each of the \eqn{n}{n} rows represents a compositional count observation.
#' @export
counts_to_proportions <- function(Y) {
    
    ## take a N times d matrix y of counts and return an N times d matrix 
    ## of proportions where each of the N rows sums to 1 
    
    if (!is.matrix(Y)) {
        stop("y must be a matrix of count observations")
    }
    N <- nrow(Y)
    d <- ncol(Y)
    Y_prop <- matrix(0, N, d)
    for (i in 1:N) {
        Y_prop[i, ] <- Y[i, ] / sum(Y[i, ])
    }
    return(Y_prop)
}