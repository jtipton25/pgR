#' Check inputs
#' 
#' this function check that the input values are properly specified
#' 
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @noRd
#' 
check_input <- function(Y, X) {
    ## check the mcmc inputs
    if (!is.matrix(Y)) 
        stop("Y must be a matrix")
    if (!is.matrix(X)) 
        stop("X must be a matrix")
    if (nrow(Y) != nrow(X))
        stop("X and Y must have the same number of rows")
}
