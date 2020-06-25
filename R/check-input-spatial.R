#' Check inputs
#' 
#' this function check that the input values are properly specified
#' 
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of spatial locations.
#' @keywords internal

check_input_spatial <- function(Y, X, locs) {
    ## check the mcmc inputs
    if (!is.matrix(Y)) 
        stop("Y must be a matrix")
    if (!is.matrix(X)) 
        stop("X must be a matrix")
    if (!is.matrix(locs)) 
        stop("locs must be a matrix")
    if (nrow(Y) != nrow(X))
        stop("X and Y must have the same number of rows")
    if (nrow(Y) != nrow(locs))
        stop("locs and Y must have the same number of rows")
}


#' Check inputs for pgSTLM
#' 
#' this function check that the input values are properly specified
#' 
#' @param Y is a \eqn{n \times d x T}{n x d X T} array of compositional count data through time.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of spatial locations.
#' @keywords internal

check_input_pgSTLM <- function(Y, X, locs) {
    ## check the mcmc inputs
    if (!is.array(Y)) 
        stop("Y must be an array")
    if (!is.matrix(X)) 
        stop("X must be a matrix")
    if (!is.matrix(locs)) 
        stop("locs must be a matrix")
    if (nrow(Y) != nrow(X))
        stop("X and Y must have the same number of rows")
    if (nrow(Y) != nrow(locs))
        stop("locs and Y must have the same number of rows")
}
