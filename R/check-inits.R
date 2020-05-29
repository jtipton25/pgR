#' Check intial values
#'
#'  this function check that the initial values for pgLM are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param inits is the list of current initial values

check_inits_pgLM <- function(Y, X, inits) {
    ## check initial conditions for regression parameters
    if(!is.null(inits)) {
        if (!is.null(inits$beta))
            if (!is_numeric_matrix(inits$beta, ncol(Y) - 1, ncol(X)))
                ## check if mu_beta is a vector of the correct dimension
                stop("if beta is specified in inits, it must be a matrix with rows equal to the number of columns of Y-1 and columns equal to the number of columns of X")
        
        
    }
}
