#' Check initial values
#'
#'  this function check that the initial values for pgLM are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param inits is the list of current initial values
#' @keywords internal

check_inits_pgLM <- function(Y, X, inits) {
    stop("The function check_inits_pgLM() has been deprecated. Please use check_inits_pg_lm() instead.")
    ## check initial conditions for regression parameters
}


#' Check initial values
#'
#'  this function check that the initial values for `pg_lm()` are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param inits is the list of current initial values
#' @keywords internal

check_inits_pg_lm <- function(Y, X, inits) {

    check_input_pg_lm(Y, X)
    
    ## check initial conditions for regression parameters
    if(!is.null(inits)) {
        if (!is.null(inits$beta))
            if (!is_numeric_matrix(inits$beta, ncol(X), ncol(Y) - 1))
                ## check if mu_beta is a vector of the correct dimension
                stop("If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
        
        
    }
}
