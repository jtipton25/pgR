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
        if (!is.null(inits[['beta']]))
            if (!is_numeric_matrix(inits[['beta']], ncol(X), ncol(Y) - 1))
                ## check if beta is a matrix of the correct dimension
                stop("If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
        if (!is.null(inits[['omega']]))
            if (!is_positive_numeric_matrix(inits[['omega']], nrow(Y), ncol(Y) - 1))
                ## check if omega is a matrix of the correct dimension
                stop("If omega is specified in inits, it must be a numeric matrix of positive values with rows equal to the number of rows of Y and columns equal to the number of columns of Y - 1.")
    }
}


#' Check initial values
#'
#'  this function check that the initial values for pgSPLM are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param inits is the list of current initial values
#' @keywords internal

check_inits_pgSPLM <- function(Y, X, inits) {
    stop("The function check_inits_pgSPLM() has been deprecated. Please use check_inits_pg_splm() instead.")
    ## check initial conditions for regression parameters
}


#' Check initial values
#'
#'  this function check that the initial values for `pg_splm()` are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param inits is the list of current initial values
#' @keywords internal

check_inits_pg_splm <- function(Y, X, locs, inits) {
    
    check_input_pg_splm(Y, X, locs)
    
    ## check initial conditions for regression parameters
    if(!is.null(inits)) {
        if (!is.null(inits[['beta']]))
            if (!is_numeric_matrix(inits[['beta']], ncol(X), ncol(Y) - 1))
                ## check if beta is a matrix of the correct dimension
                stop("If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    }
}
