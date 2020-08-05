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
    }
}


#' Check initial values
#'
#'  this function check that the initial values for pgSPLM are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param inits is the list of current initial values
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @keywords internal

check_inits_pgSPLM <- function(Y, X, locs, inits, corr_fun = "exponential") {
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
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logical input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @keywords internal

check_inits_pg_splm <- function(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE) {
    
    check_corr_fun(corr_fun)
    check_input_pg_splm(Y, X, locs)
    
    ## check initial conditions for regression parameters
    if(!is.null(inits)) {
        if (!is.null(inits[['beta']]))
            if (!is_numeric_matrix(inits[['beta']], ncol(X), ncol(Y) - 1))
                ## check if beta is a matrix of the correct dimension
                stop("If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
        
        ## check covariance function parameters
        if (shared_covariance_params) {
            if (corr_fun == "exponential") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_vector(inits[['theta']], 1))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric value when corr_fun = "exponential" and shared_covariance_params = TRUE.')
            } else if (corr_fun == "matern") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_vector(inits[['theta']], 2))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = "matern" and shared_covariance_params = TRUE.')
            }
            
            if (!is.null(inits[['tau2']]))
                if (!is_positive_numeric(inits[['tau2']], 1))
                    ## check if tau2 is the correct dimension
                    stop('If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.')
            
        } else {
            if (corr_fun == "exponential") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_vector(inits[['theta']], ncol(Y) - 1))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = "exponential" and shared_covariance_params = FALSE.')
            } else if (corr_fun == "matern") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_matrix(inits[['theta']], ncol(Y) - 1, 2))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = "matern" and shared_covariance_params = FALSE.')
            }
            
            if (!is.null(inits[['tau2']]))
                if (!is_positive_numeric(inits[['tau2']], ncol(Y) - 1))
                    ## check if tau2 is the correct dimension
                    stop('If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.')
        }
        
    }
}


#' Check initial values
#'
#'  this function check that the initial values for pgSTLM are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param inits is the list of current initial values
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @keywords internal
check_inits_pgSTLM <- function(Y, X, inits) {
    stop("The function check_inits_pgSTLM() has been deprecated. Please use check_inits_pg_stlm() instead.")
    ## check initial conditions for regression parameters
}


#' Check initial values
#'
#'  this function check that the initial values for `pg_stlm()` are properly specified
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param inits is the list of current initial values
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logical input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @keywords internal

check_inits_pg_stlm <- function(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE) {
    
    check_corr_fun(corr_fun)
    check_input_pg_stlm(Y, X, locs)
    
    ## check initial conditions for regression parameters
    if(!is.null(inits)) {
        if (!is.null(inits[['beta']]))
            if (!is_numeric_matrix(inits[['beta']], ncol(X), ncol(Y) - 1))
                ## check if beta is a matrix of the correct dimension
                stop("If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
        
        ## check covariance function parameters
        if (shared_covariance_params) {
            if (corr_fun == "exponential") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_vector(inits[['theta']], 1))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric value when corr_fun = "exponential" and shared_covariance_params = TRUE.')
            } else if (corr_fun == "matern") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_vector(inits[['theta']], 2))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = "matern" and shared_covariance_params = TRUE.')
            }
            
            if (!is.null(inits[['tau2']]))
                if (!is_positive_numeric(inits[['tau2']], 1))
                    ## check if theta is the correct dimension
                    stop('If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.')
            
        } else {
            if (corr_fun == "exponential") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_vector(inits[['theta']], ncol(Y) - 1))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = "exponential" and shared_covariance_params = FALSE.')
            } else if (corr_fun == "matern") {
                if (!is.null(inits[['theta']]))
                    if (!is_numeric_matrix(inits[['theta']], ncol(Y) - 1, 2))
                        ## check if theta is the correct dimension
                        stop('If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = "matern" and shared_covariance_params = FALSE.')
            }
            
            if (!is.null(inits[['tau2']]))
                if (!is_positive_numeric(inits[['tau2']], ncol(Y) - 1))
                    ## check if theta is the correct dimension
                    stop('If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.')
        }
        
        if (!is.null(inits[['rho']]))
            if (!is_numeric(inits[['rho']], 1) | inits[['rho']] < -1 | inits[['rho']] > 1)
                ## check if rho is of the correct dimension and value
                stop("If rho is specified in inits, it must be a numeric value between -1 and 1.")
        
        
        if (!is.null(inits[['eta']]))
            if (!is_numeric(inits[['eta']], nrow(Y) * (ncol(Y) - 1) * dim(Y)[3]) | !is.array(Y) | length(dim(Y)) != 3)
                ## check if eta is an array of the correct dimension
                stop("If eta is specified in inits, it must be a numeric array with dimension equal to the dimensions of Y.")
        
    }
}

