#' Check priors
#' 
#' this function check that the prior values are properly specified for the `pg_lm()` model 
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is the list of current prior values
#' @keywords internal

check_priors_pg_lm <- function(Y, X, priors) {
    ## check mu_beta prior
    if (!is.null(priors[['mu_beta']])) {
        if (any(is.na(priors[['mu_beta']]))) 
            stop("If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
        if(!is_numeric_vector(priors[['mu_beta']], ncol(X)))
            stop("If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    }
    
    ## check Sigma_beta prior
    if (!is.null(priors[['Sigma_beta']])) {
        if (any(is.na(priors[['Sigma_beta']]))) 
            stop("If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
        if(!is_sympd_matrix(priors[['Sigma_beta']], ncol(X)))
            stop("If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    }
}


#' Check priors
#' 
#' this function check that the prior values are properly specified for the `pg_splm()` model.
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is the list of current prior values
#' @param corr_fun is the correlation function
#' @keywords internal

check_priors_pg_splm <- function(Y, X, priors, corr_fun = "exponential") {
    
    ## check mu_beta prior
    if (!is.null(priors[['mu_beta']])) {
        if (any(is.na(priors[['mu_beta']]))) 
            stop("If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
        if(!is_numeric_vector(priors[['mu_beta']], ncol(X)))
            stop("If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    }
    
    ## check Sigma_beta prior
    if (!is.null(priors[['Sigma_beta']])) {
        if (any(is.na(priors[['Sigma_beta']]))) 
            stop("If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
        if(!is_sympd_matrix(priors[['Sigma_beta']], ncol(X)))
            stop("If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    }
    
    ## check alpha_tau prior
    if (!is.null(priors[['alpha_tau']])) {
        if (any(is.na(priors[['alpha_tau']]))) 
            stop("If priors contains a value for alpha_tau, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['alpha_tau']], 1))
            stop("If priors contains a value for alpha_tau, it must be a positive numeric value.")
    }
    ## check beta_tau prior
    if (!is.null(priors[['beta_tau']])) {
        if (any(is.na(priors[['beta_tau']]))) 
            stop("If priors contains a value for beta_tau, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['beta_tau']], 1))
            stop("If priors contains a value for beta_tau, it must be a positive numeric value.")
    }
    
    ## check mean_range prior
    if (!is.null(priors[['mean_range']])) {
        if (any(is.na(priors[['mean_range']]))) 
            stop("If priors contains a value for mean_range, it must be a numeric value.")
        if(!is_numeric(priors[['mean_range']], 1))
            stop("If priors contains a value for mean_range, it must be a numeric value.")
    }
    ## check sd_range prior
    if (!is.null(priors[['sd_range']])) {
        if (any(is.na(priors[['sd_range']]))) 
            stop("If priors contains a value for sd_range, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['sd_range']], 1))
            stop("If priors contains a value for sd_range, it must be a positive numeric value.")
    }
    
    ## check mean_nu prior
    if (!is.null(priors[['mean_nu']])) {
        if (any(is.na(priors[['mean_nu']]))) 
            stop("If priors contains a value for mean_nu, it must be a numeric value.")
        if(!is_numeric(priors[['mean_nu']], 1))
            stop("If priors contains a value for mean_nu, it must be a numeric value.")
    }
    ## check sd_range prior
    if (!is.null(priors[['sd_nu']])) {
        if (any(is.na(priors[['sd_nu']]))) 
            stop("If priors contains a value for sd_nu, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['sd_nu']], 1))
            stop("If priors contains a value for sd_nu, it must be a positive numeric value.")
    }
}

#' Check priors
#' 
#' this function check that the prior values are properly specified for the `pg_splm()` model.
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is the list of current prior values
#' @param corr_fun is the correlation function
#' @keywords internal

check_priors_pg_stlm <- function(Y, X, priors, corr_fun = "exponential") {
    
    ## check mu_beta prior
    if (!is.null(priors[['mu_beta']])) {
        if (any(is.na(priors[['mu_beta']]))) 
            stop("If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
        if(!is_numeric_vector(priors[['mu_beta']], ncol(X)))
            stop("If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    }
    
    ## check Sigma_beta prior
    if (!is.null(priors[['Sigma_beta']])) {
        if (any(is.na(priors[['Sigma_beta']]))) 
            stop("If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
        if(!is_sympd_matrix(priors[['Sigma_beta']], ncol(X)))
            stop("If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    }
    
    ## check alpha_tau prior
    if (!is.null(priors[['alpha_tau']])) {
        if (any(is.na(priors[['alpha_tau']]))) 
            stop("If priors contains a value for alpha_tau, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['alpha_tau']], 1))
            stop("If priors contains a value for alpha_tau, it must be a positive numeric value.")
    }
    ## check beta_tau prior
    if (!is.null(priors[['beta_tau']])) {
        if (any(is.na(priors[['beta_tau']]))) 
            stop("If priors contains a value for beta_tau, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['beta_tau']], 1))
            stop("If priors contains a value for beta_tau, it must be a positive numeric value.")
    }
    
    ## check mean_range prior
    if (!is.null(priors[['mean_range']])) {
        if (any(is.na(priors[['mean_range']]))) 
            stop("If priors contains a value for mean_range, it must be a numeric value.")
        if(!is_numeric(priors[['mean_range']], 1))
            stop("If priors contains a value for mean_range, it must be a numeric value.")
    }
    ## check sd_range prior
    if (!is.null(priors[['sd_range']])) {
        if (any(is.na(priors[['sd_range']]))) 
            stop("If priors contains a value for sd_range, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['sd_range']], 1))
            stop("If priors contains a value for sd_range, it must be a positive numeric value.")
    }
    
    ## check mean_nu prior
    if (!is.null(priors[['mean_nu']])) {
        if (any(is.na(priors[['mean_nu']]))) 
            stop("If priors contains a value for mean_nu, it must be a numeric value.")
        if(!is_numeric(priors[['mean_nu']], 1))
            stop("If priors contains a value for mean_nu, it must be a numeric value.")
    }
    ## check sd_range prior
    if (!is.null(priors[['sd_nu']])) {
        if (any(is.na(priors[['sd_nu']]))) 
            stop("If priors contains a value for sd_nu, it must be a positive numeric value.")
        if(!is_positive_numeric(priors[['sd_nu']], 1))
            stop("If priors contains a value for sd_nu, it must be a positive numeric value.")
    }
}
