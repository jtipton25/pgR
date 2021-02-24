#' Initialized default values of priors
#' 
#' A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @export 

default_priors_pgLM <- function(Y, X) {
    
    stop("The function default_priors_pgLM() has been deprecated. Please use default_priors_pg_lm() instead.")
    # priors <- list(
    #     mu_beta        = rep(0, ncol(X)),
    #     Sigma_beta     = diag(ncol(X))
    # )
    # return(priors)
}

#' Initialized default values of priors for `pg_lm()`
#' 
#' A function for setting up the default `pg_lm()` priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @export 

default_priors_pg_lm <- function(Y, X) {

    check_input_pg_lm(Y, X)
    priors <- list(
        mu_beta        = rep(0, ncol(X)),
        Sigma_beta     = 25 * diag(ncol(X))
    )
    return(priors)
}


#' Initialized default values of priors
#' 
#' A function for setting up the default pgSPLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @export

default_priors_pgSPLM <- function(Y, X, corr_fun = "exponential") {
    stop("The function default_priors_pgSPLM() has been deprecated. Please use default_priors_pg_splm() instead.")
    # check_corr_fun(corr_fun)
    # priors <- NULL
    # if (corr_fun == "matern") {
    #     priors <- list(
    #         mu_beta        = rep(0, ncol(X)),
    #         Sigma_beta     = diag(ncol(X)),
    #         alpha_tau      = 0.1,
    #         beta_tau       = 0.1,
    #         mean_nu        = -1,
    #         sd_nu          = 1,
    #         mean_range     = 0,
    #         sd_range       = 10
    #     )
    # } else if (corr_fun == "exponential") {
    #     priors <- list(
    #         mu_beta        = rep(0, ncol(X)),
    #         Sigma_beta     = diag(ncol(X)),
    #         alpha_tau      = 0.1,
    #         beta_tau       = 0.1,
    #         mean_range     = 0,
    #         sd_range       = 10
    #     )
    # }
    # return(priors)
}



#' Initialized default values of priors for `pg_splm()`
#' 
#' A function for setting up the default `pg_splm()` priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @export

default_priors_pg_splm <- function(Y, X, corr_fun = "exponential") {
    check_corr_fun(corr_fun)
    check_input_pg_lm(Y, X)
    priors <- NULL
    if (corr_fun == "matern") {
        priors <- list(
            mu_beta        = rep(0, ncol(X)),
            Sigma_beta     = 25 * diag(ncol(X)),
            alpha_tau      = 0.1,
            beta_tau       = 0.1,
            mean_nu        = -1,
            sd_nu          = 1,
            mean_range     = 0,
            sd_range       = 10
        )
    } else if (corr_fun == "exponential") {
        priors <- list(
            mu_beta        = rep(0, ncol(X)),
            Sigma_beta     = 25 * diag(ncol(X)),
            alpha_tau      = 0.1,
            beta_tau       = 0.1,
            mean_range     = 0,
            sd_range       = 10
        )
    }
    return(priors)
}


#' Initialized default values of priors for `pg_splm_mra()`
#' 
#' A function for setting up the default `pg_splm_mra()` priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @export

default_priors_pg_splm_mra <- function(Y, X) {
    check_input_pg_lm(Y, X)
    priors <- list(
        mu_beta        = rep(0, ncol(X)),
        Sigma_beta     = 25 * diag(ncol(X)),
        alpha_tau2     = 0.1,
        beta_tau2      = 0.1
    )
    return(priors)
}


#' Initialized default values of priors
#' 
#' A function for setting up the default pgSTLM priors
#'
#' @param Y is a \eqn{n \times d \times T}{n x d x T} array of compositional count data through time.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @export

default_priors_pgSTLM <- function(Y, X, corr_fun = "exponential") {
    stop("The function default_priors_pgSTLM() has been deprecated. Please use default_priors_pg_stlm() instead.")
    # check_corr_fun(corr_fun)
    # priors <- NULL
    # if (corr_fun == "matern") {
    #     priors <- list(
    #         mu_beta        = rep(0, ncol(X)),
    #         Sigma_beta     = diag(ncol(X)),
    #         alpha_tau      = 0.1,
    #         beta_tau       = 0.1,
    #         mean_nu        = -1,
    #         sd_nu          = 1,
    #         mean_range     = 0,
    #         sd_range       = 10
    #     )
    # } else if (corr_fun == "exponential") {
    #     priors <- list(
    #         mu_beta        = rep(0, ncol(X)),
    #         Sigma_beta     = diag(ncol(X)),
    #         alpha_tau      = 0.1,
    #         beta_tau       = 0.1,
    #         mean_range     = 0,
    #         sd_range       = 10
    #     )
    # }
    # return(priors)
}

#' Initialized default values of priors for `pg_stlm()`
#' 
#' A function for setting up the default `pg_stlm()` priors
#'
#' @param Y is a \eqn{n \times d \times T}{n x d x T} array of compositional count data through time.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @export

default_priors_pg_stlm <- function(Y, X, corr_fun = "exponential") {

    ## check the mcmc inputs
    if (!is.array(Y)) 
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    if (length(dim(Y)) != 3)
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    missing_idx <- matrix(FALSE, dim(Y)[1], dim(Y)[3])
    for (i in 1:dim(Y)[1]) {
        for (tt in 1:dim(Y)[3]) {
            if(!any(is.na(Y[i, , tt])))
                if(sum(Y[i, , tt]) == 0)
                    stop ("There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
        }
    }
    na_idx <- which(is.na(Y))
    if(length(na_idx) == 0) { 
        if (!is_integer(Y, length(Y)))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    } else {
        if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    }    
    if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
        stop("X must be a numeric matrix.")
    if (nrow(Y) != nrow(X))
        stop("Y and X must have the same number of rows.")
    
    check_corr_fun(corr_fun)
    priors <- NULL
    if (corr_fun == "matern") {
        priors <- list(
            mu_beta        = rep(0, ncol(X)),
            Sigma_beta     = 25 * diag(ncol(X)),
            alpha_tau      = 0.1,
            beta_tau       = 0.1,
            mean_nu        = -1,
            sd_nu          = 1,
            mean_range     = 0,
            sd_range       = 10
        )
    } else if (corr_fun == "exponential") {
        priors <- list(
            mu_beta        = rep(0, ncol(X)),
            Sigma_beta     = 25 * diag(ncol(X)),
            alpha_tau      = 0.1,
            beta_tau       = 0.1,
            mean_range     = 0,
            sd_range       = 10
        )
    }
    return(priors)
}