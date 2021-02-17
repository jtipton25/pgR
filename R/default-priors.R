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
    
    priors <- list(
        mu_beta        = rep(0, ncol(X)),
        Sigma_beta     = diag(ncol(X))
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
    priors <- NULL
    if (corr_fun == "matern") {
        priors <- list(
            mu_beta        = rep(0, ncol(X)),
            Sigma_beta     = diag(ncol(X)),
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
            Sigma_beta     = diag(ncol(X)),
            alpha_tau      = 0.1,
            beta_tau       = 0.1,
            mean_range     = 0,
            sd_range       = 10
        )
    }
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
    check_corr_fun(corr_fun)
    priors <- NULL
    if (corr_fun == "matern") {
        priors <- list(
            mu_beta        = rep(0, ncol(X)),
            Sigma_beta     = diag(ncol(X)),
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
            Sigma_beta     = diag(ncol(X)),
            alpha_tau      = 0.1,
            beta_tau       = 0.1,
            mean_range     = 0,
            sd_range       = 10
        )
    }
    return(priors)
}