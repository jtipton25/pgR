#' Initialized default initial conditions
#'
#'  A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is a list of prior settings. 
#' @export

default_inits_pgLM <- function(Y, X, priors) {
    stop("The function default_inits_pgLM() has been deprecated. Please use default_inits_pg_lm() instead.")
    
    # inits <- list(
    #     beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
    # )
    # return(inits)
}


#' Initialized default initial conditions
#'
#'  A function for setting up the default `pg_lm()` initial values
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is a list of prior settings. 
#' @export

default_inits_pg_lm <- function(Y, X, priors) {
    
    inits <- list(
        beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
    )
    return(inits)
}


#' Initialized default initial conditions
#'
#'  A function for setting up the default pgSPLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is a list of prior settings. 
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logical input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @export

default_inits_pgSPLM <- function(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE) {
    stop("The function default_inits_pgSPLM() has been deprecated. Please use default_inits_pg_splm() instead.")
    # J <- ncol(Y)
    # theta_mean <- NULL
    # theta_var  <- NULL
    # if (corr_fun == "matern") {
    #     theta_mean <- c(priors$mean_range, priors$mean_nu)
    #     theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    # } else if (corr_fun == "exponential") {
    #     theta_mean <- priors$mean_range
    #     theta_var  <- priors$sd_range^2
    # }
    # 
    # inits <- list(
    #     beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
    #     tau2  = if (shared_covariance_params) {
    #         min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
    #     } else {
    #         pmin(1 / stats::rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
    #     },
    #     theta = if (shared_covariance_params) {
    #         if (corr_fun == "matern") {
    #             as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -1), 0.1))
    #         } else if (corr_fun == "exponential") {
    #             pmin(pmax(stats::rnorm(1, theta_mean, sqrt(theta_var)), -1), 0.1)
    #         }
    #     } else {
    #         if (corr_fun == "matern") {
    #             pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -1), 0.1)
    #         } else if (corr_fun == "exponential") {
    #             pmin(pmax(stats::rnorm(J-1, theta_mean, sqrt(theta_var)), -1), 0.1)
    #         }
    #     }
    # )
    # return(inits)
}

#' Initialized default initial conditions
#'
#'  A function for setting up the default `pg_splm()` initial values
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is a list of prior settings. 
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logical input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @export

default_inits_pg_splm <- function(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE) {
    
    J <- ncol(Y)
    theta_mean <- NULL
    theta_var  <- NULL
    if (corr_fun == "matern") {
        theta_mean <- c(priors$mean_range, priors$mean_nu)
        theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    } else if (corr_fun == "exponential") {
        theta_mean <- priors$mean_range
        theta_var  <- priors$sd_range^2
    }
    
    inits <- list(
        beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
        tau2  = if (shared_covariance_params) {
            min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
        } else {
            pmin(1 / stats::rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
        },
        theta = if (shared_covariance_params) {
            if (corr_fun == "matern") {
                as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -1), 0.1))
            } else if (corr_fun == "exponential") {
                pmin(pmax(stats::rnorm(1, theta_mean, sqrt(theta_var)), -1), 0.1)
            }
        } else {
            if (corr_fun == "matern") {
                pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -1), 0.1)
            } else if (corr_fun == "exponential") {
                pmin(pmax(stats::rnorm(J-1, theta_mean, sqrt(theta_var)), -1), 0.1)
            }
        }
    )
    return(inits)
}


#' Initialized default initial conditions
#'
#'  A function for setting up the default `pg_stlm()` initial values
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors is a list of prior settings. 
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logical input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @export

default_inits_pg_stlm <- function(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE) {
    
    J <- ncol(Y)
    theta_mean <- NULL
    theta_var  <- NULL
    if (corr_fun == "matern") {
        theta_mean <- c(priors$mean_range, priors$mean_nu)
        theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    } else if (corr_fun == "exponential") {
        theta_mean <- priors$mean_range
        theta_var  <- priors$sd_range^2
    }
    
    inits <- list(
        beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
        tau2  = if (shared_covariance_params) {
            min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
        } else {
            pmin(1 / stats::rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
        },
        theta = if (shared_covariance_params) {
            if (corr_fun == "matern") {
                as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -1), 0.1))
            } else if (corr_fun == "exponential") {
                pmin(pmax(stats::rnorm(1, theta_mean, sqrt(theta_var)), -1), 0.1)
            }
        } else {
            if (corr_fun == "matern") {
                pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -1), 0.1)
            } else if (corr_fun == "exponential") {
                pmin(pmax(stats::rnorm(J-1, theta_mean, sqrt(theta_var)), -1), 0.1)
            }
        },
        rho = runif(1, 0.5, 1)
    )
    return(inits)
}


