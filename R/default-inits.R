#' Initialized default intial conditions
#'
#'  A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @noRd
default_inits_pgLM <- function(Y, X, priors) {
    
    inits <- list(
        beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
    )
    return(inits)
}


#' Initialized default intial conditions
#'
#'  A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors
#' @param shared_covariance_params
#' @noRd
default_inits_pgSPLM <- function(Y, X, priors, shared_covariance_params) {

    J <- ncol(Y)
    theta_mean <- c(priors$mean_range, priors$mean_nu)
    theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    
    inits <- list(
        beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
        tau2  = if (shared_covariance_params) {
            min(1 / rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
        } else {
            pmin(1 / rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
        },
        theta = if (shared_covariance_params) {
            as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), 0), 0.1))
        } else {
            pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), 0), 0.1)
        }
    )
    return(inits)
}


