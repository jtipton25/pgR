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
#' @noRd
default_inits_pgSPLM <- function(Y, X, priors) {
    
    inits <- list(
        beta  = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)),
        tau2  = min(1 / rgamma(1, priors$alpha_tau, priors$beta_tau), 10),
        theta = mvnfast::rmvn(1, c(priors$mean_nu, priors$mean_range), diag(c(priors$sd_nu^2, priors$sd_range^2)))
    )
    return(inits)
}


