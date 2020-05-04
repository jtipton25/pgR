#' Initialized default values of priors
#' 
#' A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @noRd

default_priors_pgLM <- function(Y, X) {
    
    priors <- list(
        mu_beta        = rep(0, ncol(X)),
        Sigma_beta     = diag(ncol(X))
    )
    return(priors)
}




#' Initialized default values of priors
#' 
#' A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @noRd

default_priors_pgSPLM <- function(Y, X) {
    
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
    return(priors)
}
