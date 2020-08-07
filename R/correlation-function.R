#' Generate Matern covariance matrix
#'
#' A function for generatingg the Matern covariance matrix
#'
#' @param D is an \eqn{n \times n}{n x n} matrix of pairwise distances between locations
#' @param theta is the correlation function parameters. If 
#' \code{corr_fun} is "matern", then \code{theta} is a 
#' \eqn{2}{2}-dimensional vector of Matern covariance parameter inputs
#' on the log-scale with the first element the Matern range parameter
#' \code{range} on the log-scale and second element the Matern smoothness parameter 
#' \code{nu} on the log-scale. If \code{corr_fun} is "exponential" then \code{theta} is 
#' a positive number representing the exponential correlation function
#' range parameter on the log-scale.
#' @param corr_fun is the correlation function
#' @export
 
correlation_function <- function(D, theta, corr_fun = "exponential") {
    R <- NULL
    
    if (!is_dist_matrix(D, nrow(D), ncol(D)))
        stop("D must be a distance matrix with only non-negative values.")
    
    check_corr_fun(corr_fun)
    if (corr_fun == "matern") {
        if (!is_numeric_vector(theta, 2))
            stop('theta must be a numeric vector of length 2 for the matern correlation function.')
        ## calculate the Matern correlation using parameters theta on the log scale 
        R <- geoR::matern(D, exp(theta[1]), exp(theta[2]))
    } else if (corr_fun == "exponential") {
        if (!is_numeric(theta, 1))
            stop('theta must be a numeric value for the exponential correlation function.')
        R <- exp( - D / exp(theta))
    }
    return(R)
}

#' Check if the correlation function type is valid
#'
#' this function checks if the correlation function type is valid
#' @param corr_fun is the correlation function
#' @keywords internal

check_corr_fun <- function(corr_fun) {
    if (length(corr_fun) != 1) 
        stop('corr_fun must be either "matern" or "exponential"')
    if (!(corr_fun %in% c("matern", "exponential"))) 
        stop('corr_fun must be either "matern" or "exponential"')
}



