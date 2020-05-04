#' Generate Matern covariance matrix
#'
#' A function for generatingg the Matern covariance matrix
#'
#' @param D is an \eqn{n \times n}{n x n} matrix of pairwise distances between locations
#' @param theta is a \eqn{2}{2}-dimensional vector of Matern covariance parameter inputs on the log-scale with the first element the Matern smoothness parameter \code{nu} and second element the Matern \code{range} parameter
#' @export
correlation_function <- function(D, theta) {
    ## calculate the Matern correlation using parameters theta on the log scale 
    geoR::matern(D, exp(theta[1]), exp(theta[2]))
}
