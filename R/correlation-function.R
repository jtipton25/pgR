#' Generate Matern covariance matrix
#'
#' A function for generating the Matern covariance matrix
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
        # R <- geoR::matern(D, exp(theta[1]), exp(theta[2]))
        R <- fields::Matern(D, range = exp(theta[1]), smoothness = exp(theta[2]))
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





#' Plot the correlation function
#'
#' @param d is an \eqn{n \times 1}{n x 1} vector of distances 
#' @param theta is the correlation function parameters. If 
#' \code{corr_fun} is "matern", then \code{theta} is a 
#' \eqn{2}{2}-dimensional vector of Matern covariance parameter inputs
#' on the log-scale with the first element the Matern range parameter
#' \code{range} on the log-scale and second element the Matern smoothness parameter 
#' \code{nu} on the log-scale. If \code{corr_fun} is "exponential" then \code{theta} is 
#' a positive number representing the exponential correlation function
#' range parameter on the log-scale.
#' @param corr_fun is the correlation function
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)

#'
#' @return Either a ggplot object of the model trace plots (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#' @export
#'
#' @import tidyverse
#'  

plot_correlation_function <- function(d, theta, corr_fun = "exponential", base_size = 12, file = NULL, width = 16, height = 9) {
    if (any(d < 0))
        stop("The distances d must be nonzero")
    # check if d is a matrix-like
    d <- as.matrix(d)
    R <- correlation_function(d, theta, corr_fun)
    
    p <- ggplot(data.frame(correlation = R, distance = d)) +
        geom_line(aes(x = distance, y = correlation))
    
    if (is.null(file)) {
        return(p)
    } else {
        ggsave(filename = file,
               plot = p,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}


#' Plot the covariance function
#'
#' @param d is an \eqn{n \times 1}{n x 1} vector of distances 
#' @param sigma is the nugget standard deviation
#' @param sigma is the partial sill standard deviation
#' @param theta is the correlation function parameters. If 
#' \code{corr_fun} is "matern", then \code{theta} is a 
#' \eqn{2}{2}-dimensional vector of Matern covariance parameter inputs
#' on the log-scale with the first element the Matern range parameter
#' \code{range} on the log-scale and second element the Matern smoothness parameter 
#' \code{nu} on the log-scale. If \code{corr_fun} is "exponential" then \code{theta} is 
#' a positive number representing the exponential correlation function
#' range parameter on the log-scale.
#' @param corr_fun is the correlation function
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)

#'
#' @return Either a ggplot object of the model trace plots (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#' @export
#'
#' @import tidyverse
#'  

plot_covariance_function <- function(d, tau, theta, sigma = NULL, corr_fun = "exponential", base_size = 12, file = NULL, width = 16, height = 9) {
    if (any(d < 0))
        stop("The distances d must be nonzero")
    # check if d is a matrix-like
    d <- as.matrix(d)
    R <- tau^2 * correlation_function(d, theta, corr_fun)
    if (!is.null(sigma))
        R[d == 0] <- R[d == 0] + sigma^2
    
    
    p <- ggplot(data.frame(correlation = R, distance = d)) +
        geom_line(aes(x = distance, y = correlation))
    
    if (is.null(file)) {
        return(p)
    } else {
        ggsave(filename = file,
               plot = p,
               device = "png",
               width = width,
               height = height,
               units = "in")
    }
}
