#' Check parameters
#'
#'  this function check that the parameter values are properly specified
#' @param params is the list of current parameter settings
#' @noRd

check_params <- function(params) {
    
    ## not currently checking if integers, leaving this up to the user's discretion
    ## check top level model options (MCMC configuration, model and link type)
    if (!is_positive_integer(params$n_adapt, 1))
        stop("params must contain a positive integer n_adapt")
    if (!is_positive_integer(params$n_mcmc, 1))
        stop("params must contain a positive integer n_mcmc")
    if (!is_positive_integer(params$n_thin, 1))
        stop("params must contain a positive integer n_thin")
    if (!is_positive_integer(params$n_message, 1))
        stop("params must contain a positive integer n_message")
    
    # }
}







