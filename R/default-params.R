#' Initialized default values of parameters
#' 
#' A function for setting up the default model parameters
#' @export

default_params <- function() {
    return(
        list(
            n_mcmc    = 100L,
            n_adapt   = 100L,
            n_thin    = 1L,
            n_message = 500L
        )
    )
}

