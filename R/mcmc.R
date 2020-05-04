#' MCMC wrapper for models in package \code{pollenReconstruction}
#' 
#' this function setups the model and runs the desired MCMC model.
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param params is the list of parameter settings.
#' @param priors is the list of prior settings. 
#' @param inits is the list of intial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.

mcmc <- function(Y, X, params, priors, inits = NULL, config = NULL) {
    
    ## check the input data
    check_input(Y, X)
    
    ## check the params for valid inputs
    ## function to complete parameters but keep original values
    # params <- complete_params(params)
    check_parameters(params)
    
    ## check the priors for valid inputs
    check_priors(params, priors)
    ## function to complete priors but keep original values
    # priors <- complete_priors(priors)

    check_inits(params, inits)
    
    ## other checks here

    # check_config(params, config)
    
    ## run MCMC algorithms
    # out <- list()
    # if (params$type == "bspline") {
    #     out <- mcmc_DM_bspline(y, X, params, priors)
    # }    
    # 
    # return(out)
}
    
