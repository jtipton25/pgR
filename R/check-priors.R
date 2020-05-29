#' Check priors
#' 
#' this function check that the prior values are properly specified
#' @param params is the list of current parameter settings
#' @param priors is the list of current prior values

check_priors <- function(params, priors) {

    ## not currently checking if integers, leaving this up to the user's discretion
    
    if (params$type == "linear") {
        
    }
    
    if (params$type == "bspline") {
        ## check if bspline priors are specified
        if(!is_numeric(priors$mu_mu, 1)) 
            stop("priors must contain a numeric scalar mu_mu")
        if(!is_positive_numeric(priors$sigma_mu, 1))
           stop("priors must contain a numeric scalar sigma_mu > 0")
        if (!is_numeric_vector(priors$mu_beta_0, params$df))
            stop("priors must contain a numeric vector mu_beta_0 of length df")
        if (!is_sympd_matrix(priors$Sigma_beta_0, params$df)) 
            stop("priors must contain a symmetric positive definite matrix Sigma_beta_0 of dimension df by df")
        # if (any(is.na(priors$Sigma_beta_0_inv)) || is.null(priors$Sigma_beta_0_inv) || !is.matrix(priors$Sigma_beta_0_inv) || dim(priors$Sigma_beta_0_inv) != params$df * c(1, 1) || any(eigen(priors$Sigma_beta_0)$values <= 0)  || !isSymmetric(priors$Sigma_beta_0_inv)) 
        #     stop("priors must contain a symmetric positive definite matrix Sigma_beta_0_inv of dimension df by df")
        if (!is_numeric(priors$nu, 1))
            stop("priors must contain an integer nu > 0")
        if (!is_sympd_matrix(priors$S, params$df))
            stop("priors must contain a numeric positive definite matrix S of dimension df by df")
        # if (any(is.na(priors$S_inv)) || is.null(priors$S_inv) || !is.matrix(priors$S_inv) || dim(priors$S_inv) != params$df * c(1, 1)) 
        #     stop("priors must contain a numeric positive definite matrix S_inv of dimension df by df")
    
        }
    if (params$type == "GP") {
        
    }

    if (params$type == "BUMMER") {
        
    }

    # } else  if (params$type == "GP") {
    ## add in error checks here
    # } else if (param$type == "BUMMER") {
    ## add in error checks here
    # }
}







