#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from pgLM
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of climate variables.
#' @param n_cores is the number of cores for parallel computation using openMP.
#' 
#' @export 

predict_pgLM <- function(
    out,
    X_pred,
    n_cores = 1L
) {
   
    ##
    ## check the inputs 
    ##
    if (class(out) != "pgLM")
        stop("THe MCMC object out must be of class pgLM which is the output of the pgLM() function.")
    
    ## 
    ## extract the parameters 
    ##
    
    beta      <- out$beta
    n_samples <- nrow(beta)    
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    eta_pred <- array(0, dim = c(n_samples, n_pred, J - 1))
   
    eta_pred <- sapply(1:n_samples, function(i) X_pred %*% beta[i, , ], simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)
    eta_pred <- aperm(eta_pred, c(3, 1, 2))
    
    ## convert from eta to pi
    pi_pred <- sapply(1:n_samples, function(i) eta_to_pi(eta_pred[i, , ]), simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)

    ## check in J=2
    if (length(dim(pi_pred)) == 2) {
        pi_pred <- aperm(pi_pred, c(2, 1))   
    } else {
        pi_pred <- aperm(pi_pred, c(3, 1, 2))
    }    

    return(
        list(
            eta = eta_pred, 
            pi  = pi_pred
        )
    )
}

