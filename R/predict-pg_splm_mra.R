#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from `pg_splm()`
#' @param X is a \eqn{n \times p}{n x p} matrix of covariates at the observed locations.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of covariates at the locations where predictions are to be made. 
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of locations where observations were taken.
#' @param locs_pred is a \eqn{n_pred \times 2}{n_pred x 2} matrix of locations where predictions are to be made.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logicial input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @param force is a logicial input that determines whether to allow for predictions at more locations than 10000. The default is FALSE
#' 
#' @export 

predict_pg_splm_mra <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    corr_fun,
    shared_covariance_params,
    progress = TRUE, 
    verbose = FALSE,
    force = FALSE
) {
    
    ##
    ## check the inputs
    ##
    
    if (!inherits(out, "pg_splm_mra"))
        stop("The MCMC object out must be of class pg_splm_mra which is the output of the pg_splm_mra() function.")
    
    ## 
    ## extract the parameters 
    ##
    
    beta      <- out$beta
    MRA       <- out$MRA
    tau2      <- out$tau2
    eta       <- out$eta
    alpha     <- out$alpha
    sigma2    <- out$sigma2
    n_samples <- nrow(beta)  
    N         <- nrow(X)
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    
    if (!force) {
        if (n_pred > 10000) {
            stop("Number of prediction points must be less than 10000 -- if this is the MRA model, this might be worth increasing")
        }
    }    
    
    ## generate the MRA spatial basis
    W_pred      <- mra_wendland_2d_pred(locs_pred, MRA = MRA)$W
    
    eta_pred <- array(0, dim = c(n_samples, n_pred, J-1))
    
    if (progress) {
        message("Beginning Kriging estimates")
        progressBar <- utils::txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * n_samples)   
    

    for (k in 1:n_samples) {
        for (j in 1:(J - 1)) {
            eta_pred[k, , j] <- X_pred %*% beta[k, , j] + W_pred %*% alpha[k, , j] + rnorm(n_pred, 0, sqrt(sigma2[j]))
        } 
        if (k %in% percentage_points && progress) {
            utils::setTxtProgressBar(progressBar, k / n_samples)
        }
    }
    
    if (progress) {
        close(progressBar)
    }
    
    ## convert from eta to pi
    pi_pred <- sapply(1:n_samples, function(i) eta_to_pi(eta_pred[i, , ]), simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)
    pi_pred <- aperm(pi_pred, c(3, 1, 2))
    
    return(
        list(
            eta = eta_pred, 
            pi  = pi_pred
        )
    )
}

