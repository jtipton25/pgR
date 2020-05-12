#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from pgSPLM
#' @param X is a \eqn{n \times p}{n x p} matrix of covariates at the observed locations.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of covariates at the locations where predictions are to be made. 
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of locations where observations were taken.
#' @param locs_pred is a \eqn{n_pred \times 2}{n_pred x 2} matrix of locations where predictions are to be made.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of input variables.
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param verbose is a logicial input that determines whether to print output including a progress bar.
#' 
#' @export 

predict_pgSPLM <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    n_cores = 1L,
    verbose = TRUE
) {

    beta      <- out$beta
    theta     <- out$theta
    tau2      <- out$tau2
    eta       <- out$eta
    n_samples <- nrow(beta)    
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    
    if (n_pred > 10000) {
        stop("Number of prediction points must be less than 10000")
    }
    
    D_obs      <- fields::rdist(locs)
    D_pred     <- fields::rdist(locs_pred)
    D_pred_obs <- fields::rdist(locs_pred, locs)
    
    eta_pred <- array(0, dim = c(n_samples, n_pred, J-1))
    
    if (verbose) {
        message("Beginning Kriging estimates")
        progressBar <- txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * n_samples)   
    
    ## parallelize this later
    for (i in 1:n_samples) {
        Sigma           <- tau2[i] * correlation_function(D_obs, theta[i, ])
        Sigma_unobs     <- tau2[i] * correlation_function(D_pred, theta[i, ])
        Sigma_unobs_obs <- tau2[i] * correlation_function(D_pred_obs, theta[i, ])
        Sigma_inv       <- chol2inv(chol(Sigma))        
        for (j in 1:(J - 1)) {
            pred_mean <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[i, , j] - X %*% beta[i, , j])) + X_pred %*% beta[i, , j]
            pred_var  <- Sigma_unobs - (Sigma_unobs_obs %*% Sigma_inv) %*% t(Sigma_unobs_obs)
            eta_pred[i, , j] <- mvnfast::rmvn(1, pred_mean, pred_var)
        } 
        if (i %in% percentage_points && verbose) {
            setTxtProgressBar(progressBar, i / n_samples)
        }
    }
    

    if (verbose) {
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

