#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from pgSPLM
#' @param X is a \eqn{n \times p}{n x p} matrix of covariates at the observed locations.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of covariates at the locations where predictions are to be made. 
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of locations where observations were taken.
#' @param locs_pred is a \eqn{n_pred \times 2}{n_pred x 2} matrix of locations where predictions are to be made.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @param posterior_mean_only is a logical input that flags whether to generate the full posterior predictive distribution (`posterior_mean_only = FALSE`) or just the posterior predictive distribution of the mean response (`posterior_mean_only = TRUE`). For large dataset, the full posterior predictive distribution can be expensive to compute and the posterior distribution of the mean response is much faster to calculte.
#' @importFrom stats toeplitz
#' @importFrom fields rdist
#' 
#' @export 

predict_pg_mvgp_univariate <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    corr_fun = "exponential",
    n_cores = 1L,
    progress = TRUE, 
    verbose = FALSE,
    posterior_mean_only = TRUE
) {
    
    if (posterior_mean_only) {
        message("For now, this function generates a posterior predictive draws for the posterior mean only, not the full posterior predictive distribution")
    }
    ## check the inputs
    check_corr_fun(corr_fun)
    
    ## for now, the only supported corr_fun is exponential
    if (corr_fun != "exponential") 
        stop('The only currently valid option for corr_fun is "exponential"')
    
    if (!inherits(out, "pg_mvgp_univariate"))
        stop("THe MCMC object out must be of class pg_mvgp_univariate which is the output of the pg_mvgp_univariate() function.")

    
    ## 
    ## extract the parameters 
    ##
    
    ## polya-gamma variables -- not needed
    beta      <- out$beta
    sigma2    <- out$sigma2
    eta       <- out$eta
    
    gamma     <- out$gamma
    theta     <- out$theta
    tau2      <- out$tau2
    Z         <- out$Z
    rho       <- out$rho
    n_samples <- nrow(beta)  
    N         <- nrow(X)
    n_time    <- dim(Z)[3]
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    
    if (n_pred > 15000) {
        stop("Number of prediction points must be less than 15000")
    }
    
    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0
    
    D_obs      <- rdist(locs)
    D_pred     <- rdist(locs_pred)
    D_pred_obs <- rdist(locs_pred, locs)
    
    ## intialize the covariance matric
    Sigma           <- matrix(NA, nrow(D_obs), ncol(D_obs))
    Sigma_unobs     <- matrix(NA, nrow(D_pred), ncol(D_pred))
    Sigma_unobs_obs <- matrix(NA, nrow(D_pred_obs), ncol(D_pred_obs))
    
    pred_var_chol_time <- matrix(NA, n_time, n_time)
    pred_var_chol_spac <- matrix(NA, n_pred, n_pred)
    Z_pred <- array(0, dim = c(n_samples, n_pred, n_time))
    
    if (progress) {
        message("Beginning Kriging estimates")
        progressBar <- utils::txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * n_samples)   
    
    ## parallelize this later
    
    ## the comments below are to verify that the faster calculations are equivalent
    ## to the slower but simpler mathematical representations
    ## \begin{pmatrix} \boldsymbol{\eta}_o \\ \boldsymbol{\eta}_{oos} \end{pmatrix} & \sim \operatorname{N} \left( \begin{pmatrix} \mathbf{X}_o \\ \mathbf{X}_{oos}  \end{pmatrix} \boldsymbol{\beta}, \boldsymbol{\Sigma}_{time} \otimes \begin{pmatrix} \boldsymbol{\Sigma}_o & \boldsymbol{\Sigma}_{o, oos} \\ \boldsymbol{\Sigma}_{oos, o} & \boldsymbol{\Sigma}_{oos} \end{pmatrix} \right)
    
    for (k in 1:n_samples) {

        if (corr_fun == "matern") {
            Sigma           <- tau2[k] * correlation_function(D_obs, theta[k, ], corr_fun = corr_fun)
            if (!posterior_mean_only) {
                Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k, ], corr_fun = corr_fun)
            }
            Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k, ], corr_fun = corr_fun)
        } else if (corr_fun == "exponential") {
            Sigma           <- tau2[k] * correlation_function(D_obs, theta[k], corr_fun = corr_fun)
            if (!posterior_mean_only) {
                Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k], corr_fun = corr_fun)
            }
            Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k], corr_fun = corr_fun)
        }           
        Sigma_chol <- tryCatch(
            chol(Sigma),
            error = function(e) {
                if (verbose)
                    message("The Cholesky decomposition of the observed covariance Sigma was ill-conditioned and mildy regularized.")
                num_chol_failures <- num_chol_failures + 1
                chol(Sigma + 1e-8 * diag(N))                    
            }
        )
        Sigma_inv <- chol2inv(Sigma_chol)      
        
        ## time covariance matrix
        W_time <- toeplitz(c(0, 1, rep(0, n_time - 2)))
        D_time <- rowSums(W_time)
        # Q_time <- diag(D_time) - rho[k] * W_time
        Q_time <- diag(c(1, rep(1 + rho[k]^2, n_time - 2), 1)) - rho[k] * W_time
        Sigma_time <- solve(Q_time)

        if (!posterior_mean_only) {
            pred_var_chol_time <- tryCatch(
                chol(Sigma_time),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
                    num_chol_failures <- num_chol_failures + 1
                    chol(Sigma_time + 1e-8 * diag(n_pred))                    
                }
            )     
            pred_var_chol_space <- tryCatch(
                chol(Sigma_unobs -  Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs))),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
                    num_chol_failures <- num_chol_failures + 1
                    chol(Sigma_unobs -  Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs)) + 1e-8 * diag(n_pred))                    
                }
            )     
        }
        
        pred_mean <- as.vector(Sigma_unobs_obs %*% Sigma_inv %*% (Z[k, , ] - as.vector(X %*% gamma[k, ]))) + as.vector(X_pred %*% gamma[k, ])
        
        if (posterior_mean_only) { 
            Z_pred[k, , ] <- matrix(pred_mean, n_pred, n_time)
        } else {
            Z_pred[k, , ] <- matrix(pred_mean, n_pred, n_time) + pred_var_chol_space %*% matrix(rnorm(n_pred * n_time), n_pred, n_time) %*% t(pred_var_chol_time)
        }
        
        if (k %in% percentage_points && progress) {
            utils::setTxtProgressBar(progressBar, k / n_samples)
        }
    }
    
    if (progress) {
        close(progressBar)
    }
    
    if (num_chol_failures > 0)
        warning("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore. To better aid in diagnosing the problem, run with vebose = TRUE")
    
    return(list(Z = Z_pred))
}

