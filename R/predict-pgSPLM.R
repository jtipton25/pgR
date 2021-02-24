#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from pgSPLM
#' @param X is a \eqn{n \times p}{n x p} matrix of covariates at the observed locations.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of covariates at the locations where predictions are to be made. 
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of locations where observations were taken.
#' @param locs_pred is a \eqn{n_pred \times 2}{n_pred x 2} matrix of locations where predictions are to be made.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params is a logicial input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' 
#' @export 

predict_pgSPLM <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    corr_fun,
    shared_covariance_params,
    n_cores = 1L,
    progress = TRUE, 
    verbose = FALSE
) {
    stop("predict_pgSPLM() is deprecated. Please use predict_pg_splm() instead.")
    
    ##
    ## check the inputs
    ##
    # if (class(out) != "pgSPLM")
    #     stop("THe MCMC object out must be of class pgSPLM which is the output of the pgSPLM() function.")
    
    # check_corr_fun(corr_fun)
    # 
    # ## 
    # ## extract the parameters 
    # ##
    # 
    # beta      <- out$beta
    # theta     <- out$theta
    # tau2      <- out$tau2
    # eta       <- out$eta
    # n_samples <- nrow(beta)  
    # N         <- nrow(X)
    # n_pred    <- nrow(X_pred)
    # J         <- dim(beta)[3] + 1
    # 
    # if (n_pred > 10000) {
    #     stop("Number of prediction points must be less than 10000")
    # }
    # 
    # ## add in a counter for the number of regularized Cholesky
    # num_chol_failures <- 0
    # 
    # D_obs      <- fields::rdist(locs)
    # D_pred     <- fields::rdist(locs_pred)
    # D_pred_obs <- fields::rdist(locs_pred, locs)
    # 
    # eta_pred <- array(0, dim = c(n_samples, n_pred, J-1))
    # 
    # if (progress) {
    #     message("Beginning Kriging estimates")
    #     progressBar <- utils::txtProgressBar(style = 3)
    # }
    # percentage_points <- round((1:100 / 100) * n_samples)   
    # 
    # ## parallelize this later
    # for (k in 1:n_samples) {
    #     if (shared_covariance_params) {
    #         if (corr_fun == "matern") {
    #             Sigma           <- tau2[k] * correlation_function(D_obs, theta[k, ], corr_fun = corr_fun)
    #             Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k, ], corr_fun = corr_fun)
    #             Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k, ], corr_fun = corr_fun)
    #         } else if (corr_fun == "exponential") {
    #             Sigma           <- tau2[k] * correlation_function(D_obs, theta[k], corr_fun = corr_fun)
    #             Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k], corr_fun = corr_fun)
    #             Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k], corr_fun = corr_fun)
    #         }           
    #         Sigma_chol <- tryCatch(
    #             chol(Sigma),
    #             error = function(e) {
    #                 if (verbose)
    #                     message("The Cholesky decomposition of the observed covariance Sigma was ill-conditioned and mildy regularized.")
    #                 num_chol_failures <- num_chol_failures + 1
    #                 chol(Sigma + 1e-8 * diag(N))                    
    #             }
    #         )
    #         Sigma_inv       <- chol2inv(Sigma_chol)        
    #         for (j in 1:(J - 1)) {
    #             pred_mean     <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[k, , j] - X %*% beta[k, , j])) + X_pred %*% beta[k, , j]
    #             pred_var      <- Sigma_unobs - (Sigma_unobs_obs %*% Sigma_inv) %*% t(Sigma_unobs_obs)
    #             pred_var_chol <- tryCatch(
    #                 chol(pred_var),
    #                 error = function(e) {
    #                     if (verbose)
    #                         message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
    #                     num_chol_failures <- num_chol_failures + 1
    #                     chol(pred_var + 1e-8 * diag(n_pred))                    
    #                 }
    #             )
    #             eta_pred[k, , j] <- mvnfast::rmvn(1, pred_mean, pred_var_chol, isChol = TRUE)
    #         } 
    #     } else {
    #         for (j in 1:(J - 1)) {
    #             if (corr_fun == "matern") {
    #                 Sigma           <- tau2[k, j] * correlation_function(D_obs, theta[k, j, ], corr_fun = corr_fun)
    #                 Sigma_unobs     <- tau2[k, j] * correlation_function(D_pred, theta[k, j, ], corr_fun = corr_fun)
    #                 Sigma_unobs_obs <- tau2[k, j] * correlation_function(D_pred_obs, theta[k, j, ], corr_fun = corr_fun)
    #             } else if (corr_fun == "exponential") {
    #                 Sigma           <- tau2[k, j] * correlation_function(D_obs, theta[k, j], corr_fun = corr_fun)
    #                 Sigma_unobs     <- tau2[k, j] * correlation_function(D_pred, theta[k, j], corr_fun = corr_fun)
    #                 Sigma_unobs_obs <- tau2[k, j] * correlation_function(D_pred_obs, theta[k, j], corr_fun = corr_fun)
    #             }
    #             
    #             Sigma_chol <- tryCatch(
    #                 chol(Sigma),
    #                 error = function(e) {
    #                     if (verbose)
    #                         message("The Cholesky decomposition of the observed covariance Sigma was ill-conditioned and mildy regularized.")
    #                     num_chol_failures <- num_chol_failures + 1
    #                     chol(Sigma + 1e-8 * diag(N))                    
    #                 }
    #             )
    #             Sigma_inv       <- chol2inv(Sigma_chol)        
    #             
    #             pred_mean <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[k, , j] - X %*% beta[k, , j])) + X_pred %*% beta[k, , j]
    #             pred_var  <- Sigma_unobs - (Sigma_unobs_obs %*% Sigma_inv) %*% t(Sigma_unobs_obs)
    #             pred_var_chol <- tryCatch(
    #                 chol(pred_var),
    #                 error = function(e) {
    #                     if (verbose)
    #                         message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
    #                     num_chol_failures <- num_chol_failures + 1
    #                     chol(pred_var + 1e-8 * diag(n_pred))                    
    #                 }
    #             )
    #             eta_pred[k, , j] <- mvnfast::rmvn(1, pred_mean, pred_var_chol, isChol = TRUE)
    #         } 
    #     }
    #     if (k %in% percentage_points && progress) {
    #         utils::setTxtProgressBar(progressBar, k / n_samples)
    #     }
    # }
    # 
    # if (progress) {
    #     close(progressBar)
    # }
    # 
    # ## convert from eta to pi
    # pi_pred <- sapply(1:n_samples, function(i) eta_to_pi(eta_pred[i, , ]), simplify = "array")
    # ## permute to be in order of MCMC samples (rows), 
    # ##    observations (columns), components (slices)
    # pi_pred <- aperm(pi_pred, c(3, 1, 2))
    # 
    # if (num_chol_failures > 0)
    #     warning("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore. To better aid in diagnosing the problem, run with vebose = TRUE")
    # 
    # return(
    #     list(
    #         eta = eta_pred, 
    #         pi  = pi_pred
    #     )
    # )
}

