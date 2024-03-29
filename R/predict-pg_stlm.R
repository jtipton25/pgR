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
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @param posterior_mean_only is a logical input that flags whether to generate the full posterior predictive distribution (`posterior_mean_only = FALSE`) or just the posterior predictive distribution of the mean response (`posterior_mean_only = TRUE`). For large dataset, the full posterior predictive distribution can be expensive to compute and the posterior distribution of the mean response is much faster to calculte.
#' @importFrom stats toeplitz
#' 
#' @export 

predict_pg_stlm <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    corr_fun,
    shared_covariance_params,
    progress            = TRUE, 
    verbose             = FALSE,
    posterior_mean_only = TRUE
) {

    ## check the inputs
    check_corr_fun(corr_fun)
    
    if (!inherits(out, "pg_stlm"))
        stop("The MCMC object out must be of class pg_stlm which is the output of the pg_stlm() function.")

    ## 
    ## extract the parameters 
    ##
    
    beta      <- out$beta
    theta     <- out$theta
    tau2      <- out$tau2
    eta       <- out$eta
    rho       <- out$rho
    n_samples <- nrow(beta)  
    N         <- nrow(X)
    n_time    <- dim(eta)[4]
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    
    if (n_pred > 10000 & posterior_mean_only == FALSE) {
        stop("Number of prediction points must be less than 10000 if posterior_mean_only = FALSE")
    }
    
    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0
    
    D_obs      <- fields::rdist(locs)
    D_pred <- NULL
    if (!posterior_mean_only) {
        ## only calculate if we are estimating the full posterior distribution
        D_pred     <- fields::rdist(locs_pred)
    }
    D_pred_obs <- fields::rdist(locs_pred, locs)
    
    eta_pred <- array(0, dim = c(n_samples, n_pred, J-1, n_time))
    
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
        if (shared_covariance_params) {
            Sigma           <- NULL
            Sigma_unobs     <- NULL
            Sigma_unobs_obs <- NULL
            if (corr_fun == "matern") {
                Sigma           <- tau2[k] * correlation_function(D_obs, theta[k, ], corr_fun = corr_fun)
                if (!posterior_mean_only) {
                    ## only calculate if we are estimating the full posterior distribution
                    Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k, ], corr_fun = corr_fun)
                }
                Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k, ], corr_fun = corr_fun)
            } else if (corr_fun == "exponential") {
                Sigma           <- tau2[k] * correlation_function(D_obs, theta[k], corr_fun = corr_fun)
                if (!posterior_mean_only) {
                    ## only calculate if we are estimating the full posterior distribution
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
            
            
            pred_var_chol_time  <- NULL
            pred_var_chol_space <- NULL
            if (!posterior_mean_only) {
                ## only calculate if we are estimating the full posterior distribution
                
                ## time covariance matrix
                W_time <- toeplitz(c(0, 1, rep(0, n_time - 2)))
                D_time <- rowSums(W_time)
                # Q_time <- diag(D_time) - rho[k] * W_time
                Q_time <- diag(c(1, rep(1 + rho[k]^2, n_time - 2), 1)) - rho[k] * W_time
                Sigma_time <- solve(Q_time)
                
                # pred_var <-  kronecker(Sigma_time, Sigma_unobs) - kronecker(Sigma_time %*% Q_time %*% Sigma_time, Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs)))
                # pred_var2 <- kronecker(Sigma_time, Sigma_unobs -  Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs)))
                # all.equal(pred_var, pred_var2)
                # microbenchmark::microbenchmark(
                #     kronecker(Sigma_time, Sigma_unobs) - kronecker(Sigma_time %*% Q_time %*% Sigma_time, Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs))),
                #     kronecker(Sigma_time, Sigma_unobs -  Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs))),
                #     times = 5
                # )
                
                
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
            # pred_var_chol <- kronecker(pred_var_chol_time, pred_var_chol_space)
            
            # pred_var_chol <- tryCatch(
            #     chol(pred_var),
            #     error = function(e) {
            #         if (verbose)
            #             message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
            #         num_chol_failures <- num_chol_failures + 1
            #         chol(pred_var + 1e-8 * diag(n_pred))                    
            #     }
            # )               
            # pred_var <- kronecker(Sigma_time, Sigma_unobs) - (kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv)) %*% t(kronecker(Sigma_time, Sigma_unobs_obs))
            
            # all.equal(
            #     kronecker(Sigma_time, Sigma_unobs) - (kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv)) %*% t(kronecker(Sigma_time, Sigma_unobs_obs)),
            #     kronecker(Sigma_time, Sigma_unobs) - kronecker(Sigma_time %*% Q_time %*% Sigma_time, Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs)))
            # )
            # 
            # microbenchmark::microbenchmark(
            #     kronecker(Sigma_time, Sigma_unobs) - (kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv)) %*% t(kronecker(Sigma_time, Sigma_unobs_obs)),
            #     kronecker(Sigma_time, Sigma_unobs) - kronecker(Sigma_time %*% Q_time %*% Sigma_time, Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs))),
            #     chol(kronecker(Sigma_time, Sigma_unobs) - (kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv)) %*% t(kronecker(Sigma_time, Sigma_unobs_obs))),
            #     chol(kronecker(Sigma_time, Sigma_unobs) - kronecker(Sigma_time %*% Q_time %*% Sigma_time, Sigma_unobs_obs %*% (Sigma_inv %*% t(Sigma_unobs_obs)))),
            #     times  = 1
            # )
            
            for (j in 1:(J - 1)) {
                # pred_mean     <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[k, , j, ] - X %*% beta[k, , j])) + X_pred %*% beta[k, , j]
                # pred_mean <- kronecker(Sigma_time, Sigma_unobs_obs) %*% (kronecker(Q_time, Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])
                # pred_mean <- kronecker(Sigma_time %*% Q_time, Sigma_unobs_obs %*% Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j])) + as.vector(X_pred %*% beta[k, , j])
                pred_mean <- as.vector(Sigma_unobs_obs %*% Sigma_inv %*% (eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])
                
                # all.equal(
                #     as.vector(kronecker(Sigma_time, Sigma_unobs_obs) %*% (kronecker(Q_time, Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])),
                #     as.vector(Sigma_unobs_obs %*% Sigma_inv %*% (eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])
                # )
                # all.equal(
                #     as.vector(kronecker(Sigma_time %*% Q_time, Sigma_unobs_obs %*% Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j])) + as.vector(X_pred %*% beta[k, , j])),
                #     as.vector(Sigma_unobs_obs %*% Sigma_inv %*% (eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])
                # )
                # 
                # microbenchmark::microbenchmark(
                #     as.vector(kronecker(Sigma_time, Sigma_unobs_obs) %*% (kronecker(Q_time, Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])),
                #     as.vector(kronecker(Sigma_time %*% Q_time, Sigma_unobs_obs %*% Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j])) + as.vector(X_pred %*% beta[k, , j])),
                #     as.vector(Sigma_unobs_obs %*% Sigma_inv %*% (eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j]),
                #     times = 10
                # )
                
                # all.equal(
                # )
                # microbenchmark::microbenchmark(     
                #     kronecker(Sigma_time, Sigma_unobs_obs) %*% (kronecker(Q_time, Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j]),
                #     kronecker(Sigma_time %*% Q_time, Sigma_unobs_obs %*% Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j])) + as.vector(X_pred %*% beta[k, , j]),
                #     times = 5
                # )
                #     
                # all.equal(kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv), kronecker(Sigma_time %*% Q_time, Sigma_unobs_obs %*% Sigma_inv))
                # microbenchmark::microbenchmark(
                #     kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv),
                #     kronecker(Sigma_time %*% Q_time, Sigma_unobs_obs %*% Sigma_inv), 
                #     times = 5)
                

                # eta_pred[k, , j, ] <- matrix(mvnfast::rmvn(1, pred_mean, pred_var_chol, isChol = TRUE), n_pred, n_time)
                if (posterior_mean_only) {
                    eta_pred[k, , j, ] <- matrix(pred_mean, n_pred, n_time) 
                } else {
                    eta_pred[k, , j, ] <- matrix(pred_mean, n_pred, n_time) + t(pred_var_chol_space) %*% matrix(rnorm(n_pred * n_time), n_pred, n_time) %*% pred_var_chol_time
                }
                # microbenchmark::microbenchmark(
                #     matrix(mvnfast::rmvn(1, pred_mean, pred_var_chol, isChol = TRUE), n_pred, n_time),
                #     matrix(pred_mean, n_pred, n_time) + pred_var_chol_space %*% matrix(rnorm(n_pred * n_time), n_pred, n_time) %*% t(pred_var_chol_time),
                #     times = 10
                # )
            } 
        } else {
            
            pred_var_chol_time <- NULL
            if (!posterior_mean_only) {
                ## only calculate if we are estimating the full posterior distribution
                
                ## time covariance matrix
                W_time <- toeplitz(c(0, 1, rep(0, n_time - 2)))
                D_time <- rowSums(W_time)
                # Q_time <- diag(D_time) - rho[k] * W_time
                Q_time <- diag(c(1, rep(1 + rho[k]^2, n_time - 2), 1)) - rho[k] * W_time
                Sigma_time <- solve(Q_time)
                
                pred_var_chol_time <- tryCatch(
                    chol(Sigma_time),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
                        num_chol_failures <- num_chol_failures + 1
                        chol(Sigma_time + 1e-8 * diag(n_pred))                    
                    }
                )     
            }
            
            for (j in 1:(J - 1)) {
                Sigma           <- NULL
                Sigma_unobs     <- NULL
                Sigma_unobs_obs <- NULL
                if (corr_fun == "matern") {
                    Sigma           <- tau2[k, j] * correlation_function(D_obs, theta[k, j, ], corr_fun = corr_fun)
                    if (!posterior_mean_only) {
                        ## only calculate if we are estimating the full posterior distribution
                        Sigma_unobs     <- tau2[k, j] * correlation_function(D_pred, theta[k, j, ], corr_fun = corr_fun)
                    }
                    Sigma_unobs_obs <- tau2[k, j] * correlation_function(D_pred_obs, theta[k, j, ], corr_fun = corr_fun)
                } else if (corr_fun == "exponential") {
                    Sigma           <- tau2[k, j] * correlation_function(D_obs, theta[k, j], corr_fun = corr_fun)
                    if (!posterior_mean_only) {
                        ## only calculate if we are estimating the full posterior distribution
                        Sigma_unobs     <- tau2[k, j] * correlation_function(D_pred, theta[k, j], corr_fun = corr_fun)
                    }
                    Sigma_unobs_obs <- tau2[k, j] * correlation_function(D_pred_obs, theta[k, j], corr_fun = corr_fun)
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
                Sigma_inv       <- chol2inv(Sigma_chol)
                
                pred_var_chol_space <- NULL
                if (!posterior_mean_only) {
                    ## only calculate if we are estimating the full posterior distribution
                    
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
                
                pred_mean <- as.vector(Sigma_unobs_obs %*% Sigma_inv %*% (eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])
                # pred_mean <- kronecker(Sigma_time, Sigma_unobs_obs) %*% (kronecker(Q_time, Sigma_inv) %*% as.vector(eta[k, , j, ] - as.vector(X %*% beta[k, , j]))) + as.vector(X_pred %*% beta[k, , j])
                # pred_var  <- kronecker(Sigma_time, Sigma_unobs) - (kronecker(Sigma_time, Sigma_unobs_obs) %*% kronecker(Q_time, Sigma_inv)) %*% t(kronecker(Sigma_time, Sigma_unobs_obs))
                # pred_var_chol <- tryCatch(
                #     chol(pred_var),
                #     error = function(e) {
                #         if (verbose)
                #             message("The Cholesky decomposition of the prediction covariance Sigma was ill-conditioned and mildy regularized.")
                #         num_chol_failures <- num_chol_failures + 1
                #         chol(pred_var + 1e-8 * diag(n_pred * n_time))                    
                #     }
                # )
                # eta_pred[k, , j, ] <- matrix(mvnfast::rmvn(1, pred_mean, pred_var_chol, isChol = TRUE), n_pred, n_time)
                if (posterior_mean_only) {
                    eta_pred[k, , j, ] <- matrix(pred_mean, n_pred, n_time) 
                } else {
                    eta_pred[k, , j, ] <- matrix(pred_mean, n_pred, n_time) + t(pred_var_chol_space) %*% matrix(rnorm(n_pred * n_time), n_pred, n_time) %*% pred_var_chol_time
                }
            } 
        }
        if (k %in% percentage_points && progress) {
            utils::setTxtProgressBar(progressBar, k / n_samples)
        }
    }
    
    if (progress) {
        close(progressBar)
    }
    
    ## convert from eta to pi
    pi_pred <- array(0, dim = c(n_samples, n_pred, J, n_time))
    for (k in 1:n_samples) {
        for (tt in 1:n_time) {
            pi_pred[k, , , tt] <- eta_to_pi(eta_pred[k, , , tt])
        }
    }

    if (num_chol_failures > 0)
        warning("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore. To better aid in diagnosing the problem, run with vebose = TRUE")
    
    return(
        list(
            eta = eta_pred, 
            pi  = pi_pred
        )
    )
}

