#' Bayesian Polya-gamma regression
#' 
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param Y is a \eqn{n \times J}{n x J} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of observation locations.
#' @param params is a list of parameter settings. The list
#' \code{params} must contain the following values:
#' * \code{n_adapt}: A positive integer number of adaptive MCMC iterations.
#' * \code{n_mcmc}: A positive integer number of total MCMC iterations
#' post adaptation.
#' * \code{n_thin}: A positive integer number of MCMC iterations per saved
#' sample.
#' * \code{n_message}: A positive integer number of frequency of iterations
#'  to output a progress message. For example, \code{n_message = 50}
#'  outputs progress messages every 50 iterations.
#' @param priors is the list of prior settings. 
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param shared_covariance_params is a logical input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @param n_chain is the MCMC chain id. The default is 1.
#' @param progress is a logical input that determines whether to print a progress bar.
#' @param verbose is a logical input that determines whether to print more detailed messages.
#' @importFrom stats rnorm rgamma runif dnorm
#' @importFrom fields rdist
#' @importFrom mvnfast rmvn dmvn
#' 
#' @export

pgSPLM <- function(
    Y, 
    X,
    locs, 
    params,
    priors,
    corr_fun = "exponential",
    shared_covariance_params = TRUE,
    n_cores = 1L,
    inits = NULL,
    config = NULL,
    n_chain       = 1,
    progress = FALSE,
    verbose = FALSE
    # pool_s2_tau2  = true,
    # file_name     = "DM-fit",
    # corr_function = "exponential"
) {
    
    stop("The function pgSPLM() has been deprecated. Please use pg_splm() instead.")
    # ##
    # ## Run error checks
    # ## 
    # 
    # check_input_spatial(Y, X, locs)
    # check_params(params)
    # check_corr_fun(corr_fun)
    # # check_inits_pgLM(params, inits)
    # # check_config(params, config)
    # 
    # ## add in faster parallel cholesky as needed
    # 
    # ## add in a counter for the number of regularized Cholesky
    # num_chol_failures <- 0
    # 
    # 
    # N  <- nrow(Y)
    # J  <- ncol(Y)
    # p <- ncol(X)
    # D <- rdist(locs)
    # 
    # ## Calculate Mi
    # Mi <- matrix(0, N, J-1)
    # for (i in 1:N){
    #     Mi[i,] <- sum(Y[i, ]) - c(0, cumsum(Y[i, ][1:(J - 2)]))
    # }
    # 
    # ## create an index for nonzero values
    # nonzero_idx <- Mi != 0
    # 
    # ## initialize kappa
    # kappa <- matrix(0, N, J-1)
    # for (i in 1:N) {
    #     kappa[i, ] <- Y[i, 1:(J - 1)] - Mi[i, ] / 2
    # }
    # 
    # ##
    # ## initial values
    # ##
    # 
    # ## currently using default priors
    # 
    # mu_beta        <- rep(0, p)
    # 
    # ## do I want to change this to be a penalized spline?
    # # Q_beta <- make_Q(params$p, 1) 
    # Sigma_beta     <- 10 * diag(p)
    # ## clean up this check
    # if (!is.null(priors[['mu_beta']])) {
    #     if (all(!is.na(priors[['mu_beta']]))) {
    #         mu_beta <- priors[['mu_beta']]
    #     }
    # }
    # 
    # ## clean up this check
    # if (!is.null(priors[['Sigma_beta']])) {
    #     if (all(!is.na(priors[['Sigma_beta']]))) {
    #         Sigma_beta <- priors[['Sigma_beta']]
    #     }
    # }
    # Sigma_beta_chol <- tryCatch(
    #     chol(Sigma_beta),
    #     error = function(e) {
    #         if (verbose)
    #             message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
    #         chol(Sigma_beta + 1e-8 * diag(N))                    
    #     }
    # )
    # Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    # 
    # ##
    # ## initialize beta
    # ##
    # 
    # beta <- t(rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    # ## clean up this check
    # if (!is.null(inits[['beta']])) {
    #     if (all(!is.na(inits[['beta']]))) {
    #         beta <- inits[['beta']]
    #     }
    # }
    # Xbeta <- X %*% beta
    # 
    # ##
    # ## initialize spatial Gaussian process -- share parameters across the different components
    # ##    can generalize to each component getting its own covariance
    # ##
    # 
    # theta_mean <- NULL
    # theta_var  <- NULL
    # if (corr_fun == "matern") {
    #     theta_mean <- c(priors$mean_range, priors$mean_nu)
    #     theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    # } else if (corr_fun == "exponential") {
    #     theta_mean <- priors$mean_range
    #     theta_var  <- priors$sd_range^2
    # }
    # 
    # theta <- NULL
    # if (shared_covariance_params) {
    #     if (corr_fun == "matern") {
    #         theta <- as.vector(pmin(pmax(rmvn(1, theta_mean, theta_var), -2), 0.1))            
    #     } else if (corr_fun == "exponential") {
    #         theta <- pmin(pmax(rnorm(1, theta_mean, sqrt(theta_var)), -2), 0.1)
    #     }
    # } else {
    #     if (corr_fun == "matern") {
    #         theta <- pmin(pmax(rmvn(J-1, theta_mean, theta_var), -2), 0.1)
    #     } else if (corr_fun == "exponential") {
    #         theta <- pmin(pmax(rnorm(J-1, theta_mean, sqrt(theta_var)), -2), 0.1)
    #     }
    #     
    # }
    # 
    # if (!is.null(inits[['theta']])) {
    #     if (all(!is.na(inits[['theta']]))) {
    #         theta <- inits[['theta']]
    #     }
    # }
    # ## check dimensions of theta
    # if (shared_covariance_params) {
    #     if (corr_fun == "matern")
    #         if (!is_numeric_vector(theta, 2)) 
    #             stop('If shared_covariance_params is TRUE, theta must be a numeric vector of length 2 when corr_fun is "matern"')
    #     if (corr_fun == "exponential")
    #         if (!is_numeric_vector(theta, 1)) 
    #             stop('If shared_covariance_params is TRUE, theta must be a numeric of length 1 when corr_fun is "exponential"')
    # } else {
    #     if (corr_fun == "matern")
    #         if (!is_numeric_matrix(theta, J-1, 2))
    #             stop('If shared_covariance_params is FALSE, theta must be a J-1 by 2 numeric matrix when corr_fun is "matern"')
    #     if (corr_fun == "exponential")
    #         if (!is_numeric_vector(theta, J-1)) 
    #             stop('If shared_covariance_params is FALSE, theta must be a numeric vector of length J-1 when corr_fun is "exponential"')
    # }
    # 
    # tau2 <- NULL
    # if (shared_covariance_params) {
    #     tau2 <- min(1 / rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
    # } else {
    #     tau2 <- pmin(1 / rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
    # }
    # if (!is.null(inits[['tau2']])) {
    #     if (all(!is.na(inits[['tau2']]))) {
    #         ## if tau2 passes error checks
    #         tau2 <- inits[['tau2']]
    #     }
    # }
    # 
    # ## check dimensions of tau2
    # if (shared_covariance_params) {
    #     if (!is_positive_numeric(tau2, 1))
    #         stop("If shared_covariance_params is FALSE, tau2 must be a numeric scalar")
    # } else {
    #     if (!is_positive_numeric(tau2, J-1))
    #         stop("If shared_covariance_params is TRUE, tau2 must be a J-1 positive numeric vector ")
    # }
    # 
    # 
    # Sigma <- NULL
    # if (shared_covariance_params) {
    #     Sigma <- tau2 * correlation_function(D, theta, corr_fun = corr_fun)
    # } else {
    #     Sigma <- array(0, dim = c(J-1, N, N))
    #     for (j in 1:(J-1)) {
    #         if (corr_fun == "matern") {
    #             Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = corr_fun)
    #         } else if (corr_fun == "exponential") {
    #             Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j], corr_fun = corr_fun)
    #         }
    #     }
    # }
    # 
    # 
    # 
    # Sigma_chol <- NULL
    # if (shared_covariance_params) {
    #     Sigma_chol <- tryCatch(
    #         chol(Sigma[j, , ]),
    #         error = function(e) {
    #             if (verbose)
    #                 message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #             num_chol_failures <- num_chol_failures + 1
    #             chol(Sigma + 1e-8 * diag(N))                    
    #         }
    #     )
    # } else {
    #     Sigma_chol <- array(0, dim = c(J-1, N, N))
    #     for (j in 1:(J-1)) {
    #         ## add a warning for the Cholesky function
    #         Sigma_chol[j, , ] <- tryCatch(
    #             chol(Sigma[j, , ]),
    #             error = function(e) {
    #                 if (verbose)
    #                     message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #                 num_chol_failures <- num_chol_failures + 1
    #                 chol(Sigma[j, , ] + 1e-8 * diag(N))                    
    #             }
    #         )
    #     }
    # }
    # 
    # Sigma_inv  <- NULL
    # if (shared_covariance_params) {
    #     Sigma_inv <- chol2inv(Sigma_chol)   
    # } else {
    #     Sigma_inv <- array(0, dim = c(J-1, N, N))
    #     for (j in 1:(J-1)) {
    #         Sigma_inv[j, , ] <- chol2inv(Sigma_chol[j, , ])
    #     }
    # }
    # 
    # 
    # eta  <- NULL
    # if (shared_covariance_params) {
    #     eta <- Xbeta + t(rmvn(J-1, rep(0, N), Sigma_chol, isChol = TRUE))
    # } else{
    #     eta <- matrix(0, N, J-1)
    #     for (j in 1:(J-1)) {
    #         eta[, j] <- Xbeta[, j] + t(rmvn(1, rep(0, N), Sigma_chol[j, , ], isChol = TRUE))
    #     }
    # }
    # 
    # ##
    # ## sampler config options -- to be added later
    # ## 
    # #
    # # bool sample_beta = true;
    # # if (params.containsElementNamed("sample_beta")) {
    # #     sample_beta = as<bool>(params["sample_beta"]);
    # # }
    # # 
    # 
    # ##
    # ## initialize omega
    # ##
    # 
    # omega <- matrix(0, N, J-1)
    # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    # 
    # if (!is.null(inits[['omega']])) {
    #     if (!is.na(inits[['omega']])) {
    #         omega <- inits[['omega']]
    #     }
    # }
    # 
    # Omega <- vector(mode = "list", length = J-1)
    # for (j in 1:(J - 1)) {
    #     Omega[[j]] <- diag(omega[, j])
    # }
    # 
    # ##
    # ## setup save variables
    # ##
    # 
    # n_save     <- params$n_mcmc / params$n_thin
    # beta_save  <- array(0, dim = c(n_save, p, J-1))
    # tau2_save  <- NULL
    # if (shared_covariance_params) {
    #     tau2_save <- rep(0, n_save)
    # } else {
    #     tau2_save <- matrix(0, n_save, J-1)
    # }
    # theta_save <- NULL
    # if (shared_covariance_params) {
    #     if (corr_fun == "matern") {
    #         theta_save <- matrix(0, n_save, 2)
    #     } else if (corr_fun == "exponential") {
    #         theta_save <- rep(0, n_save)
    #     }
    # } else {
    #     if (corr_fun == "matern") {
    #         theta_save <- array(0, dim = c(n_save, J-1, 2))
    #     } else if (corr_fun == "exponential") {
    #         theta_save <- matrix(0, n_save, J-1)
    #     }
    # }
    # eta_save   <- array(0, dim = c(n_save, N, J-1))
    # 
    # ## 
    # ## initialize tuning 
    # ##
    # 
    # ##
    # ## tuning variables for adaptive MCMC
    # ##
    # 
    # theta_batch           <- NULL
    # theta_accept          <- NULL
    # theta_accept_batch    <- NULL
    # lambda_theta          <- NULL
    # Sigma_theta_tune      <- NULL
    # Sigma_theta_tune_chol <- NULL
    # theta_tune            <- NULL
    # 
    # if (shared_covariance_params) {
    #     
    #     theta_accept       <- 0
    #     theta_accept_batch <- 0
    #     
    #     if (corr_fun == "matern") {
    #         theta_batch <- matrix(0, 50, 2) 
    #         lambda_theta          <- 0.05
    #         Sigma_theta_tune      <- 1.8 * diag(2) - .8
    #         Sigma_theta_tune_chol <- tryCatch(
    #             chol(Sigma_theta_tune),
    #             error = function(e) {
    #                 if (verbose)
    #                     message("The Cholesky decomposition of the Metroplois-Hastings adaptive tuning matrix for Matern parameters theta was ill-conditioned and mildy regularized.")
    #                 chol(Sigma_theta_tune + 1e-8 * diag(2))                    
    #             }
    #         )
    #     } else if (corr_fun == "exponential") {
    #         theta_tune <- mean(D) / 2
    #     }
    #     
    # } else {
    #     
    #     theta_accept       <- rep(0, J-1)
    #     theta_accept_batch <- rep(0, J-1)
    #     
    #     if (corr_fun == "matern") {
    #         theta_batch <- array(0, dim = c(50, 2, J-1))      
    #         lambda_theta     <- rep(0.05, J-1)
    #         Sigma_theta_tune <- array(0, dim = c(2, 2, J-1))
    #         for (j in 1:(J-1)) {
    #             Sigma_theta_tune[, , j] <- 1.8 * diag(2) - .8
    #         }
    #         Sigma_theta_tune_chol <- array(0, dim = c(2, 2, J-1))
    #         for (j in 1:(J-1)) {
    #             Sigma_theta_tune_chol[, , j] <- tryCatch(
    #                 chol(Sigma_theta_tune[,,j]),
    #                 error = function(e) {
    #                     if (verbose)
    #                         message("The Cholesky decomposition of the Metroplois-Hastings adaptive tuning matrix for Matern parameters theta was ill-conditioned and mildy regularized.")
    #                     chol(Sigma_theta_tune[, , j] + 1e-8 * diag(2))                    
    #                 }
    #             )
    #         }
    #     } else if (corr_fun == "exponential") {
    #         theta_tune <- rep(mean(D) / 2, J-1)
    #     }
    # }
    # 
    # ##
    # ## Starting MCMC chain
    # ##
    # 
    # message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")
    # if (progress) {
    #     progressBar <- utils::txtProgressBar(style = 3)
    # }
    # percentage_points <- round((1:100 / 100) * (params$n_adapt + params$n_mcmc))
    # 
    # 
    # for (k in 1:(params$n_adapt + params$n_mcmc)) {
    #     if (k == params$n_adapt + 1) {
    #         message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations \n")
    #     }
    #     if (k %% params$n_message == 0) {
    #         if (k <= params$n_adapt) {
    #             message("MCMC adaptation iteration ", k, " for chain ", n_chain)
    #         } else {
    #             message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
    #         }
    #     }
    #     
    #     ##
    #     ## sample Omega
    #     ##
    #     if (verbose)
    #         message("sample omega")
    #     
    #     omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    #     
    #     # for (i in 1:N) {
    #     #     for (j in 1:(J-1)) {
    #     #         if (Mi[i, j] != 0){
    #     #             omega[i, j] <- pgdraw(Mi[i, j], eta[i, j])
    #     #         }
    #     #         else {
    #     #             omega[i, j] <- 0
    #     #         }
    #     #     }
    #     # }
    #     
    #     for (j in 1:(J-1)) {
    #         Omega[[j]] <- diag(omega[, j])
    #     }
    #     
    #     ##
    #     ## sample beta -- double check these values
    #     ##
    #     
    #     ## modify this for the spatial process eta
    #     
    #     ## parallelize this update -- each group of parameteres is 
    #     ## conditionally independent given omega and kappa(y)
    #     if (verbose)
    #         message("sample beta")
    #     
    #     if (shared_covariance_params) {
    #         for (j in 1:(J-1)) {
    #             ## can make this much more efficient
    #             # Sigma_tilde <- chol2inv(chol(Sigma_beta_inv + t(X) %*% (Omega[[j]] %*% X))) 
    #             # mu_tilde    <- c(Sigma_tilde %*% (Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j]))
    #             # beta[, j]   <- rmvn(1, mu_tilde, Sigma_tilde)
    #             tXSigma_inv <- t(X) %*% Sigma_inv
    #             A <- tXSigma_inv %*% X + Sigma_beta_inv
    #             ## guarantee a symmetric matrix
    #             A <- (A + t(A)) / 2
    #             b <- tXSigma_inv %*% eta[, j] + Sigma_beta_inv %*% mu_beta
    #             beta[, j]   <- rmvn_arma(A, b)
    #         }
    #     } else {
    #         for (j in 1:(J-1)) {
    #             ## can make this much more efficient
    #             # Sigma_tilde <- chol2inv(chol(Sigma_beta_inv + t(X) %*% (Omega[[j]] %*% X))) 
    #             # mu_tilde    <- c(Sigma_tilde %*% (Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j]))
    #             # beta[, j]   <- rmvn(1, mu_tilde, Sigma_tilde)
    #             tXSigma_inv <- t(X) %*% Sigma_inv[j, , ]
    #             A <- tXSigma_inv %*% X + Sigma_beta_inv
    #             ## guarantee a symmetric matrix
    #             A <- (A + t(A)) / 2
    #             b <- tXSigma_inv %*% eta[, j] + Sigma_beta_inv %*% mu_beta
    #             beta[, j]   <- rmvn_arma(A, b)
    #         }
    #     }        
    #     Xbeta <- X %*% beta
    #     
    #     ##
    #     ## sample spatial correlation parameters theta
    #     ##
    #     
    #     if (verbose)
    #         message("sample theta")
    #     
    #     if (shared_covariance_params) {
    #         ## update a common theta for all processes
    #         theta_star <- NULL
    #         if (corr_fun == "matern") {
    #             theta_star <- as.vector(
    #                 rmvn( 
    #                     n      = 1,
    #                     mu     = theta,
    #                     sigma  = lambda_theta * Sigma_theta_tune_chol,
    #                     isChol = TRUE
    #                 )
    #             )
    #         } else if (corr_fun == "exponential") {
    #             theta_star <- rnorm(1, theta, theta_tune)
    #         }
    #         Sigma_star       <- tau2 * correlation_function(D, theta_star, corr_fun = corr_fun)
    #         ## add in faster parallel cholesky as needed
    #         Sigma_chol_star <- tryCatch(
    #             chol(Sigma_star),
    #             error = function(e) {
    #                 if (verbose)
    #                     message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #                 num_chol_failures <- num_chol_failures + 1
    #                 chol(Sigma_star + 1e-8 * diag(N))                    
    #             }
    #         )
    #         Sigma_inv_star  <- chol2inv(Sigma_chol_star)
    #         ## parallelize this
    #         mh1 <- sum(
    #             sapply(1:(J-1), function(j) {
    #                 dmvn(eta[, j], Xbeta[, j], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
    #             })
    #         ) +
    #             ## prior
    #             dmvn(theta_star, theta_mean, theta_var, log = TRUE)
    #         ## parallelize this        
    #         mh2 <- sum(
    #             sapply(1:(J-1), function(j) {
    #                 dmvn(eta[, j], Xbeta[, j], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)
    #             })
    #         ) +
    #             ## prior
    #             dmvn(theta, theta_mean, theta_var, log = TRUE)
    #         
    #         mh <- exp(mh1 - mh2)
    #         if (mh > runif(1, 0, 1)) {
    #             theta      <- theta_star
    #             Sigma      <- Sigma_star
    #             Sigma_chol <- Sigma_chol_star
    #             Sigma_inv  <- Sigma_inv_star 
    #             if (k <= params$n_adapt) {
    #                 theta_accept_batch <- theta_accept_batch + 1 / 50
    #             } else {
    #                 theta_accept <- theta_accept + 1 / params$n_mcmc
    #             }
    #         }
    #         ## adapt the tuning
    #         if (k <= params$n_adapt) {
    #             if (corr_fun == "matern") {
    #                 save_idx <- k %% 50
    #                 if ((k %% 50) == 0) {
    #                     save_idx <- 50
    #                 } 
    #                 theta_batch[save_idx, ] <- theta 
    #             }
    #             if (k %% 50 == 0) {
    #                 if (corr_fun == "matern") {
    #                     out_tuning <- update_tuning_mv(
    #                         k,
    #                         theta_accept_batch,
    #                         lambda_theta,
    #                         theta_batch,
    #                         Sigma_theta_tune,
    #                         Sigma_theta_tune_chol
    #                     )
    #                     theta_batch           <- out_tuning$batch_samples
    #                     Sigma_theta_tune      <- out_tuning$Sigma_tune
    #                     Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
    #                     lambda_theta          <- out_tuning$lambda
    #                     theta_accept_batch    <- out_tuning$accept
    #                 } else if (corr_fun == "exponential") {
    #                     out_tuning <- update_tuning(k, theta_accept_batch, theta_tune)
    #                     theta_tune         <- out_tuning$tune
    #                     theta_accept_batch <- out_tuning$accept
    #                 }
    #             }
    #         }
    #     } else {
    #         ## 
    #         ## theta varies for each component
    #         ##
    #         for (j in 1:(J-1)) {
    #             theta_star <- NULL
    #             if (corr_fun == "matern") {
    #                 theta_star <- as.vector(
    #                     rmvn( 
    #                         n      = 1,
    #                         mu     = theta[j, ],
    #                         sigma  = lambda_theta[j] * Sigma_theta_tune_chol[, , j],
    #                         isChol = TRUE
    #                     )
    #                 )
    #             } else if (corr_fun == "exponential") {
    #                 theta_star <- rnorm(1, theta[j], theta_tune)
    #             }
    #             
    #             Sigma_star      <- tau2[j] * correlation_function(D, theta_star, corr_fun = corr_fun)
    #             ## add in faster parallel cholesky as needed
    #             Sigma_chol_star <- tryCatch(
    #                 chol(Sigma_star),
    #                 error = function(e) {
    #                     if (verbose)
    #                         message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #                     num_chol_failures <- num_chol_failures + 1
    #                     chol(Sigma_star + 1e-8 * diag(N))                    
    #                 }
    #             )
    #             Sigma_inv_star  <- chol2inv(Sigma_chol_star)
    #             ## parallelize this
    #             mh1 <- dmvn(eta[, j], Xbeta[, j], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) +
    #                 ## prior
    #                 dmvn(theta_star, theta_mean, theta_var, log = TRUE)
    #             ## parallelize this        
    #             mh2 <- dmvn(eta[, j], Xbeta[, j], Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores) +
    #                 ## prior
    #                 if (corr_fun == "matern") {
    #                     dmvn(theta[j, ], theta_mean, theta_var, log = TRUE)
    #                 } else if (corr_fun == "exponential") {
    #                     dnorm(theta[j], theta_mean, sqrt(theta_var), log = TRUE)
    #                 }
    #             
    #             mh <- exp(mh1 - mh2)
    #             if (mh > runif(1, 0, 1)) {
    #                 if (corr_fun == "matern") {
    #                     theta[j, ] <- theta_star
    #                 } else if (corr_fun == "exponential") {
    #                     theta[j] <- theta_star
    #                 }
    #                 Sigma[j, , ]      <- Sigma_star
    #                 Sigma_chol[j, , ] <- Sigma_chol_star
    #                 Sigma_inv[j, , ]  <- Sigma_inv_star 
    #                 if (k <= params$n_adapt) {
    #                     theta_accept_batch[j] <- theta_accept_batch[j] + 1 / 50
    #                 } else {
    #                     theta_accept[j] <- theta_accept[j] + 1 / params$n_mcmc
    #                 }
    #             }
    #         }
    #         ## adapt the tuning
    #         if (k <= params$n_adapt) {
    #             if (corr_fun == "matern") {
    #                 save_idx <- k %% 50
    #                 if ((k %% 50) == 0) {
    #                     save_idx <- 50
    #                 } 
    #                 theta_batch[save_idx, , ] <- theta 
    #             }
    #             if (k %% 50 == 0) {
    #                 if (corr_fun == "matern") {
    #                     out_tuning <- update_tuning_mv_mat(
    #                         k,
    #                         theta_accept_batch,
    #                         lambda_theta,
    #                         theta_batch,
    #                         Sigma_theta_tune,
    #                         Sigma_theta_tune_chol
    #                     )
    #                     theta_batch           <- out_tuning$batch_samples
    #                     Sigma_theta_tune      <- out_tuning$Sigma_tune
    #                     Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
    #                     lambda_theta          <- out_tuning$lambda
    #                     theta_accept_batch    <- out_tuning$accept
    #                 } else if (corr_fun == "exponential") {
    #                     out_tuning <- update_tuning_vec(k, theta_accept_batch, theta_tune)
    #                     theta_tune         <- out_tuning$tune
    #                     theta_accept_batch <- out_tuning$accept
    #                 } 
    #             }   
    #         }        
    #     }        
    #     
    #     ##
    #     ## sample spatial process variance tau2
    #     ##
    #     
    #     if (verbose)
    #         message("sample tau2")
    #     
    #     if (shared_covariance_params) {
    #         devs       <- eta - Xbeta
    #         SS         <- sum(devs * (tau2 * Sigma_inv %*% devs))
    #         tau2       <- 1 / rgamma(1, N * (J-1) / 2 + priors$alpha_tau, SS / 2 + priors$beta_tau) 
    #         Sigma      <- tau2 * correlation_function(D, theta, corr_fun = corr_fun) 
    #         ## add in faster parallel cholesky as needed
    #         ## see https://github.com/RfastOfficial/Rfast/blob/master/src/cholesky.cpp
    #         Sigma_chol <- tryCatch(
    #             chol(Sigma),
    #             error = function(e) {
    #                 if (verbose)
    #                     message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #                 num_chol_failures <- num_chol_failures + 1
    #                 chol(Sigma + 1e-8 * diag(N))                    
    #             }
    #         )
    #         Sigma_inv  <- chol2inv(Sigma_chol) 
    #     } else {
    #         for (j in 1:(J-1)) {
    #             devs       <- eta[, j] - Xbeta[, j]
    #             SS         <- sum(devs * (tau2[j] * Sigma_inv[j, , ] %*% devs))
    #             tau2[j]    <- 1 / rgamma(1, N / 2 + priors$alpha_tau, SS / 2 + priors$beta_tau) 
    #             if (corr_fun == "matern") {
    #                 Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = corr_fun) 
    #             } else {
    #                 Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j], corr_fun = corr_fun) 
    #             }
    #             ## add in faster parallel cholesky as needed
    #             ## see https://github.com/RfastOfficial/Rfast/blob/master/src/cholesky.cpp
    #             Sigma_chol[j, , ] <- tryCatch(
    #                 chol(Sigma[j, , ]),
    #                 error = function(e) {
    #                     if (verbose)
    #                         message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #                     num_chol_failures <- num_chol_failures + 1
    #                     chol(Sigma[j, , ] + 1e-8 * diag(N))                    
    #                 }
    #             )
    #             Sigma_inv[j, , ]  <- chol2inv(Sigma_chol[j, , ])
    #         }
    #     }
    #     
    #     ##
    #     ## sample eta
    #     ##
    #     
    #     if (verbose)
    #         message("sample eta")
    #     
    #     if (shared_covariance_params) {
    #         ## double check this and add in fixed effects X %*% beta
    #         for (j in 1:(J-1)) {
    #             ## can make this much more efficient
    #             ## can this be parallelized? seems like it
    #             A        <- Sigma_inv + Omega[[j]]
    #             b        <- Sigma_inv %*% Xbeta[, j] + kappa[, j]
    #             eta[, j] <- rmvn_arma(A, b)
    #             # A <- Sigma_inv + Omega[[j]]
    #             # A_chol <- tryCatch(
    #             #     chol(A),
    #             #     error = function(e) {
    #             #         if (verbose)
    #             #             message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #             #         num_chol_failures <- num_chol_failures + 1
    #             #         chol(A + 1e-8 * diag(N))                    
    #             #     }
    #             # )
    #             # A_inv <- chol2inv(A_chol)
    #             # # A_inv <- chol2inv(chol(Sigma_inv + Omega[[j]]))
    #             # b     <- Sigma_inv %*% Xbeta[, j] + kappa[, j]
    #             # eta[, j]   <- rmvn(1, A_inv %*% b, A_inv)                
    #         }
    #     } else {      
    #         for (j in 1:(J-1)) {
    #             ## can make this much more efficient
    #             ## can this be parallelized? seems like it
    #             
    #             A        <- Sigma_inv[j, , ] + Omega[[j]]
    #             b        <- Sigma_inv[j, , ] %*% Xbeta[, j] + kappa[, j]
    #             eta[, j] <- rmvn_arma(A, b)
    #             # A <- Sigma_inv[j, , ] + Omega[[j]]
    #             # A_chol <- tryCatch(
    #             #     chol(A),
    #             #     error = function(e) {
    #             #         if (verbose)
    #             #             message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
    #             #         num_chol_failures <- num_chol_failures + 1
    #             #         chol(A + 1e-8 * diag(N))                    
    #             #     }
    #             # )
    #             # A_inv <- chol2inv(A_chol)
    #             # b     <- Sigma_inv[j, , ] %*% Xbeta[, j] + kappa[, j]
    #             # eta[, j]   <- rmvn(1, A_inv %*% b, A_inv)
    #             
    #         }
    #     }
    #     
    #     ##
    #     ## save variables
    #     ##
    #     if (k >= params$n_adapt) {
    #         if (k %% params$n_thin == 0) {
    #             save_idx                <- (k - params$n_adapt) / params$n_thin
    #             beta_save[save_idx, , ] <- beta
    #             if (shared_covariance_params) {
    #                 if (corr_fun == "matern") {
    #                     theta_save[save_idx, ]  <- theta
    #                 } else if (corr_fun == "exponential") {
    #                     theta_save[save_idx]  <- theta
    #                 }
    #                 tau2_save[save_idx]     <- tau2
    #             } else {
    #                 if (corr_fun == "matern") {
    #                     theta_save[save_idx, , ]  <- theta
    #                 } else if (corr_fun == "exponential") {
    #                     theta_save[save_idx, ]  <- theta
    #                 }
    #                 tau2_save[save_idx, ]     <- tau2
    #             }
    #             eta_save[save_idx, , ]  <- eta
    #         }
    #         
    #     }
    #     
    #     ##
    #     ## End of MCMC loop
    #     ##
    #     
    #     if (k %in% percentage_points && progress) {
    #         utils::setTxtProgressBar(progressBar, k / (params$n_adapt + params$n_mcmc))
    #     }
    # }
    # 
    # ## print out acceptance rates -- no tuning in this model
    # 
    # if (num_chol_failures > 0)
    #     warning("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")
    # 
    # ## eventually create a model class and include this as a variable in the class
    # message("Acceptance rate for theta is ", mean(theta_accept))
    # 
    # if (progress) {
    #     close(progressBar)
    # }
    # 
    # ##
    # ## return the MCMC output -- think about a better way to make this a class
    # ## 
    # 
    # out <- list(
    #     beta  = beta_save,
    #     theta = theta_save,
    #     tau2  = tau2_save,
    #     eta   = eta_save
    # )
    # class(out) <- "pgSPLM"
    # 
    # return(out)
    
}


