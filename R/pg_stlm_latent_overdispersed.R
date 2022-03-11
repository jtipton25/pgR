#' Bayesian Polya-gamma regression
#' 
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param Y is a \eqn{N \times J \times T}{N x J x T} array of compositional count data.
#' @param X is a \eqn{N \times p}{n_sites x p} matrix of climate variables.
#' @param locs is a \eqn{n_sites \times 2}{n_sites x 2} matrix of observation locations.
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
#' @param priors is a list of prior settings. 
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param shared_covariance_params is a logicial input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @param inits is the list of intial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @importFrom stats rmultinom
#' @importFrom hms as_hms
#' @importFrom BayesLogit rpg
#' @export

## polya-gamma spatial linear regression model
pg_stlm_latent_overdispersed <- function(
    Y, 
    X,
    locs, 
    params,
    priors,
    corr_fun                 = "exponential",
    n_cores                  = 1L,
    shared_covariance_params = TRUE,
    inits                    = NULL,
    config                   = NULL,
    n_chain                  = 1,
    progress                 = FALSE,
    verbose                  = FALSE
) {
    
    start <- Sys.time()
    
    ##
    ## Run error checks
    ## 
    
    check_input_pg_stlm(Y, X, locs)
    check_params(params)
    check_corr_fun(corr_fun)
    
    if (!is_positive_integer(n_cores, 1))
        stop("n_cores must be a positive integer")
    
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    ## add in faster parallel cholesky as needed
    
    ## 
    ## setup config
    ##
    
    ## do we sample the functional relationship parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
        }
    }
    
    ## do we sample the climate autocorrelation parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_rho <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_rho']])) {
            sample_rho <- config[['sample_rho']]
        }
    }
    
    ## do we sample the spatial variance parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
        }
    }
    
    
    ## do we sample the ovdiserpsion variance parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2']])) {
            sample_sigma2 <- config[['sample_sigma2']]
        }
    }
    
    sample_psi <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_psi']])) {
            sample_psi <- config[['sample_psi']]
        }
    }
    
    # fix these later
    priors$alpha_sigma <- 0.1
    priors$beta_sigma <- 0.1
    
    ## do we sample the climate spatial covariance parameters? This is
    ## primarily used to troubleshoot model fitting using simulated data
    sample_theta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_theta']])) {
            sample_theta <- config[['sample_theta']]
        }
    }
    
    ## do we sample the latent intensity parameter eta
    sample_eta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_eta']])) {
            sample_eta <- config[['sample_eta']]
        }
    }
    
    
    ##
    ## setup constats
    ##
    
    N      <- nrow(Y)
    J      <- ncol(Y)
    n_time <- dim(Y)[3]
    p      <- ncol(X)
    D      <- fields::rdist(locs)
    I      <- diag(N)
    
    ## Add in a counter for the number of regularized Cholesky factors.
    ## This is useful in correcting for numerical errors resulting in 
    ## covariance matrices that are not full rank
    num_chol_failures <- 0

    
    ## We assume a partially missing observation is the same as 
    ## fully missing. The index allows for fast accessing of missing
    ## observations
    missing_idx <- matrix(FALSE, N, n_time)
    for (i in 1:N) {
        for (tt in 1:n_time) {
            missing_idx[i, tt] <- any(is.na(Y[i, , tt]))
        }
    }
    
    message("There are ", ifelse(any(missing_idx), sum(missing_idx), "no"), " observations with missing count vectors")
    
    ## Calculate Mi and kappa
    Mi    <- array(0, dim = c(N, J - 1, n_time))
    kappa <- array(0, dim = c(N, J - 1, n_time))
    
    for (tt in 1:n_time) {
        Mi[, , tt]    <- calc_Mi(Y[, , tt])
        kappa[, , tt] <- calc_kappa(Y[, , tt], Mi[, , tt])
    }
    
    ## create an index for nonzero values
    nonzero_idx <- Mi != 0
    n_nonzero   <- sum(nonzero_idx)
    
    ##
    ## initial values
    ##
    
    ## default priors
    
    mu_beta        <- rep(0, p)
    Sigma_beta     <- 10 * diag(p)
    
    ## check if priors for mu_beta are specified
    if (!is.null(priors[['mu_beta']])) {
        if (all(!is.na(priors[['mu_beta']]))) {
            mu_beta <- priors[['mu_beta']]
        }
    }
    
    ## check if priors for Sigma_beta are specified
    if (!is.null(priors[['Sigma_beta']])) {
        if (all(!is.na(priors[['Sigma_beta']]))) {
            Sigma_beta <- priors[['Sigma_beta']]
        }
    }
    Sigma_beta_chol <- tryCatch(
        chol(Sigma_beta),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
            chol(Sigma_beta + 1e-8 * diag(N))                    
        }
    )
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##
    
    beta <- t(mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    ## check if initial value for beta is given
    if (!is.null(inits[['beta']])) {
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
        }
    }
    Xbeta <- X %*% beta
    
    ##
    ## initialize sigma2
    ##
    
    alpha_sigma2 <- 1
    beta_sigma2  <- 1
    
    if (shared_covariance_params) {
        sigma2 <- pmin(rgamma(1, alpha_sigma2, beta_sigma2), 5) 
    } else {
        sigma2 <- pmin(rgamma(J-1, alpha_sigma2, beta_sigma2), 5) 
    }
    
    ##
    ## initialize spatial Gaussian process -- share parameters across the different components
    ##    can generalize to each component getting its own covariance
    ##
    ## assume the GP parameters don't change through time
    
    theta_mean <- NULL
    theta_var  <- NULL
    if (corr_fun == "matern") {
        theta_mean <- c(priors$mean_range, priors$mean_nu)
        theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    } else if (corr_fun == "exponential") {
        theta_mean <- priors$mean_range
        theta_var  <- priors$sd_range^2
    }
    
    theta <- NULL
    if (shared_covariance_params) {
        if (corr_fun == "matern") {
            theta <- as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -2), 0.1))            
        } else if (corr_fun == "exponential") {
            theta <- pmin(pmax(rnorm(1, theta_mean, sqrt(theta_var)), -2), 0.1)
        }
    } else {
        if (corr_fun == "matern") {
            theta <- pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -2), 0.1)
        } else if (corr_fun == "exponential") {
            theta <- pmin(pmax(rnorm(J-1, theta_mean, sqrt(theta_var)), -2), 0.1)
        }
        
    }
    
    if (!is.null(inits[['theta']])) {
        if (all(!is.na(inits[['theta']]))) {
            theta <- inits[['theta']]
        }
    }
    
    ## check dimensions of theta
    if (shared_covariance_params) {
        if (corr_fun == "matern")
            if (!is_numeric_vector(theta, 2)) 
                stop('If shared_covariance_params is TRUE, theta must be a numeric vector of length 2 when corr_fun is "matern"')
        if (corr_fun == "exponential")
            if (!is_numeric_vector(theta, 1)) 
                stop('If shared_covariance_params is TRUE, theta must be a numeric of length 1 when corr_fun is "exponential"')
    } else {
        if (corr_fun == "matern")
            if (!is_numeric_matrix(theta, J-1, 2))
                stop('If shared_covariance_params is FALSE, theta must be a J-1 by 2 numeric matrix when corr_fun is "matern"')
        if (corr_fun == "exponential")
            if (!is_numeric_vector(theta, J-1)) 
                stop('If shared_covariance_params is FALSE, theta must be a numeric vector of length J-1 when corr_fun is "exponential"')
    }
    
    tau2 <- NULL
    if (shared_covariance_params) {
        tau2 <- min(1 / rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
    } else {
        tau2 <- pmin(1 / rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
    }
    if (!is.null(inits[['tau2']])) {
        if (all(!is.na(inits[['tau2']]))) {
            ## if tau2 passes error checks
            tau2 <- inits[['tau2']]
        }
    }
    
    ## check dimensions of tau2
    if (shared_covariance_params) {
        if (!is_positive_numeric(tau2, 1))
            stop("If shared_covariance_params is FALSE, tau2 must be a numeric scalar")
    } else {
        if (!is_positive_numeric(tau2, J-1))
            stop("If shared_covariance_params is TRUE, tau2 must be a J-1 positive numeric vector ")
    }
    
    
    Sigma <- NULL
    if (shared_covariance_params) {
        Sigma <- tau2 * correlation_function(D, theta, corr_fun = corr_fun)
    } else {
        Sigma <- array(0, dim = c(J - 1, N, N))
        for (j in 1:(J-1)) {
            if (corr_fun == "matern") {
                Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = corr_fun)
            } else if (corr_fun == "exponential") {
                Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j], corr_fun = corr_fun)
            }
        }
    }
    
    Sigma_chol <- NULL
    if (shared_covariance_params) {
        Sigma_chol <- tryCatch(
            chol(Sigma),
            error = function(e) {
                if (verbose)
                    message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                num_chol_failures <- num_chol_failures + 1
                chol(Sigma + 1e-8 * diag(N))                    
            }
        )
    } else {
        Sigma_chol <- array(0, dim = c(J-1, N, N))
        for (j in 1:(J-1)) {
            ## add a warning for the Cholesky function
            Sigma_chol[j, , ] <- tryCatch(
                chol(Sigma[j, , ]),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                    num_chol_failures <- num_chol_failures + 1
                    chol(Sigma[j, , ] + 1e-8 * diag(N))                    
                }
            )
        }
    }
    
    Sigma_inv  <- NULL
    if (shared_covariance_params) {
        Sigma_inv <- chol2inv(Sigma_chol)   
    } else {
        Sigma_inv <- array(0, dim = c(J-1, N, N))
        for (j in 1:(J-1)) {
            Sigma_inv[j, , ] <- chol2inv(Sigma_chol[j, , ])
        }
    }
    
    ## temporal autocorrelation
    rho <- runif(1, 0, 1)
    if (!is.null(inits[['rho']])) {
        if (all(!is.na(inits[['rho']]))) {
            ## if rho passes error checks
            rho <- inits[['rho']]
        }
    }
    
    psi  <- array(0, dim = c(N, J-1, n_time))
    for (tt in 1:n_time) {
        if (tt == 1) {
            for (j in 1:(J-1)) {
                if (shared_covariance_params) {
                    psi[, j, 1] <- t(mvnfast::rmvn(1, rep(0, N), Sigma_chol, isChol = TRUE))
                } else{
                    psi[, j, 1] <- t(mvnfast::rmvn(1, rep(0, N), Sigma_chol[j, , ], isChol = TRUE))
                }
            }
        } else {
            for (j in 1:(J-1)) {
                if (shared_covariance_params) {
                    psi[, j, tt] <- mvnfast::rmvn(1, rho * psi[, j, tt - 1], Sigma_chol, isChol = TRUE)
                } else {
                    psi[, j, tt] <- t(mvnfast::rmvn(1, rho * psi[, j, tt - 1], Sigma_chol[j, , ], isChol = TRUE))
                }
            }
        }
    }    
    
    if (!is.null(inits[['psi']])) {
        if (all(!is.na(inits[['psi']]))) {
            psi <- inits[['psi']]
        }
    }
    
    eta  <- array(0, dim = c(N, J-1, n_time))
    for (tt in 1:n_time) {
        for (j in 1:(J-1)) {
            if (shared_covariance_params) {
                eta[, j, tt] <- rnorm(N, Xbeta[, j] + psi[, j, tt], sqrt(sigma2))
            } else {
                eta[, j, tt] <- rnorm(N, Xbeta[, j] + psi[, j, tt], sqrt(sigma2[j]))
            }
        }
    }

    if (!is.null(inits[['eta']])) {
        if (all(!is.na(inits[['eta']]))) {
            eta <- inits[['eta']]
        }
    }
    
    ##
    ## sampler config options -- to be added later
    ## 
    #
    
    ##
    ## initialize omega
    ##
    
    omega <- array(0, dim = c(N, J-1, n_time))
    # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
    
    save_omega <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['save_omega']])) {
            save_omega <- config[['save_omega']]
        }
    }
    
    ##
    ## setup save variables
    ##
    
    n_save     <- params$n_mcmc / params$n_thin
    beta_save  <- array(0, dim = c(n_save, p, J-1))
    tau2_save  <- NULL
    if (shared_covariance_params) {
        tau2_save <- rep(0, n_save)
    } else {
        tau2_save <- matrix(0, n_save, J-1)
    }
    sigma2_save  <- NULL
    if (shared_covariance_params) {
        sigma2_save <- rep(0, n_save)
    } else {
        sigma2_save <- matrix(0, n_save, J-1)
    }
    theta_save <- NULL
    if (shared_covariance_params) {
        if (corr_fun == "matern") {
            theta_save <- matrix(0, n_save, 2)
        } else if (corr_fun == "exponential") {
            theta_save <- rep(0, n_save)
        }
    } else {
        if (corr_fun == "matern") {
            theta_save <- array(0, dim = c(n_save, J-1, 2))
        } else if (corr_fun == "exponential") {
            theta_save <- matrix(0, n_save, J-1)
        }
    }
    eta_save   <- array(0, dim = c(n_save, N, J-1, n_time))
    psi_save   <- array(0, dim = c(n_save, N, J-1, n_time))
    pi_save    <- array(0, dim = c(n_save, N, J, n_time))
    rho_save   <- rep(0, n_save)
    omega_save <- NULL
    if (save_omega) {
        omega_save   <- array(0, dim = c(n_save, N, J-1, n_time))
    }
    
    ## 
    ## initialize tuning 
    ##
    
    ##
    ## tuning variables for adaptive MCMC
    ##
    
    theta_batch           <- NULL
    theta_accept          <- NULL
    theta_accept_batch    <- NULL
    lambda_theta          <- NULL
    Sigma_theta_tune      <- NULL
    Sigma_theta_tune_chol <- NULL
    theta_tune            <- NULL
    
    if (shared_covariance_params) {
        
        theta_accept       <- 0
        theta_accept_batch <- 0
        
        if (corr_fun == "matern") {
            theta_batch <- matrix(0, 50, 2) 
            lambda_theta          <- 0.05
            Sigma_theta_tune      <- 1.8 * diag(2) - .8
            Sigma_theta_tune_chol <- tryCatch(
                chol(Sigma_theta_tune),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the Metroplois-Hastings adaptive tuning matrix for Matern parameters theta was ill-conditioned and mildy regularized.")
                    chol(Sigma_theta_tune + 1e-8 * diag(2))                    
                }
            )
        } else if (corr_fun == "exponential") {
            theta_tune <- mean(D) / 2
        }
        
    } else {
        
        theta_accept       <- rep(0, J-1)
        theta_accept_batch <- rep(0, J-1)
        
        if (corr_fun == "matern") {
            theta_batch <- array(0, dim = c(50, 2, J-1))      
            lambda_theta     <- rep(0.05, J-1)
            Sigma_theta_tune <- array(0, dim = c(2, 2, J-1))
            for (j in 1:(J-1)) {
                Sigma_theta_tune[, , j] <- 1.8 * diag(2) - .8
            }
            Sigma_theta_tune_chol <- array(0, dim = c(2, 2, J-1))
            for (j in 1:(J-1)) {
                Sigma_theta_tune_chol[, , j] <- tryCatch(
                    chol(Sigma_theta_tune[, , j]),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the Metroplois-Hastings adaptive tuning matrix for Matern parameters theta was ill-conditioned and mildy regularized.")
                        chol(Sigma_theta_tune[, , j] + 1e-8 * diag(2))                    
                    }
                )
            }
        } else if (corr_fun == "exponential") {
            theta_tune <- rep(mean(D) / 2, J-1)
        }
    }
    
    # tuning for tau2
    tau2_accept          <- NULL
    tau2_accept_batch    <- NULL
    tau2_tune            <- NULL
    
    if (shared_covariance_params) {
        tau2_accept       <- 0
        tau2_accept_batch <- 0
        tau2_tune         <- 0.5
    } else {
        tau2_accept       <- rep(0, J-1)
        tau2_accept_batch <- rep(0, J-1)
        tau2_tune         <- rep(0.5, J-1)
    }

    # tuning for sigma2   
    sigma2_accept          <- NULL
    sigma2_accept_batch    <- NULL
    sigma2_tune            <- NULL
    
    if (shared_covariance_params) {
        sigma2_accept       <- 0
        sigma2_accept_batch <- 0
        sigma2_tune         <- 0.5
    } else {
        sigma2_accept       <- rep(0, J-1)
        sigma2_accept_batch <- rep(0, J-1)
        sigma2_tune         <- rep(0.5, J-1)
    }
    
    
    ## tuning for rho
    rho_accept       <- 0
    rho_accept_batch <- 0
    rho_tune         <- 0.025
    
    ##
    ## Starting MCMC chain ----
    ##
    
    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")
    if (progress) {
        progressBar <- utils::txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * (params$n_adapt + params$n_mcmc))
    
    
    for (k in 1:(params$n_adapt + params$n_mcmc)) {
        if (k == params$n_adapt + 1) {
            message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations \n")
        }
        if (k %% params$n_message == 0) {
            if (k <= params$n_adapt) {
                message("MCMC adaptation iteration ", k, " for chain ", n_chain)
            } else {
                message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
            }
        }
        
        ##
        ## sample Omega ----
        ##
        
        if (verbose)
            message("sample omega")
        
        # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
        
        ##
        ## sample beta ----
        ##
        
        ## can parallelize this update -- each group of parameters is 
        ## conditionally independent given omega and kappa(y)
        
        if (sample_beta){
            if (verbose)
                message("sample beta")
            
            if (shared_covariance_params) {
                for (j in 1:(J-1)) {
                    A <- n_time * t(X) %*% X / sigma2 + Sigma_beta_inv
                    ## guarantee a symmetric matrix
                    A <- (A + t(A)) / 2
                    b <- rowSums(t(X) %*% (eta[, j, ] - psi[, j, ])) / sigma2 + Sigma_beta_inv %*% mu_beta
                    beta[, j]   <- rmvn_arma(A, b)
                }
            } else {
                for (j in 1:(J-1)) {
                    A <- n_time * t(X) %*% X / sigma2[j] + Sigma_beta_inv
                    ## guarantee a symmetric matrix
                    A <- (A + t(A)) / 2
                    b <- rowSums(t(X) %*% (eta[, j, ] - psi[, j, ])) / sigma2[j] + Sigma_beta_inv %*% mu_beta
                    beta[, j]   <- rmvn_arma(A, b)
                }
            }        
            Xbeta <- X %*% beta
        }
        
        ##
        ## sample spatial correlation parameters theta ----
        ##
        
        if(sample_theta) {      
            if (verbose)
                message("sample theta")
            
            if (shared_covariance_params) {
                ## update a common theta for all processes
                theta_star <- NULL
                if (corr_fun == "matern") {
                    theta_star <- as.vector(
                        mvnfast::rmvn( 
                            n      = 1,
                            mu     = theta,
                            sigma  = lambda_theta * Sigma_theta_tune_chol,
                            isChol = TRUE
                        )
                    )
                } else if (corr_fun == "exponential") {
                    theta_star <- rnorm(1, theta, theta_tune)
                }
                Sigma_star       <- tau2 * correlation_function(D, theta_star, corr_fun = corr_fun)
                ## add in faster parallel cholesky as needed
                Sigma_chol_star <- tryCatch(
                    chol(Sigma_star),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                        num_chol_failures <- num_chol_failures + 1
                        chol(Sigma_star + 1e-8 * diag(N))                    
                    }
                )
                Sigma_inv_star  <- chol2inv(Sigma_chol_star)
                ## parallelize this
                mh1 <- sum(
                    sapply(1:(J-1), function(j) {
                        mvnfast::dmvn(psi[, j, 1], rep(0, N), Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                    })
                ) +
                    sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(psi[, j, tt], rho * psi[, j, tt - 1], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                            })
                        }) 
                        
                    ) +
                    ## prior
                    mvnfast::dmvn(theta_star, theta_mean, theta_var, log = TRUE)
                ## parallelize this        
                mh2 <- sum(
                    sapply(1:(J-1), function(j) {
                        mvnfast::dmvn(psi[, j, 1], rep(0, N), Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
                    })
                ) +
                    sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(psi[, j, tt], rho * psi[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)
                            })
                        })
                    ) +
                    ## prior
                    mvnfast::dmvn(theta, theta_mean, theta_var, log = TRUE)
                
                mh <- exp(mh1 - mh2)
                if (mh > runif(1, 0, 1)) {
                    theta      <- theta_star
                    Sigma      <- Sigma_star
                    Sigma_chol <- Sigma_chol_star
                    Sigma_inv  <- Sigma_inv_star 
                    if (k <= params$n_adapt) {
                        theta_accept_batch <- theta_accept_batch + 1 / 50
                    } else {
                        theta_accept <- theta_accept + 1 / params$n_mcmc
                    }
                }
                ## adapt the tuning
                if (k <= params$n_adapt) {
                    if (corr_fun == "matern") {
                        save_idx <- k %% 50
                        if ((k %% 50) == 0) {
                            save_idx <- 50
                        } 
                        theta_batch[save_idx, ] <- theta 
                    }
                    if (k %% 50 == 0) {
                        if (corr_fun == "matern") {
                            out_tuning <- update_tuning_mv(
                                k,
                                theta_accept_batch,
                                lambda_theta,
                                theta_batch,
                                Sigma_theta_tune,
                                Sigma_theta_tune_chol
                            )
                            theta_batch           <- out_tuning$batch_samples
                            Sigma_theta_tune      <- out_tuning$Sigma_tune
                            Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
                            lambda_theta          <- out_tuning$lambda
                            theta_accept_batch    <- out_tuning$accept
                        } else if (corr_fun == "exponential") {
                            out_tuning         <- update_tuning(k, theta_accept_batch, theta_tune)
                            theta_tune         <- out_tuning$tune
                            theta_accept_batch <- out_tuning$accept
                        }
                    }
                }
            } else {
                ## 
                ## theta varies for each component
                ##
                for (j in 1:(J-1)) {
                    theta_star <- NULL
                    if (corr_fun == "matern") {
                        theta_star <- as.vector(
                            mvnfast::rmvn( 
                                n      = 1,
                                mu     = theta[j, ],
                                sigma  = lambda_theta[j] * Sigma_theta_tune_chol[, , j],
                                isChol = TRUE
                            )
                        )
                    } else if (corr_fun == "exponential") {
                        theta_star <- rnorm(1, theta[j], theta_tune)
                    }
                    
                    Sigma_star      <- tau2[j] * correlation_function(D, theta_star, corr_fun = corr_fun)
                    ## add in faster parallel cholesky as needed
                    Sigma_chol_star <- tryCatch(
                        chol(Sigma_star),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                            num_chol_failures <- num_chol_failures + 1
                            chol(Sigma_star + 1e-8 * diag(N))                    
                        }
                    )
                    Sigma_inv_star  <- chol2inv(Sigma_chol_star)
                    
                    ## parallelize this
                    mh1 <- mvnfast::dmvn(psi[, j, 1], rep(0, N), Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) +
                      sum(
                        sapply(2:n_time, function(tt) {
                            mvnfast::dmvn(psi[, j, tt], rho * psi[, j, tt - 1], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                        })) +
                      ## prior
                      mvnfast::dmvn(theta_star, theta_mean, theta_var, log = TRUE)

                    ## parallelize this        
                    theta_mh <- NULL
                    if (corr_fun == "matern") {
                        theta_mh <- theta[j, ]
                    } else if (corr_fun == "exponential") {
                        theta_mh <- theta[j]                        
                    }
                    mh2 <- mvnfast::dmvn(psi[, j, 1], rep(0, N), Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores) +
                        sum(
                            sapply(2:n_time, function(tt) {
                                mvnfast::dmvn(psi[, j, tt], rho * psi[, j, tt - 1], Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores)
                        })) +
                        ## prior
                        mvnfast::dmvn(theta_mh, theta_mean, theta_var, log = TRUE)
                    # mvnfast::dmvn(theta[j, , drop = FALSE], theta_mean, theta_var, log = TRUE)

                    mh <- exp(mh1 - mh2)
                    if (mh > runif(1, 0, 1)) {
                        if (corr_fun == "matern") {
                            theta[j, ] <- theta_star
                        } else if (corr_fun == "exponential") {
                            theta[j] <- theta_star
                        }
                        Sigma[j, , ]      <- Sigma_star
                        Sigma_chol[j, , ] <- Sigma_chol_star
                        Sigma_inv[j, , ]  <- Sigma_inv_star 
                        if (k <= params$n_adapt) {
                            theta_accept_batch[j] <- theta_accept_batch[j] + 1 / 50
                        } else {
                            theta_accept[j] <- theta_accept[j] + 1 / params$n_mcmc
                        }
                    }
                }
                ## adapt the tuning
                if (k <= params$n_adapt) {
                    if (corr_fun == "matern") {
                        save_idx <- k %% 50
                        if ((k %% 50) == 0) {
                            save_idx <- 50
                        } 
                        theta_batch[save_idx, , ] <- theta 
                    }
                    if (k %% 50 == 0) {
                        if (corr_fun == "matern") {
                            out_tuning <- update_tuning_mv_mat(
                                k,
                                theta_accept_batch,
                                lambda_theta,
                                theta_batch,
                                Sigma_theta_tune,
                                Sigma_theta_tune_chol
                            )
                            theta_batch           <- out_tuning$batch_samples
                            Sigma_theta_tune      <- out_tuning$Sigma_tune
                            Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
                            lambda_theta          <- out_tuning$lambda
                            theta_accept_batch    <- out_tuning$accept
                        } else if (corr_fun == "exponential") {
                            out_tuning         <- update_tuning_vec(k, theta_accept_batch, theta_tune)
                            theta_tune         <- out_tuning$tune
                            theta_accept_batch <- out_tuning$accept
                        }
                    }   
                }        
            }
        }        
        
        ##
        ## sample spatial process variance tau2 ----
        ##
        
        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            if (shared_covariance_params) {
                # ## update a common tau2 for all processes
                SS <- rep(0, J-1)
                for (j in 1:(J-1)) {
                    devs <- cbind(psi[, j, 1],
                                  psi[, j, -1] - rho * psi[, j, -n_time])
                    SS[j] <- sum(sapply(1:n_time, function(tt) devs[, tt] %*% (tau2 * Sigma_inv %*% devs[, tt])))
                }
                tau2 <- 1 / rgamma(1, priors$alpha_tau + 0.5 * N * n_time * (J-1), priors$beta_tau + 0.5 * sum(SS))
                # update Sigma
                Sigma <- tau2 * correlation_function(D, theta, corr_fun = corr_fun)
                Sigma_chol <- tryCatch(
                    chol(Sigma),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                        num_chol_failures <- num_chol_failures + 1
                        chol(Sigma + 1e-8 * diag(N))                    
                    }
                )
                Sigma_inv  <- chol2inv(Sigma_chol)
            } else {
                for (j in 1:(J-1)) {
                    devs <- cbind(psi[, j, 1],
                                  psi[, j, -1] - rho * psi[, j, -n_time])
                    SS <- sum(sapply(1:n_time, function(tt) devs[, tt] %*% (tau2[j] * Sigma_inv[j, , ] %*% devs[, tt])))
                    tau2[j] <- 1 / rgamma(1, priors$alpha_tau + 0.5 * N * n_time, priors$beta_tau + 0.5 * SS)
                    # update Sigma
                    if (corr_fun == "matern") {
                        Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = corr_fun)
                    } else if (corr_fun == "exponential") {
                        Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j], corr_fun = corr_fun)
                    }
                    Sigma_chol[j, , ] <- tryCatch(
                        chol(Sigma[j, , ]),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                            num_chol_failures <- num_chol_failures + 1
                            chol(Sigma[j, , ] + 1e-8 * diag(N))                    
                        }
                    )
                    Sigma_inv[j, , ]  <- chol2inv(Sigma_chol[j, , ])
                }
            }
        }
        
        ##
        ## sample overdispersion variance sigma2 ----
        ##
        
        if (sample_sigma2) {
            if (verbose)
                message("sample sigma2")
            
            if (shared_covariance_params) {
                SS <- sum(sapply(1:n_time, function(tt) (eta[, , tt] - Xbeta - psi[, , tt])^2))
                sigma2 <- 1 / rgamma(1, alpha_sigma2 + N * n_time * (J-1) / 2, beta_sigma2 + SS / 2)
            } else {
                for (j in 1:(J-1)) {
                    SS <- sum(sapply(1:n_time, function(tt) (eta[, j, tt] - Xbeta[, j] - psi[, j, tt])^2))
                    sigma2[j] <- 1 / rgamma(1, alpha_sigma2 + N * n_time / 2, beta_sigma2 + SS / 2)
                }
            }
        }
        
        ##
        ## sample eta ----
        ##
        
        if (sample_eta) {
            if (verbose)
                message("sample eta")
            
            ## need to add in the autocorrelation process (how does this affect the kappas?)
            for (tt in 1:n_time) {
                ## double check this full conditional
                if (shared_covariance_params) {
                    eta[, , tt] <- sapply(1:(J-1), function(j) {
                        sigma2_tilde <- 1 / (1 / sigma2 + omega[, j, tt])
                        mu_tilde     <- 1 / sigma2 * (Xbeta[, j] + psi[, j, tt]) + kappa[, j, tt]
                        return(
                            rnorm(
                                N, 
                                sigma2_tilde * mu_tilde,
                                sqrt(sigma2_tilde)
                            )
                        )
                    })
                } else {
                    eta[, , tt] <- sapply(1:(J-1), function(j) {
                        sigma2_tilde <- 1 / (1 / sigma2[j] + omega[, j, tt])
                        mu_tilde     <- 1 / sigma2[j] * (Xbeta[, j] + psi[, j, tt]) + kappa[, j, tt]
                        return(
                            rnorm(
                                N, 
                                sigma2_tilde * mu_tilde,
                                sqrt(sigma2_tilde)
                            )
                        )
                    })
                }
            }
        }
        
        ##
        ## sample psi ---- 
        ##
        
        if (sample_psi) {
            if (verbose)
                message("sample psi")
            
            ## need to add in the autocorrelation process (how does this affect the kappas?)
            for (tt in 1:n_time) {
                for (j in 1:(J-1)) {
                    A <- NULL
                    b <- NULL
                    if (tt == 1) {
                        if (shared_covariance_params) {
                            A <- (1 + rho^2) * Sigma_inv + 1 / sigma2 * I
                            b <- Sigma_inv %*% (rho * psi[, j, tt + 1]) + 1 / sigma2 * (eta[, j, tt] - Xbeta[, j])
                        } else {
                            A <- (1 + rho^2) * Sigma_inv[j, , ] + 1 / sigma2[j] * I
                            b <- Sigma_inv[j, , ] %*% (rho * psi[, j, tt + 1]) + 1 / sigma2[j] * (eta[, j, tt] - Xbeta[, j])
                        }
                    } else if (tt == n_time) {
                        if (shared_covariance_params) {
                            A <- Sigma_inv + 1 / sigma2 * I
                            b <- Sigma_inv %*% (rho * psi[, j, tt - 1]) + 1 / sigma2 * (eta[, j, tt] - Xbeta[, j])
                        } else {
                            A <- Sigma_inv[j, , ] + 1 / sigma2[j] * I
                            b <- Sigma_inv[j, , ] %*% (rho * psi[, j, tt - 1]) + 1 / sigma2[j] * (eta[, j, tt] - Xbeta[, j])
                        }
                    } else {
                        if (shared_covariance_params) {
                            A <- (1 + rho^2) * Sigma_inv + 1 / sigma2 * I
                            b <- Sigma_inv %*% (rho * (psi[, j, tt - 1] + psi[, j, tt + 1])) + 1 / sigma2 * (eta[, j, tt] - Xbeta[, j])
                        } else {
                            A <- (1 + rho^2) * Sigma_inv[j, , ] + 1 / sigma2[j] * I
                            b <- Sigma_inv[j, , ] %*% (rho * (psi[, j, tt - 1] +  psi[, j, tt + 1])) + 1 / sigma2[j] * (eta[, j, tt] - Xbeta[, j])
                        }
                    }
                    ## guarantee that A is symmetric
                    A            <- (A + t(A)) / 2
                    psi[, j, tt] <- rmvn_arma(A, b)
                }
            }
        }
        
        ##
        ## sample rho
        ##
        
        if (sample_rho) {
            if (verbose)
                message("sample rho")
            
            rho_star <- rnorm(1, rho, rho_tune)
            if (rho_star < 1 & rho_star > -1) {
                mh1 <- NULL
                mh2 <- NULL
                if (shared_covariance_params){
                    mh1 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(psi[, j, tt], rho_star * psi[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)
                            })
                        })

                    )
                    ## parallelize this
                    mh2 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(psi[, j, tt], rho * psi[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)
                            })
                        })
                    )
                } else {
                    mh1 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(psi[, j, tt], rho_star * psi[, j, tt - 1], Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores)
                            })
                        })

                    )
                    ## parallelize this
                    mh2 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(psi[, j, tt], rho * psi[, j, tt - 1], Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores)
                            })
                        })
                    )
                }

                mh <- exp(mh1 - mh2)
                if (length(mh) > 1)
                    stop("error in mh for rho")
                if (mh > runif(1, 0.0, 1.0)) {
                    rho   <- rho_star
                    if (k <= params$n_adapt) {
                        rho_accept_batch <- rho_accept_batch + 1.0 / 50.0
                    } else {
                        rho_accept <- rho_accept + 1.0 / params$n_mcmc
                    }
                }

                ## update tuning
                if (k <= params$n_adapt) {
                    if (k %% 50 == 0){
                        out_tuning <- update_tuning(
                            k,
                            rho_accept_batch,
                            rho_tune
                        )
                        rho_tune         <- out_tuning$tune
                        rho_accept_batch <- out_tuning$accept
                    }
                }
            }
        }
        
        ##
        ## save variables ----
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ] <- beta
                if (shared_covariance_params) {
                    if (corr_fun == "matern") {
                        theta_save[save_idx, ]  <- theta
                    } else if (corr_fun == "exponential") {
                        theta_save[save_idx]  <- theta
                    }
                    tau2_save[save_idx]     <- tau2
                    sigma2_save[save_idx]   <- sigma2
                } else {
                    if (corr_fun == "matern") {
                        theta_save[save_idx, , ]  <- theta
                    } else if (corr_fun == "exponential") {
                        theta_save[save_idx, ]  <- theta
                    }
                    tau2_save[save_idx, ]     <- tau2
                    sigma2_save[save_idx, ]   <- sigma2
                }
                eta_save[save_idx, , , ] <- eta
                psi_save[save_idx, , , ] <- psi
                for (tt in 1:n_time) {
                    pi_save[save_idx, , , tt]  <- eta_to_pi(eta[, , tt])
                }
                rho_save[save_idx]       <- rho
                if (save_omega) {
                    omega_save[save_idx, , , ] <- omega
                }
            }
            
        }
        
        ##
        ## End of MCMC loop
        ##
        
        if (k %in% percentage_points && progress) {
            utils::setTxtProgressBar(progressBar, k / (params$n_adapt + params$n_mcmc))
        }
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    if (num_chol_failures > 0)
        warning("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")
    
    ## eventually create a model class and include this as a variable in the class
    message("Acceptance rate for theta is ", mean(theta_accept))
    message("Acceptance rate for rho is ", mean(rho_accept))
    message("Acceptance rate for tau2 is ", mean(tau2_accept))
    message("Acceptance rate for sigma2 is ", mean(sigma2_accept))
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    if (progress) {
        close(progressBar)
    }
    
    stop    <- Sys.time()
    runtime <- stop - start
    
    message("MCMC took ", hms::as_hms(runtime))
    
    
    out <- NULL
    if (save_omega) {
        out <- list(
            beta   = beta_save,
            theta  = theta_save,
            tau2   = tau2_save,
            sigma2 = sigma2_save,
            eta    = eta_save,
            psi    = psi_save,
            pi     = pi_save,
            rho    = rho_save,
            omega  = omega_save)        
    } else {
        out <- list(
            beta   = beta_save,
            theta  = theta_save,
            tau2   = tau2_save,
            sigma2 = sigma2_save,
            eta    = eta_save,
            psi    = psi_save,
            pi     = pi_save,
            rho    = rho_save)
    }

    class(out) <- "pg_stlm_latent_overdispersed"
    
    return(out)
}


