#' Bayesian Polya-gamma regression
#' 
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param Y is a \eqn{n \times J}{n x J} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of fixed effects (like latitude, elevation, etc)
#' @param Z0 is a \eqn{n \times Q}{n x Q} matrix of observed climate variables
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
#' @param model is the form of the  polya-gamma model. Currently, this option is not active the only model is the "iid error" model. This option allows for independent species-specific overdispersion variance terms.
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#'
#' @export
#' 
#' @importFrom LaplacesDemon rinvwishart rtrunc
#' @importFrom stats rmultinom
#' @importFrom hms as_hms
#' @importFrom fields rdist
#' @importFrom BayesLogit rpg

## polya-gamma spatial linear regression model
pg_mvgp_univariate <- function(
    Y, 
    X,
    Z0,
    locs, 
    params,
    priors,
    corr_fun = "exponential",
    model    = "iid error",
    n_cores  = 1L,
    inits    = NULL,
    config   = NULL,
    verbose  = FALSE,
    n_chain  = 1
) {
    

    ##
    ## Run error checks
    ## 
    
    if(!is.vector(Z0))
        stop("Z must be a vector of climate variable inputs")
    if (!is_positive_integer(n_cores, 1))
        stop("n_cores must be a positive integer")
    
    
    
    
    
    check_corr_fun(corr_fun)
    # check_input_spatial(Y, X, locs)
    # check_params(params)
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    ## 
    ## setup config
    ##
    
    ## do we sample the functional relationship parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_beta)) {
            sample_beta <- config$sample_beta
        }
    }
    
    ## do we only use the modern climate data to estimate the functional 
    ## relationship or do we use the estimated climate as well
    sample_beta_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_beta_modern)) {
            sample_beta_modern <- config$sample_beta_modern
        }
    }
    
    ## do we sample the climate fixed effect parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_gamma <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_gamma)) {
            sample_gamma <- config$sample_gamma
        }
    }
    

    ## do we only use the modern climate data to estimate the fixed effects 
    ## or do we use the estimated climate as well
    sample_gamma_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_gamma_modern)) {
            sample_gamma_modern <- config$sample_gamma_modern
        }
    }
    
    ## do we sample the climate autocorrelation parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_rho <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_rho)) {
            sample_rho <- config$sample_rho
        }
    }
    
    ## do we sample the climate variance parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_tau2)) {
            sample_tau2 <- config$sample_tau2
        }
    }
    
    ## do we sample the climate variance parameter using only the 
    ## modern data? 
    sample_tau2_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_tau2_modern)) {
            sample_tau2_modern <- config$sample_tau2_modern
        }
    }
    
    ## do we sample the climate spatial covariance parameters? This is
    ## primarily used to troubleshoot model fitting using simulated data
    sample_theta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_theta)) {
            sample_theta <- config$sample_theta
        }
    }
    
    ## do we only use the modern climate data to estimate the spatial 
    ## covariance parameters or do we use the estimated climate as well
    sample_theta_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_theta_modern)) {
            sample_theta_modern <- config$sample_theta_modern
        }
    }
    
    ## do we sample the climate latent random effects? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_Z <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_Z)) {
            sample_Z <- config$sample_Z
        }
    }
    

    ## do we sample the overdispersion parameter
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_sigma2)) {
            sample_sigma2 <- config$sample_sigma2
        }
    }
    
    ## do we sample the overdispersion parameter only for the modern data?
    sample_sigma2_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_sigma2_modern)) {
            sample_sigma2_modern <- config$sample_sigma2_modern
        }
    }
    
    ## do we sample the latent intensity parameter eta
    sample_eta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config$sample_eta)) {
            sample_eta <- config$sample_eta
        }
    }
    

    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0
    
    N      <- nrow(Y)
    J      <- ncol(Y)
    n_time <- dim(Y)[3]
    p      <- ncol(X)
    Q      <- 1
    D      <- rdist(locs)
    
    ## we assume a partially missing observation is the same as fully missing
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
    for (i in 1:N){
        for (tt in 1:n_time) {
            if (missing_idx[i, tt]) {
                Mi[i, , tt]    <- 0
                kappa[i, , tt] <- 0
            } else {
                Mi[i, , tt] <- sum(Y[i, , tt]) - c(0, cumsum(Y[i, , tt][1:(J - 2)]))
                kappa[i, , tt] <- Y[i, 1:(J - 1), tt] - Mi[i, , tt] / 2
            }
        }
    }
    
    ## create an index for nonzero values
    nonzero_idx <- Mi != 0
    n_nonzero   <- sum(nonzero_idx)
    
    ##
    ## initial values
    ##
    
    ##
    ## priors for beta
    ##
    
    mu_beta        <- rep(0, Q + 1)
    Sigma_beta     <- 10 * diag(Q + 1)
    
    ## check if priors for mu_beta are specified
    if (!is.null(priors$mu_beta)) {
        if (all(!is.na(priors$mu_beta))) {
            mu_beta <- priors$mu_beta
        }
    }
    
    ## check if priors for Sigma_beta are specified
    if (!is.null(priors$Sigma_beta)) {
        if (all(!is.na(priors$Sigma_beta))) {
            Sigma_beta <- priors$Sigma_beta
        }
    }
    Sigma_beta_chol <- tryCatch(
        chol(Sigma_beta),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
            chol(Sigma_beta + 1e-8 * diag(Q + 1))                    
        }
    )
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##
    
    beta <- t(mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    ## check if initial value for beta is given
    if (!is.null(inits$beta)) {
        if (all(!is.na(inits$beta))) {
            beta <- inits$beta
        }
    }
    
    ##
    ## priors for gamma
    ##
    
    mu_gamma        <- rep(0, p)
    Sigma_gamma     <- 10 * diag(p)
    
    ## check if priors for mu_gamma are specified
    if (!is.null(priors$mu_gamma)) {
        if (all(!is.na(priors$mu_gamma))) {
            mu_gamma <- priors$mu_gamma
        }
    }
    
    ## check if priors for Sigma_gamma are specified
    if (!is.null(priors$Sigma_gamma)) {
        if (all(!is.na(priors$Sigma_gamma))) {
            Sigma_gamma <- priors$Sigma_gamma
        }
    }
    Sigma_gamma_chol <- tryCatch(
        chol(Sigma_gamma),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the prior covariance Sigma_gamma was ill-conditioned and mildy regularized.")
            chol(Sigma_gamma + 1e-8 * diag(p))                    
        }
    )
    Sigma_gamma_inv  <- chol2inv(Sigma_gamma_chol)
    
    ##
    ## initialize gamma
    ##
    
    gamma <- t(mvnfast::rmvn(Q, mu_gamma, Sigma_gamma_chol, isChol = TRUE))
    ## check if initial value for gamma is given
    if (!is.null(inits$gamma)) {
        if (all(!is.na(inits$gamma))) {
            gamma <- inits$gamma
        }
    }
    
    Xgamma <- as.vector(X %*% gamma)
    
    ##
    ## initialize spatial Gaussian process -- share parameters across the different components
    ##    can generalize to each component getting its own covariance
    ##
    
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
    if (corr_fun == "matern") {
        theta <- as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -2), 0.1))            
    } else if (corr_fun == "exponential") {
        theta <- pmin(pmax(rnorm(1, theta_mean, sqrt(theta_var)), -2), 0.1)
    }
    
    if (!is.null(inits$theta)) {
        if (all(!is.na(inits$theta))) {
            theta <- inits$theta
        }
    } 
    
    ## check dimensions of theta
    if (corr_fun == "matern")
        if (!is_numeric_vector(theta, 2)) 
            stop('If shared_covariance_params is TRUE, theta must be a numeric vector of length 2 when corr_fun is "matern"')
    if (corr_fun == "exponential")
        if (!is_numeric_vector(theta, 1)) 
            stop('If shared_covariance_params is TRUE, theta must be a numeric of length 1 when corr_fun is "exponential"')
    
    # tau2       <- min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
    R      <- correlation_function(D, theta, corr_fun = corr_fun)
    ## add in faster parallel cholesky as needed
    R_chol <- tryCatch(
        chol(R),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the spatial correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
            num_chol_failures <- num_chol_failures + 1
            chol(R + 1e-8 * diag(N))                    
        }
    )
    R_inv  <- chol2inv(R_chol)
    
    ## initialize tau2
    tau2    <- rgamma(1, priors$alpha_tau2, priors$beta_tau2) 
    
    if (!is.null(inits$tau2)) {
        if (all(!is.na(inits$tau2))) {
            tau2 <- inits$tau2
        }
    }
    
    tau     <- sqrt(tau2)
    
    Sigma <- tau2 * R
    Sigma_chol <- tau * R_chol
    Sigma_inv <- 1 / tau2 * R_inv
    
    ## initialize rho
    rho      <- runif(1, -1, 1)
    
    if (!is.null(inits$rho)) {
        if (all(!is.na(inits$rho))) {
            rho <- inits$rho
        }
    }
    
    ## initialize sigma2
    sigma2  <- rgamma(J-1, priors$alpha_sigma2, priors$beta_sigma2)
    if (!is.null(inits$sigma2)) {
        if (all(!is.na(inits$sigma2))) {
            sigma2 <- inits$sigma2
        }
    }
    sigma   <- sqrt(sigma2)
    
    ## initialize Z and eta
    Z   <- matrix(0, N, n_time)
    eta <- array(0, dim = c(N, J-1, n_time))
    
    for (tt in 1:n_time) {
        if (tt == 1) {
            Z[, tt]     <- Z0
        } else {
            Z[, tt]   <- Xgamma + t(Sigma_chol) %*%  rnorm(N)
        }
        eta[, , tt] <- cbind(1, Z[, tt]) %*% beta +
            sapply(1:(J-1), function(j) rnorm(N, 0, sigma[j]))
    }

    if (!is.null(inits$Z)) {
        if (all(!is.na(inits$Z))) {
            Z <- inits$Z
        }
    }
    if (!is.null(inits$eta)) {
        if (all(!is.na(inits$eta))) {
            eta <- inits$eta
        }
    }
    
    
    ##
    ## sampler config options -- to be added later
    ## 
    #
    # bool sample_beta = true;
    # if (params.containsElementNamed("sample_beta")) {
    #     sample_beta = as<bool>(params["sample_beta"]);
    # }
    # 
    
    ##
    ## initialize omega
    ##
    
    omega <- array(0, dim = c(N, J-1, n_time))
    # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
    
    if (!is.null(inits$omega)) {
        if (!is.na(inits$omega)) {
            omega <- inits$omega
        }
    }
    
    ##
    ## setup save variables
    ##
    
    n_save       <- params$n_mcmc / params$n_thin
    beta_save    <- array(0, dim = c(n_save, Q+1, J-1))
    gamma_save   <- matrix(0, n_save, p)
    rho_save     <- rep(0, n_save)
    tau2_save    <- rep(0, n_save)
    theta_save <- NULL
    if (corr_fun == "matern") {
        theta_save <- matrix(0, n_save, 2)
    } else if (corr_fun == "exponential") {
        theta_save <- rep(0, n_save)
    }
    sigma2_save  <- matrix(0, n_save, J-1)
    Z_save       <- array(0, dim = c(n_save, N, n_time))
    eta_save     <- array(0, dim = c(n_save, N, J-1, n_time))
    
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
    
    ## tuning for rho
    rho_accept       <- 0
    rho_accept_batch <- 0
    rho_tune         <- 0.025
    
    ##
    ## Starting MCMC chain
    ##
    
    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")
    
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
        ## sample Omega
        ##
        
        if (verbose)
            message("sample omega")
        
        # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
        
        ##
        ## sample eta
        ##
        
        if (sample_eta) {
            if (verbose)
                message("sample eta")
            
            for (tt in 1:n_time) {
                eta[, , tt] <- sapply(1:(J-1), function(j) {
                    sigma2_tilde <- 1 / (1 / sigma2[j] + omega[, j, tt])
                    mu_tilde     <- 1 / sigma2[j] * cbind(1, Z[, tt]) %*% beta[, j] + kappa[, j, tt]
                    rnorm(
                        N, 
                        sigma2_tilde * mu_tilde,
                        sqrt(sigma2_tilde)
                    )
                })
            }
        }
        
        ##
        ## sample sigma2
        ##
        
        if (sample_sigma2) {
            if (verbose)
                message("sample sigma2")
            
            for (j in 1:(J-1)) {
                devs      <- sapply(1:n_time, function(tt) eta[, j, tt] - cbind(1, Z[, tt]) %*% beta[, j])
                SS        <- sum(devs^2)
                sigma2[j] <- 1 / rgamma(1, N * n_time / 2 + priors$alpha_sigma2, SS / 2 + priors$beta_sigma2) 
            }
            sigma     <- sqrt(sigma2)
        }
        
        ##
        ## sample beta -- double check these values
        ##
        
        if (sample_beta) {
            if (verbose)
                message("sample beta")
            
            for (j in 1:(J-1)) {
                ## only use the modern climate state to update beta
                A <- 1 / sigma2[j] * t(cbind(1, Z[, 1])) %*% cbind(1, Z[, 1]) + Sigma_beta_inv
                b <- 1 / sigma2[j] *t(cbind(1, Z[, 1])) %*% eta[, j, 1] + Sigma_beta_inv %*% mu_beta
                if (!sample_beta_modern) {
                    A <- A + 1 / sigma2[j] * apply(sapply(2:n_time, function(tt) t(cbind(1, Z[, tt])) %*% cbind(1, Z[, tt]), simplify = "array"), c(1, 2), sum)
                    b <- b + 1 / sigma2[j] * rowSums(sapply(2:n_time, function(tt) t(cbind(1, Z[, tt])) %*% eta[, j, tt]))
                }
                ## guarantee a symmetric matrix
                A <- (A + t(A)) / 2
                beta[, j]   <- rmvn_arma(A, b)
            }
        }
        
        ##
        ## sample gamma
        ##
        
        if (sample_gamma) {
            if (verbose)
                message("sample gamma")
            
            ## check these in very fine detail -- these are very preliminary guesses
            tX_Sigma_inv <- t(X) %*% Sigma_inv
            tX_Sigma_inv_X <- tX_Sigma_inv %*% X
            A <- tX_Sigma_inv_X + Sigma_gamma_inv 
            b <- tX_Sigma_inv %*% Z[, 1] + Sigma_gamma_inv %*% mu_gamma
            if (!sample_gamma_modern) {
                A <- A + (n_time - 1) * (1 - rho)^2 * tX_Sigma_inv_X
                b <- b + (1 - rho) * rowSums(tX_Sigma_inv %*% (Z[, 2:n_time] - rho * Z[,1:(n_time - 1)]))
            }
            ## guarantee a symmetric matrix
            A      <- (A + t(A)) / 2
            gamma  <- rmvn_arma(A, b)
            Xgamma <- as.vector(X %*% gamma)
        }
        
        ##
        ## sample Z
        ##
        
        if (sample_Z) {
            if (verbose)
                message("sample Z")
            
            for(tt in 2:n_time) {
                if (tt == n_time) {
                    A_Z <- sum(beta[2, ]^2 / sigma2) * diag(N) + Sigma_inv
                    b_Z <- rowSums(sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j]))) +
                        as.vector(Sigma_inv %*% (rho * (Z[, tt - 1] - Xgamma) + Xgamma))
                    Z[, tt] <- rmvn_arma(A_Z, b_Z)
                } else {
                    A_Z <- sum(beta[2, ]^2 / sigma2) * diag(N) + (1 + rho^2) * Sigma_inv
                    b_Z <- rowSums(sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j]))) +
                        as.vector(rho * Sigma_inv %*% (Z[, tt + 1] + Z[, tt - 1] - 2 * Xgamma)) +
                        as.vector((1 + rho^2) * Sigma_inv %*% Xgamma)
                    Z[, tt] <- rmvn_arma(A_Z, b_Z)
                }       
            }
        }
        
        ##
        ## sample rho
        ##
        
        if (sample_rho) {
            if (verbose)
                message("sample rho")
            
            # rho_star <- rnorm(1, rho, rho_tune)
            # if (rho_star < 1 & rho_star > -1) {
            #     mh1 <- sum(
            #         sapply(
            #             2:n_time, 
            #             function(tt) {
            #                 mvnfast::dmvn(Z[, tt], (1 - rho_star) * Xgamma + rho_star * Z[, tt-1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
            #             }
            #         )
            #     ) 
            #     ## parallelize this        
            #     mh2 <- sum(
            #         sapply(
            #             2:n_time, 
            #             function(tt) {
            #                 mvnfast::dmvn(Z[, tt], (1 - rho) * Xgamma + rho * Z[, tt-1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
            #             }
            #         )
            #     ) 
            #     mh <- exp(mh1 - mh2)
            #     if (length(mh) > 1)
            #         stop("error in mh for rho")
            #     if (mh > runif(1, 0.0, 1.0)) {
            #         rho   <- rho_star
            #         if (k <= params$n_adapt) {
            #             rho_accept_batch <- rho_accept_batch + 1.0 / 50.0
            #         } else {
            #             rho_accept <- rho_accept + 1.0 / params$n_mcmc
            #         }
            #     }
            #     
            #     ## update tuning
            #     if (k <= params$n_adapt) {
            #         if (k %% 50 == 0){
            #             out_tuning <- update_tuning(
            #                 k,
            #                 rho_accept_batch, 
            #                 rho_tune
            #             )
            #             rho_tune         <- out_tuning$tune
            #             rho_accept_batch <- out_tuning$accept
            #         }
            #     }
            # }

            rho_vals <- rowSums(
                sapply(2:n_time, function(tt) {
                    tZ_Xgamma_Sigma_inv <- t(Z[, tt-1] - Xgamma) %*% Sigma_inv
                    c(
                        tZ_Xgamma_Sigma_inv %*% (Z[, tt-1] - Xgamma),
                        tZ_Xgamma_Sigma_inv %*% (Z[, tt] - Xgamma)
                    )
                })
            )

            a_rho <- rho_vals[1]
            b_rho <- rho_vals[2]
            rho   <- rtrunc(1, "norm", a = -1, b = 1, mean = b_rho / a_rho, sd = sqrt(1 / a_rho))
        }
        
        ##
        ## sample spatial correlation parameters theta
        ##
        
        if (sample_theta) {
            if (verbose)
                message("sample theta")
            
            theta_star <- NULL
            if (corr_fun == "matern") {
                theta_star <- as.vector(
                    rmvn( 
                        n      = 1,
                        mu     = theta,
                        sigma  = lambda_theta * Sigma_theta_tune_chol,
                        isChol = TRUE
                    )
                )
            } else if (corr_fun == "exponential") {
                theta_star <- rnorm(1, theta, theta_tune)
            }
            R_star          <- correlation_function(D, theta_star, corr_fun = corr_fun)
            ## add in faster parallel cholesky as needed
            R_chol_star    <- tryCatch(
                chol(R_star),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                    num_chol_failures <- num_chol_failures + 1
                    chol(R_star + 1e-8 * diag(N))                    
                }
            )
            R_inv_star      <- chol2inv(R_star)
            Sigma_star      <- tau2 * R_star
            Sigma_chol_star <- tau * R_chol_star
            Sigma_inv_star  <- 1 / tau2 * R_inv_star
            
            ## parallelize this
            mh1 <- mvnfast::dmvn(Z[, 1], Xgamma, Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) +
                ## prior
                mvnfast::dmvn(theta_star, theta_mean, theta_var, log = TRUE)
            if (!sample_theta_modern) {
                mh1 <- mh1 + sum(
                    sapply(
                        2:n_time, 
                        function(tt) {
                            mvnfast::dmvn(Z[, tt], (1 - rho) * Xgamma + rho * Z[, tt-1], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                        }
                    )
                ) 
            }
            ## parallelize this        
            mh2 <- mvnfast::dmvn(Z[, 1], Xgamma, Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)  +
                ## prior
                mvnfast::dmvn(theta, theta_mean, theta_var, log = TRUE)
            if(!sample_theta_modern) {
                mh2 <- mh2 + sum(
                    sapply(
                        2:n_time, 
                        function(tt) {
                            mvnfast::dmvn(Z[, tt], (1 - rho) * Xgamma + rho * Z[, tt-1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
                        }
                    )
                )
            }
            
            mh <- exp(mh1 - mh2)
            if(mh > stats::runif(1, 0, 1)) {
                theta      <- theta_star
                R          <- R_star
                R_chol     <- R_chol_star
                R_inv      <- R_inv_star 
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
        }
        
        ##
        ## sample tau2
        ##
        
        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            if (sample_tau2_modern) {
                devs <- Z[, 1] - Xgamma
                SS   <- sum(devs %*% R_inv %*% devs)
                tau2 <- 1 / rgamma(1, N / 2 + priors$alpha_tau2, SS / 2 + priors$beta_tau2) 
            } else {
                ## can we make this more efficient?
                devs      <- matrix(0, N, n_time)
                devs[, 1] <- Z[, 1] - Xgamma 
                for (tt in 2:n_time) {
                    devs[, tt] <- Z[, tt] - Xgamma - rho * (Z[, tt-1] - Xgamma)
                }
                ## double check this math later -- seems right for now
                SS <- sum(
                    sapply(1:n_time, function(tt) {
                        devs[, tt] * (R_inv %*% devs[, tt])
                    })
                )
                tau2       <- 1 / rgamma(1, N * n_time / 2 + priors$alpha_tau2, SS / 2 + priors$beta_tau2) 
            }
            tau        <- sqrt(tau2)
            Sigma      <- tau2 * R
            Sigma_chol <- tau * R_chol
            Sigma_inv  <- 1 / tau2 * R_inv 
        }

        ##
        ## save variables
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ] <- beta
                gamma_save[save_idx, ]  <- gamma
                rho_save[save_idx]      <- rho
                if (corr_fun == "matern") {
                    theta_save[save_idx, ] <- theta
                } else if (corr_fun == "exponential") {
                    theta_save[save_idx]   <- theta
                }
                tau2_save[save_idx]     <- tau2
                sigma2_save[save_idx, ] <- sigma2
                Z_save[save_idx, , ]    <- Z
                eta_save[save_idx, , , ]  <- eta
            }
        }
        
        ##
        ## End of MCMC loop
        ##
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    if (sample_theta)
        message("Acceptance rate for theta is ", theta_accept)

    # if (sample_rho)
    #     message("Acceptance rate for rho is ", rho_accept)
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    out <- list(
        beta   = beta_save,
        gamma  = gamma_save,
        rho    = rho_save,
        theta  = theta_save,
        tau2   = tau2_save,
        sigma2 = sigma2_save,
        Z      = Z_save,
        eta    = eta_save
    )
    
    class(out) <- "pg_mvgp_univariate" 
    return(out)
}