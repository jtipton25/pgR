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
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of intial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.

## polya-gamma spatially varying linear regression model
pgSPVLM <- function(
    Y, 
    X,
    locs, 
    params,
    priors,
    n_cores = 1L,
    inits = NULL,
    config = NULL,
    n_chain       = 1
    # pool_s2_tau2  = true,
    # file_name     = "DM-fit",
    # corr_function = "exponential"
) {
    
    ##
    ## Run error checks
    ## 
    
    check_input_spatial(Y, X, locs)
    check_params(params)
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    N  <- nrow(Y)
    J  <- ncol(Y)
    p <- ncol(X)
    D <- fields::rdist(locs)
    
    ## Calculate Mi
    Mi <- matrix(0, N, J-1)
    for(i in 1: N){
        Mi[i,] <- sum(Y[i, ]) - c(0, cumsum(Y[i,][1:(J-2)]))
    }
    
    ## create an index for nonzero values
    nonzero_idx <- Mi != 0
    
    ## initialize kappa
    kappa <- matrix(0, N, J-1)
    for (i in 1: N) {
        kappa[i,] <- Y[i, 1:(J - 1)]- Mi[i, ] / 2
    }
    
    ##
    ## initial values
    ##
    
    ## currently using default priors
    
    mu_beta        <- rep(0, p)
    
    ## do I want to change this to be a penalized spline?
    # Q_beta <- make_Q(params$p, 1) 
    Sigma_beta     <- 10 * diag(p)
    ## clean up this check
    if (!is.null(priors$mu_beta)) {
        if (all(!is.na(priors$mu_beta))) {
            mu_beta <- priors$mu_beta
        }
    }
    
    ## clean up this check
    if (!is.null(priors$Sigma_beta)) {
        if (all(!is.na(priors$Sigma_beta))) {
            Sigma_beta <- priors$Sigma_beta
        }
    }
    Sigma_beta_chol <- chol(Sigma_beta)
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##
    
    beta <- t(mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    ## clean up this check
    if (!is.null(inits$beta)) {
        if (all(!is.na(inits$beta))) {
            beta <- inits$beta
        }
    }
    Xbeta <- X %*% beta
    
    ##
    ## initialize omega
    ##
    
    omega <- matrix(0, N, J-1)
    omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    
    if (!is.null(inits$omega)) {
        if (!is.na(inits$omega)) {
            omega <- inits$omega
        }
    }
    
    Omega <- vector(mode = "list", length = J-1)
    for (j in 1:(J - 1)) {
        Omega[[j]] <- diag(omega[, j])
    }
    
    ##
    ## initialize spatial Gaussian process -- share parameters across the different components
    ##    can generalize to each component getting its own covariance
    ##
    
    theta_mean <- c(priors$mean_nu, priors$mean_range)
    theta_var  <- diag(c(priors$sd_nu, priors$sd_range)^2)
    theta <- as.vector(mvnfast::rmvn(1, theta_mean, theta_var))
    if (!is.null(inits$theta)) {
        if (!is.na(inits$theta)) {
            theta <- inits$theta
        }
    }
    
    tau2       <- min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
    Sigma      <- tau2 * correlation_function(D, theta)
    ## add in faster parallel cholesky as needed
    Sigma_chol <- chol(Sigma)
    Sigma_inv  <- chol2inv(Sigma_chol)
    
    eta  <- X %*% beta + t(mvnfast::rmvn(J-1, rep(0, N), Sigma_chol, isChol = TRUE))
    
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
    ## setup save variables
    ##
    
    n_save     <- params$n_mcmc / params$n_thin
    beta_save  <- array(0, dim = c(n_save, p, J-1))
    tau2_save  <- rep(0, n_save)
    theta_save <- matrix(0, n_save, 2)
    eta_save   <- array(0, dim = c(n_save, N, J-1))
    
    ## 
    ## initialize tuning 
    ##
    
    ##
    ## tuning variables for adaptive MCMC
    ##
    theta_batch           <- matrix(0, 50, 2)
    theta_accept          <- 0
    theta_accept_batch    <- 0
    lambda_theta          <- 0.05
    Sigma_theta_tune      <- 1.8 * diag(2) - .8
    Sigma_theta_tune_chol <- chol(Sigma_theta_tune)
    
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
        
        omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        
        # for (i in 1:N) {
        #     for (j in 1:(J-1)) {
        #         if(Mi[i, j] != 0){
        #             omega[i, j] <- pgdraw(Mi[i, j], eta[i, j])
        #         }
        #         else {
        #             omega[i, j] <- 0
        #         }
        #     }
        # }
        
        for (j in 1:(J-1)) {
            Omega[[j]] <- diag(omega[, j])
        }
        
        ##
        ## sample beta -- double check these values
        ##
        
        ## modify this for the spatial process eta
        
        ## parallelize this update -- each group of parameteres is 
        ## conditionally independent given omega and kappa(y)
        for (j in 1:(J-1)) {
            ## can make this much more efficient
            # Sigma_tilde <- chol2inv(chol(Sigma_beta_inv + t(X) %*% (Omega[[j]] %*% X))) 
            # mu_tilde    <- c(Sigma_tilde %*% (Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j]))
            # beta[, j]   <- mvnfast::rmvn(1, mu_tilde, Sigma_tilde)
            tXSigma_inv <- t(X) %*% Sigma_inv
            A <- tXSigma_inv %*% X + Sigma_beta_inv
            b <- tXSigma_inv %*% eta[, j] + Sigma_beta_inv %*% mu_beta
            beta[, j]   <- rmvn_arma(A, b)
        }
        
        Xbeta <- X %*% beta
        
        ##
        ## sample spatial correlation parameters theta
        ##
        
        theta_star <- mvnfast::rmvn( 
            n      = 1,
            mu     = theta,
            sigma  = lambda_theta * Sigma_theta_tune_chol,
            isChol = TRUE
        )
        Sigma_star       <- tau2 * correlation_function(D, theta_star)
        ## add in faster parallel cholesky as needed
        Sigma_chol_star <- chol(Sigma_star)
        Sigma_inv_star  <- chol2inv(Sigma_star)
        ## parallelize this
        mh1 <- sum(
            sapply(
                1:(J-1), 
                function(j) {
                    mvnfast::dmvn(eta[, j], Xbeta[, j], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) }
            )
        ) +
            ## prior
            mvnfast::dmvn(theta_star, theta_mean, theta_var, log = TRUE)
        ## parallelize this        
        mh2 <- sum(
            sapply( 
                1:(J-1),
                function(j) {
                    mvnfast::dmvn(eta[, j], Xbeta[, j], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)
                }
            )
        ) +
            ## prior
            mvnfast::dmvn(theta, theta_mean, theta_var, log = TRUE)
        
        mh <- exp(mh1 - mh2)
        if(mh > stats::runif(1, 0, 1)) {
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
            save_idx <- k %% 50
            if ((k %% 50) == 0) {
                save_idx <- 50
            } 
            theta_batch[save_idx, ] <- theta 
            if (k %% 50 == 0) {
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
                lambda_theta_tune     <- out_tuning$lambda
                theta_accept_batch    <- out_tuning$accept
            } 
        }   
        
        ##
        ## sample spatial process variance tau2
        ##
        
        devs       <- eta - Xbeta
        SS         <- sum(devs * (tau2 * Sigma_inv %*% devs))
        tau2       <- 1 / stats::rgamma(1, N * (J - 1) / 2 + priors$alpha_tau, SS / 2 + priors$beta_tau) 
        Sigma      <- tau2 * correlation_function(D, theta) 
        ## add in faster parallel cholesky as needed
        Sigma_chol <- chol(Sigma)
        Sigma_inv  <- chol2inv(Sigma_chol)
        
        
        ##
        ## sample eta
        ##
        
        ## double check this and add in fixed effects X %*% beta
        for (j in 1:(J-1)) {
            ## can make this much more efficient
            ## can this be parallelized? seems like it
            A        <- Sigma_inv + Omega[[j]]
            b        <- Sigma_inv %*% Xbeta[, j] + kappa[, j]
            eta[, j] <- rmvn_arma(A, b) 
        }
        
        ##
        ## save variables
        ##
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ] <- beta
                theta_save[save_idx, ]  <- theta
                tau2_save[save_idx]     <- tau2
                eta_save[save_idx, , ]  <- eta
            }
        }
        
        ##
        ## End of MCMC loop
        ##
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    message("Acceptance rate for theta is ", theta_accept)
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    out <- list(
        beta  = beta_save,
        theta = theta_save,
        tau2  = tau2_save,
        eta   = eta_save
    )
    class(out) <- "pgSPVLM"
    
    return(out)
}


