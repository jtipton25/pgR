#' Bayesian Polya-gamma regression
#' 
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param Y is a \eqn{n \times J}{n x J} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
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
#' @param sample_rmvn is an indicator whether the faster multivariate normal sampler is used. 
#' @importFrom BayesLogit rpg
#' @export

pg_lm <- function(
    Y, 
    X,
    params,
    priors,
    n_cores = 1L,
    inits = NULL,
    config = NULL,
    n_chain       = 1) {
    
    ##
    ## Run error checks
    ## 
    
    check_input_pg_lm(Y, X)
    check_params(params)
    
    if (!is_positive_integer(n_cores, 1))
        stop("n_cores must be a positive integer")
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    N  <- nrow(Y)
    J  <- ncol(Y)
    p  <- ncol(X)
    
    ## We assume a partially missing observation is the same as 
    ## fully missing. The index allows for fast accessing of missing
    ## observations
    missing_idx <- rep(FALSE, N)
    for (i in 1:N) {
        missing_idx[i] <- any(is.na(Y[i, ]))
    }
    
    message("There are ", ifelse(any(missing_idx), sum(missing_idx), "no"), " observations with missing count vectors")
    
    # Calculate Mi
    Mi <- calc_Mi(Y)
    
    ## create an index for nonzero values
    nonzero_idx <- Mi != 0
    n_nonzero   <- sum(nonzero_idx)
    
    # Calculate kappa
    kappa <- calc_kappa(Y, Mi)
    
    ##
    ## initial values
    ##
    
    ## currently using default priors
    
    mu_beta        <- rep(0, p)
    Sigma_beta     <- 10 * diag(p)
    ## clean up this check
    if (!is.null(priors[['mu_beta']])) {
        if (all(!is.na(priors[['mu_beta']]))) {
            mu_beta <- priors[['mu_beta']]
        }
    }
    
    ## clean up this check
    if (!is.null(priors[['Sigma_beta']])) {
        if (all(!is.na(priors[['Sigma_beta']]))) {
            Sigma_beta <- priors[['Sigma_beta']]
        }
    }
    Sigma_beta_chol <- chol(Sigma_beta)
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##

    beta <- t(mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    ## initial values for beta
    if (!is.null(inits[['beta']])) {
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
        }
    }
    
    eta  <- X %*% beta 
    
    
    ##
    ## initialize omega
    ##
    
    omega <- matrix(0, N, J-1)
    # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])

    ## initial values for omega
    if (!is.null(inits[['omega']])) {
        if (all(!is.na(inits[['omega']]))) {
            omega <- inits[['omega']]
        }
    }
    
    
    ##
    ## setup config
    ##
    
    ## do we sample the regression parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
        }
    }
    
    ## do we sample the Polya-gamma parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_omega <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_omega']])) {
            sample_omega <- config[['sample_omega']]
        }
    }
    
    ## do we save the omega random variables
    save_omega <- FALSE
    if (!is.null(config)) {
        if (!is.null(config[['save_omega']])) {
            save_omega <- config[['save_omega']]
        }
    }
    
    ##
    ## setup save variables
    ##
    
    n_save    <- params$n_mcmc / params$n_thin
    beta_save <- array(0, dim = c(n_save, p, J-1))
    eta_save  <- array(0, dim = c(n_save, N, J-1))
    omega_save <- NULL
    if (save_omega) {
        omega_save <- array(0, dim = c(n_save, N, J-1)) 
    }
    
    ## 
    ## initialize tuning - no tuning in this model
    ##
    
    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")
    
    ##
    ## Starting MCMC chain
    ##
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
        ## sample omega
        ##

        if (sample_omega) {
            # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
            omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
        }
                
        ##
        ## sample beta 
        ##
        
        if (sample_beta) {
            for (j in 1:(J-1)) {
                ## use the efficient Cholesky sampler
                A <- Sigma_beta_inv + t(X) %*% (omega[, j] * X)
                b <- as.vector(Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j])
                beta[, j]   <- rmvn_arma(A, b)
                
            }
            # update eta
            eta <- X %*% beta
        }
        
        ##
        ## save variables
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ] <- beta
                eta_save[save_idx, , ] <- eta
                if (save_omega) {
                    omega_save[save_idx, , ] <- omega
                }
            }
        }
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    ##
    ## return the MCMC output 
    ## 
    
    out <- NULL
    if (save_omega) {
        out <- list(beta = beta_save, 
                    eta  = eta_save, 
                    omega = omega_save)        
    } else {
        out <- list(beta = beta_save, 
                    eta  = eta_save)
    }
    
    class(out) <- "pg_lm"
    
    return(out)
}