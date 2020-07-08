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
#' @export

pgLM <- function(
    Y, 
    X,
    params,
    priors,
    n_cores = 1L,
    inits = NULL,
    config = NULL,
    n_chain       = 1,
    sample_rmvn = FALSE
    # pool_s2_tau2  = true,
    # file_name     = "DM-fit",
    # corr_function = "exponential"
) {

    ##
    ## Run error checks
    ## 
    
    check_input(Y, X)
    check_params(params)
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    N  <- nrow(Y)
    J  <- ncol(Y)
    p  <- ncol(X)
    
    ## Calculate Mi
    Mi <- matrix(0, N, J-1)
    for(i in 1: N){
        if (J == 2) {
            Mi[i,] <- sum(Y[i, ])
        } else {
            Mi[i,] <- sum(Y[i, ]) - c(0, cumsum(Y[i,][1:(J-2)]))
        }
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
    Sigma_beta     <- diag(p)
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
    
    eta  <- X %*% beta 
    
    
    ##
    ## initialize omega
    ##
    
    omega <- matrix(0, N, J-1)
    omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)

    
    # ## can parallelize this, see example
    # for (i in 1:N) {
    #     for (j in 1:(J-1)) {
    #         if (Mi[i, j] != 0) {
    #             omega[i, j] <- pgdraw(Mi[i, j], eta[i, j])
    #         }
    #     }
    # }
    if (!is.null(inits$omega)) {
        if (!is.na(inits$omega)) {
            omega <- inits$omega
        }
    }
    
    ## don't need a diagonal matrix form
    # Omega <- vector(mode = "list", length = J-1)
    # for (j in 1:(J - 1)) {
    #     Omega[[j]] <- diag(omega[, j])
    # }
    
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
    
    n_save    <- params$n_mcmc / params$n_thin
    beta_save <- array(0, dim = c(n_save, p, J-1))
    eta_save  <- array(0, dim = c(n_save, N, J-1))
    
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

        omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)

        
        ## don't need a diagonal matrix form
        # for (j in 1:(J-1)) {
        #     Omega[[j]] <- diag(omega[, j])
        # }
        
        ##
        ## sample beta -- double check these values
        ##
        
        ## parallelize this update -- each group of parameteres is 
        ## conditionally independent given omega and kappa(y)
        
        # beta <- 1:(J-1) %>%
        #     future_map(
        #         sample_beta,
        #         .options = future_options(
        #             globals = TRUE,#c(
        #             #     "Sigma_beta_inv", 
        #             #     "mu_beta", 
        #             #     "X", 
        #             #     "omega",
        #             #     "kappa"
        #             # ),
        #             packages = "mvnfast"
        #         )
        #     ) %>%
        #     unlist() %>%
        #     matrix(., p, J-1)
        
        if (sample_rmvn) {
            for (j in 1:(J-1)) {

                ## use the efficient Cholesky sampler
                
                ## there is an issue in rmvn_arma(A, b) -- I don't know why as
                ##    they are samples from the same distribution...
                A <- Sigma_beta_inv + t(X) %*% (omega[, j] * X)
                b <- as.vector(Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j])
                beta[, j]   <- rmvn_arma(A, b)
                # beta[, j]   <- rmvn_R(A, b)
            }
        } else {
            ## parallelization of this using furrr is not faster
            for (j in 1:(J-1)) {
                ## can make this much more efficient
                Sigma_tilde <- chol2inv(chol(Sigma_beta_inv + t(X) %*% (omega[, j] * X)))
                # Sigma_tilde <- chol2inv(chol(Sigma_beta_inv + t(X) %*% (Omega[[j]] %*% X)))
                mu_tilde    <- c(Sigma_tilde %*% (Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j]))
                beta[, j]   <- mvnfast::rmvn(1, mu_tilde, Sigma_tilde)
            }
        }

        # beta <- matrix(0, 4, 9)
        eta <- X %*% beta
        
        # message(
        #     "mean beta = ", round(mean(beta), digits = 2), 
        #     "    sd beta = ", round(sd(beta), digits = 2),
        #     "    mean eta = ", round(mean(eta), digits = 2),
        #     "    sd eta = ", round(sd(eta), digits = 2)
        # )
        
        
        ##
        ## save variables
        ##
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ] <- beta
            }
        }
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    out <- list(beta = beta_save)
    class(out) <- "pgLM"
    
    return(out)
}