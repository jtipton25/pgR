#' Bayesian Polya-gamma regression
#' 
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param Y is a \code{n \times d} matrix of compositional count data.
#' @param X is a \code{n \times p} matrix of climate variables.
#' @param params is the list of parameter settings.
#' @param priors is the list of prior settings. 
#' @param inits is the list of intial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.

## polya-gamma linear regression model
pgLM <- function(
    Y, 
    X,
    params,
    priors,
    cores = 1,
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
    
    check_input(Y, X)
    check_params(params)
    check_inits(params, inits)
    # check_config(params, config)
    
    N  <- nrow(y)
    J  <- ncol(y)
    df <- ncol(X)
    
    ## Calculate Mi
    Mi <- matrix(0, N, J-1)
    for(i in 1: N){
        Mi[i,] <- sum(y[i, ]) - c(0, cumsum(y[i,][1:(J-2)]))
    }
    
    ## initialize kappa
    kappa <- matrix(0, N, J-1)
    for (i in 1: N) {
        kappa[i,] <- y[i, 1:(J - 1)]- Mi[i, ] / 2
    }
    
    ##
    ## initial values
    ##
    
    ## currently using default priors
    
    mu_beta        <- rep(0, df)
    
    ## do I want to change this to be a penalized spline?
    # Q_beta <- make_Q(params$df, 1) 
    Sigma_beta     <- 10 * diag(x = df)
    Sigma_beta_inv <- chol2inv(chol(Sigma_beta))
    ## clean up this check
    if (!is.null(priors$mu_beta)) {
        if (!is.na(priors$mu_beta)) {
            mu_beta <- priors$mu_beta
        }
    }
    
    
    ## clean up this check
    if (!is.null(priors$Sigma_beta)) {
        if (!is.na(priors$Sigma_beta)) {
            Sigma_beta <- priors$Sigma_beta
        }
    }
    Sigma_beta_chol <- chol(Sigma_beta)
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##
    
    beta <- mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE)
    ## clean up this check
    if (!is.null(inits$beta)) {
        if (!is.na(inits$beta)) {
            beta <- inits$beta
        }
    }
    
    eta  <- X %*% beta  
    
    ##
    ## initialize omega
    ##
    
    omega <- matrix(0, N, J-1)
    
    ## can parallelize this, see example
    for (i in 1:N) {
        for (j in 1:(J-1)) {
            if (Mi[i, j] != 0) {
                omega[i, j] <- pgdraw(Mi[i, j], eta[i, j])
            }
        }
    }
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
    beta_save <- array(0, dim = c(n_save, df, J-1))
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
        ## sample Omega
        ##
        
        for (i in 1:N) {
            for (j in 1:(J-1)) {
                if(Mi[i, j] != 0){
                    omega[i, j] <- pgdraw(Mi[i, j], eta[i, j])
                }
                else {
                    omega[i, j] <- 0
                }
            }
        }
        
        for (j in 1:(J-1)) {
            Omega[[j]] <- diag(omega[, j])
        }
        
        ##
        ## sample beta -- double check these values
        ##
        
        for (j in 1:(J-1)) {
            ## can make this much more efficient
            Sigma_tilde <- chol2inv(chol(Sigma_beta_inv + t(X) %*% (Omega[[j]] %*% X))) 
            mu_tilde <- c(Sigma_tilde %*% (Sigma_beta_inv %*% mu_beta + t(X) %*% kappa[, j]))
            beta[, j] <- rmvn(1, mu_tilde, Sigma_tilde)
        }
        eta <- X %*% beta
        
        ##
        ## save variables
        ##
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ]  <- beta
            }
        }
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    return(
        list(
            beta       = beta_save
        )
    )
}