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
#' @param M The number of resolutions.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.
#' @param progress is a logical input that determines whether to print a progress bar.
#' @param verbose is a logical input that determines whether to print more detailed messages.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#' @importFrom stats rnorm rgamma runif dnorm
#' @importFrom fields rdist
#' @importFrom mvnfast rmvn dmvn
#' @importFrom BayesLogit rpg
#' @import BayesMRA
#' @import spam  
#'  
#' @export

pg_splm_mra <- function(
    Y, 
    X,
    locs, 
    params,
    priors,
    n_cores       = 1L,
    M             = 4,
    n_coarse_grid = 10,
    inits         = NULL,
    config        = NULL,
    n_chain       = 1,
    progress      = FALSE,
    verbose       = FALSE,
    use_spam      = TRUE
) {
    
    ##
    ## Run error checks
    ## 
    
    check_input_pg_splm(Y, X, locs)
    check_params(params)
    
    if (!use_spam)
        stop("The only sparse matrix pacakage available is spam")
    if (!is_positive_integer(n_cores, 1))
        stop("n_cores must be a positive integer")
    
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0
    
    
    N  <- nrow(Y)
    J  <- ncol(Y)
    p <- ncol(X)
    tX <- t(X)
    tXX <- tX %*% X
    
    ## We assume a partially missing observation is the same as 
    ## fully missing. The index allows for fast accessing of missing
    ## observations
    missing_idx <- rep(FALSE, N)
    for (i in 1:N) {
        missing_idx[i] <- any(is.na(Y[i, ]))
    }
    
    message("There are ", ifelse(any(missing_idx), sum(missing_idx), "no"), " observations with missing count vectors")
    
    ## Calculate Mi
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
    
    ## do I want to change this to be a penalized spline?
    # Q_beta <- make_Q(params$p, 1) 
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
    
    beta <- t(rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))

    ## initialize with mean 0
    beta[1, ] <- 0
    ## clean up this check
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
    
    sigma2 <- pmin(rgamma(J-1, alpha_sigma2, beta_sigma2), 5)
    
    ## initial values for sigma2
    ## add in checking later
    if (!is.null(inits[['sigma2']])) {
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
        }
    }
    
    ##
    ## setup MRA spatial basis
    ##
    
    MRA      <- mra_wendland_2d(locs, M, n_coarse_grid = n_coarse_grid, use_spam = use_spam)
    W        <- MRA$W
    n_dims   <- MRA$n_dims
    dims_idx <- MRA$dims_idx
    
    # if (RSR) {
    #     W <- IMPX %*% W
    # }
    tW <- NULL
    if (use_spam) {
        tW <- t(W)
    } else {
        tW <- Matrix::t(W)
    }
    
    tWW <- tW %*% W    
    
    ##
    ## priors for tau2
    ##
    
    alpha_tau2 <- 0.01
    beta_tau2  <- 0.01
    
    ## check if priors for alpha_tau2 are specified
    if (!is.null(priors[['alpha_tau2']])) {
        alpha_tau2 <- priors[['alpha_tau2']]
    }
    
    ## check if priors for beta_tau2 are specified
    if (!is.null(priors[['beta_tau2']])) {
        beta_tau2 <- priors[['beta_tau2']]
    }
    
    ##
    ## intialize a proper CAR structure to initialize the parameter alpha
    ##
    
    Q_alpha <- make_Q_alpha_2d(sqrt(n_dims), rep(0.999, length(n_dims)), use_spam = use_spam)
    tau2 <- matrix(0, M, J-1)
    for (j in 1:(J-1)) {
        tau2[, j] <- 100 * 2^(1:M) * pmax(rgamma(M, alpha_tau2, beta_tau2), 1)
    }
    
    if (!is.null(inits[['tau2']])) {
        if (all(!is.na(inits[['tau2']]))) {
            if (!is_positive_numeric_matrix(inits[['tau2']], M, J-1))
                stop ("If specified, inits$tau2 must be a M x J-1 matrix of positive values")
            ## if tau2 passes error checks
            tau2 <- inits[['tau2']]
        }
    }

    Q_alpha_tau2 <- vector(mode = "list", length = J-1)
    for (j in 1:(J-1)) {
        Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j], use_spam = use_spam)
    }
    
    ##
    ## initialize alpha
    ##
    
    ## define the sum-to-0 constraint for alpha 
    # eventually modify this so the options for constraint and joint are allowed
    constraints <- make_constraint(MRA, constraint = "resolution", joint = TRUE)
    A_constraint <- constraints$A_constraint
    a_constraint <- constraints$a_constraint
    
    alpha <- matrix(0, sum(n_dims), J-1)
    eta <- kappa   ## default initial value based on data Y
    
    if (use_spam) {
        for (j in 1:(J-1)) {
            ## double check this full conditional
            A_alpha    <- 1 / sigma2[j] * tWW + Q_alpha_tau2[[j]]
            b_alpha    <- 1 / sigma2[j] * tW %*% (eta[, j] - Xbeta[, j]) 
            alpha[, j] <- rmvnorm.canonical.const(1, b_alpha, A_alpha, 
                                                  A = A_constraint, a = a_constraint)
        }
    } else {
        stop("The only sparse matrix pacakage available is spam")
        # alpha[, 1] <- as.vector(chol2inv(chol(W %*% Q_alpha_tau2 %*% tW)) %*% (tW %*% (Q_alpha_tau2 %*% (Z0 - Xgamma))))
        # # alpha[, 1] <- as.vector(ginv(as.matrix(tWW)) %*% tW %*% (Z0 - Xgamma))
        # # alpha[, 1] <- as.vector(rmvn.sparse(1, rep(0, sum(n_dims)), CH = Cholesky(Q_alpha_tau2), prec = TRUE))     
        # for (tt in 2:n_time) {
        #     ## initialize with the current climate
        #     alpha[, tt] <- alpha[, 1]
        #     ## intialize with the prior process
        #     # alpha[, tt] <- as.vector(rmvn.sparse(1, rho * alpha[, tt - 1], CH = Cholesky(Q_alpha_tau2), prec = TRUE))
        # }
    }
    
    
    ## initialize an ICAR structure for fitting alpha
    Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(1, length(n_dims)), use_spam = use_spam)
    for (j in 1:(J-1)) {
        Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j], use_spam = use_spam)
    }
    
    ## precalculate the sparse cholesky structure for faster Gibbs updates
    Rstruct        <- NULL
    if (use_spam) {
        ## double check this full conditional
        A <- 1 / sigma2[j]  * tWW + Q_alpha_tau2[[j]]
        Rstruct <- chol(A)
    }
    
    
    ## initial values for alpha
    ## add in checking later
    if (!is.null(inits[['alpha']])) {
        if (all(!is.na(inits[['alpha']]))) {
            alpha <- inits[['alpha']]
        }
    }
    
    W_alpha <- W %*% alpha
    
    
    eta <- Xbeta + W_alpha + sapply(1:(J-1), function(j) rnorm(N, 0, sqrt(sigma2[j])))

    ## initial values for eta
    ## add in checking later
    if (!is.null(inits[['eta']])) {
        if (all(!is.na(inits[['eta']]))) {
            eta <- inits[['eta']]
        }
    }
    
    
    ##
    ## sampler config options 
    ## 
    
    ## do we sample the functional relationship parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
        }
    }
    
    ## do we sample the latent spatial parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_alpha <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_alpha']])) {
            sample_alpha <- config[['sample_alpha']]
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
    
    ## do we sample the climate variance parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
        }
    }
    
    ## do we sample the overdispersion parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2']])) {
            sample_sigma2 <- config[['sample_sigma2']]
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
    ## initialize omega
    ##
    
    omega <- matrix(0, N, J-1)
    # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
    
    Omega <- vector(mode = "list", length = J-1)
    for (j in 1:(J - 1)) {
        Omega[[j]] <- diag(omega[, j])
    }
    
    save_omega <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['save_omega']])) {
            save_omega <- config[['save_omega']]
        }
    }
    
    ##
    ## setup save variables
    ##
    
    n_save      <- params$n_mcmc / params$n_thin
    beta_save   <- array(0, dim = c(n_save, p, J-1))
    tau2_save   <- array(0, dim = c(n_save, M, J-1))
    alpha_save  <- array(0, dim = c(n_save, sum(n_dims), J-1))
    eta_save    <- array(0, dim = c(n_save, N, J-1))
    sigma2_save <- matrix(0, n_save, J-1)
    omega_save <- NULL
    if (save_omega) {
        omega_save   <- array(0, dim = c(n_save, N, J-1))
    }
    
    
    ## 
    ## initialize tuning 
    ##

    ## no tuning necessary for this model    
    
    ##
    ## Starting MCMC chain
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
        ## sample Omega
        ##
        
        if (verbose)
            message("sample omega")
        
        # omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        omega[nonzero_idx] <- rpg(n_nonzero, Mi[nonzero_idx], eta[nonzero_idx])
        
        for (j in 1:(J-1)) {
            Omega[[j]] <- diag(omega[, j])
        }
        
        ##
        ## sample beta 
        ##
        
        if (sample_beta) {
            if (verbose)
                message("sample beta")
            
            for (j in 1:(J-1)) {
                A <- tXX / sigma2[j] + Sigma_beta_inv
                ## guarantee a symmetric matrix
                A         <- (A + t(A)) / 2
                b         <- tX %*% (eta[, j] - W_alpha[, j]) / sigma2[j] + Sigma_beta_inv %*% mu_beta
                beta[, j] <- rmvn_arma(A, b)
            }
            Xbeta <- X %*% beta
        }
        
        ##
        ## sample spatial random effects alpha
        ##
        
        if(sample_alpha) {
            if (verbose)
                message("sample alpha")
            
            ## double check this full conditional
            for (j in 1:(J-1)) {      
                A_alpha    <- 1 / sigma2[j] * tWW + Q_alpha_tau2[[j]]
                b_alpha    <- 1 / sigma2[j] * tW %*% (eta[, j] - Xbeta[, j]) 
                alpha[, j] <- rmvnorm.canonical.const(1, b_alpha, A_alpha, 
                                                      Rstruct = Rstruct,
                                                      A = A_constraint,
                                                      a = a_constraint)
            }        
            W_alpha <- W %*% alpha
        }
        
        ##
        ## sample spatial process variance tau2
        ##
        
        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            
            ## double check this full conditional
            for (j in 1:(J-1)) {
                for (m in 1:M) {
                    devs       <- alpha[dims_idx == m, j]
                    SS         <- as.numeric(devs %*% (Q_alpha[[m]] %*% devs))
                    tau2[m, j] <- rgamma(1, alpha_tau2 + n_dims[m] / 2, beta_tau2 + SS / 2)
                }
            }
            for (j in 1:(J-1)) {
                Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j], use_spam = use_spam)
            }
        }
        
        ##
        ## sample eta
        ##
        
        if (sample_eta) {
            if (verbose)
                message("sample eta")
            
            ## double check this full conditional
            eta <- sapply(1:(J-1), function(j) {
                sigma2_tilde <- 1 / (1 / sigma2[j] + omega[, j])
                mu_tilde     <- 1 / sigma2[j] * (Xbeta[, j] + W_alpha[, j]) + kappa[, j]
                
                return(
                    rnorm(
                        N, 
                        sigma2_tilde * mu_tilde,
                        sqrt(sigma2_tilde)
                    )
                )
            })
        }
        
        ## 
        ## sample sigma2
        ## 
        
        if (sample_sigma2) {
            if (verbose)
                message("sample sigma2")
            
            for (j in 1:(J-1)) {
                SS <- sum((eta[, j] - Xbeta[, j] - W_alpha[, j])^2)
                sigma2[j] <- 1 / rgamma(1, alpha_sigma2 + N / 2, beta_sigma2 + SS / 2)
            }
        }      
        
        ##
        ## save variables
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                 <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ]  <- beta
                tau2_save[save_idx, , ]  <- tau2
                alpha_save[save_idx, , ] <- alpha
                eta_save[save_idx, , ]   <- eta
                sigma2_save[save_idx, ]  <- sigma2
                if (save_omega) {
                    omega_save[save_idx, , ] <- omega
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
    
    if (progress) {
        close(progressBar)
    }
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    out <- NULL
    if (save_omega) {
        out <- list(
            beta   = beta_save,
            tau2   = tau2_save,
            alpha  = alpha_save,
            eta    = eta_save,
            sigma2 = sigma2_save,
            MRA    = MRA,
            omega = omega_save)        
    } else {
        out <- list(
            beta   = beta_save,
            tau2   = tau2_save,
            alpha  = alpha_save,
            eta    = eta_save,
            sigma2 = sigma2_save,
            MRA    = MRA)
    }

    class(out) <- "pg_splm_mra"
    
    return(out)
    
}


