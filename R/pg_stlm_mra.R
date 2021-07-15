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
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param M The number of resolutions.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param inits is the list of intial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#' @importFrom stats rmultinom
#' @importFrom truncnorm rtruncnorm
#' @importFrom hms as_hms
#' @import BayesMRA
#' @import spam  
#' @export

## polya-gamma spatial linear regression model
pg_stlm_mra <- function(
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
    use_spam      = TRUE ## use spam or Matrix for sparse matrix operations
) {
    
    start <- Sys.time()
    
    ##
    ## Run error checks
    ## 
    
    check_input_pg_stlm(Y, X, locs)
    check_params(params)
    
    if (!use_spam)
        stop("The only sparse matrix pacakage available is spam")
    if (!is_positive_integer(n_cores, 1))
        stop("n_cores must be a positive integer")

    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    

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
    
    ## do we sample the climate variance parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
        }
    }
    
    ## do we sample the latent spatial random effect alpha
    sample_alpha <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_alpha']])) {
            sample_alpha <- config[['sample_alpha']]
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
    tX <- t(X)
    tXX <- tX %*% X
    
    ## Add in a counter for the number of regularized Cholesky factors.
    ## This is useful in correcting for numerical errors resulting in 
    ## covariance matrices that are not full rank
    num_chol_failures <- 0

    
    ## We assume a partially missing observation is the same as 
    ## fully missing. The index allows for fast accesing of missing
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
    
    sigma2 <- pmin(rgamma(J-1, alpha_sigma2, beta_sigma2), 5)

    ##
    ## initialized temporal autocorrelation
    ##
    rho <- runif(J-1, 0, 1)
    if (!is.null(inits[['rho']])) {
        if (all(!is.na(inits[['rho']]))) {
            if (!is_numeric_vector(inits[['rho']], J-1))
                stop ("If specified, inits$rho must be a vector of length J-1 with values between -1 and 1.")
            if (any(rho > 1) | any(rho < -1))
                stop ("If specified, inits$rho must be a vector of length J-1 with values between -1 and 1.")
            ## if rho passes error checks
            rho <- inits[['rho']]
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
    
    alpha_tau2 <- 1
    beta_tau2  <- 1
    
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
            if (is_positive_numeric_matrix(inits[['tau2']], M, J-1))
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
    
    ## define the sum-to-0 constraint for alpha_1
    # eventually modify this so the options for constraint and joint are allowed
    constraints <- make_constraint(MRA, constraint = "resolution", joint = TRUE)
    A_constraint <- constraints$A_constraint
    a_constraint <- constraints$a_constraint
    
    alpha <- array(0, dim = c(sum(n_dims), J-1, n_time))
    eta <- kappa   ## default initial value based on data Y to get started
    if (use_spam) {
        for (j in 1:(J-1)) {
            A_alpha_1  <- 1 / sigma2[j] * tWW + (1 + rho[j]^2) * Q_alpha_tau2[[j]]
            b_alpha    <- 1 / sigma2[j] * tW %*% (eta[, j, 1] - Xbeta[, j]) 
            alpha[, j, 1] <- rmvnorm.canonical.const(1, b_alpha, A_alpha_1, 
                                                  A = A_constraint, a = a_constraint)
        }
        for (tt in 2:n_time) {
            for (j in 1:(J-1)) {
                ## initialize with the current random effect
                alpha[, j, tt] <- alpha[, j, 1]
            }
        }
    } else {
        stop("The only sparse matrix pacakage available is spam")
        # # alpha[, 1] <- as.vector(solve(W %*% Q_alpha_tau2 %*% tW) %*% (tW %*% (Q_alpha_tau2 %*% (Z0 - Xgamma))))
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
    
    
    ## intialize an ICAR structure for fitting alpha
    Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(1, length(n_dims)), use_spam = use_spam)
    for (j in 1:(J-1)) {
        Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j], use_spam = use_spam)
    }
    
    
    ## precalculate the sparse cholesky structure for faster Gibbs updates
    Rstruct        <- NULL
    Rstruct_1      <- NULL
    Rstruct_n_time <- NULL
    if (use_spam) {
        A1        <- 1 / sigma2[1] * tWW + (1 + rho[1]^2) * Q_alpha_tau2[[1]]
        Rstruct_1 <- chol(A1)
        
        A <- 1 / sigma2[1] * tWW + (1 + rho[1]^2) * Q_alpha_tau2[[1]]
        Rstruct <- chol(A)
        
        A_n_time <- 1 / sigma2[1] * tWW + Q_alpha_tau2[[1]]
        Rstruct_n_time <- chol(A_n_time)
    }
    
    W_alpha <- array(0, dim = c(N, J-1, n_time))
    for (tt in 1:n_time) {
        W_alpha[, , tt] <- W %*% alpha[, , tt]
    }
    
 
    ##
    ## initialize the latent random process eta
    ##
    
    eta  <- array(0, dim = c(N, J-1, n_time))
    for (tt in 1:n_time) {
        eta[, , tt] <- Xbeta + W_alpha[, , tt] + sapply(1:(J-1), function(j) rnorm(N, 0, sqrt(sigma2[j])))
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
    omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    
    ##
    ## setup save variables
    ##
    
    n_save      <- params$n_mcmc / params$n_thin
    beta_save   <- array(0, dim = c(n_save, p, J-1))
    beta_save   <- array(0, dim = c(n_save, p, J-1))
    tau2_save   <- array(0, dim = c(n_save, M, J-1))
    alpha_save  <- array(0, dim = c(n_save, sum(n_dims), J-1, n_time))
    ## check memory requirements for saved objects
    eta_save      <- NULL
    eta_save_mean <- FALSE
    n_save        <- params$n_mcmc / params$n_thin
    if (!is.null(config)) {
        if (!is.null(config[['eta_save_mean']])) {
            eta_save_mean <- config[['eta_save_mean']]
        }
    }
    if (eta_save_mean) {
        eta_save      <- tryCatch(
            array(0, dim = c(N, J-1, n_time)),
            error = function(e) {
                stop('The memory needed to save the mean of the eta parameter is too large. Please contact the package maintainer for solutions.')
            }
        )
        
    } else {
        eta_save      <- tryCatch(
            array(0, dim = c(n_save, N, J-1, n_time)),
            error = function(e) {
                stop('The memory needed to save the eta parameter is too large. Either set config$save_eta_mean = TRUE or increase "params$n_thin')
            }
        )
    }
    
    pi_save     <- array(0, dim = c(n_save, N, J, n_time))
    sigma2_save <- matrix(0, n_save, J-1)
    rho_save    <- matrix(0, n_save, J-1)

    
    ##
    ## initialize tuning variables for adaptive MCMC
    ##
    
    ## tuning for rho
    rho_accept       <- rep(0, J-1)
    rho_accept_batch <- rep(0, J-1)
    rho_tune         <- rep(0.025, J-1)
    
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
        
        omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        
        ##
        ## sample beta
        ##
        
        ## can parallelize this update -- each group of parameters is 
        ## conditionally independent given omega and kappa(y)
        
        if (sample_beta){
            if (verbose)
                message("sample beta")
            
            for (j in 1:(J-1)) {
                A         <- n_time * tXX / sigma2[j] + Sigma_beta_inv
                A         <- (A + t(A)) / 2    ## guarantee a symmetric matrix
                b         <- rowSums(tX %*% (eta[, j, ] - W_alpha[, j, ])) / sigma2[j] + Sigma_beta_inv %*% mu_beta
                beta[, j] <- rmvn_arma(A, b)
            }
            Xbeta <- X %*% beta
        }
        
        ##
        ## sample spatial random effects alpha
        ##

        if (sample_alpha) {
            if (verbose)
                message("sample alpha")
            
            ## double check this full conditional
            A_alpha_1      <- vector(mode = "list", length = J-1)
            A_alpha        <- vector(mode = "list", length = J-1)
            A_alpha_n_time <- vector(mode = "list", length = J-1)
            for (j in 1:(J-1)) {
                A_alpha_1[[j]]      <- 1 / sigma2[j] * tWW + (1 + rho[j]^2) * Q_alpha_tau2[[j]]
                A_alpha[[j]]        <- 1 / sigma2[j] * tWW + (1 + rho[j]^2) * Q_alpha_tau2[[j]]
                A_alpha_n_time[[j]] <- 1 / sigma2[j] * tWW + Q_alpha_tau2[[j]]
            }
            
            ## parallelize this
            for(tt in 1:n_time) {
                for (j in 1:(J-1)) {
                    if (tt == 1) {
                        b_alpha <- 1 / sigma2[j] * tW %*% (eta[, j, 1] - Xbeta[, j]) +
                            Q_alpha_tau2[[j]] %*% as.vector(rho[j] * alpha[, j, 2])
                        alpha[, j, 1] <- tryCatch(
                            as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_1[[j]], Rstruct = Rstruct_1, A = A_constraint, a = a_constraint)),
                            error = function(e) {
                                if (verbose)
                                    message("The Cholesky decomposition conditional precision for alpha_1 was ill-conditioned and mildy regularized.")
                                num_chol_failures <- num_chol_failures + 1
                                A_alpha_1[[j]] <<- A_alpha_1[[j]] + 1e-8 * diag(sum(n_dims))
                                return(as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_1[[j]], Rstruct = Rstruct_1, A = A_constraint, a = a_constraint)))
                            })
                        
                    } else if (tt == n_time) {
                        b_alpha <- 1 / sigma2[j] * tW %*% (eta[, j, n_time] - Xbeta[, j]) +
                            Q_alpha_tau2[[j]] %*% as.vector(rho[j] * alpha[, j, n_time - 1])
                        alpha[, j, tt] <- tryCatch(
                            # as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_n_time[[j]], Rstruct = Rstruct_n_time, A = A_constraint, a = a_constraint)),
                            as.vector(rmvnorm.canonical(1, b_alpha, A_alpha_n_time[[j]], Rstruct = Rstruct_n_time)),
                            error = function(e) {
                                if (verbose)
                                    message("The Cholesky decomposition conditional precision for alpha_n_time was ill-conditioned and mildy regularized.")
                                num_chol_failures <- num_chol_failures + 1
                                A_alpha_n_time[[j]] <<- A_alpha_n_time[[j]] + 1e-8 * diag(sum(n_dims))
                                # return(as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_n_time[[j]], Rstruct = Rstruct_n_time, A = A_constraint, a = a_constraint)))
                                return(as.vector(rmvnorm.canonical(1, b_alpha, A_alpha_n_time[[j]], Rstruct = Rstruct_n_time)))
                            })
                        
                    } else {
                        b_alpha <- 1 / sigma2[j] * tW %*% (eta[, j, tt] - Xbeta[, j]) +
                            Q_alpha_tau2[[j]] %*% as.vector(rho[j] * (alpha[, j, tt - 1] + alpha[, j, tt + 1]))
                        alpha[, j, tt] <- tryCatch(
                            # as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha[[j]], Rstruct = Rstruct, A = A_constraint, a = a_constraint)),
                            as.vector(rmvnorm.canonical(1, b_alpha, A_alpha[[j]], Rstruct = Rstruct)),
                            error = function(e) {
                                if (verbose)
                                    message("The Cholesky decomposition conditional precision for alpha_t was ill-conditioned and mildy regularized.")
                                num_chol_failures <- num_chol_failures + 1
                                A_alpha[[j]] <<- A_alpha[[j]] + 1e-8 * diag(sum(n_dims))
                                # return(as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha[[j]], Rstruct = Rstruct, A = A_constraint, a = a_constraint)))
                                return(as.vector(rmvnorm.canonical(1, b_alpha, A_alpha[[j]], Rstruct = Rstruct)))
                            })
                    }       
                }
                ## update the latent process variable
                W_alpha[, , tt] <- W %*% alpha[, , tt]
            }
        }
        
        ##
        ## sample spatial process variance tau2
        ##
        
        ## double check this full conditional
        
        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            for (j in 1:(J-1)) {
                for (m in 1:M) {
                    devs <- cbind(
                        alpha[dims_idx == m, j, 1],
                        alpha[dims_idx == m, j, -1] - rho[j] * alpha[dims_idx == m, j, -n_time]
                    )
                    SS       <- sum(sapply(1:n_time, function(tt) devs[, tt] %*% (Q_alpha[[m]] %*% devs[, tt])))
                    tau2[m, j]  <- rgamma(1, alpha_tau2 + n_dims[m] * n_time / 2, beta_tau2 + SS / 2)
                }
                ## moved this inside the loop to update this parameter
                Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j], use_spam = use_spam)
            }
            ##  this was outside the loop for tau2
            # Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j], use_spam = use_spam)
        }        

        ##
        ## sample eta
        ##
        
        if (sample_eta) {
            if (verbose)
                message("sample eta")
            for (tt in 1:n_time) {
                ## double check this full conditional
                eta[, , tt] <- sapply(1:(J-1), function(j) {
                    sigma2_tilde <- 1 / (1 / sigma2[j] + omega[, j, tt])
                    mu_tilde     <- 1 / sigma2[j] * (Xbeta[, j] + W_alpha[, j, tt]) + kappa[, j, tt]
                    
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
        
        ##
        ## sample rho
        ##
        
        if (sample_rho) {
            if (verbose)
                message("sample rho")
            
            # rho_star <- rnorm(1, rho, rho_tune)
            # if (rho_star < 1 & rho_star > -1) {
            #     mh1 <- NULL
            #     mh2 <- NULL
            #     if (shared_covariance_params){
            #         mh1 <- sum(
            #             sapply(2:n_time, function(tt) {
            #                 sapply(1:(J-1), function(j) {
            #                     mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho_star * eta[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
            #                 })
            #             }) 
            #             
            #         ) 
            #         ## parallelize this        
            #         mh2 <- sum(
            #             sapply(2:n_time, function(tt) {
            #                 sapply(1:(J-1), function(j) {
            #                     mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
            #                 })
            #             }) 
            #         ) 
            #     } else {
            #         mh1 <- sum(
            #             sapply(2:n_time, function(tt) {
            #                 sapply(1:(J-1), function(j) {
            #                     mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho_star * eta[, j, tt - 1], Sigma_chol[j,,], isChol = TRUE, log = TRUE, ncores = n_cores) 
            #                 })
            #             }) 
            #             
            #         ) 
            #         ## parallelize this        
            #         mh2 <- sum(
            #             sapply(2:n_time, function(tt) {
            #                 sapply(1:(J-1), function(j) {
            #                     mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol[j,,], isChol = TRUE, log = TRUE, ncores = n_cores) 
            #                 })
            #             }) 
            #         )
            #     }
            #     
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
            
            for (j in 1:(J-1)) {
                rho_vals <- rowSums(
                    sapply(2:n_time, function(tt) {
                        t_alpha_Q <- t(alpha[, j, tt-1]) %*% Q_alpha_tau2[[j]]
                        c(
                            t_alpha_Q %*% alpha[, j, tt-1],
                            t_alpha_Q %*% alpha[, j, tt]
                        )
                    })
                )
                
                a_rho  <- rho_vals[1]
                b_rho  <- rho_vals[2]
                rho[j] <- rtruncnorm(1, a = -1, b = 1, mean = b_rho / a_rho, sd = sqrt(1 / a_rho))
            }
        }
        
        ## 
        ## sample sigma2
        ## 
        
        if (verbose)
            message("sample sigma2")
        
        for (j in 1:(J-1)) {
            SS <- sum(sapply(1:n_time, function(tt) (eta[, j, tt] - Xbeta[, j] - W_alpha[, j, tt])^2))
            sigma2[j] <- 1 / rgamma(1, alpha_sigma2 + N * n_time / 2, beta_sigma2 + SS / 2)
        }
        
        ##
        ## save variables
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                   <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ]    <- beta
                tau2_save[save_idx, , ]    <- tau2
                alpha_save[save_idx, , , ] <- alpha
                if (eta_save_mean) {
                    eta_save[, , , ] <- 1 / n_save * eta
                } else {
                    eta_save[save_idx, , , ] <- eta
                }
                sigma2_save[save_idx, ]    <- sigma2
                for (tt in 1:n_time) {
                    pi_save[save_idx, , , tt]  <- eta_to_pi(eta[, , tt])
                }
                rho_save[save_idx, ]       <- rho
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
    # message("Acceptance rate for rho is ", mean(rho_accept))
    
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    if (progress) {
        close(progressBar)
    }
    
    stop    <- Sys.time()
    runtime <- stop - start
    
    message("MCMC took ", hms::as_hms(runtime))
    
    
    out <- list(
        beta   = beta_save,
        alpha  = alpha_save,
        tau2   = tau2_save,
        eta    = eta_save,
        pi     = pi_save,
        sigma2 = sigma2_save,
        rho    = rho_save,
        MRA    = MRA
    )
    class(out) <- "pg_stlm_mra"
    
    return(out)
}


