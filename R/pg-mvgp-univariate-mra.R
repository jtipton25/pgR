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
#' @param M The number of resolutions.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. \code{n_coarse_grid = 10} results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by \code{n_padding}.
#' @param model is the form of the  polya-gamma model. Currently, this option is not active the only model is the "iid error" model. This option allows for independent species-specific overdispersion variance terms.
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param inits is the list of initial values if the user wishes to specify initial values. If these values are not specified, then the initial values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#' @param n_chain is the MCMC chain id. The default is 1.
#' 
#' @export
#' 
#' @importFrom LaplacesDemon rinvwishart rtrunc
#' @importFrom stats rmultinom lm
#' @importFrom hms as_hms
#' @importFrom fields rdist
#' @importFrom sparseMVN rmvn.sparse
#' @importFrom Matrix Cholesky
#' @importFrom MASS ginv
#' @import BayesMRA
#' @import spam  

## polya-gamma spatial linear regression model
pg_mvgp_univariate_mra <- function(
    Y, 
    X,
    Z0,
    locs, 
    params,
    priors,
    M = 4,
    n_coarse_grid = 10,
    # n_max_fine_grid = 2^6,
    model    = "iid error",
    n_cores = 1L,
    inits   = NULL,
    config  = NULL,
    verbose = FALSE,
    use_spam = TRUE, ## use spam or Matrix for sparse matrix operations
    n_chain       = 1
) {
    
    
    # RSR <- FALSE
    
    ##
    ## Run error checks
    ## 
    
    # check_input_spatial(Y, X, locs)
    # check_params(param)
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    if(!is.vector(Z0))
        stop("Z must be a vector of climate variable inputs")
    if (!is_positive_integer(n_cores, 1))
        stop("n_cores must be a positive integer")
    
    N      <- nrow(Y)
    J      <- ncol(Y)
    n_time <- dim(Y)[3]
    p      <- ncol(X)
    Q      <- 1
    D      <- rdist(locs)
    
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
    
    ## do we only use the modern climate data to estimate the functional 
    ## relationship or do we use the estimated climate as well
    sample_beta_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta_modern']])) {
            sample_beta_modern <- config[['sample_beta_modern']]
        }
    }
    
    ## do we sample the climate fixed effect parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_gamma <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_gamma']])) {
            sample_gamma <- config[['sample_gamma']]
        }
    }
    
    
    ## do we only use the modern climate data to estimate the fixed effects 
    ## or do we use the estimated climate as well
    sample_gamma_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_gamma_modern']])) {
            sample_gamma_modern <- config[['sample_gamma_modern']]
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
    
    ## do we sample the climate variance parameter using only the 
    ## modern data? 
    sample_tau2_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2_modern']])) {
            sample_tau2_modern <- config[['sample_tau2_modern']]
        }
    }
    
    ## do we sample the overdispersion parameter
    sample_sigma2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2']])) {
            sample_sigma2 <- config[['sample_sigma2']]
        }
    }
    
    ## do we sample the overdispersion parameter only for the modern data?
    sample_sigma2_modern <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2_modern']])) {
            sample_sigma2_modern <- config[['sample_sigma2_modern']]
        }
    }
    
    ## do we sample the latent intensity parameter alpha
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
    
    ## do we sample the climate observation error parameter
    sample_sigma2_0 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_sigma2_0']])) {
            sample_sigma2_0 <- config[['sample_sigma2_0']]
        }
    }
    
    ## do we sample the climate time varying mean parameter
    sample_mu <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_mu']])) {
            sample_mu <- config[['sample_mu']]
        }
    }
    
    
    ## add in a counter for the number of regularized Cholesky
    num_chol_failures <- 0
    

    ## Restricted spatial regression?
    # IMPX <- diag(N)
    # if (RSR) {
    #     IMPX <- diag(N) -  X %*% chol2inv(chol(t(X) %*% X)) %*% t(X)
    # }
    
    
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
    ## initial values
    ##
    
    tX  <- t(X)
    tXX <- tX %*% X
    
    ##
    ## priors for beta
    ##
    
    mu_beta        <- rep(0, Q + 1)
    Sigma_beta     <- 10 * diag(Q + 1)
    
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
            chol(Sigma_beta + 1e-8 * diag(Q + 1))                    
        }
    )
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##
    
    beta <- t(mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    
    ##
    ## priors for gamma
    ##
    
    mu_gamma        <- rep(0, p)
    Sigma_gamma     <- 10 * diag(p)
    
    ## check if priors for mu_gamma are specified
    if (!is.null(priors[['mu_gamma']])) {
        if (all(!is.na(priors[['mu_gamma']]))) {
            mu_gamma <- priors[['mu_gamma']]
        }
    }
    
    ## check if priors for Sigma_gamma are specified
    if (!is.null(priors[['Sigma_gamma']])) {
        if (all(!is.na(priors[['Sigma_gamma']]))) {
            Sigma_gamma <- priors[['Sigma_gamma']]
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
    tau2 <- 100 * pmax(rgamma(M, alpha_tau2, beta_tau2), 1)
    
    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)
    
    
    ##
    ## initialize gamma
    ##
    
    ## using sparse Cholesky routines -- only used for initialization
    # Sigma_init_inv <- 1 / sigma2_0 * (
    #     diag(N) - W %*% chol2inv(chol(tWW + sigma2_0 * Q_alpha_tau2)) %*% t(W)
    # )
    
    # gamma <- t(mvnfast::rmvn(Q, mu_gamma, Sigma_gamma_chol, isChol = TRUE))
    # gamma <- solve(t(X) %*% Sigma_init_inv %*% X) %*% t(X) %*% Sigma_init_inv %*% Z0
    gamma <- lm(Z0 ~ X-1)$coeff
    
    Xgamma <- as.vector(X %*% gamma)
    
    
    
    ##
    ## initialize rho
    ##
    
    rho      <- runif(1, -1, 1)
    
    ##
    ## priors for sigma2
    ##
    
    alpha_sigma2 <- 1
    beta_sigma2  <- 1
    
    ## check if priors for alpha_sigma2 are specified
    if (!is.null(priors[['alpha_sigma2']])) {
        alpha_sigma2 <- priors[['alpha_sigma2']]
    }
    
    ## check if priors for beta_sigma2 are specified
    if (!is.null(priors[['beta_sigma2']])) {
        beta_sigma2 <- priors[['beta_sigma2']]
    }
    
    ##
    ## initialize sigma2
    ##
    
    sigma2  <- rgamma(J-1, alpha_sigma2, beta_sigma2)
    sigma   <- sqrt(sigma2)
    
    ##
    ## priors for sigma2_0
    ##
    
    ## default strong shrinkage prior towards 0
    alpha_sigma2_0 <- 200
    beta_sigma2_0  <- 1
    
    ## check if priors for alpha_sigma2 are specified
    if (!is.null(priors[['alpha_sigma2']])) {
        alpha_sigma2 <- priors[['alpha_sigma2']]
    }
    
    ## check if priors for beta_sigma2 are specified
    if (!is.null(priors[['beta_sigma2']])) {
        beta_sigma2 <- priors[['beta_sigma2']]
    }
    
    ##
    ## intialize sigma2_0
    ##
    
    sigma2_0 <- pmin(1 / rgamma(1, alpha_sigma2_0, beta_sigma2_0), 0.1)
    
    ##
    ## initialize alpha
    ##
    
    ## define the sum-to-0 constraint for alpha_1
    A_constraint <- t(
        sapply(1:M, function(m) {
            tmp <- rep(0, sum(n_dims))
            tmp[dims_idx == m] <- 1
            return(tmp)
        })
    )
    a_constraint <- rep(0, M)
    
    alpha <- matrix(0, sum(n_dims), n_time)
    if (use_spam) {
        A_alpha_1 <- 1 / sigma2_0 * tWW + (1 + rho^2) * Q_alpha_tau2
        b_alpha   <- 1 / sigma2_0 * tW %*% (Z0 - Xgamma) 
        alpha[, 1] <- rmvnorm.canonical.const(1, b_alpha, A_alpha_1, 
                                              A = A_constraint, a = a_constraint)
        for (tt in 2:n_time) {
            ## initialize with the current climate
            alpha[, tt] <- alpha[, 1]
            ## intialize with the prior process
            # alpha[, tt] <- as.vector(rmvnorm.prec(1, rho * alpha[, tt-1], Q_alpha_tau2))
        }
    } else {
        stop("The only sparse matrix pacakage available is spam")
        # alpha[, 1] <- as.vector(solve(W %*% Q_alpha_tau2 %*% tW) %*% (tW %*% (Q_alpha_tau2 %*% (Z0 - Xgamma))))
        alpha[, 1] <- as.vector(chol2inv(chol(W %*% Q_alpha_tau2 %*% tW)) %*% (tW %*% (Q_alpha_tau2 %*% (Z0 - Xgamma))))
        # alpha[, 1] <- as.vector(ginv(as.matrix(tWW)) %*% tW %*% (Z0 - Xgamma))
        # alpha[, 1] <- as.vector(rmvn.sparse(1, rep(0, sum(n_dims)), CH = Cholesky(Q_alpha_tau2), prec = TRUE))     
        for (tt in 2:n_time) {
            ## initialize with the current climate
            alpha[, tt] <- alpha[, 1]
            ## intialize with the prior process
            # alpha[, tt] <- as.vector(rmvn.sparse(1, rho * alpha[, tt - 1], CH = Cholesky(Q_alpha_tau2), prec = TRUE))
        }
    }
    
    
    ## intialize an ICAR structure for fitting alpha
    Q_alpha      <- make_Q_alpha_2d(sqrt(n_dims), rep(1, length(n_dims)), use_spam = use_spam)
    Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)
    
    ## precalculate the sparse cholesky structure for faster Gibbs updates
    Rstruct        <- NULL
    Rstruct_1      <- NULL
    Rstruct_n_time <- NULL
    if (use_spam) {
        A1        <- 1 / sigma2_0 * tWW + sum(beta[2, ]^2 / sigma2) * tWW + (1 + rho^2) * Q_alpha_tau2
        Rstruct_1 <- chol(A1)
        
        A <- sum(beta[2, ]^2 / sigma2)  * tWW + (1 + rho^2) * Q_alpha_tau2
        Rstruct <- chol(A)
        
        A_n_time <- sum(beta[2, ]^2 / sigma2)  * tWW + Q_alpha_tau2
        Rstruct_n_time <- chol(A_n_time)
    }
    
    W_alpha <- NULL
    W_alpha <- W %*% alpha
    
    ## initialize the time varying mean
    ## note: we assume the intial state has mu[1] = 0
    mu        <- rep(0, n_time)
    sigma2_mu <- 0.25
    
    ## check if prior value of sigma2_mu is specified
    if (!is.null(priors[['sigma2_mu']])) {
        sigma2_mu <- priors[['sigma2_mu']]
    }
    
    
    ## initialize Z and eta
    Z   <- matrix(0, N, n_time)
    eta <- array(0, dim = c(N, J-1, n_time))
    
    for (tt in 1:n_time) {
        if (tt == 1) {
            Z[, tt]     <- Z0
        } else {
            Z[, tt]   <- Xgamma + W_alpha[, tt] + mu[tt] * rep(1, N)
        }
        eta[, , tt] <- cbind(1, Z[, tt]) %*% beta +
            sapply(1:(J-1), function(j) rnorm(N, 0, sigma[j]))
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
    omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    
    
    ##
    ## check for initial values
    ##
    
    ## initial values for beta 
    if (!is.null(inits[['beta']])) {
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
            if(!is_numeric_matrix(beta, Q + 1, J - 1))
                stop("initial value for beta must be a Q+1 by J-1 numeric matrix")
        }
    }
    
    ## intial values for sigma2
    if (!is.null(inits[['sigma2']])) {
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
            if(!is_positive_numeric(sigma2, J - 1))
                stop("initial value for sigma2 must be a J-1 vector of positive numeric values")
        }
    }
    
    ## initial values for gamma
    if (!is.null(inits[['gamma']])) {
        if (all(!is.na(inits[['gamma']]))) {
            gamma <- inits[['gamma']]
        }
    }
    
    Xgamma <- as.vector(X %*% gamma)
    
    ## initial values for alpha
    if (!is.null(inits[['alpha']])) {
        if (all(!is.na(inits[['alpha']]))) {
            alpha <- inits[['alpha']]
        }
    }
    
    
    ## initial values for omega
    if (!is.null(inits[['omega']])) {
        if (!is.na(inits[['omega']])) {
            omega <- inits[['omega']]
        }
    }
    
    ## initial values for Z
    if (!is.null(inits[['Z']])) {
        if (all(!is.na(inits[['Z']]))) {
            Z <- inits[['Z']]
        }
    }
    
    ## initial values for eta
    if (!is.null(inits[['eta']])) {
        if (all(!is.na(inits[['eta']]))) {
            eta <- inits[['eta']]
        }
    }
    
    ## initial values for tau2
    if (!is.null(inits[['tau2']])) {
        if (all(!is.na(inits[['tau2']]))) {
            tau2 <- inits[['tau2']]
        }
    }
    
    ## initial values for rho
    if (!is.null(inits[['rho']])) {
        if (all(!is.na(inits[['rho']]))) {
            rho <- inits[['rho']]
        }
    }
    ## initial values for sigma2
    if (!is.null(inits[['sigma2']])) {
        if (all(!is.na(inits[['sigma2']]))) {
            sigma2 <- inits[['sigma2']]
        }
    }
    
    ## initial values for sigma2_0
    if (!is.null(inits[['sigma2_0']])) {
        if (all(!is.na(inits[['sigma2_0']]))) {
            sigma2_0 <- inits[['sigma2_0']]
        }
    }
    
    ## initial values for mu
    if (!is.null(inits[['mu']])) {
        if (all(!is.na(inits[['mu']]))) {
            mu <- inits[['mu']]
        }
    }
    ## make sure the first value is zero
    if (mu[1] != 0) {
        message("The value for the first value mu[1] must be 0. The parameter mu[1] has been set to 0.")
        mu[1] <- 0
    }
    
    
    ##
    ## setup save variables
    ##
    
    n_save        <- params$n_mcmc / params$n_thin
    beta_save     <- array(0, dim = c(n_save, Q+1, J-1))
    gamma_save    <- matrix(0, n_save, p)
    rho_save      <- rep(0, n_save)
    tau2_save     <- matrix(0, n_save, M)
    sigma2_save   <- matrix(0, n_save, J-1)
    alpha_save    <- array(0, dim = c(n_save, sum(n_dims), n_time))
    Z_save        <- array(0, dim = c(n_save, N, n_time))
    sigma2_0_save <- rep(0, n_save)
    mu_save       <- matrix(0, n_save, n_time)
    
    
    ## 
    ## initialize tuning 
    ##
    
    ##
    ## tuning variables for adaptive MCMC
    ##
    
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
        
        omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        
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
        ## sample sigma2
        ##
        
        if (sample_sigma2) {
            if (verbose)
                message("sample sigma2")
            
            for (j in 1:(J-1)) {
                if (sample_sigma2_modern) {
                    devs      <- eta[, j, 1] - cbind(1, Z[, 1]) %*% beta[, j]
                    SS        <- sum(devs^2)
                    sigma2[j] <- 1 / rgamma(1, N / 2 + alpha_sigma2, SS / 2 + beta_sigma2) 
                } else {
                    devs      <- sapply(1:n_time, function(tt) eta[, j, tt] - cbind(1, Z[, tt]) %*% beta[, j])
                    SS        <- sum(devs^2)
                    sigma2[j] <- 1 / rgamma(1, N * n_time / 2 + alpha_sigma2, SS / 2 + beta_sigma2) 
                }
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
                b <- 1 / sigma2[j] * t(cbind(1, Z[, 1])) %*% eta[, j, 1] + Sigma_beta_inv %*% mu_beta
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
            A <- tXX / sigma2_0 + tXX * sum(beta[2, ]^2 / sigma2) + Sigma_gamma_inv
            b <- tX %*% (Z[, 1] - W_alpha[, 1]) / sigma2_0 +
                tX %*% rowSums(
                    sapply(1:(J-1), function (j) {
                        beta[2, j] / sigma2[j] * (eta[, j, 1] - beta[1, j] * rep(1, N) - beta[2, j] * W_alpha[, 1])
                        # beta[1, j] / sigma2[j] * (eta[, j, 1] - beta[1, j] * rep(1, N) - W_alpha[, 1])
                    }) 
                ) + Sigma_gamma_inv %*% mu_gamma
            if (!sample_gamma_modern) {
                A <- A + (n_time - 1) *  tXX * sum(beta[2, ]^2 / sigma2)
                b <- b + rowSums(
                    sapply(2:n_time, function(tt) {
                        tX %*% rowSums(
                            sapply(1:(J-1), function (j) {
                                # beta[1, j] / sigma2[j] * (eta[, j, tt] - beta[1, j] * rep(1, N) - W_alpha[, tt])
                                beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j] * rep(1, N) - beta[2, j] * W_alpha[, tt] - beta[2, j] * mu[tt] * rep(1, N))
                            }) 
                        )                
                    })
                )
            }
            ## guarantee a symmetric matrix
            A      <- (A + t(A)) / 2
            gamma  <- rmvn_arma(A, b)
            Xgamma <- as.vector(X %*% gamma)
        }
        
        ## update Z
        for (tt in 1:n_time) {
            if (tt == 1) {
                Z[, tt]     <- Z0
            } else {
                Z[, tt]   <- Xgamma + W_alpha[, tt] + mu[tt] * rep(1, N)
            }
        }
        
        ##
        ## sample alpha
        ##
        
        if (sample_alpha) {
            if (verbose)
                message("sample alpha")
            
            
            A_alpha_1      <- 1 / sigma2_0 * tWW + sum(beta[2, ]^2 / sigma2) * tWW + (1 + rho^2) * Q_alpha_tau2
            A_alpha_n_time <- sum(beta[2, ]^2 / sigma2) * tWW + Q_alpha_tau2
            A_alpha        <- sum(beta[2, ]^2 / sigma2) * tWW + (1 + rho^2) * Q_alpha_tau2
            
            ## parallelize this
            for(tt in 1:n_time) {
                if (tt == 1) {
                    b_alpha <- 1 / sigma2_0 * tW %*% (Z[, 1] - Xgamma - mu[tt] * rep(1, N)) +
                        tW %*% rowSums(sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, 1] - beta[1, j] * rep(1, N) - beta[2, j] * Xgamma - beta[2, j] * mu[tt] * rep(1, N)))) +
                        Q_alpha_tau2 %*% as.vector(rho * alpha[, 2])
                    alpha[, 1] <- tryCatch(
                        as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_1, Rstruct = Rstruct_1, A = A_constraint, a = a_constraint)),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition conditional precision for alpha_1 was ill-conditioned and mildy regularized.")
                            num_chol_failures <- num_chol_failures + 1
                            A_alpha_1 <<- A_alpha_1 + 1e-8 * diag(sum(n_dims))
                            return(as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_1, Rstruct = Rstruct_1, A = A_constraint, a = a_constraint)))
                        })
                    
                } else if (tt == n_time) {
                    b_alpha <- tW %*% rowSums(sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j] * rep(1, N) - beta[2, j] * Xgamma - beta[2, j] * mu[tt] * rep(1, N)))) +
                        Q_alpha_tau2 %*% as.vector(rho * alpha[, tt - 1])
                    alpha[, tt] <- tryCatch(
                        as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_n_time, Rstruct = Rstruct_n_time, A = A_constraint, a = a_constraint)),
                        # as.vector(rmvnorm.canonical(1, b_alpha, A_alpha_n_time, Rstruct = Rstruct_n_time)),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition conditional precision for alpha_1 was ill-conditioned and mildy regularized.")
                            num_chol_failures <- num_chol_failures + 1
                            A_alpha_n_time <<- A_alpha_n_time + 1e-8 * diag(sum(n_dims))
                            return(as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha_n_time, Rstruct = Rstruct_n_time, A = A_constraint, a = a_constraint)))
                            # return(as.vector(rmvnorm.canonical(1, b_alpha, A_alpha_n_time, Rstruct = Rstruct_n_time)))
                        })
                    
                } else {
                    b_alpha <- tW %*% rowSums(sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j] * rep(1, N) - beta[2, j] * Xgamma - beta[2, j] * mu[tt] * rep(1, N))))  +
                        Q_alpha_tau2 %*% as.vector(rho * alpha[, tt - 1] + rho * alpha[, tt + 1])
                    alpha[, tt] <- tryCatch(
                        as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, Rstruct = Rstruct, A = A_constraint, a = a_constraint)),
                        # as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct)),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition conditional precision for alpha_1 was ill-conditioned and mildy regularized.")
                            num_chol_failures <- num_chol_failures + 1
                            A_alpha <<- A_alpha + 1e-8 * diag(sum(n_dims))
                            return(as.vector(rmvnorm.canonical.const(1, b_alpha, A_alpha, Rstruct = Rstruct, A = A_constraint, a = a_constraint)))
                            # return(as.vector(rmvnorm.canonical(1, b_alpha, A_alpha, Rstruct = Rstruct)))
                        })
                }       
            }
        }
        
        ## update W_alpha
        W_alpha <- W %*% alpha
        
        ## update Z
        for (tt in 1:n_time) {
            if (tt == 1) {
                Z[, tt]   <- Z0
            } else {
                Z[, tt]   <- Xgamma + W_alpha[, tt] + mu[tt] * rep(1, N)
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
                    t_alpha_Q <- t(alpha[, tt-1]) %*% Q_alpha_tau2
                    c(
                        t_alpha_Q %*% alpha[, tt-1],
                        t_alpha_Q %*% alpha[, tt]
                    )
                })
            )
            
            a_rho <- rho_vals[1]
            b_rho <- rho_vals[2]
            rho   <- rtrunc(1, "norm", a = -1, b = 1, mean = b_rho / a_rho, sd = sqrt(1 / a_rho))
        }
        
        
        ##
        ## sample tau2
        ##
        
        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            if (sample_tau2_modern) {
                for (m in 1:M) {
                    devs <- alpha[dims_idx == m, 1]
                    SS       <- as.numeric(devs %*% (Q_alpha[[m]] %*% devs))
                    tau2[m]  <- rgamma(1, alpha_tau2 + n_dims[m] / 2, beta_tau2 + SS / 2)
                }
            } else {
                for (m in 1:M) {
                    devs <- cbind(
                        alpha[dims_idx == m, 1],
                        alpha[dims_idx == m, - 1] - rho * alpha[dims_idx == m, - n_time]
                    )
                    
                    SS       <- sum(sapply(1:n_time, function(tt) devs[, tt] %*% (Q_alpha[[m]] %*% devs[, tt])))
                    tau2[m]  <- rgamma(1, alpha_tau2 + n_dims[m] * n_time / 2, alpha_tau2 + SS / 2)
                }
            }
            
            Q_alpha_tau2 <- make_Q_alpha_tau2(Q_alpha, tau2, use_spam = use_spam)
            tau        <- sqrt(tau2)
        }
        
        ##
        ## sample sigma2_0
        ##
        
        if (sample_sigma2_0) {
            if (verbose)
                message("sample sigma2_0")
            
            SS <- sum((Z0 - Xgamma - W_alpha[, 1])^2)
            sigma2_0 <- 1 / rgamma(1, alpha_sigma2_0 + N / 2, beta_sigma2_0 + SS / 2)
        }
        
        
        ## 
        ## sample mu
        ## 
        
        if (sample_mu) {
            if (verbose)
                message("sample mu")
            
            for (tt in 2:n_time) {
                if (tt == n_time) {
                    a <- sapply(1:(J-1), function(j) beta[2, j]^2 / sigma2[j] * N)  + 1 / sigma2_mu
                    b <- sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j] * rep(1, N) - beta[2, j] * Xgamma - beta[2, j] * W_alpha[, tt])) +
                        1 / sigma2_mu * mu[tt - 1]
                    mu[tt] <- rnorm(1, b / a, sqrt(1 / a))
                } else {
                    a <- sapply(1:(J-1), function(j) beta[2, j]^2 / sigma2[j] * N)  + 2 / sigma2_mu
                    b <- sapply(1:(J-1), function(j) beta[2, j] / sigma2[j] * (eta[, j, tt] - beta[1, j] * rep(1, N) - beta[2, j] * Xgamma - beta[2, j] * W_alpha[, tt])) +
                        1 / sigma2_mu * (mu[tt - 1] + mu[tt + 1])
                    mu[tt] <- rnorm(1, b / a, sqrt(1 / a))
                }
            }
        }
        
        
        ##
        ## save variables
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                 <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ]  <- beta
                gamma_save[save_idx, ]   <- gamma
                rho_save[save_idx]       <- rho
                tau2_save[save_idx, ]    <- tau2
                sigma2_save[save_idx, ]  <- sigma2
                alpha_save[save_idx, , ] <- alpha
                Z_save[save_idx, , ]     <- Z
                if (eta_save_mean) {
                    eta_save[, , , ] <- 1 / n_save * eta
                } else {
                    eta_save[save_idx, , , ] <- eta
                }
                sigma2_0_save[save_idx]  <- sigma2_0
                mu_save[save_idx, ]      <- mu
            }
        }
        
        ##
        ## End of MCMC loop
        ##
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    
    if (num_chol_failures > 0)
        warning("The Cholesky decomposition for the update of alpha was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    out <- list(
        beta     = beta_save,
        gamma    = gamma_save,
        rho      = rho_save,
        tau2     = tau2_save,
        sigma2   = sigma2_save,
        alpha    = alpha_save,
        eta      = eta_save,
        Z        = Z_save,
        W        = W,
        sigma2_0 = sigma2_0_save,
        mu       = mu_save,
        MRA      = MRA
    )
    
    class(out) <- "pg_mvgp_univariate_mra" 
    
    return(out)
}