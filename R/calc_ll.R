#' Log-likelihood for pg_lm() model for use in model selection
#' 
#' this function generates the log-linkelihood for data fit using the pg_lm() 
#' function. The log-likelihood can be used for model selection and evaluation.
#' @param Y is a \eqn{n \times J}{n x J} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param out is a list of MCMC outputs from `pg_lm()`
#' 
#' @importFrom stats dmultinom
#' @importFrom mgcv choldrop
#' 
#' @export 

calc_ll_pg_lm <- function(Y, X, out) {
    ##
    ## check the inputs 
    ##
    if (!inherits(out, "pg_lm"))
        stop("The MCMC object out must be of class pg_lm which is the output of the pg_lm() function.")
    
    check_input_pg_lm(Y, X)
    
    N  <- nrow(Y)
    J  <- ncol(Y)
    p  <- ncol(X)
    
    beta      <- out$beta
    n_samples <- nrow(beta)    
    eta <- array(0, dim = c(n_samples, N, J - 1))
   
    eta <- sapply(1:n_samples, function(i) X %*% beta[i, , ], simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)
    eta <- aperm(eta, c(3, 1, 2))
    
    ## convert from eta to pi
    pi <- sapply(1:n_samples, function(i) eta_to_pi(eta[i, , ]), simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)

    ## check in J=2
    if (length(dim(pi)) == 2) {
        p <- aperm(pi, c(2, 1))   
    } else {
        pi <- aperm(pi, c(3, 1, 2))
    }    
    
    ll <- matrix(0, N, n_samples)
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        for (i in 1:N) {
            ll[i, k] <- dmultinom(Y[i, ], prob = pi[k, i, ], log = TRUE)
        }
    }
    
    return(list(ll = ll, pi = pi))
}

#' Log-likelihood for pg_splm() model for use in model selection
#' 
#' this function generates the log-linkelihood for data fit using the pg_splm() 
#' function. The log-likelihood can be used for model selection and evaluation.
#' @param Y is a \eqn{n \times J}{n x J} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param out is a list of MCMC outputs from `pg_splm()`
#' 
#' @importFrom stats dmultinom
#' 
#' @export 

calc_ll_pg_splm <- function(Y, X, out) {
    ##
    ## check the inputs 
    ##
    if (!inherits(out, "pg_splm"))
        stop("The MCMC object out must be of class pg_lm which is the output of the pg_lm() function.")
    
    check_input_pg_lm(Y, X)
    
    N  <- nrow(Y)
    J  <- ncol(Y)
    p  <- ncol(X)
    
    eta      <- out$eta
    n_samples <- dim(out$eta)[1]    
    ## convert from eta to pi
    pi <- sapply(1:n_samples, function(i) eta_to_pi(eta[i, , ]), simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)
    
    ## check in J=2
    if (length(dim(pi)) == 2) {
        p <- aperm(pi, c(2, 1))   
    } else {
        pi <- aperm(pi, c(3, 1, 2))
    }    
    
    ll <- matrix(0, N, n_samples)
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        for (i in 1:N) {
            ll[i, k] <- dmultinom(Y[i, ], prob = pi[k, i, ], log = TRUE)
        }
    }
    
    return(list(ll = ll, pi = pi))
}


#' Log-likelihood for pg_stlm() model for use in model selection
#' 
#' this function generates the log-linkelihood for data fit using the pg_splm() 
#' function. The log-likelihood can be used for model selection and evaluation.
#' @param Y is a \eqn{N \times J \times T}{N x J x T} array of compositional count data.
#' @param X is a \eqn{N \times p}{n_sites x p} matrix of climate variables.
#' @param out is a list of MCMC outputs from `pg_stlm()`
#' 
#' @importFrom stats dmultinom
#' 
#' @export 

calc_ll_pg_stlm <- function(Y, X, out) {
    ##
    ## check the inputs 
    ##
    if (!(inherits(out, "pg_stlm") || inherits(out, "pg_stlm_overdispersed")|| inherits(out, "pg_stlm_latent_overdispersed") || inherits(out, "pg_stlm_mra")))
        stop("The MCMC object out must be of class pg_stlm, pg_stlm_overdispersed, pg_stlm_latent_overdispersed, or pg_stlm_mra which are the output of the pg_stlm(), pg_stlm_overdispersed(), or pg_stlm_mra() functions, respectively.")
    
    if (!is.array(Y)) 
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    if (length(dim(Y)) != 3)
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    missing_idx <- matrix(FALSE, dim(Y)[1], dim(Y)[3])
    for (i in 1:dim(Y)[1]) {
        for (tt in 1:dim(Y)[3]) {
            if(!any(is.na(Y[i, , tt])))
                if(sum(Y[i, , tt]) == 0)
                    stop ("There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
        }
    }
    na_idx <- which(is.na(Y))
    if(length(na_idx) == 0) { 
        if (!is_integer(Y, length(Y)))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    } else {
        if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    }    
    if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
        stop("X must be a numeric matrix.")
    
    N  <- dim(Y)[1]
    J  <- dim(Y)[2]
    n_time <- dim(Y)[3]
    p  <- ncol(X)
    
    eta      <- out$eta
    n_samples <- dim(out$eta)[1]    
    ## convert from eta to pi
    pi <- array(0, dim = c(n_samples, N, J, n_time))
    for (k in 1:n_samples) {
        for(tt in 1:n_time) {
            pi[k, , , tt] <- eta_to_pi(eta[k, , , tt])
        }
    }
    
    # log-likelihood
    ll <- array(0, dim = c(n_samples, N, tt))
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        for (i in 1:N) {
            for (tt in 1:n_time) {
                if (any(is.na(Y[i, , tt]))) {
                    ll[k, i, tt] <- NA
                } else {
                    ll[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi[k, i, , tt], log = TRUE)
                }
            }
        }
    }
    
    return(list(ll = ll, pi = pi))
}


#' Leave-one-out log-likelihood for pg_stlm() model for use in model selection
#' 
#' this function generates the log-likelihood for data fit using the pg_splm() 
#' function. The log-likelihood can be used for model selection and evaluation.
#' @param Y is a \eqn{N \times J \times T}{N x J x T} array of compositional count data.
#' @param X is a \eqn{N \times p}{n_sites x p} matrix of climate variables.
#' @param locs is a \eqn{N \times 2}{N x 2} matrix of observation locations
#' @param out is a list of MCMC outputs from `pg_stlm()`, `pg_stlm_overdispersed`, or `pg_stlm_latent_overdispersed`
#' @param n_message A positive integer number of frequency of iterations
#'  to output a progress message. For example, \code{n_message = 50}
#'  outputs progress messages every 50 iterations.
#' @importFrom stats dmultinom
#' 
#' @export 

calc_ll_pg_stlm_loo <- function(Y, X, locs, out, n_message = 50) {
    ##
    ## check the inputs 
    ##
    if (!(inherits(out, "pg_stlm") || inherits(out, "pg_stlm_overdispersed") || inherits(out, "pg_stlm_latent_overdispersed") || inherits(out, "pg_stlm_mra")))
        stop("The MCMC object out must be of class pg_stlm, pg_stlm_overdispersed, pg_stlm_latent_overdispersed, or pg_stlm_mra which are the output of the pg_stlm(), pg_stlm_overdispersed(), or pg_stlm_mra() functions, respectively.")
    
    if (!is.array(Y)) 
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    if (length(dim(Y)) != 3)
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    missing_idx <- matrix(FALSE, dim(Y)[1], dim(Y)[3])
    for (i in 1:dim(Y)[1]) {
        for (tt in 1:dim(Y)[3]) {
            if(!any(is.na(Y[i, , tt])))
                if(sum(Y[i, , tt]) == 0)
                    stop ("There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
        }
    }
    na_idx <- which(is.na(Y))
    if(length(na_idx) == 0) { 
        if (!is_integer(Y, length(Y)))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    } else {
        if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    }    
    if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
        stop("X must be a numeric matrix.")
    
    ##
    ## TODO add checks for locs?
    ##
    

    # calculate the ll    
    ll <- NULL
    
    if (inherits(out, "pg_stlm")) {
        ll <- calc_ll_pg_stlm_matern(Y, X, locs, out, n_message)
    }
    
    if (inherits(out, "pg_stlm_overdispersed")) {
        ll <- calc_ll_pg_stlm_overdispersed(Y, X, locs, out, n_message)
    }
    
    if (inherits(out, "pg_stlm_latent_overdispersed")) {
        ll <- calc_ll_pg_stlm_latent_overdispersed(Y, X, locs, out, n_message)
    }
    
    if (inherits(out, "pg_stlm_mra")) {
        stop("the MRA loo-ll has not been completed")
        ll <- calc_ll_pg_stlm_mra(Y, X, locs, out, n_message)
    }
        
    return(ll)
}

calc_ll_pg_stlm_matern <- function(Y, X, locs, out, n_message = 50) {
    
    # TODO For now, this only supports non-shared covariance parameters
    
    if(is.null(dim(out$tau2)[2]))
        stop("shared covariance parameters are currently not supported")
    
    N  <- dim(Y)[1]
    J  <- dim(Y)[2]
    n_time <- dim(Y)[3]
    p  <- ncol(X)
    
    eta       <- out$eta
    beta      <- out$beta
    theta     <- out$theta
    tau2      <- out$tau2
    rho       <- out$rho
    pi        <- out$pi
    n_samples <- dim(out$eta)[1]    
    
    eta_loo    <- array(0, dim = dim(eta)) 
    pi_loo     <- array(0, dim = dim(pi)) 
    
    corr_fun <- NULL
    if (is.null(dim(theta)[3])) {
        corr_fun <- "exponential"
    } else {
        corr_fun <- "matern"
    }
    
    D <- fields::rdist(locs)
    
    # initialized non-shared covariance matrices
    Sigma <- array(D, dim = c(N, N, J-1))
    Sigma_chol <- array(D, dim = c(N, N, J-1))

    
    Sigma_obs_unobs <- array(D, dim = c(N-1, 1, J-1))
    Sigma_unobs_unobs <- array(D, dim = c(N-1, N-1, J-1))
    Sigma_unobs_unobs_inv <- array(D, dim = c(N-1, N-1, J-1))
    Sigma_obs_obs     <- rep(0, J-1)
    Sigma_obs_unobs_Sigma_unobs_unobs_inv <- array(D, dim = c(N-1, 1, J-1))
    
    message("Starting loo likelihood, running for ", n_samples, " iterations")
    for (k in 1:n_samples) {
        if (k %% n_message == 0) {
            message("On iteration ", k, " out of ", n_samples)
        }
        
        for (j in 1:(J-1)) {
            if (corr_fun == "exponential") {
                Sigma[, , j] <- correlation_function(D, theta[k, j], corr_fun = corr_fun)
            } else {
                Sigma[, , j] <- correlation_function(D, theta[k, j, ], corr_fun = corr_fun)
            }
            Sigma_chol[, , j] <- chol(Sigma[, , j])
        }
        
        for (i in 1:N) {
            
            # message("iteration ", i)
            for(j in 1:(J-1)) {
                Sigma_obs_unobs[, , j]  <- Sigma[i, -i, j]
                Sigma_unobs_unobs[, , j] <- Sigma[-i, -i, j]
                Sigma_obs_obs[j] <- Sigma[i, i, j]
                Sigma_unobs_unobs_chol <- choldrop(Sigma_chol[, , j], i)
                Sigma_unobs_unobs_inv[, , j] <- chol2inv(Sigma_unobs_unobs_chol)
                Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] <- Sigma_obs_unobs[, , j] %*% Sigma_unobs_unobs_inv[, , j] 
            }
            
            for(tt in 1:n_time) {
                for (j in 1:(J-1)) {
                    mu_bar <- X[i, ] %*% beta[k, , j]
                    tmp <- eta[k, -i, j, tt] - X[-i, , drop = FALSE] %*% beta[k, , j]
                    if (tt > 1) {
                        mu_bar <- mu_bar + rho[k] * (eta[k, i, j, tt-1] - X[i, ] %*% beta[k, , j])
                        tmp <- tmp - rho[k] * (eta[k, -i, j, tt-1] - X[-i, , drop = FALSE] %*% beta[k, , j])
                    }
                    if (tt < n_time) {
                        mu_bar <- mu_bar + rho[k] * (eta[k, i, j, tt+1] - X[i, ] %*% beta[k, , j])
                        tmp <- tmp - rho[k] * (eta[k, -i, j, tt+1] - X[-i, , drop = FALSE] %*% beta[k, , j])
                    }
                    mu_bar <- mu_bar + Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] %*% tmp
                    sigma_bar <- sqrt(Sigma_obs_obs[j] - Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] %*% Sigma_obs_unobs[, , j])
                    eta_loo[k, i, j, tt] <- rnorm(1, mu_bar, sigma_bar)
                }
            }
        }
    }
    
    
    for (k in 1:n_samples) {
        for(tt in 1:n_time) {
            pi[k, , , tt] <- eta_to_pi(eta[k, , , tt])
            pi_loo[k, , , tt] <- eta_to_pi(eta_loo[k, , , tt])
        }
    }
    
    # log-likelihood
    ll <- array(0, dim = c(n_samples, N, tt))
    ll_loo <- array(0, dim = c(n_samples, N, tt))
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        for (i in 1:N) {
            for (tt in 1:n_time) {
                if (any(is.na(Y[i, , tt]))) {
                    ll[k, i, tt] <- NA
                    ll_loo[k, i, tt] <- NA
                } else {
                    ll[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi[k, i, , tt], log = TRUE)
                    ll_loo[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi_loo[k, i, , tt], log = TRUE)
                }
            }
        }
    }

    return(list(ll = ll, pi = pi, ll_loo = ll_loo, pi_loo = pi_loo))
}

calc_ll_pg_stlm_overdispersed <- function(Y, X, locs, out, n_message = 50) {
    
    # TODO For now, this only supports non-shared covariance parameters
    
    if(is.null(dim(out$tau2)[2]))
        stop("shared covariance parameters are currently not supported")
    
    N  <- dim(Y)[1]
    J  <- dim(Y)[2]
    n_time <- dim(Y)[3]
    p  <- ncol(X)
    I_N <- diag(N)
    I_Nm1 <- diag(N - 1)
    
    eta       <- out$eta
    beta      <- out$beta
    theta     <- out$theta
    sigma2    <- out$sigma2
    tau2      <- out$tau2
    rho       <- out$rho
    pi        <- out$pi
    n_samples <- dim(out$eta)[1]    
    
    eta_loo    <- array(0, dim = dim(eta)) 
    pi_loo     <- array(0, dim = dim(pi)) 
    
    corr_fun <- NULL
    if (is.null(dim(theta)[3])) {
        corr_fun <- "exponential"
    } else {
        corr_fun <- "matern"
    }
    
    D <- fields::rdist(locs)
    
    # initialized non-shared covariance matrices
    Sigma <- array(D, dim = c(N, N, J-1))
    Sigma_chol <- array(D, dim = c(N, N, J-1))
    
    Sigma_obs_unobs <- array(D, dim = c(N-1, 1, J-1))
    Sigma_unobs_unobs <- array(D, dim = c(N-1, N-1, J-1))
    Sigma_unobs_unobs_inv <- array(D, dim = c(N-1, N-1, J-1))
    Sigma_obs_obs     <- rep(0, J-1)
    Sigma_obs_unobs_Sigma_unobs_unobs_inv <- array(D, dim = c(N-1, 1, J-1))
    
    message("Starting loo likelihood, running for ", n_samples, " iterations")
    for (k in 1:n_samples) {
        if (k %% n_message == 0) {
            message("On iteration ", k, " out of ", n_samples)
        }
        
        for (j in 1:(J-1)) {
            if (corr_fun == "exponential") {
                Sigma[, , j] <- correlation_function(D, theta[k, j], corr_fun = corr_fun) + sigma2[k, j] * I_N
            } else {
                Sigma[, , j] <- correlation_function(D, theta[k, j, ], corr_fun = corr_fun) + sigma2[k, j] * I_N
            }
            Sigma_chol[, , j] <- chol(Sigma[, , j])
        }

        for (i in 1:N) {
            # message("iteration ", i)
            for(j in 1:(J-1)) {
                Sigma_obs_unobs[, , j]  <- Sigma[i, -i, j]
                Sigma_unobs_unobs[, , j] <- Sigma[-i, -i, j]
                Sigma_obs_obs[j] <- Sigma[i, i, j]
                Sigma_unobs_unobs_chol <- choldrop(Sigma_chol[, , j], i)
                Sigma_unobs_unobs_inv[, , j] <- chol2inv(Sigma_unobs_unobs_chol)
                Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] <- Sigma_obs_unobs[, , j] %*% Sigma_unobs_unobs_inv[, , j] 
            } 
            
            for(tt in 1:n_time) {
                for (j in 1:(J-1)) {
                    mu_bar <- X[i, ] %*% beta[k, , j]
                    tmp <- eta[k, -i, j, tt] - X[-i, , drop = FALSE] %*% beta[k, , j]
                    if (tt > 1) {
                        mu_bar <- mu_bar + rho[k] * (eta[k, i, j, tt-1] - X[i, ] %*% beta[k, , j])
                        tmp <- tmp - rho[k] * (eta[k, -i, j, tt-1] - X[-i, , drop = FALSE] %*% beta[k, , j])
                    }
                    if (tt < n_time) {
                        mu_bar <- mu_bar + rho[k] * (eta[k, i, j, tt+1] - X[i, ] %*% beta[k, , j])
                        tmp <- tmp - rho[k] * (eta[k, -i, j, tt+1] - X[-i, , drop = FALSE] %*% beta[k, , j])
                    }
                    mu_bar <- mu_bar + Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] %*% tmp
                    sigma_bar <- sqrt(Sigma_obs_obs[j] - Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] %*% Sigma_obs_unobs[, , j])
                    eta_loo[k, i, j, tt] <- rnorm(1, mu_bar, sigma_bar)
                }
            }
            
            
        }
    }
    
    
    for (k in 1:n_samples) {
        for(tt in 1:n_time) {
            pi[k, , , tt] <- eta_to_pi(eta[k, , , tt])
            pi_loo[k, , , tt] <- eta_to_pi(eta_loo[k, , , tt])
        }
    }
    
    # log-likelihood
    ll <- array(0, dim = c(n_samples, N, tt))
    ll_loo <- array(0, dim = c(n_samples, N, tt))
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        for (i in 1:N) {
            for (tt in 1:n_time) {
                if (any(is.na(Y[i, , tt]))) {
                    ll[k, i, tt] <- NA
                    ll_loo[k, i, tt] <- NA
                } else {
                    ll[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi[k, i, , tt], log = TRUE)
                    ll_loo[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi_loo[k, i, , tt], log = TRUE)
                }
            }
        }
    }
    
    return(list(ll = ll, pi = pi, ll_loo = ll_loo, pi_loo = pi_loo))
}


calc_ll_pg_stlm_latent_overdispersed <- function(Y, X, locs, out, n_message = 50) {
    
    # TODO For now, this only supports non-shared covariance parameters
    
    if(is.null(dim(out$tau2)[2]))
        stop("shared covariance parameters are currently not supported")
    
    N  <- dim(Y)[1]
    J  <- dim(Y)[2]
    n_time <- dim(Y)[3]
    p  <- ncol(X)
    I_N <- diag(N)
    I_Nm1 <- diag(N - 1)
    
    eta       <- out$eta
    beta      <- out$beta
    theta     <- out$theta
    sigma2    <- out$sigma2
    psi       <- out$psi
    tau2      <- out$tau2
    rho       <- out$rho
    pi        <- out$pi
    n_samples <- dim(out$eta)[1]    
    
    psi_loo    <- array(0, dim = dim(psi)) 
    eta_loo    <- array(0, dim = dim(eta)) 
    pi_loo     <- array(0, dim = dim(pi)) 
    
    corr_fun <- NULL
    if (is.null(dim(theta)[3])) {
        corr_fun <- "exponential"
    } else {
        corr_fun <- "matern"
    }
    
    D <- fields::rdist(locs)
    
    # initialized non-shared covariance matrices
    Sigma <- array(D, dim = c(N, N, J-1))
    Sigma_chol <- array(D, dim = c(N, N, J-1))
    
    Sigma_obs_unobs <- array(D, dim = c(N-1, 1, J-1))
    Sigma_unobs_unobs <- array(D, dim = c(N-1, N-1, J-1))
    Sigma_unobs_unobs_inv <- array(D, dim = c(N-1, N-1, J-1))
    Sigma_obs_obs     <- rep(0, J-1)
    Sigma_obs_unobs_Sigma_unobs_unobs_inv <- array(D, dim = c(N-1, 1, J-1))
    
    message("Starting loo likelihood, running for ", n_samples, " iterations")
    for (k in 1:n_samples) {
        if (k %% n_message == 0) {
            message("On iteration ", k, " out of ", n_samples)
        }
        
        for (j in 1:(J-1)) {
            if (corr_fun == "exponential") {
                Sigma[, , j] <- correlation_function(D, theta[k, j], corr_fun = corr_fun)
            } else {
                Sigma[, , j] <- correlation_function(D, theta[k, j, ], corr_fun = corr_fun)
            }
            Sigma_chol[, , j] <- chol(Sigma[, , j])
        }
        
        
        for (i in 1:N) {
            # message("iteration ", i)
            for(j in 1:(J-1)) {
                Sigma_obs_unobs[, , j]  <- Sigma[i, -i, j]
                Sigma_unobs_unobs[, , j] <- Sigma[-i, -i, j]
                Sigma_obs_obs[j] <- Sigma[i, i, j]
                Sigma_unobs_unobs_chol <- choldrop(Sigma_chol[, , j], i)
                Sigma_unobs_unobs_inv[, , j] <- chol2inv(Sigma_unobs_unobs_chol)
                Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] <- Sigma_obs_unobs[, , j] %*% Sigma_unobs_unobs_inv[, , j] 
            } 
            
            
            for(tt in 1:n_time) {
                for (j in 1:(J-1)) {
                    mu_bar <- 0
                    tmp <- psi[k, -i, j, tt] 
                    if (tt > 1) {
                        mu_bar <- mu_bar + rho[k] * psi[k, i, j, tt-1]
                        tmp <- tmp - rho[k] * psi[k, -i, j, tt-1]
                    }
                    if (tt < n_time) {
                        mu_bar <- mu_bar + rho[k] * psi[k, i, j, tt+1]
                        tmp <- tmp - rho[k] * psi[k, -i, j, tt+1]
                    }
                    mu_bar <- mu_bar + Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] %*% tmp
                    sigma_bar <- sqrt(Sigma_obs_obs[j] - Sigma_obs_unobs_Sigma_unobs_unobs_inv[, , j] %*% Sigma_obs_unobs[, , j])
                    psi_loo[k, i, j, tt] <- rnorm(1, mu_bar, sigma_bar)
                    eta_loo[k, i, j, tt] <- X[i, ] %*% beta[k, , j] + psi_loo[k, i, j, tt] + rnorm(1, 0, sqrt(sigma2[k, j]))
                }
            }
            
            
        }
    }
    
    
    for (k in 1:n_samples) {
        for(tt in 1:n_time) {
            pi[k, , , tt] <- eta_to_pi(eta[k, , , tt])
            pi_loo[k, , , tt] <- eta_to_pi(eta_loo[k, , , tt])
        }
    }
    
    # log-likelihood
    ll <- array(0, dim = c(n_samples, N, tt))
    ll_loo <- array(0, dim = c(n_samples, N, tt))
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        for (i in 1:N) {
            for (tt in 1:n_time) {
                if (any(is.na(Y[i, , tt]))) {
                    ll[k, i, tt] <- NA
                    ll_loo[k, i, tt] <- NA
                } else {
                    ll[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi[k, i, , tt], log = TRUE)
                    ll_loo[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi_loo[k, i, , tt], log = TRUE)
                }
            }
        }
    }
    
    return(list(ll = ll, pi = pi, ll_loo = ll_loo, pi_loo = pi_loo))
}


calc_ll_pg_stlm_mra <- function(Y, X, out, locs, n_message = 50) {
    N  <- dim(Y)[1]
    J  <- dim(Y)[2]
    n_time <- dim(Y)[3]
    p  <- ncol(X)
    
    eta      <- out$eta
    n_samples <- dim(out$eta)[1]    
    ## convert from eta to pi
    pi <- array(0, dim = c(n_samples, N, J, n_time))
    for (k in 1:n_samples) {
        for(tt in 1:n_time) {
            pi[k, , , tt] <- eta_to_pi(eta[k, , , tt])
        }
    }
    
    # log-likelihood
    ll <- array(0, dim = c(n_samples, N, tt))
    ## parallelize/vectorize this later
    for (k in 1:n_samples) {
        if (k %% n_message == 0) {
            message("On iteration ", k, " out of ", n_samples)
        }
        
        for (i in 1:N) {
            for (tt in 1:n_time) {
                if (any(is.na(Y[i, , tt]))) {
                    ll[k, i, tt] <- NA
                } else {
                    ll[k, i, tt] <- dmultinom(Y[i, , tt], prob = pi[k, i, , tt], log = TRUE)
                }
            }
        }
    }
    
    return(list(ll = ll, pi = pi))
}


