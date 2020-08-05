context("test check functions")

test_that("check_corr_function", {
    expect_silent(check_corr_fun(corr_fun = "exponential"))
    expect_silent(check_corr_fun(corr_fun = "matern"))
    expect_error(check_corr_fun(corr_fun = "aaa"), 'corr_fun must be either "matern" or "exponential"')
    expect_error(check_corr_fun(corr_fun = NA), 'corr_fun must be either "matern" or "exponential"')
    expect_error(check_corr_fun(corr_fun = NULL), 'corr_fun must be either "matern" or "exponential"')
    expect_error(check_corr_fun(corr_fun = 43), 'corr_fun must be either "matern" or "exponential"')
    expect_error(check_corr_fun(corr_fun = c("exponential", "matern")), 'corr_fun must be either "matern" or "exponential"')
})

test_that("check_priors_pg_lm", {
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(20), 10, 2)
    expect_silent(check_priors_pg_lm(Y, X, priors = NULL))
    priors <- default_priors_pg_lm(Y, X)
    expect_silent(check_priors_pg_lm(Y, X, priors))
    priors$mu_beta <- c("aaa", "aaa")
    expect_error(check_priors_pg_lm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X")
    priors$mu_beta <- 1:3
    expect_error(check_priors_pg_lm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X")
    priors$mu_beta <- c(1, NA)
    expect_error(check_priors_pg_lm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X")
    priors$mu_beta <- NULL
    
    priors$Sigma_beta <- diag(3)
    expect_error(check_priors_pg_lm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X")
    priors$Sigma_beta <- 1:3
    expect_error(check_priors_pg_lm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X")
    priors$Sigma_beta <- diag(3)
    priors$Sigma_beta[1, 1] <- NA
    expect_error(check_priors_pg_lm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X")
})

test_that("check_priors_pg_splm", {
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(20), 10, 2)
    expect_silent(check_priors_pg_splm(Y, X, priors = NULL))
    ## check the exponential covariance
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    expect_silent(check_priors_pg_splm(Y, X, priors))
    priors <- default_priors_pg_splm(Y, X, corr_fun = "matern")
    expect_silent(check_priors_pg_splm(Y, X, priors))
    
    priors$mu_beta <- c("aaa", "aaa")
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- NULL
    
    priors$Sigma_beta <- diag(3)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- diag(3)
    priors$Sigma_beta[1, 1] <- NA
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- NULL
    
    priors$alpha_tau <- "aaa"
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- -1
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- NULL
    
    priors$beta_tau <- "aaa"
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- -1
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- NULL
    
    priors$mean_range <- "aaa"
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mean_range, it must be a numeric value.")
    priors$mean_range <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mean_range, it must be a numeric value.")
    priors$mean_range <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mean_range, it must be a numeric value.")
    priors$mean_range <- NULL
    
    priors$sd_range <- "aaa"
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- -1
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- NULL
    
    priors$mean_nu <- "aaa"
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mean_nu, it must be a numeric value.")
    priors$mean_nu <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mean_nu, it must be a numeric value.")
    priors$mean_nu <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for mean_nu, it must be a numeric value.")
    priors$mean_nu <- NULL
    
    priors$sd_nu <- "aaa"
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- 1:3
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- c(1, NA)
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- -1
    expect_error(check_priors_pg_splm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- NULL

})


test_that("check_priors_pg_stlm", {
    Y <- array(1:160, dim = c(10, 4, 4))
    X <- matrix(rnorm(20), 10, 2)
    expect_silent(check_priors_pg_stlm(Y, X, priors = NULL))
    ## check the exponential covariance
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "exponential")
    expect_silent(check_priors_pg_stlm(Y, X, priors))
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "matern")
    expect_silent(check_priors_pg_stlm(Y, X, priors))
    
    priors$mu_beta <- c("aaa", "aaa")
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mu_beta, it must be a numeric vector of length equal to the number of columns of X.")
    priors$mu_beta <- NULL
    
    priors$Sigma_beta <- diag(3)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- diag(3)
    priors$Sigma_beta[1, 1] <- NA
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for Sigma_beta, it must be a symmertic positive definite matrix with rows and columns equal to the number of columns of X.")
    priors$Sigma_beta <- NULL
    
    priors$alpha_tau <- "aaa"
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- -1
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for alpha_tau, it must be a positive numeric value.")
    priors$alpha_tau <- NULL
    
    priors$beta_tau <- "aaa"
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- -1
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for beta_tau, it must be a positive numeric value.")
    priors$beta_tau <- NULL
    
    priors$mean_range <- "aaa"
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mean_range, it must be a numeric value.")
    priors$mean_range <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mean_range, it must be a numeric value.")
    priors$mean_range <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mean_range, it must be a numeric value.")
    priors$mean_range <- NULL
    
    priors$sd_range <- "aaa"
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- -1
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_range, it must be a positive numeric value.")
    priors$sd_range <- NULL
    
    priors$mean_nu <- "aaa"
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mean_nu, it must be a numeric value.")
    priors$mean_nu <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mean_nu, it must be a numeric value.")
    priors$mean_nu <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for mean_nu, it must be a numeric value.")
    priors$mean_nu <- NULL
    
    priors$sd_nu <- "aaa"
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- 1:3
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- c(1, NA)
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- -1
    expect_error(check_priors_pg_stlm(Y, X, priors), "If priors contains a value for sd_nu, it must be a positive numeric value.")
    priors$sd_nu <- NULL
})



test_that("check_inits_pg_lm", {
    
    ## Check initial conditions for pg_lm
    expect_error(check_inits_pgLM(), "The function check_inits_pgLM\\(\\) has been deprecated. Please use check_inits_pg_lm\\(\\) instead.")

    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_lm(Y, X)
    inits  <- default_inits_pg_lm(Y, X, priors)
    expect_silent(check_inits_pg_lm(Y, X, inits))

    X <- matrix(rnorm(20), 10, 4)    
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 10, 5)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")

    Y <- matrix(1:50, 5, 10)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_lm(Y, X, inits), "Y and X must have the same number of rows")
    
    Y <- matrix(1:50, 5, 10)
    X <- rnorm(10)
    expect_error(check_inits_pg_lm(Y, X, inits), "X must be a numeric matrix")        
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(10), 10, 1)
    params <- default_params()
    priors <- default_priors_pg_lm(Y, X)
    inits  <- default_inits_pg_lm(Y, X, priors)
    
    inits$beta <- NA
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")

    inits$beta <- 1:3
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(1:9, 3, 3)
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(c(1:2, NA), 3, 1)
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    expect_identical(
        {
            set.seed(111)
            inits  <- default_inits_pg_lm(Y, X, priors)
        }, 
        {
            set.seed(111)
            inits <- list(
                beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
            )
            
        }
    )
})    


test_that("check_inits_pg_splm", {
    
    ## Check initial conditions for pg_lm
    expect_error(check_inits_pgSPLM(), "The function check_inits_pgSPLM\\(\\) has been deprecated. Please use check_inits_pg_splm\\(\\) instead.")
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(20), 10, 2)
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    inits  <- default_inits_pg_splm(Y, X, priors)
    expect_silent(default_inits_pg_splm(Y, X, inits))
    
    X <- matrix(rnorm(20), 10, 4)    
    expect_error(default_inits_pg_splm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 10, 5)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(default_inits_pg_splm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 5, 10)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(default_inits_pg_splm(Y, X, inits), "Y and X must have the same number of rows")
    
    Y <- matrix(1:50, 5, 10)
    X <- rnorm(10)
    expect_error(default_inits_pg_splm(Y, X, inits), "X must be a numeric matrix")        
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(10), 10, 1)
    params <- default_params()
    priors <- default_priors_pg_splm(Y, X)
    inits  <- default_inits_pg_splm(Y, X, priors)
    
    inits$beta <- NA
    expect_error(check_inits_pg_splm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- 1:3
    expect_error(check_inits_pg_splm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(1:9, 3, 3)
    expect_error(check_inits_pg_splm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(c(1:2, NA), 3, 1)
    expect_error(check_inits_pg_splm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    expect_identical(
        {
            set.seed(111)
            inits  <- default_inits_pg_splm(Y, X, priors)
        }, 
        {
            set.seed(111)
            inits <- list(
                beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
            )
            
        }
    )

    ## check inits for pg_spvlm
    
    ## check inits for pg_stlm
})    

test_that("check_inits_pg_stlm", {
    
    ## Check initial conditions for pg_lm
    expect_error(check_inits_pgSTLM(), "The function check_inits_pgSTLM\\(\\) has been deprecated. Please use check_inits_pg_stlm\\(\\) instead.")
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(20), 10, 2)
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "exponential")
    inits  <- default_inits_pg_stlm(Y, X, priors)
    expect_silent(check_inits_pg_stlm(Y, X, inits))
    
    X <- matrix(rnorm(20), 10, 4)    
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 10, 5)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 5, 10)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_lm(Y, X, inits), "Y and X must have the same number of rows")
    
    Y <- matrix(1:50, 5, 10)
    X <- rnorm(10)
    expect_error(check_inits_pg_lm(Y, X, inits), "X must be a numeric matrix")        
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(10), 10, 1)
    params <- default_params()
    priors <- default_priors_pg_lm(Y, X)
    inits  <- default_inits_pg_lm(Y, X, priors)
    
    inits$beta <- NA
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- 1:3
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(1:9, 3, 3)
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(c(1:2, NA), 3, 1)
    expect_error(check_inits_pg_lm(Y, X, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    expect_identical(
        {
            set.seed(111)
            inits  <- default_inits_pg_lm(Y, X, priors)
        }, 
        {
            set.seed(111)
            inits <- list(
                beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
            )
            
        }
    )
    
    ## check inits for pg_spvlm

})    

test_that("check_input", {
    
    ## check input for pg_lm
    
    ## check input for pg_splm
    
    ## check input for pg_spvlm
    
    ## check input for pg_stlm
    
})

test_that("check_params", {
    params <- default_params()
    expect_invisible(check_params(params))
    
    params <- default_params()
    params$n_mcmc <- 1:5
    expect_error(check_params(params), "params must contain a positive integer n_mcmc")   
    params$n_mcmc <- 22.3
    expect_error(check_params(params), "params must contain a positive integer n_mcmc")    
    params$n_mcmc <- -1
    expect_error(check_params(params), "params must contain a positive integer n_mcmc")    
    params$n_mcmc <- NA
    expect_error(check_params(params), "params must contain a positive integer n_mcmc")  
    params$n_mcmc <- NULL
    expect_error(check_params(params), "params must contain a positive integer n_mcmc")  
    params$n_mcmc <- "aaa"
    expect_error(check_params(params), "params must contain a positive integer n_mcmc")    
    
    params <- default_params()
    params$n_adapt <- 1:5
    expect_error(check_params(params), "params must contain a positive integer n_adapt")   
    params$n_adapt <- 22.3
    expect_error(check_params(params), "params must contain a positive integer n_adapt")    
    params$n_adapt <- -1
    expect_error(check_params(params), "params must contain a positive integer n_adapt")    
    params$n_adapt <- NA
    expect_error(check_params(params), "params must contain a positive integer n_adapt")  
    params$n_adapt <- NULL
    expect_error(check_params(params), "params must contain a positive integer n_adapt")  
    params$n_adapt <- "aaa"
    expect_error(check_params(params), "params must contain a positive integer n_adapt")   
    
    params <- default_params()
    params$n_thin <- 1:5
    expect_error(check_params(params), "params must contain a positive integer n_thin")   
    params$n_thin <- 22.3
    expect_error(check_params(params), "params must contain a positive integer n_thin")    
    params$n_thin <- -1
    expect_error(check_params(params), "params must contain a positive integer n_thin")    
    params$n_thin <- NA
    expect_error(check_params(params), "params must contain a positive integer n_thin")  
    params$n_thin <- NULL
    expect_error(check_params(params), "params must contain a positive integer n_thin")  
    params$n_thin <- "aaa"
    expect_error(check_params(params), "params must contain a positive integer n_thin")   
    
    params <- default_params()
    params$n_message <- 1:5
    expect_error(check_params(params), "params must contain a positive integer n_message")    
    params$n_message <- 22.3
    expect_error(check_params(params), "params must contain a positive integer n_message")    
    params$n_message <- -1
    expect_error(check_params(params), "params must contain a positive integer n_message")    
    params$n_message <- NA
    expect_error(check_params(params), "params must contain a positive integer n_message")  
    params$n_message <- NULL
    expect_error(check_params(params), "params must contain a positive integer n_message")  
    params$n_message <- "aaa"
    expect_error(check_params(params), "params must contain a positive integer n_message")    
    
    params$n_message <- 22.3
    expect_error(check_params(params), "params must contain a positive integer n_message")    
    params$n_message <- -1
    expect_error(check_params(params), "params must contain a positive integer n_message")    
    params$n_message <- NA
    expect_error(check_params(params), "params must contain a positive integer n_message")  
    params$n_message <- NULL
    expect_error(check_params(params), "params must contain a positive integer n_message")  
    params$n_message <- "aaa"
    expect_error(check_params(params), "params must contain a positive integer n_message")    
})







