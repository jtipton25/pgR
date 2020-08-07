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
    inits$beta <- NULL
    
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
    inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE)
    
    expect_error(check_inits_pgSPLM(Y, X, inits), "The function check_inits_pgSPLM\\(\\) has been deprecated. Please use check_inits_pg_splm\\(\\) instead.")
    
    for (j in c("exponential", "matern")) {
        for (k in c(TRUE, FALSE)) {
            priors <- default_priors_pg_splm(Y, X, corr_fun = j)
            inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = j, shared_covariance_params = k)
            expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = j, shared_covariance_params = k))
        }
    }
    
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = "exponential" and shared_covariance_params = FALSE.')
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.')
    inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = FALSE)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.')
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.')
    
    
    priors <- default_priors_pg_splm(Y, X, corr_fun = "matern")
    inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "matern", shared_covariance_params = TRUE)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.')
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.')
    inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "matern", shared_covariance_params = FALSE)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = \"exponential\" and shared_covariance_params = FALSE.')
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.')
    
    inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "matern", shared_covariance_params = FALSE)
    X <- matrix(rnorm(20), 10, 4)    
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 10, 5)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 5, 10)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "Y and X must have the same number of rows")
    
    Y <- matrix(1:50, 5, 10)
    X <- rnorm(10)
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "X must be a numeric matrix")        
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_splm(Y, X)
    inits  <- default_inits_pg_splm(Y, X, priors)
    
    inits$beta <- NA
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    inits$beta <- 1:3
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    inits$beta <- matrix(1:9, 3, 3)
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    inits$beta <- matrix(c(1:2, NA), 3, 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    inits$beta <- NULL
    
    inits$theta <- NA
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.")
    inits$theta <- "aaa"
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.")
    inits$theta <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.")
    inits$theta <- 1
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    inits$theta <- -1
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    
    inits$tau2 <- rep(5, ncol(Y)-1)
    inits$theta <- rep(NA, ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = \"exponential\" and shared_covariance_params = FALSE.")
    inits$theta <- rep("aaa", ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = \"exponential\" and shared_covariance_params = FALSE.")
    inits$theta <- 1:(ncol(Y) - 1)
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    inits$theta <- -(1:(ncol(Y) - 1))
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    
    inits$tau2 <- 5
    inits$theta <- rep(NA, 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.")
    inits$theta <- rep("aaa", 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.")
    inits$theta <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.")
    inits$theta <- 1:2
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    inits$theta <- -c(1:2)
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    
    inits$tau2 <- rep(5, ncol(Y)- 1)
    inits$theta <- matrix(NA, ncol(Y) - 1, 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.")
    inits$theta <- matrix("aaa", ncol(Y)-1, 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.")
    inits$theta <- matrix(1:(ncol(Y) * 2 - 2), ncol(Y) -1, 2)
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    inits$theta <- matrix(-(1:(ncol(Y) * 2 - 2)), ncol(Y) -1, 2)
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    

    inits$theta <- 1
    inits$tau2 <- NA
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- "aaa"
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- -3
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    
    inits$theta <- rep(5, ncol(Y)-1)
    inits$tau2 <- rep(NA, ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- rep("aaa", ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- -c(1:(ncol(Y) - 1))
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    
    inits$theta <- c(5, 5)
    inits$tau2 <- rep(NA, 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- rep("aaa", 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- -1
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    
    inits$theta <- matrix(5, ncol(Y) - 1, 2)
    inits$tau2 <- rep(NA, ncol(Y)- 1)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- matrix("aaa", ncol(Y)-1, 2)
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- -c(1:(ncol(Y) - 1))
    expect_error(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_silent(check_inits_pg_splm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    
    
    # expect_identical(
    #     {
    #         set.seed(111)
    #         inits  <- default_inits_pg_splm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE)
    #     }, 
    #     {
    #         set.seed(111)
    #         inits <- list(
    #             beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)),
    #             tau2 = min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10),
    #             theta = pmin(pmax(stats::rnorm(1, theta_mean, sqrt(theta_var)), -1), 0.1)
    #         )
    #         
    #     }
    # )

    ## check inits for pg_spvlm
    
    ## check inits for pg_stlm
})    



test_that("check_inits_pgSTLM", {
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(rnorm(20), 10, 2)
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_stlm(Y, X)
    inits  <- default_inits_pg_stlm(Y, X, priors)
    expect_error(check_inits_pgSTLM(Y, X, inits), "The function check_inits_pgSTLM\\(\\) has been deprecated. Please use check_inits_pg_stlm\\(\\) instead.")
})
test_that("check_inits_pg_stlm", {
    
    ## Check initial conditions for pg_lm
    # expect_error(check_inits_pgSTLM(), "The function check_inits_pgSTLM\\(\\) has been deprecated. Please use check_inits_pg_stlm\\(\\) instead.")
    
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(rnorm(20), 10, 2)
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()

    for (j in c("exponential", "matern")) {
        for (k in c(TRUE, FALSE)) {
            priors <- default_priors_pg_stlm(Y, X, corr_fun = j)
            inits  <- default_inits_pg_stlm(Y, X, priors, corr_fun = j, shared_covariance_params = k)
            expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = j, shared_covariance_params = k))
        }
    }
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "exponential")
    inits  <- default_inits_pg_stlm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = "exponential" and shared_covariance_params = FALSE.')
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.')
    inits  <- default_inits_pg_stlm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = FALSE)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.')
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.')
    
    
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "matern")
    inits  <- default_inits_pg_stlm(Y, X, priors, corr_fun = "matern", shared_covariance_params = TRUE)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.')
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.')
    inits  <- default_inits_pg_stlm(Y, X, priors, corr_fun = "matern", shared_covariance_params = FALSE)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), 'If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = \"exponential\" and shared_covariance_params = FALSE.')
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), 'If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.')

    
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "matern")
    inits  <- default_inits_pg_stlm(Y, X, priors, corr_fun = "matern")
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern"))
    
    X <- matrix(rnorm(20), 10, 4)    
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    Y <- matrix(1:50, 10, 5)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- list(1:50)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    
    Y <- matrix(1:50, 5, 10)
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")

    Y <- array(1:250, dim=c(5, 10, 5))
    X <- matrix(rnorm(10), 10, 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), " and X must have the same number of rows.")
    
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- rnorm(10)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "X must be a numeric matrix")        
    
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(rnorm(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "matern")
    inits  <- default_inits_pg_stlm(Y, X, priors)
    
    inits$beta <- NA
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- 1:3
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(1:9, 3, 3)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    
    inits$beta <- matrix(c(1:2, NA), 3, 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits), "If beta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of X and columns equal to the number of columns of Y - 1.")
    inits$beta <- NULL
    
    inits$theta <- NA
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.")
    inits$theta <- "aaa"
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.")
    inits$theta <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric value when corr_fun = \"exponential\" and shared_covariance_params = TRUE.")
    inits$theta <- 1
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    
    inits$tau2 <- rep(5, ncol(Y)-1)
    inits$theta <- rep(NA, ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = \"exponential\" and shared_covariance_params = FALSE.")
    inits$theta <- rep("aaa", ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric vector of length equal to the number of columns of Y - 1 when corr_fun = \"exponential\" and shared_covariance_params = FALSE.")
    inits$theta <- 1:(ncol(Y) - 1)
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    
    inits$tau2 <- 5
    inits$theta <- rep(NA, 2)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.")
    inits$theta <- rep("aaa", 2)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.")
    inits$theta <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If theta is specified in inits, it must be a numeric vector of length 2 when corr_fun = \"matern\" and shared_covariance_params = TRUE.")
    inits$theta <- 1:2
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    
    inits$tau2 <- rep(5, ncol(Y)- 1)
    inits$theta <- matrix(NA, ncol(Y) - 1, 2)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.")
    inits$theta <- matrix("aaa", ncol(Y)-1, 2)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If theta is specified in inits, it must be a numeric matrix with rows equal to the number of columns of Y - 1 and 2 columns when corr_fun = \"matern\" and shared_covariance_params = FALSE.")
    inits$theta <- matrix(1:(ncol(Y) * 2 - 2), ncol(Y) -1, 2)
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    
    
    inits$theta <- 1
    inits$tau2 <- NA
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- "aaa"
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- -1
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    
    inits$theta <- rep(5, ncol(Y) - 1)
    inits$tau2 <- rep(NA, ncol(Y)-1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- rep("aaa", ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- -c(1:(ncol(Y) - 1))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    
    inits$theta <- rep(5, 2)
    inits$tau2 <- NA
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- rep("aaa", 2)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- -1
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If tau2 is specified in inits, it must be a positive numeric value when shared_covariance_params = TRUE.")
    inits$tau2 <- 1
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    
    inits$theta <- matrix(5, ncol(Y) - 1, 2)
    inits$tau2 <- rep(NA, ncol(Y)- 1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- rep("aaa", ncol(Y)-1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- -c(1:(ncol(Y) - 1))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If tau2 is specified in inits, it must be a positive numeric vector of length equal to the number of columns Y - 1 when shared_covariance_params = TRUE.")
    inits$tau2 <- 1:(ncol(Y) - 1)
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))

    inits <- NULL
    inits$rho <- 0.5
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    inits$rho <- -0.5
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    
    inits$rho <- -1.1
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    
    inits$rho <- 1.1
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    
    inits$rho <- NA
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    
    inits$rho <- "aaa"
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    
    inits$rho <- 1:5
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE), "If rho is specified in inits, it must be a numeric value between -1 and 1.")
    
    inits <- NULL
    inits$eta <- array(1, dim = c(dim(Y)[1], dim(Y)[2] - 1, dim(Y)[3]))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    expect_silent(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    
    inits$eta <- c(array(1, dim = c(dim(Y)[1], dim(Y)[2] - 1, dim(Y)[3])))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))

    inits$eta <- array(1, dim = c(dim(Y)[1], dim(Y)[2] - 1, dim(Y)[3], 5))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    
    inits$eta <- list(1)
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "exponential", shared_covariance_params = FALSE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = TRUE))
    expect_error(check_inits_pg_stlm(Y, X, locs, inits, corr_fun = "matern", shared_covariance_params = FALSE))
    
    
    # expect_identical(
    #     {
    #         set.seed(111)
    #         inits  <- default_inits_pg_lm(Y, X, priors)
    #     }, 
    #     {
    #         set.seed(111)
    #         inits <- list(
    #             beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
    #         )
    #         
    #     }
    # )
    
    ## check inits for pg_spvlm

})    

test_that("check_input_pg_lm", {
    
    Y <- matrix(1:40, 10, 4)
    X <- matrix(1:20, 10, 2)
    expect_silent(check_input_pg_lm(Y, X))
    
    Y <- rnorm(10)
    expect_error(check_input_pg_lm(Y, X), "Y must be an integer matrix.")
    Y <- matrix(1:40, 10, 4)
    Y[1, 1] <- NA
    expect_error(check_input_pg_lm(Y, X), "Y must be an integer matrix.")
    Y <- matrix(1:20, 5, 4)
    expect_error(check_input_pg_lm(Y, X), "Y and X must have the same number of rows.")
    Y <- matrix(1:20, 5, 4)
    X <- matrix(1:30, 15, 2)
    expect_error(check_input_pg_lm(Y, X), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(check_input_pg_lm(Y, X), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, ] <- 0
    expect_error(check_input_pg_lm(Y, X), "There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")

})

test_that("check_input_pg_splm", {
    Y <- matrix(1:40, 10, 4)
    X <- matrix(1:20, 10, 2)
    locs <- matrix(runif(20), 10, 2)
    expect_silent(check_input_pg_splm(Y, X, locs))
    
    Y <- rnorm(10)
    expect_error(check_input_pg_splm(Y, X, locs), "Y must be an integer matrix.")
    Y <- matrix(1:40, 10, 4)
    Y[1, 1] <- NA
    expect_error(check_input_pg_splm(Y, X, locs), "Y must be an integer matrix.")
    Y <- matrix(1:20, 5, 4)
    expect_error(check_input_pg_splm(Y, X, locs), "Y and X must have the same number of rows.")
    Y <- matrix(1:20, 5, 4)
    X <- matrix(1:30, 15, 2)
    expect_error(check_input_pg_splm(Y, X, locs), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(check_input_pg_splm(Y, X, locs), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, ] <- 0
    expect_error(check_input_pg_splm(Y, X, locs), "There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    Y <- matrix(1:40, 10, 4)
    locs <- matrix(1:30, 15, 2)
    expect_error(check_input_pg_splm(Y, X, locs), "Y and locs must have the same number of rows.")
    locs <- matrix(1:30, 10, 3)
    expect_error(check_input_pg_splm(Y, X, locs), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    locs <- matrix("aaa", 10, 2)
    expect_error(check_input_pg_splm(Y, X, locs), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")

})

test_that("check_input_pg_stlm", {
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(1:20, 10, 2)
    locs <- matrix(runif(20), 10, 2)
    expect_error(check_input_pgSTLM(Y, X, locs), "The function check_input_pgSTLM\\(\\) is now deprecated. Please use check_input_pg_stlm\\(\\).")
})

test_that("check_input_pg_stlm", {
    
    
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(1:20, 10, 2)
    locs <- matrix(runif(20), 10, 2)
    expect_silent(check_input_pg_stlm(Y, X, locs))
    expect_error(check_input_pgSTLM(Y, X, locs, "The function check_input_pgSTLM() is now deprecated. Please use check_input_pg_stlm()."))
    Y <- array(1:200, dim = c(5, 4, 5, 2))
    expect_error(check_input_pg_stlm(Y, X, locs), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")

    Y <- array(1:200, dim = c(10, 4, 5))
    # Y[1, 1, 1] <- NA
    # expect_error(check_input_pg_stlm(Y, X, locs), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y[1, 1, 1] <- 0.5
    expect_error(check_input_pg_stlm(Y, X, locs), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y[1, 1, 2] <- NA
    expect_error(check_input_pg_stlm(Y, X, locs), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- matrix(1:20, 5, 4)
    expect_error(check_input_pg_stlm(Y, X, locs), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(1:30, 15, 2)
    expect_error(check_input_pg_stlm(Y, X, locs), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(check_input_pg_stlm(Y, X, locs), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, , 1] <- 0
    expect_error(check_input_pg_stlm(Y, X, locs), "There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    Y <- array(1:200, dim = c(10, 4, 5))
    locs <- matrix(1:30, 15, 2)
    expect_error(check_input_pg_stlm(Y, X, locs), "Y and locs must have the same number of rows.")
    locs <- matrix(1:30, 10, 3)
    expect_error(check_input_pg_stlm(Y, X, locs), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    locs <- matrix("aaa", 10, 2)
    expect_error(check_input_pg_stlm(Y, X, locs), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    
    
})

# test_that("check_input_pt_svlm", {
# 
#     ## check input for pg_spvlm
# 
# })

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







