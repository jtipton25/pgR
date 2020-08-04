context("test check functions")

test_that("check_inits", {
    
    ## Check initial conditions for pg_lm
    expect_error(check_inits_pgLM(), "The function check_inits_pgLM\\(\\) has been deprecated. Please use check_inits_pg_lm\\(\\) instead.")

    Y <- matrix(1:40, 10, 4)
    X <- matrix(rnorm(10), 10, 1)
    params <- default_params()
    priors <- default_priors_pg_lm(Y, X)
    inits  <- default_inits_pg_lm(Y, X, priors)
    expect_silent(check_inits_pg_lm(Y, X, inits))

    X <- matrix(rnorm(20), 10, 2)    
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
    
    ## check inits for pg_splm
    
    ## check inits for pg_spvlm
    
    ## check inits for pg_stlm
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


test_that("check_priors", {
    
})

test_that("check_corr_function", {
    
})




