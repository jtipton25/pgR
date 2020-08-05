context("mcmc functions")

test_that("pgLM", {
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    params <- default_params()
    expect_error(pgLM(Y, X, params), "The function pgLM\\(\\) has been deprecated. Please use pg_lm\\(\\) instead.")
})


test_that("pg_lm", {
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    expect_error(pg_lm(Y, X), 'argument "params" is missing, with no default')
    params <- default_params()
    expect_error(pg_lm(Y, X, params), 'argument "priors" is missing, with no default')
    priors <- default_priors_pg_lm(Y, X)
    out <- pg_lm(Y, X, params, priors)
    expect_true(class(out) == "pg_lm")
    
})


test_that("pgSPLM", {
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    expect_error(pgSPLM(Y, X, locs, params), "The function pgSPLM\\(\\) has been deprecated. Please use pg_splm\\(\\) instead.")
})


test_that("pg_splm", {
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    locs <- matrix(runif(20), 10, 2)
    expect_error(pg_splm(Y, X), 'argument "locs" is missing, with no default')
    params <- default_params()
    expect_error(pg_splm(Y, X, locs), 'argument "params" is missing, with no default')
    params <- default_params()
    expect_error(pg_splm(Y, X, locs, params), 'argument "priors" is missing, with no default')
    
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    out <- pg_splm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = TRUE)
    expect_true(class(out) == "pg_splm")
    out <- pg_splm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = FALSE)
    expect_true(class(out) == "pg_splm")
    
    priors <- default_priors_pg_splm(Y, X, corr_fun = "matern")
    out <- pg_splm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = TRUE)
    expect_true(class(out) == "pg_splm")
    out <- pg_splm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = FALSE)
    expect_true(class(out) == "pg_splm")
})

