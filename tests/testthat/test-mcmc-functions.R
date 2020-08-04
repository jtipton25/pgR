context("mcmc functions")

test_that("pgLM", {
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    expect_error(pgLM(Y, X), "The function pgLM\\(\\) has been deprecated. Please use pg_lm\\(\\) instead.")
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