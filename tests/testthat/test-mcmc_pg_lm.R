# pgLM -- deprecated -----------------------------------------------------------

test_that("pgLM", {
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    params <- default_params()
    expect_error(pgLM(Y, X, params), "The function pgLM\\(\\) has been deprecated. Please use pg_lm\\(\\) instead.")
})

# pg_lm ------------------------------------------------------------------------

test_that("pg_lm", {
    set.seed(2021)
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    expect_error(pg_lm(Y, X), 'argument "params" is missing, with no default')
    params <- default_params()
    expect_error(pg_lm(Y, X, params), 'argument "priors" is missing, with no default')
    priors <- default_priors_pg_lm(Y, X)
    suppressMessages(out <- pg_lm(Y, X, params, priors))
    capture_output(expect_snapshot_value(out, style = "serialize"))
    expect_true(class(out) == "pg_lm")
    
    ## check inputs
    expect_error(pg_lm(Y, X, params, priors, n_cores = -2), "n_cores must be a positive integer")
    Y <- rnorm(10)
    expect_error(pg_lm(Y, X, params, priors), "Y must be an integer matrix.")
    Y <- matrix(1:40, 10, 4)
    Y[1, 1] <- NA
    expect_error(pg_lm(Y, X, params, priors), "Y must be an integer matrix.")
    Y <- matrix(1:20, 5, 4)
    expect_error(pg_lm(Y, X, params, priors), "Y and X must have the same number of rows.")
    Y <- matrix(1:20, 5, 4)
    X <- matrix(1:30, 15, 2)
    expect_error(pg_lm(Y, X, params, priors), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(pg_lm(Y, X, params, priors), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, ] <- 0
    expect_error(pg_lm(Y, X, params, priors), "There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
})
