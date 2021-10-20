# pg_splm_mra ------------------------------------------------------------------

test_that("pg_splm_mra", {
    set.seed(2021)
    Y <- matrix(1:40, 10, 4)
    X <- as.matrix(1:10)
    locs <- matrix(runif(20), 10, 2)
    expect_error(pg_splm_mra(Y, X), 'argument "locs" is missing, with no default')
    params <- default_params()
    expect_error(pg_splm_mra(Y, X, locs), 'argument "params" is missing, with no default')
    params <- default_params()
    suppressMessages(expect_error(pg_splm_mra(Y, X, locs, params), 'argument "priors" is missing, with no default'))
    params$n_adapt <- 50
    params$n_mcmc <- 50
    
    priors <- default_priors_pg_splm_mra(Y, X)
    suppressMessages(out <- pg_splm_mra(Y, X, locs, params, priors, M = 2, n_coarse_grid = 4))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(class(out) == "pg_splm_mra")
    suppressMessages(out <- pg_splm_mra(Y, X, locs, params, priors, M = 3, n_coarse_grid = 6))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(class(out) == "pg_splm_mra")
    
    ## check inputs into pg_splm_mra function
    
    expect_error(pg_splm_mra(Y, X, locs, params, priors, M = -2), "the number of resolutions M must be a positive integer")
    expect_error(pg_splm_mra(Y, X, locs, params, priors, n_coarse_grid = -2), "n_coarse_grid must be a positive integer")
    expect_error(pg_splm_mra(Y, X, locs, params, priors, n_cores = -2), "n_cores must be a positive integer")
    expect_error(pg_splm_mra(Y, X, locs, params, priors, use_spam = -2), "use_spam must be either TRUE or FALSE")
    expect_error(pg_splm_mra(Y, X, locs, params, priors, use_spam = FALSE), "The only sparse matrix pacakage available is spam")
    
    Y <- rnorm(10)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "Y must be an integer matrix.")
    Y <- matrix(1:40, 10, 4)
    Y[1, 1] <- NA
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "Y must be an integer matrix.")
    Y <- matrix(1:20, 5, 4)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "Y and X must have the same number of rows.")
    Y <- matrix(1:20, 5, 4)
    X <- matrix(1:30, 15, 2)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, ] <- 0
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    Y <- matrix(1:40, 10, 4)
    locs <- matrix(1:30, 15, 2)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "Y and locs must have the same number of rows.")
    locs <- matrix(1:30, 10, 3)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    locs <- matrix("aaa", 10, 2)
    expect_error(pg_splm_mra(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
})
