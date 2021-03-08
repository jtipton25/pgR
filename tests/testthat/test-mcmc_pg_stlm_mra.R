# pg_stlm_mra ----------------------------------------------------------------------

test_that("pg_stlm_mra", {
    set.seed(2021)
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- as.matrix(1:10)
    locs <- matrix(runif(20), 10, 2)
    expect_error(pg_stlm_mra(Y, X), 'argument "locs" is missing, with no default')
    params <- default_params()
    params$n_mcmc <- 10
    params$n_adapt <- 10
    expect_error(pg_stlm_mra(Y, X, locs), 'argument "params" is missing, with no default')
    expect_error(suppressMessages(pg_stlm_mra(Y, X, locs, params)), 'argument "priors" is missing, with no default')
    
    priors <- default_priors_pg_stlm_mra(Y, X)
    suppressMessages(out <- pg_stlm_mra(Y, X, locs, params, priors, M = 2, n_coarse_grid = 5))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(inherits(out, "pg_stlm_mra"))
    suppressMessages(out <- pg_stlm_mra(Y, X, locs, params, priors, M = 3, n_coarse_grid = 4))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(inherits(out, "pg_stlm_mra"))
    
    # Check inputs
    expect_error(pg_stlm_mra(Y, X, locs, params, priors, n_cores = -2), "n_cores must be a positive integer")
    Y[1, , ] <- 0.5
    expect_error(
        pg_stlm_mra(Y, X, locs, params, priors),
        "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- array(1:200, dim = c(5, 4, 5, 2))
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    
    Y <- array(1:200, dim = c(10, 4, 5))
    # Y[1, 1, 1] <- NA
    # expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y[1, 1, 1] <- 0.5
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y[1, 1, 2] <- NA
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- matrix(1:20, 5, 4)
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(1:30, 15, 2)
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, , 1] <- 0
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    Y <- array(1:200, dim = c(10, 4, 5))
    locs <- matrix(1:30, 15, 2)
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "Y and locs must have the same number of rows.")
    locs <- matrix(1:30, 10, 3)
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    locs <- matrix("aaa", 10, 2)
    expect_error(pg_stlm_mra(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    
})
