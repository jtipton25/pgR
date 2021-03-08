

# pgSTLM -- deprecated ---------------------------------------------------------

test_that("pgSTLM", {
    Y <- array(1:400, dim = c(10, 4, 10))
    X <- as.matrix(1:10)
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "exponential")
    expect_error(pgSTLM(Y, X, locs, params, priors), "The function pgSTLM\\(\\) has been deprecated. Please use pg_stlm\\(\\) instead")
})


# pg_stlm ----------------------------------------------------------------------

test_that("pg_stlm", {
    set.seed(2021)
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- as.matrix(1:10)
    locs <- matrix(runif(20), 10, 2)
    expect_error(pg_stlm(Y, X), 'argument "locs" is missing, with no default')
    params <- default_params()
    params$n_mcmc <- 10
    params$n_adapt <- 10
    expect_error(pg_stlm(Y, X, locs), 'argument "params" is missing, with no default')
    expect_error(suppressMessages(pg_stlm(Y, X, locs, params)), 'argument "priors" is missing, with no default')
    
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "exponential")
    suppressMessages(out <- pg_stlm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = TRUE))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(inherits(out, "pg_stlm"))
    suppressMessages(out <- pg_stlm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = FALSE))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(inherits(out, "pg_stlm"))
    
    priors <- default_priors_pg_stlm(Y, X, corr_fun = "matern")
    suppressMessages(out <- pg_stlm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = TRUE))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(class(out) == "pg_stlm")
    suppressMessages(out <- pg_stlm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = FALSE))
    # suppressMessages(expect_snapshot_value(out, style = "serialize"))
    expect_true(class(out) == "pg_stlm")
    
    # Check inputs
    expect_error(pg_stlm(Y, X, locs, params, priors, n_cores = -2), "n_cores must be a positive integer")
    Y[1, , ] <- 0.5
    expect_error(
        pg_stlm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = FALSE),
        "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- array(1:200, dim = c(5, 4, 5, 2))
    expect_error(pg_stlm(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    
    Y <- array(1:200, dim = c(10, 4, 5))
    # Y[1, 1, 1] <- NA
    # expect_error(pg_stlm(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y[1, 1, 1] <- 0.5
    expect_error(pg_stlm(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y[1, 1, 2] <- NA
    expect_error(pg_stlm(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- matrix(1:20, 5, 4)
    expect_error(pg_stlm(Y, X, locs, params, priors), "Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    Y <- array(1:200, dim = c(10, 4, 5))
    X <- matrix(1:30, 15, 2)
    expect_error(pg_stlm(Y, X, locs, params, priors), "Y and X must have the same number of rows.")
    X <- matrix("aaa", 10, 2)
    expect_error(pg_stlm(Y, X, locs, params, priors), "X must be a numeric matrix.")
    X <- matrix(1:20, 10, 2)
    Y[1, , 1] <- 0
    expect_error(pg_stlm(Y, X, locs, params, priors), "There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    Y <- array(1:200, dim = c(10, 4, 5))
    locs <- matrix(1:30, 15, 2)
    expect_error(pg_stlm(Y, X, locs, params, priors), "Y and locs must have the same number of rows.")
    locs <- matrix(1:30, 10, 3)
    expect_error(pg_stlm(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    locs <- matrix("aaa", 10, 2)
    expect_error(pg_stlm(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    
})
