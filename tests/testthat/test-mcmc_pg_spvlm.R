# pgSPVLM -- deprecated --------------------------------------------------------

# test_that("pgSPVLM", {
#     Y <- matrix(1:40, 10, 4)
#     X <- as.matrix(1:10)
#     params <- default_params()
#     expect_error(pgSPVLM(Y, X, params), "The function pgSPVLM\\(\\) has been deprecated. Please use pg_spvlm\\(\\) instead.")
# })

# pg_spvlm ---------------------------------------------------------------------

# test_that("pg_spvlm", {
#     set.seed(2021)
#     Y <- matrix(1:40, 10, 4)
#     X <- as.matrix(1:10)
#     locs <- matrix(runif(20), 10, 2)
#     expect_error(pg_spvlm(Y, X), 'argument "locs" is missing, with no default')
#     params <- default_params()
#     expect_error(pg_spvlm(Y, X, locs), 'argument "params" is missing, with no default')
#     params <- default_params()
#     expect_error(pg_spvlm(Y, X, locs, params), 'argument "priors" is missing, with no default')
#     
#     priors <- default_priors_pg_spvlm(Y, X, corr_fun = "exponential")
#     suppressMessages(out <- pg_spvlm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = TRUE))
#     expect_snapshot_value(out, style = "serialize")
#     expect_true(class(out) == "pg_spvlm")
#     suppressMessages(out <- pg_spvlm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = FALSE))
#     expect_snapshot_value(out, style = "serialize")
#     expect_true(class(out) == "pg_spvlm")
#     
#     priors <- default_priors_pg_spvlm(Y, X, corr_fun = "matern")
#     suppressMessages(out <- pg_spvlm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = TRUE))
#     expect_snapshot_value(out, style = "serialize")
#     expect_true(class(out) == "pg_spvlm")
#     suppressWarnings(suppressMessages(out <- pg_spvlm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = FALSE)))
#     expect_snapshot_value(out, style = "serialize")
#     expect_true(class(out) == "pg_spvlm")
#     
#     ## check inputs into pg_spvlm function
#     expect_error(pg_spvlm(Y, X, locs, params, priors, n_cores = -2), "n_cores must be a positive integer")
#     Y <- rnorm(10)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "Y must be an integer matrix.")
#     Y <- matrix(1:40, 10, 4)
#     Y[1, 1] <- NA
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "Y must be an integer matrix.")
#     Y <- matrix(1:20, 5, 4)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "Y and X must have the same number of rows.")
#     Y <- matrix(1:20, 5, 4)
#     X <- matrix(1:30, 15, 2)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "Y and X must have the same number of rows.")
#     X <- matrix("aaa", 10, 2)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "X must be a numeric matrix.")
#     X <- matrix(1:20, 10, 2)
#     Y[1, ] <- 0
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
#     Y <- matrix(1:40, 10, 4)
#     locs <- matrix(1:30, 15, 2)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "Y and locs must have the same number of rows.")
#     locs <- matrix(1:30, 10, 3)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
#     locs <- matrix("aaa", 10, 2)
#     expect_error(pg_spvlm(Y, X, locs, params, priors), "locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
# })
