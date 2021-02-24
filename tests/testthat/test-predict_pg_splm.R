# predict_pgSPLM -- deprecated -------------------------------------------------

test_that("predict_pgSPLM", {
    set.seed(111)
    Y <- matrix(1:40, 10, 4)
    X <- cbind(1, as.matrix(1:10))
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    X_pred <- cbind(1, rnorm(5))
    locs_pred <- matrix(runif(10), 5, 2)
    
    # check the different MCMC settings
    suppressMessages(out <- pg_splm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = TRUE))
    expect_error(predict_pgSPLM(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE, progress = FALSE), "predict_pgSPLM\\(\\) is deprecated. Please use predict_pg_splm\\(\\) instead.")
})

# predict_pg_splm --------------------------------------------------------------

test_that("predict_pg_splm", {
    set.seed(111)
    Y <- matrix(1:40, 10, 4)
    X <- cbind(1, as.matrix(1:10))
    locs <- matrix(runif(20), 10, 2)
    params <- default_params()
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    X_pred <- cbind(1, rnorm(5))
    locs_pred <- matrix(runif(10), 5, 2)
    
    # check the different MCMC settings
    suppressMessages(out <- pg_splm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = TRUE))
    capture.output(expect_snapshot_value(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE, progress = FALSE), style = "serialize"))
    
    suppressMessages(out <- pg_splm(Y, X, locs, params, priors, corr_fun = "exponential", shared_covariance_params = FALSE))
    capture.output(expect_snapshot_value(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = FALSE, progress = FALSE), style = "serialize"))
    
    priors <- default_priors_pg_splm(Y, X, corr_fun = "matern")
    suppressWarnings(suppressMessages(out <- pg_splm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = TRUE)))
    capture.output(expect_snapshot_value(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "matern", shared_covariance_params = TRUE, progress = FALSE), style = "serialize"))
    
    suppressWarnings(suppressMessages(out <- pg_splm(Y, X, locs, params, priors, corr_fun = "matern", shared_covariance_params = FALSE)))
    capture.output(expect_snapshot_value(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "matern", shared_covariance_params = FALSE, progress = FALSE), style = "serialize"))
    
    # check the MCMC inputs
    expect_error(predict_pg_splm(out, cbind(X, 1), X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "The number of colums of X must be equal to the number of columns of beta in the object out")
    expect_error(predict_pg_splm(out, X, cbind(X_pred, 1), locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "The number of colums of X_pred must be equal to the number of columns of beta in the object out")
    expect_error(predict_pg_splm(out, X, X_pred, cbind(locs, 1), locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "locs must be a numeric matrix with 2 columns")
    expect_error(predict_pg_splm(out, X, X_pred, locs, cbind(locs_pred, 1), corr_fun = "exponential", shared_covariance_params = TRUE), "locs_pred must be a numeric matrix with 2 columns")
    expect_error(predict_pg_splm(out, rbind(X, 1), X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "The number of rows of X must equal the number of rows of locs")
    expect_error(predict_pg_splm(out, X, rbind(X_pred, 1), locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "X_pred and locs_pred must be numeric matrices with the same number of rows")
    expect_error(predict_pg_splm(out, X, X_pred, rbind(locs, 1), locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "The number of rows of X must equal the number of rows of locs")
    expect_error(predict_pg_splm(out, X, X_pred, locs, rbind(locs_pred, 1), corr_fun = "exponential", shared_covariance_params = TRUE), "X_pred and locs_pred must be numeric matrices with the same number of rows")
    
    X[1] <- NA
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "X must be a numeric matrix")
    X[1] <- "aa"
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "X must be a numeric matrix")
    X <- cbind(1, as.matrix(1:10))
    X_pred[1] <- NA
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "X_pred must be a numeric matrix")
    X_pred[1] <- "aa"
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "X_pred must be a numeric matrix")
    X_pred <- cbind(1, rnorm(5))
    locs[1] <- NA
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "locs must be a numeric matrix with 2 columns")
    locs[1] <- "aa"
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "locs must be a numeric matrix with 2 columns")
    locs <- cbind(1, as.matrix(1:10))
    locs_pred[1] <- NA
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "locs_pred must be a numeric matrix with 2 columns")
    locs_pred[1] <- "aa"
    expect_error(predict_pg_splm(out, X, X_pred, locs, locs_pred, corr_fun = "exponential", shared_covariance_params = TRUE), "locs_pred must be a numeric matrix with 2 columns")
    locs_pred <- matrix(rnorm(10), 5, 2)
    
    
    class(out) <- "aa"
    expect_error(predict_pg_splm(out, X), "The MCMC object out must be of class pg_splm which is the output of the pg_splm\\(\\) function.")
    
})
