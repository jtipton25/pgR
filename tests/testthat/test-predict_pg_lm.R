# predict_pgLM -- deprecated ---------------------------------------------------

test_that("predict_pgLM", {
    set.seed(111)
    Y <- matrix(1:40, 10, 4)
    X <- cbind(1, as.matrix(1:10))
    params <- default_params()
    params$n_adapt <- 5
    params$n_mcmc <- 5
    priors <- default_priors_pg_lm(Y, X)
    suppressMessages(out <- pg_lm(Y, X, params, priors))
    # check successful output
    expect_error(predict_pgLM(out, X), "predict_pgLM\\(\\) has been deprecated. Please use predict_pg_lm\\(\\) instead.")
    
})


# predict_pg_lm ----------------------------------------------------------------

test_that("predict_pg_lm", {
    set.seed(111)
    Y <- matrix(1:40, 10, 4)
    X <- cbind(1, as.matrix(1:10))
    params <- default_params()
    priors <- default_priors_pg_lm(Y, X)
    suppressMessages(out <- pg_lm(Y, X, params, priors))
    # check successful output
    capture.output(expect_snapshot_value(predict_pg_lm(out, X), style = "serialize"))
    
    expect_error(predict_pg_lm(out, matrix("aa", 10, 2)), "X_pred must be a numeric matrix")
    expect_error(predict_pg_lm(out, cbind(X, 1)), "The number of colums of X_pred must be equal to the number of columns of beta in the object out")
    class(out) <- "aa"
    expect_error(predict_pg_lm(out, X), "The MCMC object out must be of class pg_lm which is the output of the pg_lm\\(\\) function.")
})
