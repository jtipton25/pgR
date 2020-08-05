context("testing functions")


test_that("checking default settings of params -- keep updated as model list grows", {
    expect_identical(
        default_params(),
        list(
            n_mcmc    = 100L,
            n_adapt   = 100L,
            n_thin    = 1L,
            n_message = 500L
        )
    )
})


test_that("checking default settings of priors", {
    
    Y <- stats::rmultinom(100, 50, rep(1/5, 5))
    X <- as.matrix(stats::rnorm(100))
    
    expect_error(default_priors_pgLM(Y, X), "default_priors_pgLM\\(\\) has been deprecated. Please use default_priors_pg_lm\\(\\) instead.")

    expect_identical(    
        default_priors_pg_lm(Y, X),
        list(
            mu_beta    = 0,
            Sigma_beta = as.matrix(1)
        )
    )
    
    Y <- stats::rmultinom(100, 50, rep(1/5, 5))
    X <- matrix(stats::rnorm(100 * 50), 100, 50)
    
    expect_identical(    
        default_priors_pg_lm(Y, X),
        list(
            mu_beta    = rep(0, 50),
            Sigma_beta = diag(50)
        )
    )
    
    expect_error(default_priors_pgSPLM(Y, X, corr_fun = "exponential"), "default_priors_pgSPLM\\(\\) has been deprecated. Please use default_priors_pg_splm\\(\\) instead.")
    expect_equal(
        default_priors_pg_splm(Y, X, corr_fun = "exponential"),
        list(        
            mu_beta    = rep(0, 50),
            Sigma_beta = diag(50),
            alpha_tau  = 0.1,
            beta_tau   = 0.1,
            mean_range = 0,
            sd_range   = 10
        )
    )
    expect_equal(
        default_priors_pg_splm(Y, X, corr_fun = "matern"),
        list(        
            mu_beta    = rep(0, 50),
            Sigma_beta = diag(50),
            alpha_tau  = 0.1,
            beta_tau   = 0.1,
            mean_nu    = -1,
            sd_nu      = 1,
            mean_range = 0,
            sd_range   = 10
        )
    )
 
    expect_error(
        default_priors_pg_splm(Y, X, corr_fun = "adfs"),
        'corr_fun must be either "matern" or "exponential"'
    )
    
})


test_that("checking default settings of inits", {

    Y      <- stats::rmultinom(100, 50, rep(1/5, 5))
    X      <- matrix(stats::rnorm(100 * 50), 100, 50)
    priors <- default_priors_pg_lm(Y, X)
    
    expect_error(default_inits_pgLM(Y, X, priors), "default_inits_pgLM\\(\\) has been deprecated. Please use default_inits_pg_lm\\(\\) instead.")
    
    expect_identical(
        {
            set.seed(111)
            default_inits_pg_lm(Y, X, priors)
        },
        {
            set.seed(111)
            list(beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)))
        }
    )

    Y      <- stats::rmultinom(100, 50, rep(1/5, 5))
    X      <- matrix(stats::rnorm(100 * 50), 100, 50)
    priors <- default_priors_pg_splm(Y, X, corr_fun = "exponential")
    
    expect_error(default_inits_pgSPLM(Y, X, priors), "default_inits_pgSPLM\\(\\) has been deprecated. Please use default_inits_pg_splm\\(\\) instead.")
    
    expect_identical(
        {
            set.seed(111)
            default_inits_pg_splm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = TRUE)
        },
        {
            set.seed(111)
            theta_mean <- priors$mean_range
            theta_var  <- priors$sd_range^2
            
            list(
                beta  = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)),
                tau2  = min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10),
                theta = pmin(pmax(stats::rnorm(1, theta_mean, sqrt(theta_var)), -1), 0.1)
            )
        }
    )
    
    expect_identical(
        {
            set.seed(111)
            default_inits_pg_splm(Y, X, priors, corr_fun = "exponential", shared_covariance_params = FALSE)
        },
        {
            set.seed(111)
            J <- ncol(Y)
            theta_mean <- priors$mean_range
            theta_var  <- priors$sd_range^2
            
            list(
                beta  = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)),
                tau2  = pmin(1 / stats::rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10),
                theta = pmin(pmax(stats::rnorm(J-1, theta_mean, sqrt(theta_var)), -1), 0.1)
            )
        }
    )
    
    
    
    
    priors <- default_priors_pg_splm(Y, X, corr_fun = "matern")
    expect_error(default_inits_pgSPLM(Y, X, priors), "default_inits_pgSPLM\\(\\) has been deprecated. Please use default_inits_pg_splm\\(\\) instead.")
    expect_identical(
        {
            set.seed(111)
            default_inits_pg_splm(Y, X, priors, corr_fun = "matern", shared_covariance_params = TRUE)
        },
        {
            set.seed(111)
            theta_mean <- c(priors$mean_range, priors$mean_nu)
            theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
            
            list(
                beta  = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)),
                tau2  = min(1 / stats::rgamma(1, priors$alpha_tau, priors$beta_tau), 10),
                theta =                 as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -1), 0.1))
            )
        }
    )
    
    expect_identical(
        {
            set.seed(111)
            default_inits_pg_splm(Y, X, priors, corr_fun = "matern", shared_covariance_params = FALSE)
        },
        {
            set.seed(111)
            J <- ncol(Y)
            theta_mean <- c(priors$mean_range, priors$mean_nu)
            theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
            
            list(
                beta  = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta)),
                tau2  = pmin(1 / stats::rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10),
                theta = pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -1), 0.1)
            )
        }
    )
    
})



test_that("make_alpha function", {
    Xbs <- matrix(1:8, 2, 4)
    beta <- matrix(4:-7, 4, 3)
    params <- default_params()
    params$link = "log"
    expect_equal(make_alpha(Xbs, beta, params), exp(Xbs %*% beta))
    params$link = "tobit"
    expect_equal(make_alpha(Xbs, beta, params), pmax(Xbs %*% beta, params$tol))
})

# test_that("setup_splines function", {
#     params <- default_params()
#     params$type = "GP"
#     expect_error(setup_splines(params), "the model type must be bspline to run setup_splines")
# })
# expect_error(make_alpha(Y, X, params))




if (require(LaplacesDemon)) {
    test_that("dvmn_arma_mc function", {
        ## add in error checking and type checking
        N     <- 100
        d     <- 6
        mu    <- stats::rnorm(6)
        Sigma <- LaplacesDemon::rwishart(d+2, diag(d))
        X     <- mvnfast::rmvn(N, mu, Sigma)
        
        expect_equal(
            mvnfast::dmvn(X, mu, Sigma, log = TRUE),
            as.vector(dmvnrm_arma_mc(X, mu, Sigma, logd = TRUE))
        )
        expect_equal(
            mvnfast::dmvn(X, mu, Sigma, log = TRUE),
            as.vector(dmvnrm_arma_mc(X, mu, Sigma, logd = TRUE, cores = 4))
        )
        
        # expect_error(setup_splines(params), "the model type must be bspline to run setup_splines")
    })
}


if (require(LaplacesDemon)) {
    test_that("rmvn_arma function", {
        ## check that the rmvn_arma function generates draws from the
        ## same distribution -- verify Monte Carlo mean and covariance
        ## are the same up to Monte Carlo error
        set.seed(11)
        N     <- 10000
        d     <- 6
        mu    <- stats::rnorm(6)
        Sigma <- LaplacesDemon::rwishart(d + 2, diag(d))
        
        X <- mvnfast::rmvn(N, solve(Sigma) %*% mu, solve(Sigma))
        Y <- t(sapply(1:N, function(i) rmvn_arma(Sigma, mu)))
        expect_equal(
            colMeans(X), colMeans(Y), tol = 0.05
        )
        
        expect_equal(
            var(X), var(Y), tol = 0.05
        )
    })
}
