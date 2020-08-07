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
        expect_equal(
            mvnfast::dmvn(X, mu, Sigma, log = FALSE),
            as.vector(dmvnrm_arma_mc(X, mu, Sigma, logd = FALSE, cores = 4))
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

test_that("log_sum_exp", {
    x <- 1:5
    expect_equal(logsumexp(x), 5.45191439593759)
    x <- -(1:5)
    expect_equal(logsumexp(x), -0.548085604062407)
    x <- rep("aaa", 5)
    expect_error(logsumexp(x), "x must be a numeric vector")
    x <- c(1:5, NA)
    expect_error(logsumexp(x), "x must be a numeric vector")
})

test_that("log_sum_exp", {
    x <- 1:5
    expect_equal(softmax(x), c(0.0116562309560396, 0.0316849207961243, 0.0861285444362687, 
                               0.234121657252737, 0.636408646558831))
    x <- -(1:5)
    expect_equal(softmax(x), c(0.636408646558831, 0.234121657252737, 0.0861285444362687, 0.0316849207961243, 
                               0.0116562309560396))
    x <- rep("aaa", 5)
    expect_error(softmax(x), "x must be a numeric vector")
    x <- c(1:5, NA)
    expect_error(softmax(x), "x must be a numeric vector")
})


test_that("counts_to_propotions", {
    Y <- matrix(1:20, 5, 4)
    expect_equal(
        counts_to_proportions(Y), 
        structure(c(0.0294117647058824, 0.0526315789473684, 0.0714285714285714, 
                    0.0869565217391304, 0.1, 0.176470588235294, 0.184210526315789, 
                    0.19047619047619, 0.195652173913043, 0.2, 0.323529411764706, 
                    0.315789473684211, 0.30952380952381, 0.304347826086957, 0.3, 
                    0.470588235294118, 0.447368421052632, 0.428571428571429, 0.41304347826087, 
                    0.4), .Dim = 5:4)
    )
    
    Y <- rnorm(1:20)
    expect_error(counts_to_proportions(Y), "Y must be a matrix of integer count observations.")
    Y <- matrix("aaa", 5, 4)
    expect_error(counts_to_proportions(Y), "Y must be a matrix of integer count observations.")
    Y <- array(1:60, dim = c(5, 4, 3))
    expect_error(counts_to_proportions(Y), "Y must be a matrix of integer count observations.")
    Y <- matrix(1:20, 5, 4)
    Y[1, 1] <- NA
    expect_error(counts_to_proportions(Y), "Y must be a matrix of integer count observations.")
    Y[1, 1] <- 0.5
    expect_error(counts_to_proportions(Y), "Y must be a matrix of integer count observations.")
    
})


test_that("expit", {
    x <- 1:5
    expect_equal(pgR::expit(x), c(0.731058578630005, 0.880797077977882, 0.952574126822433, 0.982013790037908, 0.993307149075715))
    x <- rep(NA, 5)
    expect_error(pgR::expit(x), "x must be a numeric value")
    x <- rep("aaa", 5)
    expect_error(pgR::expit(x), "x must be a numeric value")
})

test_that("logit", {
    set.seed(111)
    x <- runif(5)
    expect_equal(pgR::logit(x), c(0.376303669311382, 0.976841852046962, -0.530406842159003, 0.0597130568639106, -0.499478564685536))
    x <- rep(NA, 5)
    expect_error(pgR::logit(x), "p must be a numeric value between 0 and 1")
    x <- rep("aaa", 5)
    expect_error(pgR::logit(x), "p must be a numeric value between 0 and 1")
    x <- rep(5, 5)
    expect_error(pgR::logit(x), "p must be a numeric value between 0 and 1")
    x <- list(runif(5))
    expect_error(pgR::logit(x), "p must be a numeric value between 0 and 1")
})


test_that("eta_to_pi", {
    eta <- 1:5
    expect_equal(
        eta_to_pi(eta), 
        structure(c(0.731058578630005, 0.880797077977882, 0.952574126822433, 
                    0.982013790037908, 0.993307149075715, 0.268941421369995, 0.119202922022118, 
                    0.0474258731775666, 0.0179862099620915, 0.00669285092428473), 
                  .Dim = c(5L, 2L))
    )
    eta <- -c(1:5)
    expect_equal(
        eta_to_pi(eta), 
        structure(c(0.268941421369995, 0.119202922022118, 0.0474258731775668, 
                    0.0179862099620916, 0.00669285092428486, 0.731058578630005, 0.880797077977882, 
                    0.952574126822433, 0.982013790037908, 0.993307149075715), 
                  .Dim = c(5L, 2L))
    )
    
    eta <- matrix(1:10, 5, 2)
    expect_equal(
        eta_to_pi(eta), 
        structure(c(0.731058578630005, 0.880797077977882, 0.952574126822433, 
                    0.982013790037908, 0.993307149075715, 0.268276430583737, 0.119094322057633, 
                    0.0474099689048091, 0.0179839905613397, 0.00669254708311723, 
                    0.000664990786257702, 0.000108599964484299, 1.59042727575776e-05, 
                    2.21940075187388e-06, 3.0384116750555e-07), 
                  .Dim = c(5L, 3L))
    )
    
    eta <- array(1:40, dim = c(5, 4, 2))
    expect_error(eta_to_pi(eta), "eta must be either a numeric vector or a numeric matrix.")
    eta <- rep(NA, 5)
    expect_error(eta_to_pi(eta), "eta must be either a numeric vector or a numeric matrix.")
    eta <- c(rep(5, 5), NA)
    expect_error(eta_to_pi(eta), "eta must be either a numeric vector or a numeric matrix.")
    eta <- matrix(NA, 5)
    expect_error(eta_to_pi(eta), "eta must be either a numeric vector or a numeric matrix.")
    eta <- list(1:5)
    expect_error(eta_to_pi(eta), "eta must be either a numeric vector or a numeric matrix.")
    eta <- rep("aaa", 5)
    expect_error(eta_to_pi(eta), "eta must be either a numeric vector or a numeric matrix.")
    
})

test_that("correlation_function", {
    set.seed(111)
    locs <- matrix(runif(20), 10, 2)
    D <- fields::rdist(locs)
    expect_equal(
        correlation_function(D, 1, "exponential"),
        structure(c(1, 0.950544242987673, 0.820757332543112, 0.827653703967322, 
                    0.846216556777341, 0.926998651622854, 0.773617821198024, 0.858345671681754, 
                    0.897763830064757, 0.831148700962243, 0.950544242987673, 1, 0.79232529415432, 
                    0.807126574185252, 0.814774460221953, 0.882412263089181, 0.737055236295147, 
                    0.85574707276131, 0.861281867842839, 0.79218156356849, 0.820757332543112, 
                    0.79232529415432, 1, 0.947767844096066, 0.967662446693328, 0.86880136646441, 
                    0.871273071661644, 0.714491308867832, 0.911720195996076, 0.798017731952593, 
                    0.827653703967322, 0.807126574185252, 0.947767844096066, 1, 0.937628654673637, 
                    0.859864551125443, 0.826110871192016, 0.713097863136091, 0.903506943974864, 
                    0.771183779468682, 0.846216556777341, 0.814774460221953, 0.967662446693328, 
                    0.937628654673637, 1, 0.897796974371921, 0.873602245297216, 0.738242057923419, 
                    0.941522426625363, 0.820097937485444, 0.926998651622854, 0.882412263089181, 
                    0.86880136646441, 0.859864551125443, 0.897796974371921, 1, 0.834514640883412, 
                    0.822114657990748, 0.951037266863276, 0.874166732257841, 0.773617821198024, 
                    0.737055236295147, 0.871273071661644, 0.826110871192016, 0.873602245297216, 
                    0.834514640883412, 1, 0.704808666036871, 0.849334796466678, 0.847201982791776, 
                    0.858345671681754, 0.85574707276131, 0.714491308867832, 0.713097863136091, 
                    0.738242057923419, 0.822114657990748, 0.704808666036871, 1, 0.783427946879552, 
                    0.813094049980482, 0.897763830064757, 0.861281867842839, 0.911720195996076, 
                    0.903506943974864, 0.941522426625363, 0.951037266863276, 0.849334796466678, 
                    0.783427946879552, 1, 0.845929338255793, 0.831148700962243, 0.79218156356849, 
                    0.798017731952593, 0.771183779468682, 0.820097937485444, 0.874166732257841, 
                    0.847201982791776, 0.813094049980482, 0.845929338255793, 1), 
                  .Dim = c(10L, 10L))
    )
        
    expect_equal(
        correlation_function(D, c(1, 2), "matern"),
        structure(c(1, 0.999899342610697, 0.998474659565074, 0.998601045452874, 
                    0.998909686788764, 0.999775187664083, 0.997425962482931, 0.999087523818759, 
                    0.99954499907465, 0.99866263075097, 0.999899342610697, 1, 0.997882316257325, 
                    0.998205336386718, 0.998359686031022, 0.999387886536294, 0.996365629557816, 
                    0.99905095972196, 0.999127846498432, 0.997879018284853, 0.998474659565074, 
                    0.997882316257325, 1, 0.999887398519482, 0.999957719146688, 0.999226382163455, 
                    0.999257306792017, 0.995589163309586, 0.999665828588167, 0.998010412436649, 
                    0.998601045452874, 0.998205336386718, 0.999887398519482, 1, 0.999837725312446, 
                    0.999108513302754, 0.998573334742808, 0.995537924613427, 0.999597201716159, 
                    0.997362470198384, 0.998909686788764, 0.998359686031022, 0.999957719146688, 
                    0.999837725312446, 1, 0.999545310444539, 0.9992857938117, 0.996403779119561, 
                    0.999857936227045, 0.998462232732003, 0.999775187664083, 0.999387886536294, 
                    0.999226382163455, 0.999108513302754, 0.999545310444539, 1, 0.998720398281737, 
                    0.998500050119059, 0.999901390108536, 0.999292602758792, 0.997425962482931, 
                    0.996365629557816, 0.999257306792017, 0.998573334742808, 0.9992857938117, 
                    0.998720398281737, 1, 0.995224894075583, 0.998957162071176, 0.998924822873654, 
                    0.999087523818759, 0.99905095972196, 0.995589163309586, 0.995537924613427, 
                    0.996403779119561, 0.998500050119059, 0.995224894075583, 1, 0.997672153868553, 
                    0.998326488035227, 0.99954499907465, 0.999127846498432, 0.999665828588167, 
                    0.999597201716159, 0.999857936227045, 0.999901390108536, 0.998957162071176, 
                    0.997672153868553, 1, 0.99890525192009, 0.99866263075097, 0.997879018284853, 
                    0.998010412436649, 0.997362470198384, 0.998462232732003, 0.999292602758792, 
                    0.998924822873654, 0.998326488035227, 0.99890525192009, 1), 
                  .Dim = c(10L, 10L), 
                  .Dimnames = list(NULL, NULL))
    )
    
    expect_error(correlation_function(D, 1:5, "exponential"), "theta must be a numeric value for the exponential correlation function.")
    expect_error(correlation_function(D, NA, "exponential"), "theta must be a numeric value for the exponential correlation function.")
    expect_error(correlation_function(D, "aaa", "exponential"), "theta must be a numeric value for the exponential correlation function.")
    expect_error(correlation_function(D, matrix(1:6, 3, 2), "exponential"), "theta must be a numeric value for the exponential correlation function.")

    expect_error(correlation_function(D, 1:5, "matern"), "theta must be a numeric vector of length 2 for the matern correlation function.")
    expect_error(correlation_function(D, NA, "matern"), "theta must be a numeric vector of length 2 for the matern correlation function.")
    expect_error(correlation_function(D, "aaa", "matern"), "theta must be a numeric vector of length 2 for the matern correlation function.")
    expect_error(correlation_function(D, matrix(1:6, 3, 2), "matern"), "theta must be a numeric vector of length 2 for the matern correlation function.")
    
    D[1, ] <- - D[1, ]
    expect_error(correlation_function(D, 1, "exponential"), "D must be a distance matrix with only non-negative values.")
    expect_error(correlation_function(D, rep(1, 2), "matern"), "D must be a distance matrix with only non-negative values.")
    
})


test_that("pgdraw", {
    b <- 1:5
    set.seed(111)
    c <- rnorm(5)
    expect_equal(pgdraw(b, c), c(0.338935612254746, 1.03322724430168, 0.262630823432546, 1.40281087351243, 1.27332846908986))
    expect_error(pgdraw(1:5, 4), "b parameter must either be of length one, or the same length as the c parameter")
    expect_error(pgdraw(-5, 4), "b parameter must contain only positive integers")
    expect_error(pgdraw(5.5, 4), "b parameter must contain only positive integers")
    expect_error(pgdraw(5, 4, cores = 3.5), "cores must be a positive integer")
})

test_that("pgdraw.moments", {
    expect_equal(pgdraw.moments(2, 3), list(mu = 0.301716084548289, var = 0.0234847516762738))
    expect_error(pgdraw.moments(1:4, 4), "b must be a positive integer value")
    expect_error(pgdraw.moments(-5, 4), "b must be a positive integer value")
    expect_error(pgdraw.moments(5.5, 4), "b must be a positive integer value")
    expect_error(pgdraw.moments(NA, 4), "b must be a positive integer value")
    expect_error(pgdraw.moments("aaa", 4), "b must be a positive integer value")

    expect_error(pgdraw.moments(5, 1:4), "c must be a numeric value")  
    expect_error(pgdraw.moments(5, NA), "c must be a numeric value")
    expect_error(pgdraw.moments(5, "aaa"), "c must be a numeric value")    
})

