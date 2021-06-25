# Simulate some space-time data ------------------------------------------------

library(pgR)
library(mvnfast)
library(tidyverse)
library(patchwork)
# library for model evaluation and comparison
library(loo)

set.seed(11)
N <- 20^2
J <- 6
p <- 2
n_time <- 4
# setup the spatial process
locs <- expand.grid(
    seq(0, 1, length = sqrt(N)),
    seq(0, 1, length = sqrt(N))
)
D <- fields::rdist(locs)
# use the same parameters for each of the J-1 components
tau2 <- 4.5
theta <- log(0.5) 
rho <- 0.9
Sigma <- tau2 * correlation_function(D, theta, corr_fun = "exponential") 
psi <- array(0, dim = c(N, J-1, n_time))
psi[, , 1] <- t(rmvn(J-1, rep(0, N), Sigma))
for (tt in 2:n_time) {
    for (j in 1:(J-1)) {
        psi[, j, tt] <- drop(rmvn(1, psi[, j, tt-1], Sigma))
    }
}

## setup the fixed effects process
X <- cbind(1, locs[, 1], locs[, 2])
beta <- matrix(0, ncol(X), J-1)
beta[-1, ] <- matrix(c(4, -4, -3, 4, 5, -2, -3, 2, 1, -2), 2, 5)
    # matrix(rnorm((J-1) * (ncol(X)-1), 0, 2.25), ncol(X)-1, (J-1))
beta[1, ] <- rnorm(J-1, 0, 0.5)
## make the intercepts smaller to reduce stochastic ordering effect
beta[1, ] <- beta[1, ] - seq(from = 2, to = 0, length.out = J-1) - 2

# make the spatial covariate more informative
# beta[3, ] <- 3 * beta[3, ]

eta <- array(0, dim = c(N, J-1, n_time))
pi <- array(0, dim = c(N, J, n_time))
for (tt in 1:n_time) {
    eta[, , tt] <- X %*% beta + psi[, , tt]
    pi[, , tt] <- eta_to_pi(eta[, , tt])
}


Y <- array(0, dim = c(N, J, n_time))
Y_prop <- array(0, dim = c(N, J, n_time))
for (tt in 1:n_time) {
    for (i in 1:N) {
        Y[i, , tt] <- rmultinom(1, 500, pi[i, , tt])
    }
    Y_prop[, , tt] <- counts_to_proportions(Y[, , tt])
}




## view the marginal effect of the covariates

dimnames(eta) <- list(sites = 1:dim(eta)[1], species = 1:dim(eta)[2], time = 1:dim(eta)[3])
dat_eta <- as.data.frame.table(eta, responseName = "eta")
Xbeta <- X %*% beta
dimnames(Xbeta) <- list(sites = 1:dim(Xbeta)[1], species = 1:dim(Xbeta)[2])
dat_Xbeta <- as.data.frame.table(Xbeta, responseName = "Xbeta")
dimnames(X) <- list(sites = 1:dim(X)[1], variable = 1:dim(X)[2])
dat_X <- as.data.frame.table(X, responseName = "X")
p_effects <- dat_eta %>%
    left_join(dat_Xbeta) %>%
    left_join(dat_X) %>%
    filter(time == 1) %>%
    filter(variable %in% c(2, 3)) %>%
    ggplot(aes(x = X, y = eta)) +
    geom_point() +
    facet_grid(variable ~ species) +
    stat_smooth(method = "lm")



## Plot the simulated data (at time 1)

## put data into data.frame for plotting
dat <- data.frame(
    Y       = c(Y_prop[, , 1]),
    lon     = locs[, 1],
    lat     = locs[, 2],
    species = factor(rep(1:J, each = N)), 
    pi      = c(pi[, , 1])
)
p_simulated <- ggplot(data = dat, aes(x = lon, y = lat, fill = pi)) +
    geom_raster() +
    scale_fill_viridis_c() +
    facet_wrap(~ species) +
    ggtitle("Simulated data") +
    theme(legend.position = "none")

dat_psi <- data.frame(
    lon     = locs[, 1],
    lat     = locs[, 2],
    species = factor(rep(1:(J-1), each = N)), 
    psi      = c(psi[, , 1])
)
dat_Xbeta <- data.frame(
    lon     = locs[, 1],
    lat     = locs[, 2],
    species = factor(rep(1:(J-1), each = N)), 
    Xbeta      = c(X %*% beta)
)

zlims <- range(range(psi[, , 1]), range(X %*% beta))
p_psi <- ggplot(data = dat_psi, aes(x = lon, y = lat, fill = psi)) +
    geom_raster() +
    scale_fill_viridis_c(limits = zlims) +
    facet_wrap(~ species) +
    ggtitle("Simulated spatio-temporal process") +
    theme(legend.position = "none")
p_Xbeta <- ggplot(data = dat_Xbeta, aes(x = lon, y = lat, fill =Xbeta)) +
    geom_raster() +
    scale_fill_viridis_c(limits = zlims) +
    facet_wrap(~ species) +
    ggtitle("Simulated fixed effects") +
    theme(legend.position = "none")
(p_effects + p_simulated) / (p_psi + p_Xbeta)


dimnames(pi) <- list(
    site      = 1:dim(pi)[1],
    species   = 1:dim(pi)[2],
    time      = 1:dim(pi)[3])

dat_locs <- data.frame(
    lon     = locs[, 1],
    lat     = locs[, 2],
    site    = 1:N )

pi_map <- as.data.frame.table(pi, responseName = "pi_truth") %>%
    mutate(site = as.numeric(site)) %>%
    left_join(dat_locs) %>%
    ggplot(aes(x = lon, y = lat, fill = pi_truth)) +
    geom_raster() +
    facet_grid(time ~ species) +
    ggtitle("Simulated spatio-temporal process") +
    theme(legend.position = "none")

pi_map

# fit model with covariates ----------------------------------------------------
params <- default_params()
params$n_adapt <- 1000
params$n_mcmc <- 1000
params$n_message <- 50
params$n_thin <- 1
priors <- default_priors_pg_stlm(Y, X)
inits  <- default_inits_pg_stlm(Y, X, priors, shared_covariance_params = TRUE)

if (file.exists(here::here("results", "pg_stlm-with-covariates.RData"))) {
    load(here::here("results", "pg_stlm-with-covariates.RData"))
} else {
    start <- Sys.time()
    out_X <- pg_stlm(Y, as.matrix(X), as.matrix(locs), params, priors, n_cores = 32L)
    stop <- Sys.time()
    runtime <- stop - start
    
    save(out_X, runtime, file = here::here("results", "pg_stlm-with-covariates.RData"))
}



dimnames(eta) <- list(
    site      = 1:dim(eta)[1],
    species   = 1:dim(eta)[2],
    time      = 1:dim(eta)[3])

eta_post <- out_X$eta
dimnames(eta_post) <- list(
    iteration = 1:dim(eta_post)[1],
    site      = 1:dim(eta_post)[2],
    species   = 1:dim(eta_post)[3],
    time      = 1:dim(eta_post)[4])


dat_eta <- as.data.frame.table(eta_post, responseName = "eta") %>%
    mutate(iteration = as.numeric(iteration))
# dat_eta_truth <- 

dat_plot <- dat_eta %>%
    group_by(site, species, time) %>%
    summarize(
        eta_mean = mean(eta), 
        upper_95_eta = quantile(eta, prob = 0.975),
        lower_95_eta = quantile(eta, prob = 0.025)) %>%
    left_join(as.data.frame.table(eta, responseName = "eta_truth"))
p_plot_covariates <- ggplot(dat_plot, aes(x = eta_truth, y = eta_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_95_eta, ymax = upper_95_eta)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated eta") +
    xlab("Simulated eta") +
    ylab("Estimated eta") +
    facet_grid(time ~ species) +
    theme_bw(base_size = 12)


p_plot_covariates

dimnames(pi) <- list(
    site      = 1:dim(pi)[1],
    species   = 1:dim(pi)[2],
    time      = 1:dim(pi)[3])

pi_post <- out_X$pi
dimnames(pi_post) <- list(
    iteration = 1:dim(pi_post)[1],
    site      = 1:dim(pi_post)[2],
    species   = 1:dim(pi_post)[3],
    time      = 1:dim(pi_post)[4])


dat_pi <- as.data.frame.table(pi_post, responseName = "pi") %>%
    mutate(iteration = as.numeric(iteration))
# dat_pi_truth <- 
    
dat_plot <- dat_pi %>%
    group_by(site, species, time) %>%
    summarize(
        pi_mean = mean(pi), 
        upper_95_pi = quantile(pi, prob = 0.975),
        lower_95_pi = quantile(pi, prob = 0.025)) %>%
    left_join(as.data.frame.table(pi, responseName = "pi_truth"))
p_pi_covariates <- ggplot(dat_plot, aes(x = pi_truth, y = pi_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_95_pi, ymax = upper_95_pi)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated pi") +
    xlab("Simulated pi") +
    ylab("Estimated pi") +
    facet_grid(time ~ species) +
    theme_bw(base_size = 12)


p_plot_covariates + p_pi_covariates

dat_locs <- data.frame(
    lon     = locs[, 1],
    lat     = locs[, 2],
    site    = 1:N 
)
pi_map_covariates <- dat_pi %>%
    group_by(site, species, time) %>%
    summarize(pi_mean = mean(pi)) %>%
    mutate(site = as.numeric(site)) %>%
    left_join(dat_locs) %>%
    ggplot(aes(x = lon, y = lat, fill = pi_mean)) +
    geom_raster() +
    facet_grid(time ~ species) +
    ggtitle("Fit with covariates") +
    theme(legend.position = "none")

pi_map_covariates


# fit model without covariates -------------------------------------------------

priors <- default_priors_pg_stlm(Y, as.matrix(rep(1, N)))
inits  <- default_inits_pg_stlm(Y, as.matrix(rep(1, N)), priors, shared_covariance_params = TRUE)

if (file.exists(here::here("results", "pg_stlm-no-covariates.RData"))) {
    load(here::here("results", "pg_stlm-no-covariates.RData"))
} else {
    start <- Sys.time()
    out_no_X <- pg_stlm(Y, as.matrix(rep(1, N)), as.matrix(locs), params, priors, n_cores = 32L)
    stop <- Sys.time()
    runtime <- stop - start
    
    save(out_no_X, runtime, file = here::here("results", "pg_stlm-no-covariates.RData"))
}

eta_post <- out_no_X$eta
dimnames(eta_post) <- list(
    iteration = 1:dim(eta_post)[1],
    site      = 1:dim(eta_post)[2],
    species   = 1:dim(eta_post)[3],
    time      = 1:dim(eta_post)[4])


dat_eta <- as.data.frame.table(eta_post, responseName = "eta") %>%
    mutate(iteration = as.numeric(iteration))
# dat_eta_truth <- 
    
    dat_plot <- dat_eta %>%
    group_by(site, species, time) %>%
    summarize(
        eta_mean = mean(eta), 
        upper_95_eta = quantile(eta, prob = 0.975),
        lower_95_eta = quantile(eta, prob = 0.025)) %>%
    left_join(as.data.frame.table(eta, responseName = "eta_truth"))
p_plot_no_covariates <- ggplot(dat_plot, aes(x = eta_truth, y = eta_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_95_eta, ymax = upper_95_eta)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated eta") +
    xlab("Simulated eta") +
    ylab("Estimated eta") +
    facet_grid(time ~ species) +
    theme_bw(base_size = 12)


dimnames(pi) <- list(
    site      = 1:dim(pi)[1],
    species   = 1:dim(pi)[2],
    time      = 1:dim(pi)[3])

pi_post <- out_no_X$pi
dimnames(pi_post) <- list(
    iteration = 1:dim(pi_post)[1],
    site      = 1:dim(pi_post)[2],
    species   = 1:dim(pi_post)[3],
    time      = 1:dim(pi_post)[4])


dat_pi <- as.data.frame.table(pi_post, responseName = "pi") %>%
    mutate(iteration = as.numeric(iteration))
# dat_pi_truth <- 

dat_plot <- dat_pi %>%
    group_by(site, species, time) %>%
    summarize(
        pi_mean = mean(pi), 
        upper_95_pi = quantile(pi, prob = 0.975),
        lower_95_pi = quantile(pi, prob = 0.025)) %>%
    left_join(as.data.frame.table(pi, responseName = "pi_truth"))
p_pi_no_covariates <- ggplot(dat_plot, aes(x = pi_truth, y = pi_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_95_pi, ymax = upper_95_pi)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated pi") +
    xlab("Simulated pi") +
    ylab("Estimated pi") +
    facet_grid(time ~ species) +
    theme_bw(base_size = 12)

p_plot_no_covariates + p_pi_no_covariates

(p_plot_covariates + p_plot_no_covariates) / (p_pi_covariates + p_pi_no_covariates)


dat_locs <- data.frame(
    lon     = locs[, 1],
    lat     = locs[, 2],
    site    = 1:N 
)
pi_map_no_covariates <- dat_pi %>%
    group_by(site, species, time) %>%
    summarize(pi_mean = mean(pi)) %>%
    mutate(site = as.numeric(site)) %>%
    left_join(dat_locs) %>%
    ggplot(aes(x = lon, y = lat, fill = pi_mean)) +
    geom_raster() +
    facet_grid(time ~ species) +
    ggtitle("Fit without covariates") +
    theme(legend.position = "none")

pi_map + pi_map_covariates + pi_map_no_covariates
# this code sends a notification to my phone 
# pushoverr::pushover(message = "Finished pgR model comparison fit")

# compare model fit using loo --------------------------------------------------

# Need to fiture out how to deal with missing observations here

# make a value NA to test how this works with calc_ll
Y[1, , 1] <- NA

# model fit using covariates
ll_X <- calc_ll_pg_stlm(Y, X, out_X)
# fit model without covariates
ll_no_X <- calc_ll_pg_stlm(Y, as.matrix(rep(1, N)), out_no_X)

# convert ll to a matrix for loo and WAIC (need to figure out missing values)
ll_mat_X <- t(apply(ll_X$ll, 1, c))
ll_mat_no_X <- t(apply(ll_no_X$ll, 1, c))

# drop the NA values

ll_mat_X <- ll_mat_X[, apply(ll_mat_X, 2, function(x) !any(is.na(x)))]
ll_mat_no_X <- ll_mat_no_X[, apply(ll_mat_no_X, 2, function(x) !any(is.na(x)))]

# WAIC using LaplacesDemon package
WAIC_LD_X <- LaplacesDemon::WAIC(ll_mat_X)
WAIC_LD_no_X <- LaplacesDemon::WAIC(ll_mat_no_X)

WAIC_LD_X$WAIC
WAIC_LD_no_X$WAIC


# I think the warnings from loo come from the autocorrelation
# I'm not really sure how to address this going forward but maybe the ideas
# presented here are a good start... https://mc-stan.org/loo/articles/loo2-non-factorizable.html

# WAIC using loo package
loo_waic_X <- waic(ll_mat_X)
loo_waic_no_X <- waic(ll_mat_no_X)
loo_compare(loo_waic_X, loo_waic_no_X)

r_eff_X <- relative_eff(ll_mat_X, chain_id = rep(1, nrow(ll_mat_X)))
r_eff_no_X <- relative_eff(ll_mat_no_X, chain_id = rep(1, nrow(ll_mat_no_X)))
loo_X <- loo(ll_mat_X, cores = 32L, r_eff = r_eff_X)
loo_no_X <- loo(ll_mat_no_X, cores = 32L, r_eff = r_eff_no_X)
comp <- loo_compare(loo_X, loo_no_X)

# the loo diagnostics don't look good -- maybe if we increase the sample size?
print(comp)
plot(loo_X)
plot(loo_no_X)





