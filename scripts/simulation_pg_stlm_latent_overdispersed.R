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
tau2 <- runif(J-1, 0.5, 4.5)
theta <- matrix(c(rnorm(J-1, log(4.5), 0.2), rnorm(J-1, log(0.5), 0.2)), J-1, 2)
rho <- 0.9
sigma <- runif(J-1, 0.5, 4)
Sigma <- array(0, dim = c(J-1, N, N))
for (j in 1:(J-1)) {
    Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = "matern")
}
psi <- array(0, dim = c(N, J-1, n_time))
for (j in 1:(J-1)) {
    psi[, j, 1] <- t(rmvn(1, rep(0, N), Sigma[j, , ]))    
}

for (tt in 2:n_time) {
    for (j in 1:(J-1)) {
        psi[, j, tt] <- drop(rmvn(1, rho * psi[, j, tt-1], Sigma[j, , ]))
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
    for (j in 1:(J-1)) {
        eta[, j, tt] <- X %*% beta[, j] + psi[, j, tt] + rnorm(N, 0, sigma[j])
    }
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
params$n_adapt <- 500#1000
params$n_mcmc <- 500#1000
params$n_message <- 50
params$n_thin <- 1
priors <- default_priors_pg_stlm(Y, X, corr_fun = "matern")
inits  <- default_inits_pg_stlm(Y, X, priors, shared_covariance_params = FALSE, corr_fun = "matern")
config <- NULL

if (file.exists(here::here("results", "pg_stlm-latent-overdispersed.RData"))) {
    load(here::here("results", "pg_stlm-latent-overdispersed.RData"))
} else {
    
    # inits <- list(
    #     beta = beta,
    #     eta = eta,
    #     psi = psi,
    #     sigma2 = sigma^2,
    #     theta = theta,
    #     tau2 = tau2,
    #     rho = rho)
    # config <- list(
    #     sample_beta = FALSE,
    #     sample_eta = FALSE,
    #     sample_psi = FALSE,
    #     sample_sigma2 = FALSE,
    #     sample_theta = FALSE,
    #     sample_tau2 = FALSE,
    #     sample_rho = FALSE)
    start <- Sys.time()
    out <- pg_stlm_latent_overdispersed(Y, as.matrix(X), as.matrix(locs), params, 
                                        corr_fun = "matern",
                                        inits = inits, config = config,
                                        priors, n_cores = 1L, shared_covariance_params = FALSE)
    stop <- Sys.time()
    runtime <- stop - start
    
    save(out, runtime, file = here::here("results", "pg_stlm-latent-overdispersed.RData"))
}

runtime


# Plot posterior estimates for beta ----
layout(matrix(1:6, 3, 2))
for (i in 1:5) {
    matplot(out$beta[, , i], type = 'l', main = paste("species", i))
    abline(h = beta[, i], col = 1:nrow(beta))
}

plot(apply(out$beta, c(2, 3), mean), beta)
abline(0, 1, col = "red")

## plot eta estimates ----
eta_post_mean <- apply(out$eta, c(2, 3, 4), mean)
dimnames(eta_post_mean) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
eta_post_lower <- apply(out$eta, c(2, 3, 4), quantile, prob = 0.025)
dimnames(eta_post_lower) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
eta_post_upper <- apply(out$eta, c(2, 3, 4), quantile, prob = 0.975)
dimnames(eta_post_upper) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
dimnames(eta) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
dat_eta           <- as.data.frame.table(eta, responseName = "truth")
dat_eta_post_mean <- as.data.frame.table(eta_post_mean, responseName = "mean")
dat_eta_post_lower <- as.data.frame.table(eta_post_lower, responseName = "lower")
dat_eta_post_upper <- as.data.frame.table(eta_post_upper, responseName = "upper")

dat_plot <- Reduce(
    function(x, y) merge(x, y, all = TRUE),
    list(dat_eta, dat_eta_post_mean, dat_eta_post_lower, dat_eta_post_upper)
)

dat_plot %>% 
    ggplot(aes(x = truth, y = mean, color = species)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    facet_grid(time ~ species) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated eta")

## plot psi estimates ----
psi_post_mean <- apply(out$psi, c(2, 3, 4), mean)
dimnames(psi_post_mean) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
psi_post_lower <- apply(out$psi, c(2, 3, 4), quantile, prob = 0.025)
dimnames(psi_post_lower) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
psi_post_upper <- apply(out$psi, c(2, 3, 4), quantile, prob = 0.975)
dimnames(psi_post_upper) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
dimnames(psi) <- list(
    site = 1:N,
    species = 1:(J-1),
    time = 1:n_time
)
dat_psi           <- as.data.frame.table(psi, responseName = "truth")
dat_psi_post_mean <- as.data.frame.table(psi_post_mean, responseName = "mean")
dat_psi_post_lower <- as.data.frame.table(psi_post_lower, responseName = "lower")
dat_psi_post_upper <- as.data.frame.table(psi_post_upper, responseName = "upper")

dat_plot <- Reduce(
    function(x, y) merge(x, y, all = TRUE),
    list(dat_psi, dat_psi_post_mean, dat_psi_post_lower, dat_psi_post_upper)
)

dat_plot %>% 
    ggplot(aes(x = truth, y = mean, color = species)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    facet_grid(time ~ species) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated psi")




## plot pi estimates ----
pi_post_mean <- apply(out$pi, c(2, 3, 4), mean)
dimnames(pi_post_mean) <- list(
    site = 1:N,
    species = 1:J,
    time = 1:n_time
)

dimnames(pi) <- list(
    site = 1:N,
    species = 1:J,
    time = 1:n_time
)
dat_pi           <- as.data.frame.table(pi, responseName = "truth")
dat_pi_post_mean <- as.data.frame.table(pi_post_mean, responseName = "mean")

# dat_plot <- Reduce(
#   function(x, y) merge(x, y, all = TRUE), 
#   list(dat_eta, dat_eta_post_mean)
# )
dat_plot <- Reduce(
    function(x, y) merge(x, y, all = TRUE),
    list(dat_pi, dat_pi_post_mean)
)

dat_plot %>% 
    ggplot(aes(x = truth, y = mean, color = species)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    facet_grid(time ~ species) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated pi")



## plot beta estimates ----
dat_plot <- data.frame(
    beta = c(
        c(apply(out$beta, c(2, 3), mean)), 
        c(beta)
    ),
    type = rep(c("estimate", "truth"), each = (J-1) * ncol(X)),
    species = factor(rep(1:(J-1), each = ncol(X))),
    variable = 1:ncol(X)
)

dat_plot %>%
    pivot_wider(names_from = type, values_from = beta) %>%
    ggplot(aes(x = estimate, y = truth)) +
    # scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    facet_wrap(~ species, nrow = 8) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated beta")


## plot pi estimates ----
pi_post_mean <- apply(out$pi, c(2, 3, 4), mean)
dimnames(pi_post_mean) <- list(
    site = 1:N,
    species = 1:J,
    time = 1:n_time
)

dimnames(pi) <- list(
    site = 1:N,
    species = 1:J,
    time = 1:n_time
)
dat_pi           <- as.data.frame.table(pi, responseName = "truth")
dat_pi_post_mean <- as.data.frame.table(pi_post_mean, responseName = "mean")

# dat_plot <- Reduce(
#   function(x, y) merge(x, y, all = TRUE), 
#   list(dat_eta, dat_eta_post_mean)
# )
dat_plot <- Reduce(
    function(x, y) merge(x, y, all = TRUE),
    list(dat_pi, dat_pi_post_mean)
)

dat_plot %>% 
    ggplot(aes(x = truth, y = mean, color = species)) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    geom_point(alpha = 0.5) +
    facet_grid(time ~ species) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Estimated vs. simulated pi")

# plot sigma2 ----
data.frame(sigma_truth = c(sigma), 
           parameter = 1:(J-1)) %>% 
    left_join(
        data.frame(sigma = sqrt(c(out$sigma2)),
                   iteration = 1:nrow(out$sigma2), 
                   parameter = rep(1:ncol(out$sigma2), each = nrow(out$sigma2)))) %>%
    ggplot(aes(x = sigma)) +
    geom_histogram() +
    facet_wrap(~parameter) +
    geom_vline(aes(xintercept = sigma_truth, color = "red"))

# plot theta ----
thetas <- out$theta
dimnames(thetas) <- list(
    iteration = 1:dim(thetas)[1],
    species = 1:dim(thetas)[2],
    parameter = 1:dim(thetas)[3]
)

data.frame(theta_truth = c(theta),
           species = 1:(J-1), 
           parameter = rep(1:2, each = J-1)) %>%
    left_join(as.data.frame.table(thetas, responseName = "theta") %>%
                  mutate(iteration = as.integer(iteration), 
                         species = as.integer(species), 
                         parameter = as.integer(parameter))) %>%
    ggplot(aes(x = theta)) +
    geom_histogram() +
    facet_grid(species~parameter) +
    geom_vline(aes(xintercept = theta_truth, color = "red"))


# plot tau2 ----
data.frame(tau2_truth = c(tau2), 
           parameter = 1:(J-1)) %>% 
    left_join(
        data.frame(tau2 = c(out$tau2),
                   iteration = 1:nrow(out$tau2), 
                   parameter = rep(1:ncol(out$tau2), each = nrow(out$tau2)))) %>%
    ggplot(aes(x = tau2)) +
    geom_histogram() +
    facet_wrap(~parameter) +
    geom_vline(aes(xintercept = tau2_truth, color = "red"))

# plot rho ----
data.frame(rho_truth = rho,
           rho = c(out$rho),
           iteration = 1:length(out$rho)) %>%
    ggplot(aes(x = rho)) +
    geom_histogram() +
    geom_vline(aes(xintercept = rho_truth, color = "red"))





