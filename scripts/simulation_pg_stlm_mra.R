# Simulate some space-time data ------------------------------------------------

library(pgR)
library(mvnfast)
library(tidyverse)
library(patchwork)
library(BayesMRA)
library(spam)
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
rho <- 0.9
sigma <- 0.7








## generate the MRA basis
M <- 3
MRA <- mra_wendland_2d(as.matrix(locs), M = M, n_coarse_grid = 8)
W <- MRA$W

## generate the MRA spatially-correlated random effects precision matrix
tau2 <- sapply(1:(J-1), function(j) 2^(1:M) * 1000)

Q_alpha <-  make_Q_alpha_2d(sqrt(MRA$n_dims), rep(0.999, length(MRA$n_dims)))
Q_alpha_tau2 <- vector(mode = "list", length = J-1)
for (j in 1:(J-1)) {
    Q_alpha_tau2[[j]] <- make_Q_alpha_tau2(Q_alpha, tau2[, j])
}

## define the sum-to-0 constraint for alpha
constraints <- make_constraint(MRA, constraint = "resolution", joint = TRUE)
A_constraint <- constraints$A_constraint
a_constraint <- constraints$a_constraint


## temporal autocorrelation parameter rho
rho <- runif(J-1, 0.8, 1)

## simulate the spatial random effects
alpha <- array(0, dim = c(sum(MRA$n_dims), J-1, n_time))
for (j in 1:(J-1)) {
    alpha[, j, 1] <- rmvnorm.prec.const(n = 1, mu = rep(0, sum(MRA$n_dims)), Q = Q_alpha_tau2[[j]], A = A_constraint, a = a_constraint)
}
for (tt in 2:n_time) {
    for (j in 1:(J-1)) {
        alpha[, j, tt] <- rmvnorm.prec.const(n = 1, mu = rho[j] * alpha[, j, tt-1], Q = Q_alpha_tau2[[j]])
    }  
}
W_alpha <- array(0, dim = c(N, J-1, n_time))
for (tt in 1:n_time) {
    W_alpha[, , tt] <- W %*% alpha[, , tt]
}

# check if first field is mean 0
apply(W_alpha, c(2, 3), sum)[, 1]

## setup the fixed effects process
X <- cbind(1, matrix(runif(N*p), N, p))
beta <- matrix(rnorm((J-1) * ncol(X), 0, 0.25), ncol(X), (J-1))
## make the intercepts smaller to reduce stochastic ordering effect
beta[1, ] <- beta[1, ] - seq(from = 2, to = 0, length.out = J-1)

## add in some residual error
sigma2 <- rgamma(J-1, 1, 5)

eta <- array(0, dim = c(N, J-1, n_time))
pi      <- array(0, dim = c(N, J, n_time), 
                 dimnames = list(site = 1:N, species = 1:J, time = 1:n_time))
Y_prop  <- array(0, dim = c(N, J, n_time), 
                 dimnames = list(site = 1:N, species = 1:J, time = 1:n_time))
Y      <- array(0, dim = c(N, J, n_time))
for (tt in 1:n_time) {
    eta[, , tt] <- X %*% beta + W_alpha[, , tt] + sapply(1:(J-1), function(j) rnorm(N, 0, sqrt(sigma2[j])))
    pi[, , tt] <- eta_to_pi(eta[, , tt])
    for (i in 1:N) {
        Y[i, , tt] <- rmultinom(1, 500, pi[i, , tt])
    }
    Y_prop[, , tt]     <- counts_to_proportions(Y[, , tt])
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
    psi      = c(W_alpha[, , 1])
)
dat_Xbeta <- data.frame(
    lon     = locs[, 1],
    lat     = locs[, 2],
    species = factor(rep(1:(J-1), each = N)), 
    Xbeta      = c(X %*% beta)
)

zlims <- range(range(W_alpha[, , 1]), range(X %*% beta))
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
inits  <- default_inits_pg_stlm(Y, X, priors, shared_covariance_params = FALSE)

if (file.exists(here::here("results", "pg_stlm-mra.RData"))) {
    load(here::here("results", "pg_stlm-mra.RData"))
} else {
    start <- Sys.time()
    out <- pg_stlm_mra(Y, as.matrix(X), as.matrix(locs), params, priors, n_cores = 1L, M = 3)
    stop <- Sys.time()
    runtime <- stop - start
    
    save(out, runtime, file = here::here("results", "pg_stlm-mra.RData"))
}

runtime

# Plot posterior estimates ----
layout(matrix(1:6, 3, 2))
for (i in 1:5) {
    matplot(out$beta[, , i], type = 'l', main = paste("species", i))
    abline(h = beta[, i], col = 1:nrow(beta))
}

plot(apply(out$beta, c(2, 3), mean), beta)
abline(0, 1, col = "red")

## plot eta estimates
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


## plot pi estimates
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



## plot beta estimates
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


## plot pi estimates
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


data.frame(sigma2 = c(out$sigma2),
           iteration = 1:nrow(out$sigma2), 
           parameter = rep(1:ncol(out$sigma2), each = nrow(out$sigma2))) %>%
    ggplot(aes(x = sqrt(sigma2))) +
    geom_histogram() +
    facet_wrap(~parameter) +
    geom_vline(aes(xintercept = sigma, color = "red"))





