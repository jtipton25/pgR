library(pgR)
library(tidyverse)
library(patchwork)
library(emdist)
## define a metric for similarity in distributions

layout(matrix(1:6, 3, 2, byrow = TRUE))


plot_samples <- function(b = 10, c = 0, n_rep = 10000, threshold = 1, n_curve = 500) {
    
    if (!pgR:::is_integer(b, 1))
        stop("b must be an integer")
    b <- as.integer(b)
    
    moments <- pgR::pgdraw.moments(b, c)
    curve_lower <- moments$mu - 5 * sqrt(moments$var)
    curve_upper <- moments$mu + 5 * sqrt(moments$var)
    time_approx <- system.time(samples_approx <- replicate(n_rep, pgR:::rcpp_pgdraw_approx(b, c, threshold = threshold)))[3]
    time_exact <- system.time(samples_exact <- replicate(n_rep, pgR:::rcpp_pgdraw(b, c)))[3]

    hist(samples_approx, freq = FALSE, main = paste("Approximation", format(time_approx, digits = 2, nsmall = 2)), breaks = 100)
    curve(dnorm(x, moments$mu, sqrt(moments$var)), from = curve_lower, to = curve_upper, add = TRUE, col = "red", n = n_curve)
    hist(samples_exact, freq = FALSE, main = paste("Exact", format(time_exact, digits = 2, nsmall = 2)), breaks = 100)
    curve(dnorm(x, moments$mu, sqrt(moments$var)), from = curve_lower, to = curve_upper, add = TRUE, col = "red", n = n_curve)
}

layout(matrix(1:6, 3, 2, byrow = TRUE))
plot_samples(5, - 80, n_curve = 1000)
plot_samples(5, - 0.8, n_curve = 1000)
plot_samples(5, - 0.08, n_curve = 1000)

layout(matrix(1:6, 3, 2, byrow = TRUE))
plot_samples(5, - 8, n_curve = 1000)
plot_samples(5, - 0.8, n_curve = 1000)
plot_samples(5, - 0.008, n_curve = 1000)

layout(matrix(1:6, 3, 2, byrow = TRUE))
plot_samples(50, - 8, n_curve = 1000)
plot_samples(50, - 0.8, n_curve = 1000)
plot_samples(50, - 0.008, n_curve = 1000)

b <- 100
c <- 7
moments=pgR::pgdraw.moments(b, c)
hist(replicate(10000, pgR:::rcpp_pgdraw_approx(b, c)), freq = FALSE, main = "Approximation", breaks = 100)
curve(dnorm(x, moments$mu, sqrt(moments$var)), from = -2, to = 20, add = TRUE, col = "red")
hist(replicate(10000, pgR:::rcpp_pgdraw(b, c)), freq = FALSE, main = "Exact", breaks = 100)
curve(dnorm(x, moments$mu, sqrt(moments$var)), from = -2, to = 20, add = TRUE, col = "red")

b <- 200
c <- -4
moments=pgR::pgdraw.moments(b, c)
hist(replicate(10000, pgR:::rcpp_pgdraw_approx(b, c, threshold = 10)), freq = FALSE, main = "Approximation", breaks = 100)
curve(dnorm(x, moments$mu, sqrt(moments$var)), from = -2, to = 20, add = TRUE)
hist(replicate(10000, pgR:::rcpp_pgdraw(b, c)), freq = FALSE, main = "Exact", breaks = 100)
curve(dnorm(x, moments$mu, sqrt(moments$var)), from = -2, to = 20, add = TRUE)


b <- 100
c = -40
microbenchmark::microbenchmark(
    pgR:::rcpp_pgdraw_approx(b, c),
    pgR:::rcpp_pgdraw(b, c)
)


ks.test(
    replicate(1000, pgR:::rcpp_pgdraw_approx(b, c)),
    replicate(1000, pgR:::rcpp_pgdraw(b, c))
)



b <- as.integer(seq(1, 50, by = 1))
c <- seq(-20, 20, by = 0.25)
grid <- as.matrix(expand.grid(b, c))
# mean_vec <- matrix(0, nrow(grid), 3)
# sd_vec <- matrix(0, nrow(grid), 3)
# timings_vec <- matrix(0, nrow(grid), 3)
n_rep <- 1000

# tst <- matrix(0, nrow(grid), n_rep)

compare_samplers <- function(grid_row, n_rep) {

    timings_exact <- system.time(samples_exact <- replicate(n_rep, pgR:::rcpp_pgdraw(b = as.integer(grid_row[1]), c = grid_row[2])))[3]
    timings_approx <- system.time(samples_approx <- replicate(n_rep, pgR:::rcpp_pgdraw_approx(grid_row[1], grid_row[2], threshold = 1)))[3]
    # moments <- pgR::pgdraw.moments(grid_row[1], grid_row[2])
    # timings_approx_rnorm <- system.time(samples_approx_rnorm <- rnorm(n_rep, moments$mu, sqrt(moments$var)))[3]
    ks_test <- ks.test(samples_exact, samples_approx)
    data.frame(
        b                    = grid_row[1],
        c                    = grid_row[2],
        mean_pg              = mean(samples_exact),
        mean_pg_approx       = mean(samples_approx),
        # mean_pg_approx_rnorm =  mean(samples_approx_rnorm),
        sd_pg                = sd(samples_exact),
        sd_pg_approx         = sd(samples_approx),
        # sd_pg_approx_rnorm   = sd(samples_approx_rnorm),
        time_pg              = timings_exact,
        time_pg_approx       = timings_approx,
        # time_pg_approx_rnorm = timings_approx_rnorm
        ks_stat              = ks_test$statistic,
        ks_pval              = ks_test$p.value
    )
}

if (file.exists(here::here("results", "pg_approx.RData"))) {
    load(file = here::here("results", "pg_approx.RData"))
} else {
    library(doParallel)
    library(foreach)
    registerDoParallel(cl = 6)
    dat_out <- foreach(
        i = 1:nrow(grid),
        .inorder = FALSE,
        .combine = rbind,
        .packages = "pgR") %dopar%
        compare_samplers(grid_row = as.vector(grid[i, ]), n_rep = n_rep)
    
    save(dat_out, file = here::here("results", "pg_approx.RData"))
}


zlims = range(c(c(dat_out$mean_pg), c(dat_out$mean_pg_approx)))
zlims_sd = range(c(c(dat_out$sd_pg), c(dat_out$sd_pg_approx)))
p1 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = mean_pg)) +
    geom_raster() +
    scale_fill_viridis_c(limits = zlims) +
    xlab("b") +
    ylab("c") +
    ggtitle("mean of exact samples")

p2 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = mean_pg_approx)) +
    geom_raster() +
    scale_fill_viridis_c(limits = zlims) +
    xlab("b") +
    ylab("c") +
    ggtitle("mean of approx samples")

p3 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = sd_pg)) +
    geom_raster() +
    scale_fill_viridis_c(limits = zlims_sd) +
    xlab("b") +
    ylab("c") +
    ggtitle("sd of exact samples")

p4 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = sd_pg_approx)) +
    geom_raster() +
    scale_fill_viridis_c(limits = zlims_sd) +
    xlab("b") +
    ylab("c") +
    ggtitle("sd of approx samples")

(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")



g1 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = abs(mean_pg - mean_pg_approx))) +
    geom_raster() +
    scale_fill_viridis_c() +
    xlab("b") +
    ylab("c") +
    ggtitle("deviation of approximations")

g2 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = abs(sd_pg - sd_pg_approx))) +
    geom_raster() +
    scale_fill_viridis_c() +
    xlab("b") +
    ylab("c") +
    ggtitle("deviation of sds approximations")

g1 + g2 + plot_layout(guides = "collect")

dat_out %>%
    ggplot(aes(x = b, y = c, fill = ks_stat)) +
    geom_raster() +
    scale_fill_viridis_c() +
    xlab("b") +
    ylab("c") +
    ggtitle("KS-test statistics")

dat_out %>%
    ggplot(aes(x = b, y = c, fill = ks_pval)) +
    geom_raster() +
    scale_fill_viridis_c(begin = 0.5, limits = c(0.05, 1), na.value = "blue") +
    xlab("b") +
    ylab("c") +
    ggtitle("KS-test p-values")

dat_out %>%
    subset(b > 15) %>%
    summarize(mean(ks_pval < 0.05))

##
## make sure the two approximate methods are consistent with each other
##

time_lims <- range(c(dat_out$time_pg, dat_out$time_pg_approx))
p1 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = time_pg)) +
    geom_raster() +
    scale_fill_viridis_c(limits = time_lims) +
    xlab("b") +
    ylab("c") +
    ggtitle("pg time")

p2 <- dat_out %>%
    ggplot(aes(x = b, y = c, fill = time_pg_approx)) +
    geom_raster() +
    scale_fill_viridis_c(limits = time_lims) +
    xlab("b") +
    ylab("c") +
    ggtitle("approx time")
p1 + p2 + plot_layout(guides = "collect")

dat_out %>%
    ggplot(aes(x = b, y = c, fill = time_pg / time_pg_approx)) +
    geom_raster() +
    colorspace::scale_fill_continuous_diverging(palette = "Blue-Red", mid = 1) +
    xlab("b") +
    ylab("c") +
    ggtitle("relative speed up")









