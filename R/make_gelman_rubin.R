#' Calculate the Gelman-Rubin statistic for custom MCMC output
#' 
#' A function to calculate Gelman-Rubin R-hat statistics for 
#' custom coded model MCMC output=
#' @param out \code{out} is a list object where each element of
#' the list is the results from a single MCMC chain. Within
#' each MCMC chain, there will be a list of objects of various
#' sizes where the first dimension of each of these objects 
#' is the number of MCMC iterations.

make_gelman_rubin <- function (out) {
  n_chains <- length(out)
  n_variables <- dim(out[[1]])[2]
  n_samples <- dim(out[[1]])[1]
  var_names <- colnames(out[[1]])
  Rhat <- vector("list", length = n_variables)
  
  Psi <- sapply(out, colMeans, simplify = TRUE)
  Psi_bar <- apply(Psi, 1, mean)
  Psi_diff <- (Psi - Psi_bar)^2
  ## Figure out how to make this faster
  W_diff <- array(0, dim = c(n_variables, n_chains, n_samples))
  for (i in 1:n_chains) {
    for (j in 1:n_samples) {
      W_diff[, i, j] <- (out[[i]][j, ] - Psi[, i])^2
    }
  }
  
  B <- n_samples / (n_chains - 1) * apply(Psi_diff, 1, sum)
  W <- 1 / n_chains * 1 / (n_samples - 1) * apply(W_diff, 1, sum)
  Rhat <- sqrt(((n_samples - 1) / n_samples * W + 
                  1 / n_samples * B) / W)
  names(Rhat) <- var_names
  return(Rhat)
}