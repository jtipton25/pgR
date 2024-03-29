# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' A function for evaluating a multivariate normal density in parallel
#'
#' @param x An \eqn{N \times d}{N x d} \code{matrix} of multivariate Gaussian observations.
#' @param mean \code{mean} A \eqn{d} \code{vector} of the Gaussian distribution mean
#' @param sigma A \eqn{d \times d}{d x d} positive definite covariance \code{matrix} 
#' @param logd A logical value indicating whether the funtion should return the density on the log scale
#' @param cores An integer that gives the number of cores for openMP parallelization
#'   
#' @export
dmvnrm_arma_mc <- function(x, mean, sigma, logd = FALSE, cores = 1L) {
    .Call('_pgR_dmvnrm_arma_mc', PACKAGE = 'pgR', x, mean, sigma, logd, cores)
}

#' A function for evaluating drawing random Polya-gamma random variables in parallel
#'
#' @param b An \eqn{N}{N} \code{vector} of Polya-gamma parameters
#' @param c An \eqn{N}{N} \code{vector} of Polya-gamma parameters
#' @param cores An integer that gives the number of cores for openMP parallelization
#' @param threshold An integer that gives the number b at which a normal approximation (central limit theorm) is used. Default is 170.
#'   
#' @export
#' @name rcpp_pgdraw
#'
#' @keywords internal
NULL

#' A function for drawing approximate random Polya-gamma random variables in parallel
#'
#' @param b An \eqn{N}{N} \code{vector} of Polya-gamma parameters
#' @param c An \eqn{N}{N} \code{vector} of Polya-gamma parameters
#' @param cores An integer that gives the number of cores for openMP parallelization
#'   
#' @export
#' @name rcpp_pgdraw_approx
#'
#' @keywords internal
NULL

rcpp_pgdraw <- function(b, c, cores = 1L, threshold = 170L) {
    .Call('_pgR_rcpp_pgdraw', PACKAGE = 'pgR', b, c, cores, threshold)
}

rcpp_pgdraw_approx <- function(b, c, cores = 1L, threshold = 170L) {
    .Call('_pgR_rcpp_pgdraw_approx', PACKAGE = 'pgR', b, c, cores, threshold)
}

#' A function for sampling from conditional multivariate normal distributions with mean A^{-1}b and covariance matrix A^{-1}.
#'
#' @param A \code{A} A \eqn{d \times d} \code{matrix} for the Gaussian full conditional distribution precision matrix.
#' @param b \code{b} A \eqn{d} \code{vector} for the Gaussian full conditional distribution mean.
#' 
#' @export
rmvn_arma <- function(A, b) {
    .Call('_pgR_rmvn_arma', PACKAGE = 'pgR', A, b)
}

#' A function for sampling from conditional multivariate normal distributions with mean A^{-1}b and covariance matrix A^{-1}.
#'
#' @param a \code{a} A scalar for the Gaussian full conditional distribution precision.
#' @param b \code{b} A \eqn{d} \code{vector} for the Gaussian full conditional distribution mean.
#' 
#' @export
rmvn_arma_scalar <- function(a, b) {
    .Call('_pgR_rmvn_arma_scalar', PACKAGE = 'pgR', a, b)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_pgR_RcppExport_registerCCallable', PACKAGE = 'pgR')
})
