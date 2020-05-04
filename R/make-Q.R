#' Generate CAR precision matrix
#'
#' A function for setting up a conditional autoregressive (CAR) precision matrix for use as a prior in Bayesian penalized splines
#'
#' @param n is the length of the coefficient vector for the penalized bspline model  
#' @param phi is a number between -1 and 1 that defines the strength of the  autoregressive process. Typically this will be set to 1 for use as a prior in Bayesian penalized splines
#' @export
make_Q <- function(n, phi) {
    if (n <= 2)
        stop("n must be an integer larger or equal to 2")
    if (phi < -1 || phi > 1)
        stop("phi must be between -1 and 1")
    W <- toeplitz(c(0, 1, rep(0, n-2)))
    D <- rowSums(W)
    Q <- diag(D) - phi * W
    return(Q)
}
