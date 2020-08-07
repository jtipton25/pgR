
##
## Define some helper functions
##

#' The inverse of the logit function
#'
#' @param x the logit-transformed probability
#' @export

expit <- function(x) {
    if (!is_numeric(x, length(x)))
        stop("x must be a numeric value.")
    
    1 / (1 + exp(-x))    
}

#' The inverse of the logit function
#'
#' @param p is the probability
#' @export 

logit <- function(p) {
    if (any(is.na(p)))
        stop("p must be a numeric value between 0 and 1.")
    if (!is_numeric(p, length(p)))
        stop("p must be a numeric value between 0 and 1.")
    if (any(p <= 0) | any(p >= 1))
        stop("p must be a numeric value between 0 and 1.")
    log(p) - log(1 - p)
}

#' Convert from eta to pi
#'
#' @param eta The latent intensity vector to be transformed into a probability
#' @export

eta_to_pi <- function(eta) {
    ## convert eta to a probability vector pi
    ## can make this more general by first checking if vector vs. matrix and then
    ## calculating the response

    if (!is_numeric(eta, length(eta)))
        stop("eta must be either a numeric vector or a numeric matrix.")
    if (!(is.matrix(eta) | is.vector(eta)))
        stop("eta must be either a numeric vector or a numeric matrix.")
    
    pi <- NULL
    if (is.vector(eta)) {
        N <- length(eta) 
        pi <- matrix(0, N, 2)
        pi[, 1] <- expit(eta)
        pi[, 2] <- 1.0 - pi[, 1]
    } else {
        N <- nrow(eta)
        J <- ncol(eta) + 1
        pi <- matrix(0, N, J)
        stick <- rep(1, N)
        for (j in 1:(J - 1)) {
            pi[, j] <- expit(eta[, j]) * stick
            stick <- stick - pi[, j]
        }
        pi[, J] <- stick
    }
    return(pi)
}

