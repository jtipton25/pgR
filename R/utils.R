
##
## Define some helper functions
##

#' The inverse of the logit function
#'
#' @param x the logit-transformed probability
#' @export

expit <- function(x) {
    1 / (1 + exp(-x))    
}

#' The inverse of the logit function
#'
#' @param p is the probability
#' @export 

logit <- function(p) {
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
    N <- nrow(eta)
    J <- ncol(eta) + 1
    pi <- matrix(0, N, J)
    stick <- rep(1, N)
    for (j in 1:(J - 1)) {
        pi[, j] <- expit(eta[, j]) * stick
        stick <- stick - pi[, j]
    }
    pi[, J] <- stick
    return(pi)
}

