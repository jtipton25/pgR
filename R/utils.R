
##
## Define some helper functions
##

expit <- function(x) {
    1 / (1 + exp(-x))    
}

logit <- function(p) {
    log(p) - log(1 - p)
}

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

