
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


#' Calculate the intermediate Polya-gamma value Mi
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @export
#' 
calc_Mi <- function(Y) {
    if (!is.matrix(Y)) 
        stop("Y must be an integer matrix.")
    na_idx <- which(is.na(Y))
    if(length(na_idx) == 0) { 
        if (!is_integer(Y, length(Y)))
            stop("Y must be an integer matrix.")
    } else {
        if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
            stop("Y must be an integer matrix")
    }    
    if (any(rowSums(Y) == 0, na.rm = TRUE))
        stop ("There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    
    N <- nrow(Y)
    J <- ncol(Y)
    Mi <- matrix(0, N, J-1)
    for(i in 1:N){
        if (any(is.na(Y[i, ]))) {
            Mi[i, ] <- 0
        } else {
            if (J == 2) {
                Mi[i, ] <- sum(Y[i, ])
            } else {
                Mi[i, ] <- sum(Y[i, ]) - c(0, cumsum(Y[i, ][1:(J-2)]))
            }
        }
    }
    return(Mi)
}


#' Calculate the intermediate Polya-gamma value kappa
#'
#' @param Y is a \eqn{n \times J}{n x J} matrix of compositional count data.
#' @param Mi is a \eqn{n \times J-1}{n x J-1} matrix of transformed observations for Polya-gamma that is the output of the `calc_Mi()` function
#' @export
#' 
calc_kappa <- function(Y, Mi) {
    if (!is.matrix(Y)) 
        stop("Y must be an integer matrix.")
    na_idx <- which(is.na(Y))
    if(length(na_idx) == 0) { 
        if (!is_integer(Y, length(Y)))
            stop("Y must be an integer matrix.")
    } else {
        if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
            stop("Y must be an integer matrix")
    }    
    if (any(rowSums(Y) == 0, na.rm = TRUE))
        stop ("There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")

    # if (!is_numeric_matrix(Mi, nrow(Mi), ncol(Mi))) 
    #     stop("Mi must be a numeric matrix with number of rows and columns equal to Y")
    if ((nrow(Y) != nrow(Mi)) || (ncol(Y) -1 != ncol(Mi)))
        stop("Mi must be a numeric matrix with number of rows equal to the number of rows of Y and number of columns equal to one less than the number of columns of Y")    
    N <- nrow(Y)
    J <- ncol(Y)
    kappa <- matrix(0, N, J-1)
    for (i in 1:N) {
        if (any(is.na(Y[i, ]))) {
            kappa[i, ] <- 0
        } else {
            kappa[i, ] <- Y[i, 1:(J-1)] - Mi[i, ] / 2
        }
    }
    return(kappa)
}



