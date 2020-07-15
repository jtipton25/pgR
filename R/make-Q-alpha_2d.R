#' Generate CAR precision matrix
#'
#' A function for setting up a conditional autoregressive (CAR) precision matrix for use as a prior in Bayesian penalized splines
#'
#' @param n is the length of the coefficient vector for the penalized bspline model  
#' @param phi is a number between -1 and 1 that defines the strength of the  autoregressive process. Typically this will be set to 1 for use as a prior in Bayesian penalized splines
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE}(see Matrix package for details)
#' @return a list of  \eqn{n \times n}{n x n} sparse spam matrices or Matrix matrices of class "dgCMatrix" (see Matrix package for details)
#' @importFrom igraph as_adjacency_matrix make_lattice
#' @importFrom Matrix Diagonal colSums
#' @importFrom spam spam
#' 
#' @export

make_Q_alpha_2d <- function(n_dims, phi, use_spam = TRUE) {
    ## n_dims is a vector of length M that contains the dimensions of each resolution of the process
    ## phi is a vector of length M that contains the CAR parameter
    if (length(n_dims) != length(phi))
        stop("n_dims and phi must both be vectors of length M")
    M <- length(n_dims)
    Q_alpha <- vector(mode = "list", length = M)
    for (m in 1:M) {
        W <- as_adjacency_matrix(make_lattice(dimvector = n_dims[[m]]), sparse = TRUE)
        D <- Diagonal(x = colSums(W))
        Q_alpha[[m]] <- D - phi[m] * W
        if (use_spam) {
            ## use the spam package for sparse matrices
            Q_alpha[[m]] <- spam(c(as.matrix(Q_alpha[[m]])), nrow = n_dims[[m]], ncol = n_dims[[m]])
        }
    }
    return(Q_alpha)
}
