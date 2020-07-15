#' Generate CAR precision matrix
#'
#' A function for setting up a conditional autoregressive (CAR) precision matrix for use as a prior in Bayesian penalized splines
#'
#' @param n is the length of the coefficient vector for the penalized bspline model  
#' @param phi is a number between -1 and 1 that defines the strength of the  autoregressive process. Typically this will be set to 1 for use as a prior in Bayesian penalized splines
#' @return an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" (see Matrix package for details)
#' @importFrom igraph as_adjacency_matrix make_lattice
#' @importFrom Matrix Diagonal colSums
#' 
#' @export

make_Q <- function(n, phi) {
    if (n <= 2)
        stop("n must be an integer larger or equal to 2")
    if (phi < -1 || phi > 1)
        stop("phi must be between -1 and 1")
    W <- as_adjacency_matrix(
        make_lattice(
            length = sqrt(n), 
            dim    = 2
        ),
        sparse = TRUE
    )
    D <- Diagonal(x = colSums(W))
    Q <- D - phi * W
    return(Q)
}
