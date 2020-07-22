#' Generate CAR precision matrix
#'
#' A function for setting up a conditional autoregressive (CAR) precision matrix for use as a prior in Bayesian penalized splines
#'
#' @param n_dims is a vector of length M (the number of MRA resolutions) for the precision matrix.
#' @param phi is a vector of length M where each element is a number between -1 and 1 that defines the strength of the  autoregressive process. Typically this will be set to 1 for use as a penalized prior in Bayesian models.
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (\code{use_spam = TRUE}) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" \code{use_spam = FALSE} (see spam and Matrix packages for details).
#' @param prec_model defines the form of precision matrix to be used. The precision matrix options are a conditional autoregressive (\code{prec_model = "CAR"}) or a simulatneous autogressive model (\code{prec_model = "SAR"}).
#' @return a list of  \eqn{n \times n}{n x n} sparse spam matrices or Matrix matrices of class "dgCMatrix" (see Matrix package for details)
#' @importFrom igraph as_adjacency_matrix make_lattice
#' @importFrom Matrix Diagonal colSums
#' @importFrom spam spam
#' 
#' @export

make_Q_alpha_2d <- function(n_dims, phi, use_spam = TRUE, prec_model = "CAR") {
    
    if (!(prec_model %in% c("CAR", "SAR")))
        stop('The only valid options for prec_model are "CAR" and "SAR"')
    ## n_dims is a vector of length M that contains the dimensions of each resolution of the process
    ## phi is a vector of length M that contains the CAR parameter
    if (length(n_dims) != length(phi))
        stop("n_dims and phi must both be vectors of length M")
    M <- length(n_dims)
    Q_alpha <- vector(mode = "list", length = M)
    for (m in 1:M) {
        W <- as_adjacency_matrix(make_lattice(length = n_dims[[m]], dim = 2), sparse = TRUE)
        D <- Diagonal(x = colSums(W))
        if (prec_model == "CAR") {
            Q_alpha[[m]] <- D - phi[m] * W
        } else if (prec_model == "SAR") {
            B <- diag(nrow(W)) - phi * W %*% Diagonal(x = 1 / colSums(W))
            Q_alpha[[m]] <- t(B) %*% B
        }
        if (use_spam) {
            ## use the spam package for sparse matrices
            # Q_alpha[[m]] <- spam(c(as.matrix(Q_alpha[[m]])), nrow = n_dims[[m]]^2, ncol = n_dims[[m]]^2)
            Q_alpha[[m]] <- as.spam.dgCMatrix(Q_alpha[[m]])
        }
    }
    return(Q_alpha)
}
