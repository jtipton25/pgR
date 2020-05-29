#' Contruct the Dirichlet-multinomial parameter matrix \code{alpha}
#'
#' A function for calculating the Dirichlet-multinomial parameters from the underlying latent parameters
#'
#' @param Xbs is an \eqn{n \times q} matrix that are inputs into the functional form
#' @param beta is an \eqn{q} vector of inputs into the functional form
#' @param params is a list of parameters for each functional type and link function. See more details for details 

make_alpha <- function(Xbs, beta, params) {
    link  <- params$link
    tol   <- params$tol
    alpha <- matrix(0, nrow(Xbs), ncol(beta))
    if (link == "log") {
        alpha <- exp(Xbs %*% beta)
    } 
    if (link == "tobit") {
        alpha <- pmax(Xbs %*% beta, tol)
    }
    return(alpha)
}