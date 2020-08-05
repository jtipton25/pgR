#' Check inputs
#' 
#' this function check that the input values are properly specified
#' 
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @keywords internal

check_input_pg_lm <- function(Y, X) {
    ## check the mcmc inputs
    if (!is.matrix(Y)) 
        stop("Y must be an integer matrix.")
    if (!is_integer_matrix(Y, nrow(Y), ncol(Y))) 
        stop("Y must be an integer matrix.")
    if (any(rowSums(Y) == 0))
        stop ("There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
        stop("X must be a numeric matrix.")
    if (nrow(Y) != nrow(X))
        stop("Y and X must have the same number of rows.")
}


#' Check inputs
#' 
#' this function check that the input values are properly specified
#' 
#' @param Y is a \eqn{n \times d}{n x d} matrix of compositional count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of spatial locations.
#' @keywords internal

check_input_pg_splm <- function(Y, X, locs) {
    ## check the mcmc inputs
    if (!is.matrix(Y)) 
        stop("Y must be an integer matrix.")
    if (!is_integer_matrix(Y, nrow(Y), ncol(Y)))
        stop("Y must be an integer matrix.")
    if (any(rowSums(Y) == 0))
        stop ("There must not be a row of counts that are all 0s. Please change any observations that have 0 total count to a vector of NAs.")
    if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
        stop("X must be a numeric matrix.")
    if (!is_numeric_matrix(locs, nrow(locs), 2)) 
        stop("locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    if (nrow(Y) != nrow(X))
        stop("Y and X must have the same number of rows.")
    if (nrow(Y) != nrow(locs))
        stop("Y and locs must have the same number of rows.")
}


#' Check inputs for pgSTLM
#' 
#' this function check that the input values are properly specified
#' 
#' @param Y is a \eqn{n \times d x T}{n x d X T} array of compositional count data through time.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of spatial locations.
#' @keywords internal

check_input_pgSTLM <- function(Y, X, locs) {
    stop("The function check_input_pgSTLM() is now deprecated. Please use check_input_pg_stlm().")
    # ## check the mcmc inputs
    # if (!is.array(Y)) 
    #     stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    # if (length(dim(Y)) != 3)
    #     stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    # na_idx <- which(is.na(Y))
    # if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
    #     stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    # if (any(apply(Y, c(2, 3), sum) == 0))
    #     stop ("There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs")
    # if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
    #     stop("X must be a numeric matrix.")
    # if (!is_numeric_matrix(locs, nrow(locs), ncol(locs)))
    #     stop("locs must be a numeric matrix.")
    # if (nrow(Y) != nrow(X))
    #     stop("Y and X must have the same number of rows.")
    # if (nrow(Y) != nrow(locs))
    #     stop("Y and locs must have the same number of rows.")
}


#' Check inputs for `pg_stlm()`
#' 
#' this function check that the input values of `pg_stlm()` are properly specified
#' 
#' @param Y is a \eqn{n \times d x T}{n x d X T} array of compositional count data through time.
#' @param X is a \eqn{n \times p}{n x p} matrix of climate variables.
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of spatial locations.
#' @keywords internal

check_input_pg_stlm <- function(Y, X, locs) {
    ## check the mcmc inputs
    if (!is.array(Y)) 
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    if (length(dim(Y)) != 3)
        stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    missing_idx <- matrix(FALSE, dim(Y)[1], dim(Y)[3])
    for (i in 1:dim(Y)[1]) {
        for (tt in 1:dim(Y)[3]) {
            if(!any(is.na(Y[i, , tt])))
                if(sum(Y[i, , tt]) == 0)
                    stop ("There must not be an observation vector that is all 0s. Please change any observations that have 0 total count to a vector of NAs.")
        }
    }
    na_idx <- which(is.na(Y))
    if(length(na_idx) == 0) { 
        if (!is_integer(Y, length(Y)))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    } else {
        if (!is_integer(Y[-na_idx], length(Y[-na_idx])))
            stop("Y must be a 3 dimensional array of integer values with rows representing the locations, columns representing the species, and the third dimension representing time.")
    }    
    if (!is_numeric_matrix(X, nrow(X), ncol(X))) 
        stop("X must be a numeric matrix.")
    if (!is_numeric_matrix(locs, nrow(locs), 2))
        stop("locs must be a numeric matrix with rows equal to the number of rows of Y and 2 columns.")
    if (nrow(Y) != nrow(X))
        stop("Y and X must have the same number of rows.")
    if (nrow(Y) != nrow(locs))
        stop("Y and locs must have the same number of rows.")
}
