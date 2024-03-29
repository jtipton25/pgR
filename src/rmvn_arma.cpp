#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace arma;

//' A function for sampling from conditional multivariate normal distributions with mean A^{-1}b and covariance matrix A^{-1}.
//'
//' @param A \code{A} A \eqn{d \times d} \code{matrix} for the Gaussian full conditional distribution precision matrix.
//' @param b \code{b} A \eqn{d} \code{vector} for the Gaussian full conditional distribution mean.
//' 
//' @export
//[[Rcpp::export]]
arma::vec rmvn_arma(arma::mat& A, arma::vec& b){
    int ncols = A.n_cols;
    arma::mat A_chol(ncols, ncols);
    bool success = false;
    int counter = 0;
    while(success == false && counter < 100){
        success = arma::chol(A_chol, A);
        if(success == false){
            counter++;
            A += arma::mat(ncols, ncols, arma::fill::eye) * 1e-6;
        }
    }
    // this option uses R's random number seed
    arma::vec devs = rnorm(ncols);
    // arma::vec devs = rnorm(ncols);
    arma::vec temp = solve(trimatl(A_chol.t()), b);
    return arma::vec(solve(trimatu(A_chol), temp + devs));
}
