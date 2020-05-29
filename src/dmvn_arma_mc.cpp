// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif


static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
    arma::uword const n = trimat.n_cols;
    
    for(unsigned j = n; j-- > 0;){
        double tmp(0.);
        for(unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

//' A function for evaluating a multivariate normal density in parallel
//'
//' @param x An \eqn{N \times d}{N x d} \code{matrix} of multivariate Gaussian observations.
//' @param mean \code{mean} A \eqn{d} \code{vector} of the Gaussian distribution mean
//' @param sigma A \eqn{d \times d}{d x d} positive definite covariance \code{matrix} 
//' @param logd A logical value indicating whether the funtion should return the density on the log scale
//' @param cores An integer that gives the number of cores for openMP parallelization
//'   
//' @export
//[[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = false,
                         int const cores = 1) {  
    using arma::uword;
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
    uword const n = x.n_rows, 
        xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
        constants = -(double)xdim/2.0 * log2pi, 
        other_terms = rootisum + constants;
    
    arma::rowvec z;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(z)
#endif
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);   
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
    
    if (logd)
        return out;
    return exp(out);
}
        
        
        
        