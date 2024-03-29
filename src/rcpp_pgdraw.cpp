/*
 * This file implements the Polya-gamma sampler PG(1,z).
 * This is a C++ implementation of Algorithm 6 in PhD thesis of Jesse 
 * Bennett Windle, 2013
 * URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
 *
 * References:
 *
 *   Jesse Bennett Windle
 *   Forecasting High-Dimensional, Time-Varying Variance-Covariance Matrices
 *   with High-Frequency Data and Sampling Polya-Gamma Random Variates for
 *   Posterior Distributions Derived from Logistic Likelihoods  
 *   PhD Thesis, 2013   
 *
 *   Damien, P. & Walker, S. G. Sampling Truncated Normal, Beta, and Gamma Densities 
 *   Journal of Computational and Graphical Statistics, 2001, 10, 206-215
 *
 *   Chung, Y.: Simulation of truncated gamma variables 
 *   Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
 *
 * (c) Copyright Enes Makalic and Daniel F Schmidt, 2018
 */

//#include <omp.h>

#include <cmath>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// FCN prototypes
double samplepg(double);
double exprnd(double);
double tinvgauss(double, double);
double truncgamma();
double randinvg(double);
double aterm(int, double, double);

//' A function for evaluating drawing random Polya-gamma random variables in parallel
//'
//' @param b An \eqn{N}{N} \code{vector} of Polya-gamma parameters
//' @param c An \eqn{N}{N} \code{vector} of Polya-gamma parameters
//' @param cores An integer that gives the number of cores for openMP parallelization
//' @param threshold An integer that gives the number b at which a normal approximation (central limit theorm) is used. Default is 170.
//'   
//' @export
//' @name rcpp_pgdraw
//'
//' @keywords internal

//[[Rcpp::export]]
NumericVector rcpp_pgdraw(NumericVector b, NumericVector c, int cores = 1, int threshold = 170)
{
    int i;    // the loop index
    int const m = b.size();
    int const n = c.size();
    NumericVector y(n);
    
    if (m > 1 && n > m) {
        stop("Inadmissible input: when b is a vector, it must be as long as c.");
    }
    
#ifdef _OPENMP
    omp_set_num_threads(cores);
    
    // int n_per_thread = n / cores;
    
#pragma omp parallel for shared(b, c, y) private(i) schedule(dynamic)
#endif
    // for (int i = 0; i < n; i++)
    for (i = 0; i < n; i++)
    {
        
        double bi = (double) ((m > 1) ? b[i] : b[0]);
        double ci = (double) ((m > 1) ? c[i] : c[0]);
        
        double pg_sample = 0.0;
        
        if (bi < threshold) {
            for (int j = 0; j < bi; j++)
            {
                pg_sample += samplepg(ci);
            }
            y[i] = pg_sample;
        } else {
            double z = fabs(0.5 * ci);
            double E_y;
            double sigma2_y;
            if (z > 1e-12) {
                E_y = 0.25 * (bi * tanh(z) / z);
                sigma2_y = 0.0625 * ((bi + 1.0) * bi * pow(tanh(z) / z, 2) + bi * ((tanh(z) - z) / pow(z, 3))) - pow(E_y, 2);
            } else {
                E_y = 0.25 * (bi * (1.0 - 1.0/3.0) * pow(z, 2) + 2.0 / 15.0 * pow(z, 4) - 17.0 / 315.0 * pow(z, 6));
                sigma2_y = 0.0625 * ((bi + 1.0) * bi * pow(1.0 - 1.0 / 3.0 * pow(z, 2) + 2.0 / 15.0 * pow(z, 4) - 17.0 / 315.0 * pow(z, 6), 2) + 
                    bi * ((-1.0 / 3.0) + 2.0 / 15.0 * pow(z, 2) - 17.0 / 315.0 * pow(z, 4))) - pow(E_y, 2);                }
            y[i] = R::rnorm(E_y, sqrt(sigma2_y));
        }
    }
    
    return y;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
    //  PG(b, z) = 0.25 * J*(b, z/2)
    z = (double)std::fabs((double)z) * 0.5;
    
    // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
    double t = MATH_2_PI;
    
    // Compute p, q and the ratio q / (q + p)
    // (derived from scratch; derivation is not in the original paper)
    double K = z*z/2.0 + MATH_PI2/8.0;
    double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
    double logK = (double)std::log(K);
    double Kt = K * t;
    double w = (double)std::sqrt(MATH_PI_2);
    
    double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
    double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
    double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
    double ratio = 1.0 / (1.0 + p_over_q); 
    
    double u, X;
    
    // Main sampling loop; page 130 of the Windle PhD thesis
    while(1) 
    {
        // Step 1: Sample X ? g(x|z)
        u = R::runif(0.0,1.0);
        if(u < ratio) {
            // truncated exponential
            X = t + exprnd(1.0)/K;
        }
        else {
            // truncated Inverse Gaussian
            X = tinvgauss(z, t);
        }
        
        if (!std::isfinite(X)) continue;
        
        // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
        int i = 1;
        double Sn = aterm(0, X, t);
        double U = R::runif(0.0,1.0) * Sn;
        int asgn = -1;
        bool even = false;
        
        while(1) 
        {
            Sn = Sn + asgn * aterm(i, X, t);
            
            // Accept if n is odd
            if(!even && (U <= Sn)) {
                X = X * 0.25;
                return X;
            }
            
            // Return to step 1 if n is even
            if(even && (U > Sn)) {
                break;
            }
            
            even = !even;
            asgn = -asgn;
            i++;
        }
    }
    return X;
}

// Generate exponential distribution random variates
double exprnd(double mu)
{
    return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
    double f = 0;
    if(x <= t) {
        f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
    }
    else {
        f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
    }    
    return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
    // sampling
    double u = R::rnorm(0.0,1.0);
    double V = u*u;
    double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
    
    if(R::runif(0.0,1.0) > mu /(mu+out)) {    
        out = mu*mu / out; 
    }    
    return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables 
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
    double c = MATH_PI_2;
    double X, gX;
    
    bool done = false;
    while(!done)
    {
        X = exprnd(1.0) * 2.0 + c;
        gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
        
        if(R::runif(0.0,1.0) <= gX) {
            done = true;
        }
    }
    
    return X;  
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
    double X, u;
    double mu = 1.0/z;
    
    // Pick sampler
    if(mu > t) {
        // Sampler based on truncated gamma 
        // Algorithm 3 in the Windle (2013) PhD thesis, page 128
        while(1) {
            u = R::runif(0.0, 1.0);
            X = 1.0 / truncgamma();
            
            if ((double)std::log(u) < (-z*z*0.5*X)) {
                break;
            }
        }
    }  
    else {
        // Rejection sampler
        X = t + 1.0;
        while(X >= t) {
            X = randinvg(mu);
        }
    }    
    return X;
}



// Approximate algorithm described in Glynn et. al. https://g-lynn.github.io/files/GlynnTokdarHowardBanks_2019.pdf and code https://github.com/G-Lynn/DLTM/blob/master/Cpp/rpgApprox.cpp
//' A function for drawing approximate random Polya-gamma random variables in parallel
//'
//' @param b An \eqn{N}{N} \code{vector} of Polya-gamma parameters
//' @param c An \eqn{N}{N} \code{vector} of Polya-gamma parameters
//' @param cores An integer that gives the number of cores for openMP parallelization
//'   
//' @export
//' @name rcpp_pgdraw_approx
//'
//' @keywords internal

//[[Rcpp::export]]
NumericVector rcpp_pgdraw_approx(NumericVector b, NumericVector c, int cores = 1, int threshold = 170)
{
    int const m = b.size();
    int const n = c.size();
    NumericVector y(n);
    
    if (m > 1 && n > m) {
        stop("Inadmissible input: when b is a vector, it must be as long as c.");
    }
    
#ifdef _OPENMP
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n; i++)
    {
        int bi = (int) ((m > 1) ? b[i] : b[0]);
        
        double pg_sample = 0.0;
        
        if (bi < threshold) {
            for (int j = 0; j < bi; j++)
            {
                pg_sample += samplepg(c[i]);
            }
            y[i] = pg_sample;
        } else {
            double z = fabs(0.5 * c[i]);
            double E_y = 1.0 / 4.0;
            double sigma2_y = 1.0 / 24.0;
            if (z > 1e-12) {
                E_y = 0.25 * (bi * tanh(z) / z);
                sigma2_y = 0.0625 * ((bi + 1.0) * bi * pow(tanh(z) / z, 2.0) + bi * ((tanh(z) - z) / pow(z, 3))) - pow(E_y, 2.0);
            }
            y[i] = R::rnorm(E_y, sqrt(sigma2_y));
        }
    }
    
    return y;
}