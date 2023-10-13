#include <Rcpp.h>
using namespace Rcpp;

#define PI 3.141592653589793238462643383279502884197169399375105820974

// Inverse gaussian
double rinvgauss(const double & mu)
{
  double U = R::rnorm(0.0, 1.0);
  double V = std::pow(U, 2.0);
  double X = mu + 0.5 * mu * (mu * V - std::sqrt(4.0 * mu * V + std::pow(mu, 2.0) * std::pow(V, 2.0)));
  if (R::runif(0.0, 1.0) > mu/(mu + X)) {    
    X = std::pow(mu, 2.0) / X; 
  }    
  return X;
}

// Truncated gamma dist
double rtgamma()
{
  double out, val;
  while(1)
  {
    out = R::rexp(1.0) * 2.0 + PI/2;
    val = std::sqrt(PI/2) / std::sqrt(out);
    if(R::runif(0.0, 1.0) <= val) {
      break;
    }
  }
  return out;  
}

// Truncated inverse-gaussian (for PG draw)
// Algorithm 4 of Windle (2013)
double rtinvgauss(const double &z,
                  const double &t)
{
  double X, U;
  double mu = 1.0 / z;
  
  if(mu > t) {
    while (1) {
      U = R::runif(0.0, 1.0);
      X = 1.0 / rtgamma();
      if (std::log(U) < (-z*z*0.5*X)) {
        break;
      } 
    }
  } else {
    X = t + 1.0;
    while (X >= t) {
      X = rinvgauss(mu);
    }
  }     
  return X;
}

// Calculate a_n(x): piecewise coefficients
// Based on Windle (2013), equation (2.14) ~ (2.15)
double a_n(const int &n,
           const double &x,
           const double &t)
{
  double out;
  if(x <= t) {
    out = std::log(PI) + std::log(n + 0.5) + 1.5 * (std::log(2) - std::log(PI * x)) - 2 * std::pow(n + 0.5, 2.0) / x;
  }
  else {
    out = std::log(PI) + std::log(n + 0.5) - x * std::pow(PI, 2.0) * std::pow(n + 0.5, 2.0) / 2;
  }    
  return std::exp(out);
}

// Generate PG random variable
double rpg_raw(const double& c) 
{
  // Reference
  // Algorithm 6 of Windle (2013; PhD dissertation) 
  // -> https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
  // Schmidt's pgdraw package (R)
  // -> https://github.com/cran/pgdraw/tree/master
  
  // Aux vars
  double z = std::fabs(c) * 0.5;
  double t = 0.636619772367581343075535053490057448137838582961825794990;
  double K = std::pow(PI, 2.0) / 8.0 + std::pow(z, 2.0) / 2.0;
  double logA = std::log(4.0) - std::log(PI) - z;
  double logK = std::log(K);
  double Kt = K * t;
  double w = std::sqrt(PI/2);
  double logf1 = logA + R::pnorm(w * (t * z - 1), 0.0, 1.0, 1, 1) + logK + Kt;
  double logf2 = logA + 2 * z + R::pnorm(-w * (t * z + 1), 0.0, 1.0, 1, 1) + logK + Kt;
  double p_over_q = std::exp(logf1) + std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  // LOOP
  double U, X;
  while (1) 
  {
    Rcpp::checkUserInterrupt();
    // Generate U ~ U(0, 1)
    U = R::runif(0.0, 1.0);
    
    // Sample X
    if (U < ratio) {
      // Truncated exp dist.
      X = t + R::rexp(1.0)/K;
    } else {
      // Truncated Inv-Gauss dist.
      X = rtinvgauss(z, t);
    }
    
    // Accept or not
    int n = 0;
    double Sn = a_n(0, X, t);
    double U = R::runif(0.0, 1.0) * Sn;
    int a_sign = -1;
    bool even = false;
    while (1) {
      n += 1;
      Sn = Sn + a_sign * a_n(n, X, t);
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even; // even -> odd or odd -> even
      a_sign = -a_sign; // sign switching
    }
  }
  return X;
}

// rpg wrapper
double rpg(const double &b,
           const double &c)
{
  double out = 0;
  for (int i = 0; i < b; i++) {
    out += rpg_raw(c);
  }
  return out;
}

// Approximation of PG draw from Glynn et al. (2019, BA).
// https://projecteuclid.org/journals/bayesian-analysis/volume-14/issue-1/Bayesian-Analysis-of-Dynamic-Linear-Topic-Models/10.1214/18-BA1100.full
double rpg_approx(const double& b,
                  const double& c) 
{
  double mean = 1 / (2 * c) * std::tanh(c / 2);
  double var = 1.0 / (4 * std::pow(c, 3.0)) * (std::sinh(c) - c) * std::pow((1 / std::cosh(c/2)), 2.0);
  return R::rnorm(b * mean, std::sqrt(b * var));
}