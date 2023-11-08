//[[Rcpp::depends(pg)]]
#include <pg.h>
#include <Rcpp.h>
using namespace Rcpp;

double sample_gamma(double alpha)
{
  int check;
  double out;
  double u, v, w, x, y, z, b, c;
  
  b = alpha - 1;
  c = 3 * alpha - 0.75;
  check = 0;
  while (check == 0) {
    u = R::runif(0.0, 1.0);
    v = R::runif(0.0, 1.0);
    w = u * (1 - u);
    y = std::sqrt(c / w) * (u - 0.5);
    x = b + y;
    if (x > 0) {
      z = 64 * std::pow(v, 2) * std::pow(w, 3);
      if (z <= (1 - (2 * std::pow (y, 2) / x))) {
        check = 1;
        out = x;
      } else if ((2 * (b * std::log(x / b) - y)) >= std::log(z)) {
        check = 1;
        out = x;
      } else {
        check = 0;
      }
    }
  }
  
  return out;
}

double rgamma(double shape, 
              double rate)
{
  double out;
  
  if (shape > 1) {
    out = sample_gamma(shape) / rate;
  } else if (shape == 1) {
    out = -std::log(R::runif(0.0, 1.0)) / rate;
  } else {
    out = sample_gamma(shape + 1) * std::pow(R::runif(0.0, 1.0), 1 / shape) / rate;
  }
  return out;
}

arma::vec rmultinom2(int &size, 
                     arma::vec &prob) 
{
  IntegerVector outcome(prob.size());
  rmultinom(size, prob.begin(), prob.size(), outcome.begin());
  return as<arma::vec>(wrap((outcome)));
}

// Approximation of PG draw from Glynn et al. (2019, BA).
// https://projecteuclid.org/journals/bayesian-analysis/volume-14/issue-1/Bayesian-Analysis-of-Dynamic-Linear-Topic-Models/10.1214/18-BA1100.full
double rpg_approx(const double &b,
                  const double &c) 
{
  double mean = 1 / (2 * c) * std::tanh(c / 2);
  double var = 1.0 / (4 * std::pow(c, 3.0)) * (std::sinh(c) - c) * std::pow((1 / std::cosh(c/2)), 2.0);
  return R::rnorm(b * mean, std::sqrt(b * var));
}

// Default rpg
double rpg(const double &b,
           const double &c)
{
  return pg::rpg_scalar_hybrid(b, c);
}