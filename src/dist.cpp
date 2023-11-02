//[[Rcpp::depends(pg)]]
#include <pg.h>
#include <Rcpp.h>
using namespace Rcpp;

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