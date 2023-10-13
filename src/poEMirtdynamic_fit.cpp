#include "poEMirtdynamic.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List poEMirtdynamic_fit(const arma::cube &Y,
                        const arma::mat &N,
                        const arma::mat &alpha_init,
                        const arma::mat &beta_init,
                        const arma::mat &theta_init,
                        const std::vector<arma::vec> &unique_categories,
                        const arma::mat &timemap,
                        const arma::vec &item_timemap,
                        const arma::vec &item_match,
                        const arma::mat &a0,
                        const arma::mat &A0,
                        const arma::mat &b0,
                        const arma::mat &B0,
                        const arma::vec &m0,
                        const arma::vec &C0,
                        const arma::vec &Delta,
                        const arma::vec &constraint,
                        const bool &std,
                        const int &maxit,
                        const int &verbose,
                        const double &tol,
                        const bool &compute_ll) 
{
  // Instance
  poEMirtdynamic Model(Y,
                       N,
                       alpha_init,
                       beta_init,
                       theta_init,
                       unique_categories,
                       timemap,
                       item_timemap,
                       item_match,
                       a0,
                       A0,
                       b0,
                       B0,
                       m0,
                       C0,
                       Delta,
                       constraint,
                       std,
                       maxit,
                       verbose,
                       tol,
                       compute_ll);
  
  // Fitting
  Model.fit();
  
  // Result
  List out = Model.output();
  return out;
}