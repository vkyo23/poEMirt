#include "poEMirtdynamic_gibbs.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List poEMirtdynamic_gibbs_fit(const arma::cube &Y,
                              const arma::mat &N,
                              arma::mat alpha,
                              arma::mat beta,
                              arma::mat theta,
                              const std::vector<arma::vec> &unique_categories,
                              const arma::mat &timemap,
                              const arma::mat &timemap2,
                              const arma::vec &item_timemap,
                              const arma::vec &item_match,
                              const arma::mat &a0,
                              const arma::mat &A0,
                              const arma::mat &b0,
                              const arma::mat &B0,
                              const arma::vec &m0,
                              const arma::vec &C0,
                              const arma::vec &Delta,
                              const bool &PG_approx,
                              const arma::vec &constraint,
                              const bool &std,
                              const int &iter,
                              const int &warmup,
                              const int &thin,
                              const bool &save_item_parameters,
                              const int &verbose) 
{
  // Instance
  poEMirtdynamic_gibbs Model(Y,
                             N,
                             alpha,
                             beta,
                             theta,
                             unique_categories,
                             timemap,
                             timemap2,
                             item_timemap,
                             item_match,
                             a0,
                             A0,
                             b0,
                             B0,
                             m0,
                             C0,
                             Delta,
                             PG_approx,
                             constraint,
                             std,
                             iter,
                             warmup,
                             thin,
                             save_item_parameters,
                             verbose);
  
  // Fitting
  Model.fit();
  
  // Result
  List out = Model.output();
  return out;
}