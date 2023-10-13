#include "poEMirtbase_gibbs.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List poEMirtbase_gibbs_fit(const arma::cube &Y,
                           const arma::mat &N,
                           arma::mat alpha,
                           arma::mat beta,
                           arma::vec theta,
                           const std::vector<arma::vec> &unique_categories,
                           const arma::mat &a0,
                           const arma::mat &A0,
                           const arma::mat &b0,
                           const arma::mat &B0,
                           const bool& PG_approx,
                           const int &constraint,
                           const bool &std,
                           const int &iter,
                           const int &warmup,
                           const int &thin,
                           const bool &save_item_parameters,
                           const int &verbose)
{
  // Instance
  poEMirtbase_gibbs Model(Y,
                          N,
                          alpha,
                          beta,
                          theta,
                          unique_categories,
                          a0,
                          A0,
                          b0,
                          B0,
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