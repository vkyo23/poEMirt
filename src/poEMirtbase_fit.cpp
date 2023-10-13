#include "poEMirtbase.h"

//[[Rcpp::export]]
List poEMirtbase_fit(const arma::cube &Y,
                     const arma::mat &N,
                     const arma::mat &alpha_init,
                     const arma::mat &beta_init,
                     const arma::vec &theta_init,
                     const std::vector<arma::vec> &unique_categories,
                     const arma::mat &a0,
                     const arma::mat &A0,
                     const arma::mat &b0,
                     const arma::mat &B0,
                     const int &constraint,
                     const bool &std,
                     const int &maxit,
                     const int &verbose,
                     const double &tol,
                     const bool &compute_ll) 
{
  // Instance
  poEMirtbase Model(Y,
                    N,
                    alpha_init,
                    beta_init,
                    theta_init,
                    unique_categories,
                    a0,
                    A0,
                    b0,
                    B0,
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