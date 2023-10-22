#include <RcppArmadillo.h>
#include "helper.h"
using namespace Rcpp;
using namespace arma;

cube calc_predY(const mat &N,
                const mat &alpha, 
                const mat &beta, 
                const vec &theta,
                const std::vector<vec> &unique_categories) 
{
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols + 1;
  cube predY(I, J, K);
  predY.fill(NA_REAL);
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      if (N(i, j) != 0) {
        double stick = 1.0;
        vec unq = unique_categories[j] - 1;
        vec pi_(unq.size());
        for (int k = 0; k < (unq.size()-1); k++) {
          double psi = alpha(j, unq[k]) + beta(j, unq[k]) * theta[i];
          pi_[k] = (std::exp(psi) / (1 + std::exp(psi))) * stick;
          stick -= pi_[k];
        }
        pi_[unq.size() - 1] = stick;
        int N_ij = N(i, j);
        vec tmp = rmultinom2(N_ij, pi_);
        for (int k = 0; k < unq.size(); k++) {
          predY(i, j, unq[k]) = tmp[k];
        }
      }
    }
  }
  return predY;
}

cube dyn_calc_predY(const mat &N,
                    const mat &alpha, 
                    const mat &beta, 
                    const mat &theta,
                    const std::vector<vec> &unique_categories,
                    const vec &item_timemap) 
{
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols + 1;
  cube predY(I, J, K);
  predY.fill(NA_REAL);
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      if (N(i, j) != 0) {
        int t = item_timemap[j];
        if (theta(i, t) != 0) {
          double stick = 1.0;
          vec unq = unique_categories[j] - 1;
          vec pi_(unq.size());
          for (int k = 0; k < (unq.size()-1); k++) {
            double psi = alpha(j, unq[k]) + beta(j, unq[k]) * theta(i, t);
            pi_[k] = (std::exp(psi) / (1 + std::exp(psi))) * stick;
            stick -= pi_[k];
          }
          pi_[unq.size() - 1] = stick;
          int N_ij = N(i, j);
          vec tmp = rmultinom2(N_ij, pi_);
          for (int k = 0; k < unq.size(); k++) {
            predY(i, j, unq[k]) = tmp[k];
          }
        }
      }
    }
  }
  return predY;
}


cube calc_probY(const mat &alpha, 
                const mat &beta, 
                const vec &theta,
                const std::vector<vec> &unique_categories) 
{
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols + 1;
  cube probY(I, J, K);
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      double stick = 1.0;
      vec unq = unique_categories[j] - 1;
      vec pi_(unq.size());
      for (int k = 0; k < (unq.size()-1); k++) {
        double psi = alpha(j, unq[k]) + beta(j, unq[k]) * theta[i];
        pi_[k] = (std::exp(psi) / (1 + std::exp(psi))) * stick;
        stick -= pi_[k];
      }
      pi_[unq.size() - 1] = stick;
      for (int k = 0; k < unq.size(); k++) {
        probY(i, j, unq[k]) = pi_[k];
      }
    }
  }
  return probY;
}

cube dyn_calc_probY(const mat &alpha, 
                    const mat &beta, 
                    const mat &theta,
                    const std::vector<vec> &unique_categories,
                    const vec &item_timemap) 
{
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols + 1;
  cube probY(I, J, K);
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      int t = item_timemap[j];
      if (theta(i, t) != 0) {
        double stick = 1.0;
        vec unq = unique_categories[j] - 1;
        vec pi_(unq.size());
        for (int k = 0; k < (unq.size()-1); k++) {
          double psi = alpha(j, unq[k]) + beta(j, unq[k]) * theta(i, t);
          pi_[k] = (std::exp(psi) / (1 + std::exp(psi))) * stick;
          stick -= pi_[k];
        }
        pi_[unq.size() - 1] = stick;
        for (int k = 0; k < unq.size(); k++) {
          probY(i, j, unq[k]) = pi_[k];
        }
      }
    }
  }
  return probY;
}

//[[Rcpp::export]]
arma::cube prediction(const arma::mat &N,
                      const arma::mat &alpha, 
                      const arma::mat &beta, 
                      const arma::mat &theta,
                      const std::vector<arma::vec> unique_categories,
                      const arma::vec &item_timemap,
                      const String &model,
                      const String &type) 
{
  cube out;
  if (type == "prob") {
    if (model == "static") {
      out = calc_probY(alpha,
                       beta,
                       theta.col(0),
                       unique_categories);
    } else {
      out = dyn_calc_probY(alpha,
                           beta,
                           theta,
                           unique_categories,
                           item_timemap);
    }
  } else {
    if (model == "static") {
      out = calc_predY(N,
                       alpha,
                       beta,
                       theta.col(0),
                       unique_categories);
    } else {
      out = dyn_calc_predY(N,
                           alpha,
                           beta,
                           theta,
                           unique_categories,
                           item_timemap);
    }
  }
  return out;
}
