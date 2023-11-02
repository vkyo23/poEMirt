#include <RcppArmadillo.h>
#include "dist.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List construct_sb_auxs(const arma::cube &Y,
                       const arma::mat &N,
                       const std::vector<arma::vec> &unique_categories) 
{
  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = Y.n_slices - 1;
  
  cube sb_check(I, J, K);
  cube Nks(I, J, K);
  cube S(I, J, K);
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j];
      double Ycum = 0;
      for (int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          Nks(i, j, unq[k]) = N(i, j) - Ycum;
          Ycum += Y(i, j, unq[k]);
          if (Nks(i, j, unq[k]) > 0) {
            sb_check(i, j, unq[k]) = 1;
            S(i, j, unq[k]) = Y(i, j, unq[k]) - Nks(i, j, unq[k]) / 2.0;
          }
        }
      }
    }
  }
  
  List L = List::create(
    Named("sb_check") = sb_check,
    Named("Nks") = Nks,
    Named("S") = S
  );
  return L;
}

//[[Rcpp::export]]
List get_dynamic_info (const arma::mat &N,
                       const arma::vec &item_timemap,
                       const arma::mat &timemap) 
{
  const int& I = timemap.n_rows;
  const int& J = item_timemap.size();
  const int& T = timemap.n_cols;
  
  std::vector<std::vector<uvec>> ITJ;
  std::vector<uvec> IT;
  
  for (int i = 0; i < I; i++) {
    int j = 0;
    std::vector<int> ind_t_index;
    std::vector<uvec> item_it_list;
    for (int t = 0; t < T; t++) {
      if (timemap(i, t) == 1) {
        ind_t_index.push_back(t);
      }
      std::vector<int> item_it;
      for (; j < J; j++) {
        if (item_timemap[j] == t) {
          if (N(i, j) != 0) {
            item_it.push_back(j);
          }
          if (item_timemap[j+1] == t+1) {
            break;
          }
        }
      }
      if (timemap(i, t) == 1) {
        item_it_list.push_back(as<uvec>(wrap(item_it)));
      }
    }
    ITJ.push_back(item_it_list);
    IT.push_back(as<uvec>(wrap(ind_t_index)));
  }
  
  List L = List::create(
    Named("ITJ") = ITJ,
    Named("IT") = IT
  );
  
  return L;
}


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
        vec unq = unique_categories[j];
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
          vec unq = unique_categories[j];
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
      vec unq = unique_categories[j];
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
        vec unq = unique_categories[j];
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