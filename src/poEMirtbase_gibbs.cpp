#include "poEMirtbase_gibbs.h"

poEMirtbase_gibbs::poEMirtbase_gibbs(const cube &Y,
                                     const cube &S,
                                     const cube &Nks,
                                     mat alpha,
                                     mat beta,
                                     vec theta,
                                     const std::vector<vec>& unique_categories,
                                     const mat &a0,
                                     const mat &A0,
                                     const mat &b0,
                                     const mat &B0,
                                     const bool& PG_approx,
                                     const int &constraint,
                                     const bool &std,
                                     const int &iter,
                                     const int &warmup,
                                     const int &thin,
                                     const bool &save_item_parameters,
                                     const int &verbose)
  :
  Y(Y),
  S(S),
  Nks(Nks),
  alpha(alpha),
  beta(beta),
  theta(theta),
  unique_categories(unique_categories),
  a0(a0),
  A0(A0),
  b0(b0),
  B0(B0),
  PG_approx(PG_approx),
  constraint(constraint),
  std(std),
  iter(iter),
  warmup(warmup),
  thin(thin),
  save_item_parameters(save_item_parameters),
  verbose(verbose),
  I(S.n_rows),
  J(S.n_cols),
  K(alpha.n_cols)
{
  Omega = cube(I, J, K);
  theta_store = mat(I, iter / thin);
}

poEMirtbase_gibbs::~poEMirtbase_gibbs()
{
  
}

void poEMirtbase_gibbs::draw_Omega()
{
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          double psi = alpha(j, unq[k]) + beta(j, unq[k]) * theta[i];
          if (Nks(i, j, unq[k]) != 0) {
            if (Nks(i, j, unq[k]) < 20 && !PG_approx) {
              Omega(i, j, unq[k]) = rpg(Nks(i, j, unq[k]), psi);
            } else {
              Omega(i, j, unq[k]) = rpg_approx(Nks(i, j, unq[k]), psi);
            }
          } else {
            break;
          }
        }
      }
    }
  }
}

void poEMirtbase_gibbs::draw_theta()
{
  for (int i = 0; i < I; i++) {
    double sig_part = 0;
    double mu_part = 0;
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]) * std::pow(beta(j, unq[k]), 2.0);
            mu_part += beta(j, unq[k]) * (S(i, j, unq[k]) - Omega(i, j, unq[k]) * alpha(j, unq[k]));
          } else {
            break;
          }
        } 
      } 
    }
    sig_part = sig_part + 1.0;
    theta[i] = R::rnorm(mu_part / sig_part, std::sqrt(1 / sig_part));
  }
  if (theta[constraint] < 0) {
    theta = -theta;
  }
  if (std) {
    theta = (theta - mean(theta)) / stddev(theta);
  }
}

void poEMirtbase_gibbs::draw_beta()
{
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j];
    for (unsigned int k = 0; k < (unq.size()-1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]) * std::pow(theta[i], 2.0);
            mu_part += theta[i] * (S(i, j, unq[k]) - Omega(i, j, unq[k]) * alpha(j, unq[k]));
          }
        } 
      }
      sig_part += 1 / B0(j, unq[k]);
      mu_part += b0(j, unq[k]) / B0(j, unq[k]);
      beta(j, unq[k]) = R::rnorm(mu_part / sig_part, std::sqrt(1 / sig_part));
    }
  }
}

void poEMirtbase_gibbs::draw_alpha() 
{
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j];
    for (unsigned int k = 0; k < (unq.size() - 1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]);
            mu_part += S(i, j, unq[k]) - Omega(i, j, unq[k]) * theta[i] * beta(j, unq[k]);
          }
        } 
      }
      sig_part += 1 / A0(j, unq[k]);
      mu_part +=  a0(j, unq[k]) / A0(j, unq[k]);
      alpha(j, unq[k]) = R::rnorm(mu_part / sig_part, std::sqrt(1 / sig_part));
    }
  }
}

void poEMirtbase_gibbs::fit() 
{
  int g = 0;
  int total_iter = iter + warmup;
  if (warmup != 0) {
    for (; g < warmup; g++) {
      checkUserInterrupt();
      draw_Omega();
      draw_theta();
      draw_beta();
      draw_alpha();
      if (g == 0 || (g+1) % verbose == 0) {
        Rcout << "* Warmup " << g + 1 << " / " << total_iter << endl;
      }
    }
  }
  
  for (; g < total_iter; g++) {
    checkUserInterrupt();
    draw_Omega();
    draw_theta();
    draw_beta();
    draw_alpha();
    
    if (g % thin == 0) {
      theta_store.col((g - warmup) / thin) = theta;
      if (save_item_parameters) {
        beta_store.push_back(beta);
        alpha_store.push_back(alpha);
      }
    }
    
    if ((g+1) % verbose == 0) {
      Rcout << "* Sampling " << g + 1 << " / " << total_iter << endl;
    }
  }
  modeloutput = List::create(
    Named("alpha") = alpha_store,
    Named("beta") = beta_store,
    Named("theta") = theta_store
  );
}

List poEMirtbase_gibbs::output() 
{
  return modeloutput;
}

//[[Rcpp::export]]
List poEMirtbase_gibbs_fit(const arma::cube &Y,
                           const arma::cube &S,
                           const arma::cube &Nks,
                           arma::mat alpha,
                           arma::mat beta,
                           arma::vec theta,
                           const std::vector<arma::vec>& unique_categories,
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
                          S,
                          Nks,
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