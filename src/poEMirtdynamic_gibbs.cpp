#include "poEMirtdynamic_gibbs.h"

poEMirtdynamic_gibbs::poEMirtdynamic_gibbs(const cube &Y,
                                           const cube &S,
                                           const cube &Nks,
                                           const cube &sb_check,
                                           mat alpha,
                                           mat beta,
                                           mat theta,
                                           const std::vector<vec> &unique_categories,
                                           const std::vector<vec> &uJ_J,
                                           const mat &timemap2,
                                           const vec &item_timemap,
                                           const std::vector<uvec> &IT,
                                           const std::vector<std::vector<uvec>> &ITJ,
                                           const mat &a0,
                                           const mat &A0,
                                           const mat &b0,
                                           const mat &B0,
                                           const vec &m0,
                                           const vec &C0,
                                           const vec &Delta,
                                           const bool &alpha_fix,
                                           const bool &PG_approx,
                                           const vec &constraint,
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
  sb_check(sb_check),
  alpha(alpha),
  beta(beta),
  theta(theta),
  unique_categories(unique_categories),
  uJ_J(uJ_J),
  timemap2(timemap2),
  item_timemap(item_timemap),
  IT(IT),
  ITJ(ITJ),
  a0(a0),
  A0(A0),
  b0(b0),
  B0(B0),
  m0(m0),
  C0(C0),
  Delta(Delta),
  alpha_fix(alpha_fix),
  PG_approx(PG_approx),
  constraint(constraint),
  std(std),
  iter(iter),
  warmup(warmup),
  thin(thin),
  save_item_parameters(save_item_parameters),
  verbose(verbose),
  I(Y.n_rows),
  J(Y.n_cols),
  K(alpha.n_cols),
  T(timemap2.n_cols)
{
  Omega = cube(I, J, K);
}

poEMirtdynamic_gibbs::~poEMirtdynamic_gibbs()
{
  
}

void poEMirtdynamic_gibbs::draw_Omega()
{
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      if (theta(i, item_timemap[j]) != 0) {
        vec unq = unique_categories[j];
        for (int k = 0; k < (unq.size() - 1); k++) {
          double psi = alpha(j, unq[k]) + beta(j, unq[k]) * theta(i, item_timemap[j]);
          if (!NumericVector::is_na(Y(i, j, unq[k]))) {
            if (Nks(i, j, unq[k]) > 0) {
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
}

void poEMirtdynamic_gibbs::draw_theta()
{
  for (int i = 0; i < I; i++) {
    Rcpp::checkUserInterrupt();
    
    //find attending session
    uvec times_i = IT[i];
    int T_i = times_i.size();
    
    mat sbc_i = sb_check.row(i);
    mat Omega_i = Omega.row(i);
    mat S_i = S.row(i);
    if (sb_check.n_slices == 1) {
      sbc_i = sbc_i.t();
      Omega_i = Omega_i.t();
      S_i = S_i.t();
    }
    
    // FFBS variables
    vec a(T_i); // Prior mean
    vec R(T_i); // Prior var
    vec m(T_i); // Posterior mean
    vec C(T_i); // Posterior var
    
    // Forward filtering
    for (int t = 0; t < T_i; t++) {
      Rcpp::checkUserInterrupt();
      int t_i = times_i[t];
      uvec Jts = ITJ[i][t];
      if (t == 0) {
        a[0] = m0[i];
        R[0] = C0[i] + Delta[i];
      } else {
        a[t] = m[t-1];
        R[t] = C[t-1] + Delta[i];
      }
      if (timemap2(i, t_i) != 0) { // for smoothing
        mat sbc_iJts = sbc_i.rows(Jts);
        mat Omega_iJts = Omega_i.rows(Jts);
        mat S_iJts = S_i.rows(Jts);
        mat alpha_Jts = alpha.rows(Jts);
        mat beta_Jts = beta.rows(Jts) % sbc_iJts;
        vec btmp = sum(beta_Jts, 1);
        
        mat f = beta_Jts * a[t]; // Forecast
        mat Q = R[t] * (btmp * btmp.t()); // Forecast var;
        mat inner = 1.0 / Omega_iJts;
        inner.elem(find_nonfinite(inner)).zeros();
        Q.diag() += sum(inner, 1);
        mat A = R[t] * sum(beta_Jts.t() * Q.i(), 0); // Kalman gain
        C[t] = R[t] - as<double>(wrap(A * Q * A.t())); // Posterior var
        mat frac = S_iJts % inner;
        mat y = (frac - alpha_Jts) % sbc_iJts; // Outcome
        vec e = sum(y - f, 1); // Forecast error
        m[t] = a[t] + as<double>(wrap(A * e)); // Posterior mean
      } else {
        m[t] = a[t];
        C[t] = R[t];
      }
    }
    // Backward sampling
    // for t = T_i
    theta(i, times_i[T_i - 1]) = std::sqrt(C[T_i-1]) * R::rnorm(0, 1) + m[T_i-1];
    for (int t = (T_i - 2); t >= 0; t--) {
      double B = C[t] * (1.0 / R[t+1]);
      double h = m[t] + B * (theta(i, times_i[t+1]) - a[t+1]);
      double H = C[t] - std::pow(B, 2.0) * R[t+1];
      theta(i, times_i[t]) = std::sqrt(H) * R::rnorm(0, 1) + h;
    }
  }
  
  for (int t = 0; t < T; t++) {
    if (theta(constraint[t], t) < 0) {
      theta.col(t) = -theta.col(t);
    }
    // if (std) {
    //   vec theta_t = theta.col(t);
    //   theta.col(t) = (theta_t - mean(theta_t.elem(find(theta_t != 0)))) / stddev(theta_t.elem(find(theta_t != 0)));
    //   theta_t.elem(find(theta_t != 0)).ones();
    //   theta.col(t) = theta.col(t) % theta_t;
    // }
  }
  if (std) {
    mat tmp = theta;
    tmp.elem(find(tmp != 0)).ones();
    double m = mean(theta.elem(find(theta != 0)));
    double sd = stddev(theta.elem(find(theta != 0)));
    theta = ((theta - m) / sd) % tmp;
  }
}

void poEMirtdynamic_gibbs::draw_beta()
{
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j];
    for (unsigned int k = 0; k < (unq.size() - 1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]) * std::pow(theta(i, item_timemap[j]), 2.0);
            mu_part += theta(i, item_timemap[j]) * (S(i, j, unq[k]) - Omega(i, j, unq[k]) * alpha(j, unq[k]));
          }
        }
      }
      sig_part += 1 / B0(j, unq[k]);
      mu_part += b0(j, unq[k]) / B0(j, unq[k]);
      beta(j, unq[k]) = R::rnorm(mu_part / sig_part, std::sqrt(1.0 / sig_part));
    }
  }
}

void poEMirtdynamic_gibbs::draw_alpha() 
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
            mu_part += S(i, j, unq[k]) - beta(j, unq[k]) * Omega(i, j, unq[k]) * theta(i, item_timemap[j]);
          }
        }  
      }
      sig_part += 1 / A0(j, unq[k]);
      mu_part += a0(j, unq[k]) / A0(j, unq[k]);
      alpha(j, unq[k]) = R::rnorm(mu_part / sig_part, std::sqrt(1.0 / sig_part));
    }
  }
}

void poEMirtdynamic_gibbs::draw_alpha_fixed()
{
  for (unsigned int uj = 0; uj < uJ_J.size(); uj++) {
    rowvec sig(K);
    rowvec mu(K);
    vec Juj = uJ_J[uj];
    for (unsigned int jj = 0; jj < Juj.size(); jj++) {
      int j = Juj[jj];
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size()-1); k++) {
        for (unsigned int i = 0; i < I; i++) {
          if (!NumericVector::is_na(Y(i, j, unq[k]))) {
            if (Nks(i, j, unq[k]) > 0) {
              sig[unq[k]] += Omega(i, j, unq[k]);
              mu[unq[k]] += S(i, j, unq[k]) - Omega(i, j, unq[k]) * (theta(i, item_timemap[j]) * beta(j, unq[k]));
            }
          }
        }
        sig[unq[k]] += 1.0 / A0(j, unq[k]);
        mu[unq[k]] += a0(j, unq[k]) / A0(j, unq[k]); 
      }
    }
    
    uvec non0 = find(sig != 0);
    rowvec draw(K);
    for (unsigned int jj = 0; jj < Juj.size(); jj++) {
      int j = Juj[jj];
      for (unsigned int k = 0; k < non0.size(); k++) {
        if (jj == 0) {
          draw[non0[k]] = R::rnorm(mu[non0[k]] / sig[non0[k]], std::sqrt(1.0 / sig[non0[k]]));
        }
      }
      alpha.row(j) = draw;
    }
  }
}

void poEMirtdynamic_gibbs::fit() 
{
  int g = 0;
  int total_iter = iter + warmup;
  if (warmup != 0) {
    for (; g < warmup; g++) {
      checkUserInterrupt();
      draw_Omega();
      draw_theta();
      draw_beta();
      if (alpha_fix) {
        draw_alpha_fixed();
      } else {
        draw_alpha();
      }
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
    if (alpha_fix) {
      draw_alpha_fixed();
    } else {
      draw_alpha();
    }
    
    if (g % thin == 0) {
      theta_store.push_back(theta);
      if (save_item_parameters) {
        beta_store.push_back(beta);
        alpha_store.push_back(alpha);
      }
    }
    
    if (g == 0 || (g+1) % verbose == 0) {
      Rcout << "* Sampling " << g + 1 << " / " << total_iter << endl;
    }
  }
  
  modeloutput = List::create(
    Named("alpha") = alpha_store,
    Named("beta") = beta_store,
    Named("theta") = theta_store
  );
}

List poEMirtdynamic_gibbs::output() 
{
  return modeloutput;
}

//[[Rcpp::export]]
List poEMirtdynamic_gibbs_fit(const arma::cube &Y,
                              const arma::cube &S,
                              const arma::cube &Nks,
                              const arma::cube &sb_check,
                              arma::mat alpha,
                              arma::mat beta,
                              arma::mat theta,
                              const std::vector<arma::vec> &unique_categories,
                              const std::vector<arma::vec> &uJ_J,
                              const arma::mat &timemap2,
                              const arma::vec &item_timemap,
                              const std::vector<arma::uvec> &IT,
                              const std::vector<std::vector<arma::uvec>> &ITJ,
                              const arma::mat &a0,
                              const arma::mat &A0,
                              const arma::mat &b0,
                              const arma::mat &B0,
                              const arma::vec &m0,
                              const arma::vec &C0,
                              const arma::vec &Delta,
                              const bool &alpha_fix,
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
                             S,
                             Nks,
                             sb_check,
                             alpha,
                             beta,
                             theta,
                             unique_categories,
                             uJ_J,
                             timemap2,
                             item_timemap,
                             IT,
                             ITJ,
                             a0,
                             A0,
                             b0,
                             B0,
                             m0,
                             C0,
                             Delta,
                             alpha_fix,
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