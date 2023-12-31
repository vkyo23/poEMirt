#include "poEMirtdynamic.h"

poEMirtdynamic::poEMirtdynamic(const cube &Y,
                               const cube &S,
                               const cube &Nks,
                               mat alpha_old,
                               mat beta_old,
                               mat theta_old,
                               double Delta,
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
                               const double &g0,
                               const double &h0,
                               const vec &constraint,
                               const bool &fix_alpha,
                               const bool &fix_beta,
                               const bool &estimate_Delta,
                               const bool &std,
                               const int &maxit,
                               const int &verbose,
                               const double &tol,
                               const bool &compute_ll)
  :
  Y(Y),
  S(S),
  Nks(Nks),
  alpha_old(alpha_old),
  beta_old(beta_old),
  theta_old(theta_old),
  Delta(Delta),
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
  g0(g0),
  h0(h0),
  constraint(constraint),
  fix_alpha(fix_alpha),
  fix_beta(fix_beta),
  estimate_Delta(estimate_Delta),
  std(std),
  maxit(maxit),
  verbose(verbose),
  tol(tol),
  compute_ll(compute_ll),
  I(S.n_rows),
  J(S.n_cols),
  K(alpha_old.n_cols),
  T(timemap2.n_cols)
{
  Omega = cube(I, J, K);
  convmat = mat(maxit, 3);
  check = false;
  iter = 0;
  converge = true;
  ll = 0.0;
}

poEMirtdynamic::~poEMirtdynamic()
{
  
}

void poEMirtdynamic::calc_ll()
{
  ll = 0.0;
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      if (theta_old(i, item_timemap[j]) != 0) {
        vec unq = unique_categories[j]; 
        for (unsigned int k = 0; k < (unq.size() - 1); k++) {
          if (!NumericVector::is_na(Y(i, j, unq[k]))) {
            double psi = alpha_old(j, unq[k]) + beta_old(j, unq[k]) * theta_old(i, item_timemap[j]);
            if (Nks(i, j, unq[k]) > 0) {
              ll += S(i, j, unq[k]) * psi - Omega(i, j, unq[k]) * std::pow(psi, 2.0) / 2.0;
            } else {
              break;
            }
          }
        } 
      }
    }
  }
  log_likelihood.push_back(ll);
}

void poEMirtdynamic::get_EOmega()
{
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      if (theta_old(i, item_timemap[j]) != 0) {
        vec unq = unique_categories[j]; 
        for (unsigned int k = 0; k < (unq.size() - 1); k++) {
          if (!NumericVector::is_na(Y(i, j, unq[k]))) {
            double psi = alpha_old(j, unq[k]) + beta_old(j, unq[k]) * theta_old(i, item_timemap[j]);
            if (Nks(i, j, unq[k]) > 0) {
              Omega(i, j, unq[k]) = (Nks(i, j, unq[k]) / (2 * psi)) * std::tanh(psi / 2);
            } else {
              break;
            }
          }
        } 
      }
    }
  }
}

mat poEMirtdynamic::update_theta()
{
  mat draw(I, T);
  for (int i = 0; i < I; i++) {
    vec X_store(T);
    mat A(T, T);
    vec B(T);
    
    //find attending session
    uvec times_i = IT[i];
    int T_i = times_i.size();
    
    for (int t = 0; t < T_i; t++) {
      int t_i = times_i[t];
      double sig_part = 0;
      double mu_part = 0;
      if (timemap2(i, t_i) == 1) {
        uvec item_index_t = ITJ[i][t];
        for (unsigned int j = 0; j < item_index_t.size(); j++) {
          int jt = item_index_t[j];
          vec unq = unique_categories[jt];
          for (unsigned int k = 0; k < (unq.size() - 1); k++) {
            if (!NumericVector::is_na(Y(i, jt, unq[k]))) {
              if (Nks(i, jt, unq[k]) > 0) {
                sig_part += Omega(i, jt, unq[k]) * std::pow(beta_old(jt, unq[k]), 2.0);
                mu_part += beta_old(jt, unq[k]) * (S(i, jt, unq[k]) - Omega(i, jt, unq[k]) * alpha_old(jt, unq[k]));
              } else {
                break;
              }
            }
          }
        }
      }
      if (t == 0) {
        double sig = sig_part + 1 / C0[i] + 1 / Delta;
        double mu = mu_part + m0[i] / C0[i];
        
        A(t, t) = sig;
        A(t, t+1) = -1 / Delta;
        A(t+1, t) = -1 / Delta;
        B(t) = mu;
      } else if (t != (T_i - 1)) {
        double sig = sig_part + 2 / Delta + 1;
        double mu = mu_part;
        
        A(t, t) = sig;
        A(t, t+1) = -1 / Delta;
        A(t+1, t) = -1 / Delta;
        B(t) = mu;
      } else if (t == (T_i - 1)) {
        double sig = sig_part + 1 / Delta + 1;
        double mu = mu_part;
        A(t, t) = sig;
        B(t) = mu;
      }
    }
    uvec diag_w = find(A.diag() != 0);
    mat AA = A.submat(diag_w, diag_w);
    vec BB = B.rows(diag_w);
    vec X = inv(AA) * BB;
    X_store.rows(times_i) = X;
    draw.row(i) = X_store.t();
  }
  
  for (int t = 0; t < T; t++) {
    if (draw(constraint[t], t) < 0) {
      draw.col(t) = -draw.col(t);
    }
    // if (std) {
    //   vec draw_t = draw.col(t);
    //   draw.col(t) = (draw_t - mean(draw_t.elem(find(draw_t != 0)))) / stddev(draw_t.elem(find(draw_t != 0)));
    //   draw_t.elem(find(draw_t != 0)).ones();
    //   draw.col(t) = draw.col(t) % draw_t;
    // }
  }
  if (std) {
    mat tmp = draw;
    tmp.elem(find(tmp != 0)).ones();
    double m = mean(draw.elem(find(draw != 0)));
    double sd = stddev(draw.elem(find(draw != 0)));
    draw = ((draw - m) / sd) % tmp;
  }
  return draw;
}

void poEMirtdynamic::update_Delta()
{
  double tmp = 0.0;
  double tmp2 = 0.0;
  for (int i = 0; i < I; i++) {
    uvec times_i = IT[i];
    int T_i = times_i.size();
    tmp2 += (double)T_i - 1;
    for (int t = 1; t < T_i; t++) {
      tmp += std::pow(theta(i, times_i[t]) - theta(i, times_i[t-1]), 2.0);
    }
  }
  Delta = (h0 + tmp) / (tmp2 + g0 + 2.0);
}

mat poEMirtdynamic::update_beta()
{
  mat draw(J, K);
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j];
    for (unsigned int k = 0; k < (unq.size() - 1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]) * std::pow(theta(i, item_timemap[j]), 2.0);
            mu_part += theta(i, item_timemap[j]) * (S(i, j, unq[k]) - Omega(i, j, unq[k]) * alpha_old(j, unq[k]));
          }
        }
      }
      sig_part += 1 / B0(j, unq[k]);
      mu_part += b0(j, unq[k]) / B0(j, unq[k]);
      draw(j, unq[k]) = mu_part / sig_part;
    }
  }
  return draw;
}

mat poEMirtdynamic::update_beta_fixed()
{
  mat draw(J, K);
  for (unsigned int uj = 0; uj < uJ_J.size(); uj++) {
    rowvec sig(K);
    rowvec mu(K);
    vec Juj = uJ_J[uj];
    for (unsigned int jj = 0; jj < Juj.size(); jj++) {
      int j = Juj[jj];
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size()-1); k++) {
        for (int i = 0; i < I; i++) {
          if (!NumericVector::is_na(Y(i, j, unq[k]))) {
            if (Nks(i, j, unq[k]) > 0) {
              sig[unq[k]] += Omega(i, j, unq[k]) * std::pow(theta(i, item_timemap[j]), 2.0);
              mu[unq[k]] += theta(i, item_timemap[j]) * (S(i, j, unq[k]) - Omega(i, j, unq[k]) * alpha_old(j, unq[k]));
            }
          }
        }
        sig[unq[k]] += 1.0 / B0(j, unq[k]);
        mu[unq[k]] += b0(j, unq[k]) / B0(j, unq[k]); 
      }
    }
    sig.elem(find(sig == 0)).ones();
    for (unsigned int jj = 0; jj < Juj.size(); jj++) {
      int j = Juj[jj];
      draw.row(j) = mu / sig;
    }
  }
  return draw;
}

mat poEMirtdynamic::update_alpha() 
{
  mat draw(J, K);
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
      draw(j, unq[k]) = mu_part / sig_part;
    }
  }
  return draw;
}

mat poEMirtdynamic::update_alpha_fixed()
{
  mat draw(J, K);
  for (unsigned int uj = 0; uj < uJ_J.size(); uj++) {
    rowvec sig(K);
    rowvec mu(K);
    vec Juj = uJ_J[uj];
    for (unsigned int jj = 0; jj < Juj.size(); jj++) {
      int j = Juj[jj];
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size()-1); k++) {
        for (int i = 0; i < I; i++) {
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
    sig.elem(find(sig == 0)).ones();
    for (unsigned int jj = 0; jj < Juj.size(); jj++) {
      int j = Juj[jj];
      draw.row(j) = mu / sig;
    }
  }
  return draw;
}

void poEMirtdynamic::convcheck(int g)
{
  vec tmp_alpha1 = alpha_old.elem(find(alpha_old != 0)); 
  vec tmp_alpha2 = alpha.elem(find(alpha != 0));
  vec tmp_beta1 = beta_old.elem(find(beta_old != 0));
  vec tmp_beta2 = beta.elem(find(beta != 0));
  vec tmp_theta1 = theta_old.elem(find(theta_old != 0));
  vec tmp_theta2 = theta.elem(find(theta != 0));
  convmat(g, 0) = cor(tmp_alpha1, tmp_alpha2).min();
  convmat(g, 1) = cor(tmp_beta1, tmp_beta2).min();
  convmat(g, 2) = cor(tmp_theta1, tmp_theta2).min();
  check = ((1 - convmat.row(g).min()) < tol);
}

// EM
void poEMirtdynamic::fit() 
{
  for (int g = 0; g < maxit; g++) {
    checkUserInterrupt();
    
    // Estep
    get_EOmega();
    
    // Mstep
    theta = update_theta();
    if (estimate_Delta) {
      update_Delta();
    }
    if (fix_beta) {
      beta = update_beta_fixed();
    } else {
      beta = update_beta();
    }
    if (fix_alpha) {
      alpha = update_alpha_fixed();
    } else {
      alpha = update_alpha();
    }
    
    if (g != 0) {
      convcheck(g);
      if (compute_ll) {
        calc_ll();
      }
      if (check) {
        convmat = convmat.rows(0, g);
        iter = g;
        break;
      } else if (g == (maxit - 1)) {
        converge = false;
        iter = maxit;
        break;
      } else if (g % verbose == 0) {
        Rcout << "  - Iteration " << g << ": eval = " << (1 - convmat.row(g).min()) << '\n';
      }
    }
    
    theta_old = theta;
    beta_old = beta;
    alpha_old = alpha;
  }
  
  modeloutput = List::create(
    Named("alpha") = alpha,
    Named("beta") = beta,
    Named("theta") = theta,
    Named("Delta") = Delta,
    Named("iter") = iter,
    Named("conv") = convmat,
    Named("converge") = converge,
    Named("log_likelihood") = log_likelihood
  );
}

List poEMirtdynamic::output() 
{
  return modeloutput;
}

//[[Rcpp::export]]
List poEMirtdynamic_fit(const arma::cube &Y,
                        const arma::cube &S,
                        const arma::cube &Nks,
                        arma::mat alpha_old,
                        arma::mat beta_old,
                        arma::mat theta_old,
                        double Delta,
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
                        const double &g0,
                        const double &h0,
                        const arma::vec &constraint,
                        const bool &fix_alpha,
                        const bool &fix_beta,
                        const bool &estimate_Delta,
                        const bool &std,
                        const int &maxit,
                        const int &verbose,
                        const double &tol,
                        const bool &compute_ll) 
{
  // Instance
  poEMirtdynamic Model(Y,
                       S,
                       Nks,
                       alpha_old,
                       beta_old,
                       theta_old,
                       Delta,
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
                       g0,
                       h0,
                       constraint,
                       fix_alpha,
                       fix_beta,
                       estimate_Delta,
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