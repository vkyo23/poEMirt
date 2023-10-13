#include "poEMirtdynamic.h"

poEMirtdynamic::poEMirtdynamic(const cube &Y,
                               const mat &N,
                               const mat &alpha_init,
                               const mat &beta_init,
                               const mat &theta_init,
                               const std::vector<vec> &unique_categories,
                               const mat &timemap,
                               const vec &item_timemap,
                               const vec &item_match,
                               const mat &a0,
                               const mat &A0,
                               const mat &b0,
                               const mat &B0,
                               const vec &m0,
                               const vec &C0,
                               const vec &Delta,
                               const vec &constraint,
                               const bool &std,
                               const int &maxit,
                               const int &verbose,
                               const double &tol,
                               const bool &compute_ll)
  :
  Y(Y),
  N(N),
  alpha_init(alpha_init),
  beta_init(beta_init),
  theta_init(theta_init),
  unique_categories(unique_categories),
  timemap(timemap),
  item_timemap(item_timemap),
  item_match(item_match),
  a0(a0),
  A0(A0),
  b0(b0),
  B0(B0),
  m0(m0),
  C0(C0),
  Delta(Delta),
  constraint(constraint),
  std(std),
  maxit(maxit),
  verbose(verbose),
  tol(tol),
  compute_ll(compute_ll),
  I(Y.n_rows),
  J(Y.n_cols),
  K(alpha_init.n_cols),
  T(timemap.n_cols)
{
  Omega = cube(I, J, K);
  convmat = mat(maxit, 3);
  log_likelihood = vec(maxit);
  check = false;
  iter = 0;
  converge = true;
  ll = 0.0;
}

poEMirtdynamic::~poEMirtdynamic()
{
  
}

void poEMirtdynamic::get_EOmega()
{
  if (compute_ll) {
    ll = 0.0;
  }
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      if (theta_old(i, item_timemap[j]) != 0) {
        vec unq = unique_categories[j]-1; 
        for (int k = 0; k < (unq.size() - 1); k++) {
          double psi = alpha_old(j, unq[k]) + beta_old(j, unq[k]) * theta_old(i, item_timemap[j]);
          if (Nks(i, j, unq[k]) > 0) {
            Omega(i, j, unq[k]) = (Nks(i, j, unq[k]) / (2 * psi)) * std::tanh(psi / 2);
            if (compute_ll) {
              ll += S(i, j, unq[k]) * psi - Omega(i, j, unq[k]) * std::pow(psi, 2.0) / 2.0;
            }
          } else {
            break;
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
    uvec times_i = time_list[i];
    int T_i = times_i.size();
    
    for (int t = 0; t < T_i; t++) {
      int t_i = times_i[t];
      uvec item_index_t = item_time_list[t_i];
      int Jt = item_index_t.size();
      double sig_part = 0;
      double mu_part = 0;
      for (int j = 0; j < Jt; j++) {
        int jt = item_index_t[j];
        vec unq = unique_categories[jt]-1;
        for (int k = 0; k < (unq.size() - 1); k++) {
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
      if (t == 0) {
        double sig = sig_part + 1 / C0[i] + 1 / Delta[i];
        double mu = mu_part + m0[i] / C0[i];
        
        A(t, t) = sig;
        A(t, t+1) = -1 / Delta[i];
        A(t+1, t) = -1 / Delta[i];
        B(t) = mu;
      } else if (t != (T_i - 1)) {
        double sig = sig_part + 2 / Delta[i] + 1;
        double mu = mu_part;
        
        A(t, t) = sig;
        A(t, t+1) = -1 / Delta[i];
        A(t+1, t) = -1 / Delta[i];
        B(t) = mu;
      } else if (t == (T_i - 1)) {
        double sig = sig_part + 1 / Delta[i] + 1;
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
    if (std) {
      vec draw_t = draw.col(t);
      draw.col(t) = (draw_t - mean(draw_t.elem(find(draw_t != 0)))) / stddev(draw_t.elem(find(draw_t != 0)));
      draw_t.elem(find(draw_t != 0)).ones();
      draw.col(t) = draw.col(t) % draw_t;
    }
  }
  return draw;
}


mat poEMirtdynamic::update_beta()
{
  mat draw(J, K);
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j] - 1;
    for (int k = 0; k < (unq.size() - 1); k++) {
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
      draw(j, unq[k]) = (1 / sig_part) * mu_part;
    }
  }
  return draw;
}

mat poEMirtdynamic::update_alpha() 
{
  mat draw(J, K);
  bool flag;
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j] - 1;
    for (int k = 0; k < (unq.size() - 1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      if (!NumericVector::is_na(item_match[j])) {
        if (draw(item_match[j], unq[k]) != 0) {
          draw(j, unq[k]) = draw(item_match[j], unq[k]);
          flag = false;
        } else {
          flag = true;
        }
      } else {
        flag = true;
      }
      if (flag) {
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
        draw(j, unq[k]) = (1 / sig_part) * mu_part;
      }
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
  convmat(g-1, 0) = cor(tmp_alpha1, tmp_alpha2).min();
  convmat(g-1, 1) = cor(tmp_beta1, tmp_beta2).min();
  convmat(g-1, 2) = cor(tmp_theta1, tmp_theta2).min();
  
  check = ((1 - convmat.row(g-1).min()) < tol);
}

// EM
void poEMirtdynamic::fit() 
{
  // Initial values
  alpha_old = alpha_init;
  beta_old = beta_init;
  theta_old = theta_init;
  
  // Auxiliary components for dynamic estimation
  List time_info = get_time_info(item_timemap, timemap);
  time_list = as<std::vector<uvec>>(wrap(time_info["time_list"]));
  item_time_list = as<std::vector<uvec>>(wrap(time_info["item_time_list"]));
  
  // Aux 2
  List tmp = construct_sb_auxs(Y, N, unique_categories);
  Nks = as<cube>(wrap(tmp["Nks"]));
  S = as<cube>(wrap(tmp["S"]));
  
  for (int g = 0; g < maxit; g++) {
    checkUserInterrupt();
    
    // Estep
    get_EOmega();
    if (g != 0) {
      if (compute_ll) {
        log_likelihood[g-1] = ll;
      }
    }
    
    // Mstep
    theta = update_theta();
    beta = update_beta();
    alpha = update_alpha();
    
    if (g != 0) {
      convcheck(g);
      if (check) {
        uvec seqq = as<uvec>(wrap(seq(0, g-1)));
        convmat = convmat.rows(seqq);
        if (compute_ll) {
          get_EOmega();
          log_likelihood[g] = ll;
        }
        log_likelihood = log_likelihood.rows(0, g);
        iter = g;
        break;
      } else if (g == (maxit - 1)) {
        converge = false;
        iter = maxit;
        break;
      } else if (g % verbose == 0) {
        Rcout << "Iteration " << g << ": eval = " << (1 - convmat.row(g-1).min()) << '\n';
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