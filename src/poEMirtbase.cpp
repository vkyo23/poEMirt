#include "poEMirtbase.h"

poEMirtbase::poEMirtbase(const cube &Y,
                         const mat &N,
                         const mat &alpha_init,
                         const mat &beta_init,
                         const vec &theta_init,
                         const std::vector<vec> &unique_categories,
                         const mat &a0,
                         const mat &A0,
                         const mat &b0,
                         const mat &B0,
                         const int &constraint,
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
  a0(a0),
  A0(A0),
  b0(b0),
  B0(B0),
  constraint(constraint),
  std(std),
  maxit(maxit),
  verbose(verbose),
  tol(tol),
  compute_ll(compute_ll),
  I(Y.n_rows),
  J(Y.n_cols),
  K(alpha_init.n_cols)
{
  Omega = cube(I, J, K);
  convmat = mat(maxit, 3);
  log_likelihood = vec(maxit);
  check = false;
  iter = 0;
  converge = true;
  ll = 0.0;
}

poEMirtbase::~poEMirtbase()
{
  
}

void poEMirtbase::get_EOmega()
{
  if (compute_ll) {
    ll = 0;
  }
  for (int i = 0; i < I; i++) {
    if (compute_ll) {
      ll -= std::pow(theta_old[i], 2.0) / 2.0;
    }
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j]-1;
      for (int k = 0; k < (unq.size() - 1); k++) {
        double psi = alpha_old(j, unq[k]) + beta_old(j, unq[k]) * theta_old[i];
        if (Nks(i, j, unq[k]) > 0) {
          Omega(i, j, unq[k]) = (Nks(i, j, unq[k]) / (2 * psi)) * std::tanh(psi / 2);
          if (compute_ll) {
            ll += S(i, j, unq[k]) * psi - Omega(i, j, unq[k]) * std::pow(psi, 2.0) / 2.0;
            if (i == 0) {
              ll = ll - std::pow(beta_old(j, unq[k]) - b0(j, unq[k]), 2.0) / (2.0 * B0(j, unq[k])) -
                std::pow(alpha_old(j, unq[k]) - a0(j, unq[k]), 2.0) / (2.0 * A0(j, unq[k]));
            }
          }
        } else {
          break;
        }
      }
    }
  }
}

vec poEMirtbase::update_theta()
{
  vec draw(I);
  for (int i = 0; i < I; i++) {
    double sig_part = 0;
    double mu_part = 0;
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j]-1;
      for (int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]) * std::pow(beta_old(j, unq[k]), 2.0);
            mu_part += S(i, j, unq[k]) * beta_old(j, unq[k]) - beta_old(j, unq[k]) * Omega(i, j, unq[k]) * alpha_old(j, unq[k]);
          } else {
            break;
          }
        } else {
          break;
        }
      } 
    }
    sig_part = sig_part + 1.0;
    draw[i] = (1 / sig_part) * mu_part;
  }
  if (draw[constraint] < 0) {
    draw = -draw;
  }
  if (std) {
    draw = (draw - mean(draw)) / stddev(draw);
  }
  return draw;
}

mat poEMirtbase::update_beta()
{
  mat draw(J, K);
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j]-1;
    for (int k = 0; k < (unq.size()-1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          sig_part += Omega(i, j, unq[k]) * std::pow(theta[i], 2.0);
          mu_part += theta[i] * (S(i, j, unq[k]) - Omega(i, j, unq[k]) * alpha_old(j, unq[k]));
        } 
      }
      sig_part += 1 / B0(j, unq[k]);
      mu_part += b0(j, unq[k]) / B0(j, unq[k]);
      draw(j, unq[k]) = (1 / sig_part) * mu_part;
    }
  }
  return draw;
}

mat poEMirtbase::update_alpha() 
{
  mat draw(J, K);
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j]-1;
    for (int k = 0; k < (unq.size() - 1); k++) {
      double sig_part = 0;
      double mu_part = 0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          sig_part += Omega(i, j, unq[k]);
          mu_part += S(i, j, unq[k]) - Omega(i, j, unq[k]) * theta[i] * beta(j, unq[k]);
        } 
      }
      sig_part += 1 / A0(j, unq[k]);
      mu_part += a0(j, unq[k]) / A0(j, unq[k]);
      draw(j, unq[k]) = (1 / sig_part) * mu_part;
    }
  }
  return draw;
}

void poEMirtbase::convcheck(int g)
{
  vec tmp_alpha1 = alpha_old.elem(find(alpha_old != 0)); 
  vec tmp_alpha2 = alpha.elem(find(alpha != 0));
  vec tmp_beta1 = beta_old.elem(find(beta_old != 0));
  vec tmp_beta2 = beta.elem(find(beta != 0));
  convmat(g-1, 0) = cor(tmp_alpha1, tmp_alpha2).min();
  convmat(g-1, 1) = cor(tmp_beta1, tmp_beta2).min();
  convmat(g-1, 2) = cor(theta, theta_old).min();
  
  check = ((1 - convmat.row(g-1).min()) < tol);
}

void poEMirtbase::fit() 
{
  
  // Assign initial values
  alpha_old = alpha_init;
  beta_old = beta_init;
  theta_old = theta_init;
  
  // Aux
  List tmp = construct_sb_auxs(Y, N, unique_categories);
  S = as<cube>(wrap(tmp["S"]));
  Nks = as<cube>(wrap(tmp["Nks"]));
  
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

List poEMirtbase::output() 
{
  return modeloutput;
}