#include "poEMirtbase.h"

poEMirtbase::poEMirtbase(const cube &Y,
                         const cube &S,
                         const cube &Nks,
                         mat alpha_old,
                         mat beta_old,
                         vec theta_old,
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
  S(S),
  Nks(Nks),
  alpha_old(alpha_old),
  beta_old(beta_old),
  theta_old(theta_old),
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
  I(S.n_rows),
  J(S.n_cols),
  K(alpha_old.n_cols)
{
  Omega = cube(I, J, K);
  convmat = mat(maxit, 3);
  check = false;
  iter = 0;
  converge = true;
  ll = 0.0;
}

poEMirtbase::~poEMirtbase()
{
  
}

void poEMirtbase::calc_ll()
{
  ll = 0.0;
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          double psi = alpha_old(j, unq[k]) + beta_old(j, unq[k]) * theta_old[i];
          if (Nks(i, j, unq[k]) > 0) {
            ll += S(i, j, unq[k]) * psi - Omega(i, j, unq[k]) * std::pow(psi, 2.0) / 2.0;
          } else {
            break;
          }
        }
      }
    }
  }
  log_likelihood.push_back(ll);
}

void poEMirtbase::get_EOmega()
{
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          double psi = alpha_old(j, unq[k]) + beta_old(j, unq[k]) * theta_old[i];
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

vec poEMirtbase::update_theta()
{
  vec draw(I);
  for (int i = 0; i < I; i++) {
    double sig_part = 0;
    double mu_part = 0;
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j];
      for (unsigned int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          if (Nks(i, j, unq[k]) > 0) {
            sig_part += Omega(i, j, unq[k]) * std::pow(beta_old(j, unq[k]), 2.0);
            mu_part += S(i, j, unq[k]) * beta_old(j, unq[k]) - beta_old(j, unq[k]) * Omega(i, j, unq[k]) * alpha_old(j, unq[k]);
          } else {
            break;
          }
        } 
      } 
    }
    sig_part = sig_part + 1.0;
    draw[i] = mu_part / sig_part;
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
    vec unq = unique_categories[j];
    for (unsigned int k = 0; k < (unq.size()-1); k++) {
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
      draw(j, unq[k]) = mu_part / sig_part;
    }
  }
  return draw;
}

mat poEMirtbase::update_alpha() 
{
  mat draw(J, K);
  for (int j = 0; j < J; j++) {
    vec unq = unique_categories[j];
    for (unsigned int k = 0; k < (unq.size() - 1); k++) {
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
      draw(j, unq[k]) = mu_part / sig_part;
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
  convmat(g, 0) = cor(tmp_alpha1, tmp_alpha2).min();
  convmat(g, 1) = cor(tmp_beta1, tmp_beta2).min();
  convmat(g, 2) = cor(theta, theta_old).min();
  check = ((1 - convmat.row(g).min()) < tol);
}

void poEMirtbase::fit() 
{
  for (int g = 0; g < maxit; g++) {
    checkUserInterrupt();
    
    // Estep
    get_EOmega();
    
    // Mstep
    theta = update_theta();
    beta = update_beta();
    alpha = update_alpha();
  
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

//[[Rcpp::export]]
List poEMirtbase_fit(const arma::cube &Y,
                     const arma::cube &S,
                     const arma::cube &Nks,
                     arma::mat alpha_old,
                     arma::mat beta_old,
                     arma::vec theta_old,
                     const std::vector<arma::vec>& unique_categories,
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
                    S,
                    Nks,
                    alpha_old,
                    beta_old,
                    theta_old,
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