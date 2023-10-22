#ifndef __POEMIRTBASE__INCLUDED__
#define __POEMIRTBASE__INCLUDED__

#include <RcppArmadillo.h>
#include "helper.h"
using namespace Rcpp;
using namespace arma;

class poEMirtbase
{
public:
  
  poEMirtbase(const cube &Y,
              const mat &N,
              const mat &alpha_init,
              const mat &beta_init,
              const vec &theta_init,
              const std::vector<vec>& unique_categories,
              const mat &a0,
              const mat &A0,
              const mat &b0,
              const mat &B0,
              const int &constraint,
              const bool &std,
              const int &maxit,
              const int &verbose,
              const double &tol,
              const bool &compute_ll);
  ~poEMirtbase();
  
  void get_EOmega();
  mat update_alpha();
  mat update_beta();
  vec update_theta();
  void convcheck(int g);
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const mat &N;
  
  mat alpha, alpha_old;
  const mat &alpha_init;
  mat beta, beta_old;
  const mat &beta_init;
  vec theta, theta_old;
  const vec &theta_init;
  cube Omega;
  
  const std::vector<vec> &unique_categories;
  
  const mat &a0;
  const mat &A0;
  const mat &b0;
  const mat &B0;
  
  cube Nks;
  cube S;
  
  const int &constraint;
  const bool &std;
  const int &maxit;
  const int &verbose;
  const double &tol;
  const bool &compute_ll;
  double ll;
  
  mat convmat;
  bool check;
  vec log_likelihood;
  
  int iter;
  bool converge;
  
  const int I, J, K;
};

#endif
