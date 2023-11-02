#ifndef __POEMIRTBASE__INCLUDED__
#define __POEMIRTBASE__INCLUDED__

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class poEMirtbase
{
public:
  
  poEMirtbase(const cube &Y,
              const cube &S,
              const cube &Nks,
              mat alpha_old,
              mat beta_old,
              vec theta_old,
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
  void calc_ll();
  void convcheck(int g);
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const cube &S;
  const cube &Nks;
  
  mat alpha, alpha_old;
  mat beta, beta_old;
  vec theta, theta_old;
  cube Omega;
  
  const std::vector<vec> &unique_categories;
  
  const mat &a0;
  const mat &A0;
  const mat &b0;
  const mat &B0;
  
  const int &constraint;
  const bool &std;
  const int &maxit;
  const int &verbose;
  const double &tol;
  const bool &compute_ll;
  double ll;
  
  mat convmat;
  bool check;
  std::vector<double> log_likelihood;
  
  int iter;
  bool converge;
  
  const int I, J, K;
};

#endif