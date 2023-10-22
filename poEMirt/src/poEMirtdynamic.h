#ifndef __poEMIRTDYNAMIC__INCLUDED__
#define __poEMIRTDYNAMIC__INCLUDED__

#include <RcppArmadillo.h>
#include "helper.h"
using namespace Rcpp;
using namespace arma;

class poEMirtdynamic
{
public:
  
  poEMirtdynamic(const cube &Y,
                 const mat &N,
                 const mat &alpha_init,
                 const mat &beta_init,
                 const mat &theta_init,
                 const std::vector<vec> &unique_categories,
                 const mat &timemap,
                 const mat &timemap2,
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
                 const bool &compute_ll);
  ~poEMirtdynamic();
  
  void get_EOmega();
  mat update_alpha();
  mat update_beta();
  mat update_theta();
  void convcheck(int g);
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const mat &N;
  
  mat alpha;
  mat alpha_old;
  const mat &alpha_init;
  mat beta;
  mat beta_old;
  const mat &beta_init;
  mat theta;
  mat theta_old;
  const mat &theta_init;
  cube Omega;
  
  const std::vector<vec> &unique_categories;
  const mat &timemap;
  const mat &timemap2;
  const vec &item_timemap;
  const vec &item_match;
  std::vector<uvec> time_list;
  std::vector<std::vector<uvec>> ind_item_time_list;
  
  const mat &a0;
  const mat &A0;
  const mat &b0;
  const mat &B0;
  const vec &m0;
  const vec &C0;
  const vec &Delta;
  
  cube S;
  cube Nks;
  
  const vec &constraint;
  const bool &std;
  const int &maxit;
  const int &verbose;
  const double &tol;
  const bool &compute_ll;
  
  mat convmat;
  bool check;
  double ll;
  vec log_likelihood;
  
  int iter;
  bool converge;
  
  const int I;
  const int J;
  const int K;
  const int T;
};

#endif
