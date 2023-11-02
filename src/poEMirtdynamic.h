#ifndef __POEMIRTDYNAMIC__INCLUDED__
#define __POEMIRTDYNAMIC__INCLUDED__

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class poEMirtdynamic
{
public:
  
  poEMirtdynamic(const cube &Y,
                 const cube &S,
                 const cube &Nks,
                 mat alpha_old,
                 mat beta_old,
                 mat theta_old,
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
                 const vec &constraint,
                 const bool &alpha_fix,
                 const bool &std,
                 const int &maxit,
                 const int &verbose,
                 const double &tol,
                 const bool &compute_ll);
  ~poEMirtdynamic();
  
  void get_EOmega();
  mat update_alpha();
  mat update_alpha_fixed();
  mat update_beta();
  mat update_theta();
  void calc_ll();
  void convcheck(int g);
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const cube &S;
  const cube &Nks;
  
  mat alpha;
  mat alpha_old;
  mat beta;
  mat beta_old;
  mat theta;
  mat theta_old;
  cube Omega;
  
  const std::vector<vec> &unique_categories;
  const std::vector<vec> &uJ_J;
  const mat &timemap2;
  const vec &item_timemap;
  const std::vector<uvec> &IT;
  const std::vector<std::vector<uvec>> &ITJ;
  
  const mat &a0;
  const mat &A0;
  const mat &b0;
  const mat &B0;
  const vec &m0;
  const vec &C0;
  const vec &Delta;
  
  const bool &alpha_fix;
  const vec &constraint;
  const bool &std;
  const int &maxit;
  const int &verbose;
  const double &tol;
  const bool &compute_ll;
  
  mat convmat;
  bool check;
  double ll;
  std::vector<double> log_likelihood;
  
  int iter;
  bool converge;
  
  const int I;
  const int J;
  const int K;
  const int T;
};

#endif