#ifndef __POEMIRTDYNAMIC_GIBBS__INCLUDED__
#define __POEMIRTDYNAMIC_GIBBS__INCLUDED__

#include <RcppArmadillo.h>
#include "dist.h"
using namespace Rcpp;
using namespace arma;

class poEMirtdynamic_gibbs
{
public:
  
  poEMirtdynamic_gibbs(const cube &Y,
                       const cube &S,
                       const cube &Nks,
                       const cube &sb_check,
                       mat alpha,
                       mat beta,
                       mat theta,
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
                       const bool &PG_approx,
                       const int &iter,
                       const int &warmup,
                       const int &thin,
                       const bool &save_item_parameters,
                       const int &verbose);
  ~poEMirtdynamic_gibbs();
  
  void draw_Omega();
  void draw_alpha();
  void draw_alpha_fixed();
  void draw_beta();
  void draw_beta_fixed();
  void draw_theta();
  void draw_Delta();
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const cube &S;
  const cube &Nks;
  const cube &sb_check;
  
  List alpha_store;
  List beta_store;
  List theta_store;
  std::vector<double> Delta_store;
  
  mat alpha;
  mat beta;
  mat theta;
  cube Omega;
  double Delta;
  
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
  const double &g0;
  const double &h0;
  
  const vec &constraint;
  const bool &fix_alpha;
  const bool &fix_beta;
  const bool &estimate_Delta;
  const bool &std;
  const bool &PG_approx;
  const int &iter;
  const int &warmup;
  const int &thin;
  const bool &save_item_parameters;
  const int &verbose;
  
  const int I;
  const int J;
  const int K;
  const int T;
};

#endif