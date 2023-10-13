#ifndef __poEMIRTDYNAMIC_GIBBS__INCLUDED__
#define __poEMIRTDYNAMIC_GIBBS__INCLUDED__

#include <RcppArmadillo.h>
#include "helper.h"
#include "pgdraw.h"
#include "dist.h"
using namespace Rcpp;
using namespace arma;

class poEMirtdynamic_gibbs
{
public:
  
  poEMirtdynamic_gibbs(const cube &Y,
                       const mat &N,
                       mat alpha,
                       mat beta,
                       mat theta,
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
                       const bool &PG_approx,
                       const vec &constraint,
                       const bool &std,
                       const int &iter,
                       const int &warmup,
                       const int &thin,
                       const bool &save_item_parameters,
                       const int &verbose);
  ~poEMirtdynamic_gibbs();
  
  void draw_Omega();
  void draw_alpha();
  void draw_beta();
  void draw_theta();
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const mat &N;
  
  List alpha_store;
  List beta_store;
  List theta_store;
  
  mat alpha;
  mat beta;
  mat theta;
  cube Omega;
  
  const std::vector<vec> &unique_categories;
  const mat &timemap;
  const mat &timemap2;
  const vec &item_timemap;
  const vec &item_match;
  List time_info;
  std::vector<uvec> time_list;
  std::vector<uvec> item_time_list;
  cube sb_check;
  cube S;
  cube Nks;
  
  const mat &a0;
  const mat &A0;
  const mat &b0;
  const mat &B0;
  const vec &m0;
  const vec &C0;
  const vec &Delta;
  
  const bool &PG_approx;
  const vec &constraint;
  const bool &std;
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