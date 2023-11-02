#ifndef __POEMIRTBASE_GIBBS__INCLUDED__
#define __POEMIRTBASE_GIBBS__INCLUDED__

#include <RcppArmadillo.h>
#include "dist.h"
using namespace Rcpp;
using namespace arma;


class poEMirtbase_gibbs
{
public:
  
  poEMirtbase_gibbs(const cube &Y,
                    const cube &S,
                    const cube &Nks,
                    mat alpha,
                    mat beta,
                    vec theta,
                    const std::vector<vec>& unique_categories,
                    const mat &a0,
                    const mat &A0,
                    const mat &b0,
                    const mat &B0,
                    const bool& PG_approx,
                    const int &constraint,
                    const bool &std,
                    const int &iter,
                    const int &warmup,
                    const int &thin,
                    const bool &save_item_parameters,
                    const int &verbose);
  ~poEMirtbase_gibbs();
  
  void draw_Omega();
  void draw_alpha();
  void draw_beta();
  void draw_theta();
  void fit();
  List output();
  
  List modeloutput;
  
private:
  
  const cube &Y;
  const cube &S;
  const cube &Nks;
  
  mat alpha;
  mat beta;
  vec theta;
  cube Omega;
  
  List alpha_store;
  List beta_store;
  mat theta_store;
  
  const std::vector<vec> &unique_categories;
  
  const mat &a0;
  const mat &A0;
  const mat &b0;
  const mat &B0;
  
  const bool& PG_approx;
  const int &constraint;
  const bool &std;
  const int &iter;
  const int &warmup;
  const int &thin;
  const bool &save_item_parameters;
  const int &verbose;
  
  const int I, J, K;
};

#endif