#ifndef __HELPER_H__
#define __HELPER_H__

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

vec rmultinom2(int &size, 
               vec &prob);

List get_time_info (const mat &N,
                    const vec &item_timemap,
                    const mat &timemap);

List construct_sb_auxs(const cube &Y,
                       const mat &N,
                       const std::vector<vec> &unique_categories);

#endif
