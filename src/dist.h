#ifndef __DIST_H__
#define __DIST_H__

#include <Rcpp.h>
using namespace Rcpp;

double rpg(const double& b,
           const double& c);

double rpg_approx(const double& b,
                  const double& c);

#endif