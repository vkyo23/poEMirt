#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

vec rmultinom2(int &size, 
               vec &prob) 
{
  IntegerVector outcome(prob.size());
  rmultinom(size, prob.begin(), prob.size(), outcome.begin());
  return as<vec>(wrap((outcome)));
}

List get_time_info (const mat &N,
                    const vec &item_timemap,
                    const mat &timemap) 
{
  const int& I = timemap.n_rows;
  const int& J = item_timemap.size();
  const int& T = timemap.n_cols;
  std::vector<std::vector<uvec>> ind_item_time_list;
  std::vector<uvec> ind_time_list;
  for (int i = 0; i < I; i++) {
    int j = 0;
    std::vector<int> ind_t_index;
    std::vector<uvec> item_it_list;
    for (int t = 0; t < T; t++) {
      if (timemap(i, t) == 1) {
        ind_t_index.push_back(t);
      }
      std::vector<int> item_it;
      for (; j < J; j++) {
        if (item_timemap[j] == t) {
          if (N(i, j) != 0) {
            item_it.push_back(j);
          }
          if (item_timemap[j+1] == t+1) {
            break;
          }
        }
      }
      if (timemap(i, t) == 1) {
        item_it_list.push_back(as<uvec>(wrap(item_it)));
      }
    }
    ind_item_time_list.push_back(item_it_list);
    ind_time_list.push_back(as<uvec>(wrap(ind_t_index)));
  }
  
  List L = List::create(Named("ind_item_time_list") = ind_item_time_list,
                        Named("time_list") = ind_time_list);
  
  return L;
}

// To do: Constructing auxiliary variables for SB calculation
List construct_sb_auxs(const cube &Y,
                       const mat &N,
                       const std::vector<vec> &unique_categories) 
{
  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = Y.n_slices - 1;
  
  cube sb_check(I, J, K);
  cube Nks(I, J, K);
  cube S(I, J, K);
  for (int i = 0; i < I; i++) {
    for (int j = 0; j < J; j++) {
      vec unq = unique_categories[j]-1;
      double Ycum = 0;
      for (int k = 0; k < (unq.size() - 1); k++) {
        if (!NumericVector::is_na(Y(i, j, unq[k]))) {
          Nks(i, j, unq[k]) = N(i, j) - Ycum;
          Ycum += Y(i, j, unq[k]);
          if (Nks(i, j, unq[k]) > 0) {
            sb_check(i, j, unq[k]) = 1;
            S(i, j, unq[k]) = Y(i, j, unq[k]) - Nks(i, j, unq[k]) / 2.0;
          }
        }
      }
    }
  }
  
  List L = List::create(
    Named("sb_check") = sb_check,
    Named("Nks") = Nks,
    Named("S") = S
  );
  return L;
}
