// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// poEMirtbase_fit
List poEMirtbase_fit(const arma::cube& Y, const arma::cube& S, const arma::cube& Nks, arma::mat alpha_old, arma::mat beta_old, arma::vec theta_old, const std::vector<arma::vec>& unique_categories, const arma::mat& a0, const arma::mat& A0, const arma::mat& b0, const arma::mat& B0, const int& constraint, const bool& std, const int& maxit, const int& verbose, const double& tol, const bool& compute_ll);
RcppExport SEXP _poEMirt_poEMirtbase_fit(SEXP YSEXP, SEXP SSEXP, SEXP NksSEXP, SEXP alpha_oldSEXP, SEXP beta_oldSEXP, SEXP theta_oldSEXP, SEXP unique_categoriesSEXP, SEXP a0SEXP, SEXP A0SEXP, SEXP b0SEXP, SEXP B0SEXP, SEXP constraintSEXP, SEXP stdSEXP, SEXP maxitSEXP, SEXP verboseSEXP, SEXP tolSEXP, SEXP compute_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Nks(NksSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_old(alpha_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_old(theta_oldSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type unique_categories(unique_categoriesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const int& >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< const bool& >::type std(stdSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const bool& >::type compute_ll(compute_llSEXP);
    rcpp_result_gen = Rcpp::wrap(poEMirtbase_fit(Y, S, Nks, alpha_old, beta_old, theta_old, unique_categories, a0, A0, b0, B0, constraint, std, maxit, verbose, tol, compute_ll));
    return rcpp_result_gen;
END_RCPP
}
// poEMirtbase_gibbs_fit
List poEMirtbase_gibbs_fit(const arma::cube& Y, const arma::cube& S, const arma::cube& Nks, arma::mat alpha, arma::mat beta, arma::vec theta, const std::vector<arma::vec>& unique_categories, const arma::mat& a0, const arma::mat& A0, const arma::mat& b0, const arma::mat& B0, const bool& PG_approx, const int& constraint, const bool& std, const int& iter, const int& warmup, const int& thin, const bool& save_item_parameters, const int& verbose);
RcppExport SEXP _poEMirt_poEMirtbase_gibbs_fit(SEXP YSEXP, SEXP SSEXP, SEXP NksSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP thetaSEXP, SEXP unique_categoriesSEXP, SEXP a0SEXP, SEXP A0SEXP, SEXP b0SEXP, SEXP B0SEXP, SEXP PG_approxSEXP, SEXP constraintSEXP, SEXP stdSEXP, SEXP iterSEXP, SEXP warmupSEXP, SEXP thinSEXP, SEXP save_item_parametersSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Nks(NksSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type unique_categories(unique_categoriesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const bool& >::type PG_approx(PG_approxSEXP);
    Rcpp::traits::input_parameter< const int& >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< const bool& >::type std(stdSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool& >::type save_item_parameters(save_item_parametersSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(poEMirtbase_gibbs_fit(Y, S, Nks, alpha, beta, theta, unique_categories, a0, A0, b0, B0, PG_approx, constraint, std, iter, warmup, thin, save_item_parameters, verbose));
    return rcpp_result_gen;
END_RCPP
}
// poEMirtdynamic_fit
List poEMirtdynamic_fit(const arma::cube& Y, const arma::cube& S, const arma::cube& Nks, arma::mat alpha_old, arma::mat beta_old, arma::mat theta_old, const std::vector<arma::vec>& unique_categories, const std::vector<arma::vec>& uJ_J, const arma::mat& timemap2, const arma::vec& item_timemap, const std::vector<arma::uvec>& IT, const std::vector<std::vector<arma::uvec>>& ITJ, const arma::mat& a0, const arma::mat& A0, const arma::mat& b0, const arma::mat& B0, const arma::vec& m0, const arma::vec& C0, const arma::vec& Delta, const arma::vec& constraint, const bool& alpha_fix, const bool& std, const int& maxit, const int& verbose, const double& tol, const bool& compute_ll);
RcppExport SEXP _poEMirt_poEMirtdynamic_fit(SEXP YSEXP, SEXP SSEXP, SEXP NksSEXP, SEXP alpha_oldSEXP, SEXP beta_oldSEXP, SEXP theta_oldSEXP, SEXP unique_categoriesSEXP, SEXP uJ_JSEXP, SEXP timemap2SEXP, SEXP item_timemapSEXP, SEXP ITSEXP, SEXP ITJSEXP, SEXP a0SEXP, SEXP A0SEXP, SEXP b0SEXP, SEXP B0SEXP, SEXP m0SEXP, SEXP C0SEXP, SEXP DeltaSEXP, SEXP constraintSEXP, SEXP alpha_fixSEXP, SEXP stdSEXP, SEXP maxitSEXP, SEXP verboseSEXP, SEXP tolSEXP, SEXP compute_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Nks(NksSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_old(alpha_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta_old(theta_oldSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type unique_categories(unique_categoriesSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type uJ_J(uJ_JSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type timemap2(timemap2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type item_timemap(item_timemapSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type IT(ITSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<arma::uvec>>& >::type ITJ(ITJSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< const bool& >::type alpha_fix(alpha_fixSEXP);
    Rcpp::traits::input_parameter< const bool& >::type std(stdSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const bool& >::type compute_ll(compute_llSEXP);
    rcpp_result_gen = Rcpp::wrap(poEMirtdynamic_fit(Y, S, Nks, alpha_old, beta_old, theta_old, unique_categories, uJ_J, timemap2, item_timemap, IT, ITJ, a0, A0, b0, B0, m0, C0, Delta, constraint, alpha_fix, std, maxit, verbose, tol, compute_ll));
    return rcpp_result_gen;
END_RCPP
}
// poEMirtdynamic_gibbs_fit
List poEMirtdynamic_gibbs_fit(const arma::cube& Y, const arma::cube& S, const arma::cube& Nks, const arma::cube& sb_check, arma::mat alpha, arma::mat beta, arma::mat theta, const std::vector<arma::vec>& unique_categories, const std::vector<arma::vec>& uJ_J, const arma::mat& timemap2, const arma::vec& item_timemap, const std::vector<arma::uvec>& IT, const std::vector<std::vector<arma::uvec>>& ITJ, const arma::mat& a0, const arma::mat& A0, const arma::mat& b0, const arma::mat& B0, const arma::vec& m0, const arma::vec& C0, const arma::vec& Delta, const bool& alpha_fix, const bool& PG_approx, const arma::vec& constraint, const bool& std, const int& iter, const int& warmup, const int& thin, const bool& save_item_parameters, const int& verbose);
RcppExport SEXP _poEMirt_poEMirtdynamic_gibbs_fit(SEXP YSEXP, SEXP SSEXP, SEXP NksSEXP, SEXP sb_checkSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP thetaSEXP, SEXP unique_categoriesSEXP, SEXP uJ_JSEXP, SEXP timemap2SEXP, SEXP item_timemapSEXP, SEXP ITSEXP, SEXP ITJSEXP, SEXP a0SEXP, SEXP A0SEXP, SEXP b0SEXP, SEXP B0SEXP, SEXP m0SEXP, SEXP C0SEXP, SEXP DeltaSEXP, SEXP alpha_fixSEXP, SEXP PG_approxSEXP, SEXP constraintSEXP, SEXP stdSEXP, SEXP iterSEXP, SEXP warmupSEXP, SEXP thinSEXP, SEXP save_item_parametersSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Nks(NksSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type sb_check(sb_checkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type unique_categories(unique_categoriesSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type uJ_J(uJ_JSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type timemap2(timemap2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type item_timemap(item_timemapSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type IT(ITSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<arma::uvec>>& >::type ITJ(ITJSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type alpha_fix(alpha_fixSEXP);
    Rcpp::traits::input_parameter< const bool& >::type PG_approx(PG_approxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< const bool& >::type std(stdSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool& >::type save_item_parameters(save_item_parametersSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(poEMirtdynamic_gibbs_fit(Y, S, Nks, sb_check, alpha, beta, theta, unique_categories, uJ_J, timemap2, item_timemap, IT, ITJ, a0, A0, b0, B0, m0, C0, Delta, alpha_fix, PG_approx, constraint, std, iter, warmup, thin, save_item_parameters, verbose));
    return rcpp_result_gen;
END_RCPP
}
// construct_sb_auxs
List construct_sb_auxs(const arma::cube& Y, const arma::mat& N, const std::vector<arma::vec>& unique_categories);
RcppExport SEXP _poEMirt_construct_sb_auxs(SEXP YSEXP, SEXP NSEXP, SEXP unique_categoriesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type unique_categories(unique_categoriesSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_sb_auxs(Y, N, unique_categories));
    return rcpp_result_gen;
END_RCPP
}
// get_dynamic_info
List get_dynamic_info(const arma::mat& N, const arma::vec& item_timemap, const arma::mat& timemap);
RcppExport SEXP _poEMirt_get_dynamic_info(SEXP NSEXP, SEXP item_timemapSEXP, SEXP timemapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type item_timemap(item_timemapSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type timemap(timemapSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dynamic_info(N, item_timemap, timemap));
    return rcpp_result_gen;
END_RCPP
}
// prediction
arma::cube prediction(const arma::mat& N, const arma::mat& alpha, const arma::mat& beta, const arma::mat& theta, const std::vector<arma::vec> unique_categories, const arma::vec& item_timemap, const String& model, const String& type);
RcppExport SEXP _poEMirt_prediction(SEXP NSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP thetaSEXP, SEXP unique_categoriesSEXP, SEXP item_timemapSEXP, SEXP modelSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec> >::type unique_categories(unique_categoriesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type item_timemap(item_timemapSEXP);
    Rcpp::traits::input_parameter< const String& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const String& >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(prediction(N, alpha, beta, theta, unique_categories, item_timemap, model, type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_poEMirt_poEMirtbase_fit", (DL_FUNC) &_poEMirt_poEMirtbase_fit, 17},
    {"_poEMirt_poEMirtbase_gibbs_fit", (DL_FUNC) &_poEMirt_poEMirtbase_gibbs_fit, 19},
    {"_poEMirt_poEMirtdynamic_fit", (DL_FUNC) &_poEMirt_poEMirtdynamic_fit, 26},
    {"_poEMirt_poEMirtdynamic_gibbs_fit", (DL_FUNC) &_poEMirt_poEMirtdynamic_gibbs_fit, 29},
    {"_poEMirt_construct_sb_auxs", (DL_FUNC) &_poEMirt_construct_sb_auxs, 3},
    {"_poEMirt_get_dynamic_info", (DL_FUNC) &_poEMirt_get_dynamic_info, 3},
    {"_poEMirt_prediction", (DL_FUNC) &_poEMirt_prediction, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_poEMirt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
