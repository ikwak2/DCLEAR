// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dist_w_ham
arma::mat dist_w_ham(const arma::mat X, const arma::vec InfoW);
RcppExport SEXP _DCLEAR_dist_w_ham(SEXP XSEXP, SEXP InfoWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type InfoW(InfoWSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_w_ham(X, InfoW));
    return rcpp_result_gen;
END_RCPP
}
// dist_w_ham2
arma::mat dist_w_ham2(const arma::mat X, const arma::vec InfoW);
RcppExport SEXP _DCLEAR_dist_w_ham2(SEXP XSEXP, SEXP InfoWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type InfoW(InfoWSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_w_ham2(X, InfoW));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DCLEAR_dist_w_ham", (DL_FUNC) &_DCLEAR_dist_w_ham, 2},
    {"_DCLEAR_dist_w_ham2", (DL_FUNC) &_DCLEAR_dist_w_ham2, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DCLEAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}