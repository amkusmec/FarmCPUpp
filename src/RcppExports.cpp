// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// QuickLM
Rcpp::List QuickLM(Rcpp::NumericVector ys, Rcpp::NumericMatrix Xs, SEXP pBigMat, IntegerVector observations, int npcs, int nqtn);
RcppExport SEXP _FarmCPUpp_QuickLM(SEXP ysSEXP, SEXP XsSEXP, SEXP pBigMatSEXP, SEXP observationsSEXP, SEXP npcsSEXP, SEXP nqtnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type npcs(npcsSEXP);
    Rcpp::traits::input_parameter< int >::type nqtn(nqtnSEXP);
    rcpp_result_gen = Rcpp::wrap(QuickLM(ys, Xs, pBigMat, observations, npcs, nqtn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FarmCPUpp_QuickLM", (DL_FUNC) &_FarmCPUpp_QuickLM, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_FarmCPUpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
