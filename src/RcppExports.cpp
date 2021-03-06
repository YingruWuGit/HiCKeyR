// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// segment
void segment(std::string argv);
RcppExport SEXP _HiCKeyR_segment(SEXP argvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type argv(argvSEXP);
    segment(argv);
    return R_NilValue;
END_RCPP
}
// segHeatMap
NumericMatrix segHeatMap(std::string argv, int s, int e);
RcppExport SEXP _HiCKeyR_segHeatMap(SEXP argvSEXP, SEXP sSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type argv(argvSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(segHeatMap(argv, s, e));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HiCKeyR_segment", (DL_FUNC) &_HiCKeyR_segment, 1},
    {"_HiCKeyR_segHeatMap", (DL_FUNC) &_HiCKeyR_segHeatMap, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_HiCKeyR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
