// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppSMC.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bpf_loglike_sv
double bpf_loglike_sv(arma::vec measurements, unsigned long lNumber, arma::vec starting_vals);
RcppExport SEXP _SVmodelRcppSMC_bpf_loglike_sv(SEXP measurementsSEXP, SEXP lNumberSEXP, SEXP starting_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type measurements(measurementsSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type lNumber(lNumberSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type starting_vals(starting_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(bpf_loglike_sv(measurements, lNumber, starting_vals));
    return rcpp_result_gen;
END_RCPP
}
// sv_model_pmmh_cpp
Rcpp::List sv_model_pmmh_cpp(arma::vec measurements, unsigned long lNumber, unsigned long lMCMCits, arma::vec starting_vals, arma::vec rw_mh_var, const int num_progress_outputs);
RcppExport SEXP _SVmodelRcppSMC_sv_model_pmmh_cpp(SEXP measurementsSEXP, SEXP lNumberSEXP, SEXP lMCMCitsSEXP, SEXP starting_valsSEXP, SEXP rw_mh_varSEXP, SEXP num_progress_outputsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type measurements(measurementsSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type lNumber(lNumberSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type lMCMCits(lMCMCitsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type starting_vals(starting_valsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rw_mh_var(rw_mh_varSEXP);
    Rcpp::traits::input_parameter< const int >::type num_progress_outputs(num_progress_outputsSEXP);
    rcpp_result_gen = Rcpp::wrap(sv_model_pmmh_cpp(measurements, lNumber, lMCMCits, starting_vals, rw_mh_var, num_progress_outputs));
    return rcpp_result_gen;
END_RCPP
}
// sv_model_al_tracking_impl
Rcpp::List sv_model_al_tracking_impl(arma::vec measurements, arma::vec starting_vals, unsigned long lNumber, const double resampleFreq);
RcppExport SEXP _SVmodelRcppSMC_sv_model_al_tracking_impl(SEXP measurementsSEXP, SEXP starting_valsSEXP, SEXP lNumberSEXP, SEXP resampleFreqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type measurements(measurementsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type starting_vals(starting_valsSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type lNumber(lNumberSEXP);
    Rcpp::traits::input_parameter< const double >::type resampleFreq(resampleFreqSEXP);
    rcpp_result_gen = Rcpp::wrap(sv_model_al_tracking_impl(measurements, starting_vals, lNumber, resampleFreq));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SVmodelRcppSMC_bpf_loglike_sv", (DL_FUNC) &_SVmodelRcppSMC_bpf_loglike_sv, 3},
    {"_SVmodelRcppSMC_sv_model_pmmh_cpp", (DL_FUNC) &_SVmodelRcppSMC_sv_model_pmmh_cpp, 6},
    {"_SVmodelRcppSMC_sv_model_al_tracking_impl", (DL_FUNC) &_SVmodelRcppSMC_sv_model_al_tracking_impl, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SVmodelRcppSMC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
