// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// generate_frequency_table
Rcpp::DataFrame generate_frequency_table(std::string head_seq, std::string seq_id, int seq_size, int ws, float cut_off, Rcpp::StringVector motifs, int k);
RcppExport SEXP _DNAmotif_generate_frequency_table(SEXP head_seqSEXP, SEXP seq_idSEXP, SEXP seq_sizeSEXP, SEXP wsSEXP, SEXP cut_offSEXP, SEXP motifsSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type head_seq(head_seqSEXP);
    Rcpp::traits::input_parameter< std::string >::type seq_id(seq_idSEXP);
    Rcpp::traits::input_parameter< int >::type seq_size(seq_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type ws(wsSEXP);
    Rcpp::traits::input_parameter< float >::type cut_off(cut_offSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type motifs(motifsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_frequency_table(head_seq, seq_id, seq_size, ws, cut_off, motifs, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DNAmotif_generate_frequency_table", (DL_FUNC) &_DNAmotif_generate_frequency_table, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_DNAmotif(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
