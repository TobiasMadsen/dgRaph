// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// evalDFG
double evalDFG(NumericVector x);
RcppExport SEXP PGMscore_evalDFG(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        double __result = evalDFG(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mgfDFG
double mgfDFG(IntegerVector varDimensions, List facPotentials, List facNeighbors);
RcppExport SEXP PGMscore_mgfDFG(SEXP varDimensionsSEXP, SEXP facPotentialsSEXP, SEXP facNeighborsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type varDimensions(varDimensionsSEXP );
        Rcpp::traits::input_parameter< List >::type facPotentials(facPotentialsSEXP );
        Rcpp::traits::input_parameter< List >::type facNeighbors(facNeighborsSEXP );
        double __result = mgfDFG(varDimensions, facPotentials, facNeighbors);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// PGMExpectCpp
NumericVector PGMExpectCpp(IntegerVector varDimensions, List facPotentials1, List facPotentials2, List facNeighbors);
RcppExport SEXP PGMscore_PGMExpectCpp(SEXP varDimensionsSEXP, SEXP facPotentials1SEXP, SEXP facPotentials2SEXP, SEXP facNeighborsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type varDimensions(varDimensionsSEXP );
        Rcpp::traits::input_parameter< List >::type facPotentials1(facPotentials1SEXP );
        Rcpp::traits::input_parameter< List >::type facPotentials2(facPotentials2SEXP );
        Rcpp::traits::input_parameter< List >::type facNeighbors(facNeighborsSEXP );
        NumericVector __result = PGMExpectCpp(varDimensions, facPotentials1, facPotentials2, facNeighbors);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sampleISCpp
DataFrame sampleISCpp(int N, double alpha, IntegerVector varDimensions, List facPotentialsNull, List facNeighbors, List facPotentialsFg);
RcppExport SEXP PGMscore_sampleISCpp(SEXP NSEXP, SEXP alphaSEXP, SEXP varDimensionsSEXP, SEXP facPotentialsNullSEXP, SEXP facNeighborsSEXP, SEXP facPotentialsFgSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type N(NSEXP );
        Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type varDimensions(varDimensionsSEXP );
        Rcpp::traits::input_parameter< List >::type facPotentialsNull(facPotentialsNullSEXP );
        Rcpp::traits::input_parameter< List >::type facNeighbors(facNeighborsSEXP );
        Rcpp::traits::input_parameter< List >::type facPotentialsFg(facPotentialsFgSEXP );
        DataFrame __result = sampleISCpp(N, alpha, varDimensions, facPotentialsNull, facNeighbors, facPotentialsFg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
