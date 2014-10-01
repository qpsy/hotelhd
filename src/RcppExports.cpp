// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// caltrSigmaSqC
double caltrSigmaSqC(NumericMatrix X_in, int n, int p);
RcppExport SEXP hotelhd_caltrSigmaSqC(SEXP X_inSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X_in(X_inSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        double __result = caltrSigmaSqC(X_in, n, p);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// caltrSigma12C
double caltrSigma12C(NumericMatrix X1_in, NumericMatrix X2_in, int n1, int n2, int p);
RcppExport SEXP hotelhd_caltrSigma12C(SEXP X1_inSEXP, SEXP X2_inSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP pSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X1_in(X1_inSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type X2_in(X2_inSEXP );
        Rcpp::traits::input_parameter< int >::type n1(n1SEXP );
        Rcpp::traits::input_parameter< int >::type n2(n2SEXP );
        Rcpp::traits::input_parameter< int >::type p(pSEXP );
        double __result = caltrSigma12C(X1_in, X2_in, n1, n2, p);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
