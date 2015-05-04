#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double caltrSigmaSqC(NumericMatrix X_in, int n, int p) {
	mat X(X_in.begin(), n, p, false), Xbar_jk(1, p), SigmaSq(p, p);
	rowvec Xj(p), Xk(p);
	IntegerVector v = seq_len(n) - 1;
	IntegerVector idxR(2);
	uvec idx(2);
 	int j, k;

	SigmaSq.zeros();

	for (j = 0; j < n-1; ++j) {
    Xj = X.row(j);
		//Rcout << "j " << j << std::endl;
	  idxR(0) = j;

    for (k = j+1; k < n; ++k) {
      Xk = X.row(k);
			//Rcout << " k " << k << std::endl;
			//Rcout << "j " << j+1 << " " << "k " << k+1 << std::endl;
			idxR(1) = k;
			idx = as<uvec>(setdiff(v, idxR));
			Xbar_jk = mean(X.rows(idx));
			SigmaSq += (Xj - Xbar_jk).t() * Xj * (Xk - Xbar_jk).t() * Xk;
			//Rcout << SigmaSq  << std::endl;
    }
  }

	//Rcout << 2 * sum(SigmaSq.diag())  << std::endl;
	double trSigmaSq =  (2 * trace(SigmaSq)) / (n*(n-1));
	//return List::create(Named("Sig") = trSigmaSq,
	//										Named("k") = Xk);
	return trSigmaSq;
}


// [[Rcpp::export]]
double caltrSigma12C(NumericMatrix X1_in, NumericMatrix X2_in,
										 int n1, int n2, int p) {
	mat X1(X1_in.begin(), n1, p, false), X2(X2_in.begin(), n2, p, false),
		X1bar_j(1, p), X2bar_k(1, p), Sigma12(p, p);
	rowvec X1j(p), X2k(p);
	IntegerVector v1 = seq_len(n1) - 1, v2 = seq_len(n2) - 1,
		idxR1(1), idxR2(1);
	uvec idx1(1), idx2(1);
 	int j, k;

	Sigma12.zeros();

	for (j = 0; j < n1; ++j) {
    X1j = X1.row(j);
	  idxR1(0) = j;
		idx1 = as<uvec>(setdiff(v1, idxR1));
		X1bar_j = mean(X1.rows(idx1));

    for (k = 0; k < n2; ++k) {
      X2k = X2.row(k);
			idxR2(0) = k;
			idx2 = as<uvec>(setdiff(v2, idxR2));
			X2bar_k = mean(X2.rows(idx2));
			Sigma12 += (X1j - X1bar_j).t() * X1j * (X2k - X2bar_k).t() * X2k;
    }
  }

	double trSigma12 =  trace(Sigma12) / (n1*n2);
	return trSigma12;
}

// [[Rcpp::export]]
double calcM2C(NumericVector Z, NumericMatrix Omega, int ndim) {
  NumericVector M = Z*Z / diag(Omega);
  std::sort(M.begin(), M.end(), std::greater<int>());
  return sum(M[seq(0, ndim-1)]);
}
