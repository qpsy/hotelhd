#' Adaptive Thresholding for Sparse Covariance Matrix Estimation
#'
#' test
#'
#' @param X an n by p data matrix.
#' @param delta a numeric value. A tuning parameter.
#'   The default value of delta is 2, and is taken from Cai and Liu (2011).
#'   Also, it can be chosen empirically through cross-validation.
#'
#' @return Test results
#'
#' @references
#' will be added
#'
#' @export
adacov <- function(X, delta=2) {
  ## the adaptive lasso thresholding
  ##  : a type of thresholing functions being used in Cai and Liu (2011)
  s_lambda <- function(S, lambda, eta=4) {
    S * (1 - (abs(lambda / S))^eta)
    # S * (1 - (abs(S / lambda))^eta)
  }

  n <- NROW(X)
  p <- NCOL(X)
  S <- var(X)

  ## theta: after a little algebra
  ss <- (sweep(X, 2, colMeans(X)))^2 # n by p
  theta <- (t(ss) %*% ss + (2-n)*(var(X))^2) / n # p by p

  lambda <- delta * sqrt(log(p)*theta / n) # p by p

  s_lambda(S, lambda)
}
