#' Simulation data for BS96
#'
#' TEST
#'
#' @param n1 nobs of X1
#' @param n2 nobs of X2
#' @param p number of variables
#' @param rho rho
#' @param eta eta
#' @param dist gamma or norm
#'
#' @author Dongjun You
#'
#' @return Test results
#'#'
#' @references
#' will be added
#'
#' @examples
#' \dontrun{
#' res <- genBS96(n1=25, n2=20, p=40, rho=0.5, eta=0, dist="norm")
#' res[["phi"]]
#' res[["lambda"]]
#' d1 <- res[["d1"]]
#' d2 <- res[["d2"]]
#' hotelhd(d1, d2, method="H")
#' hotelhd(d1, d2, method="D")
#' hotelhd(d1, d2, method="BS")
#' hotelhd(d1, d2, method="CQ")
#' }
#'
#' @export
genBS96 <- function(n1, n2, p, rho, eta, dist=c("gamma", "norm")) {
  gen_U <- function(n) {
    if (dist == "gamma") {
      matrix(rgamma(n*(p+1), shape=4, scale=1), ncol=p+1)
    } else {
      matrix(rnorm(n*(p+1), mean=0, sd=1), ncol=p+1)
    }
  }

  dist <- match.arg(dist)

  if (dist == "gamma") {
    Sigma <- diag(rep(4*(1 + rho*rho), p))
    Sigma[(row(Sigma) + 1) == col(Sigma)] <- 4*rho
    Sigma[(row(Sigma) - 1) == col(Sigma)] <- 4*rho
  } else {
    Sigma <- (1 - rho) * diag(nrow=p, ncol=p) +
        rho * matrix(1, nrow=p, ncol=p)
  }

  lambda <- eigen(Sigma)$values
  phi <- lambda[1] / lambda[p]

  # assume mu2 = 0
  trSigma2 <- sum(diag(Sigma %*% Sigma))
  # assume the p elements of mu1 have the same value mu1i
  mu1i <- sqrt(eta * sqrt(trSigma2) / p)

  U1 <- gen_U(n1)
  d1 <- U1[, 1:p] + rho * U1[, 2:(p+1)] + mu1i
  U2 <- gen_U(n2)
  d2 <- U2[, 1:p] + rho * U2[, 2:(p+1)]

  invisible(list(d1=d1, d2=d2, phi=phi, lambda=lambda))
}
