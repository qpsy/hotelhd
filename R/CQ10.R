#' Simulation data for CQ10
#'
#' TEST
#'
#' @param n nobs of X1 and X2
#' @param p number of variables
#' @param perc percentages
#' @param dist gamma or norm
#' @param dependence two or full
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
#' res <- genCQ10(n=124, p=500, perc=50)
#' X1 <- res[["X1"]]
#' X2 <- res[["X2"]]
#' hotelhd(X1, X2, method="H")
#' hotelhd(X1, X2, method="D")
#' hotelhd(X1, X2, method="BS")
#' hotelhd(X1, X2, method="CQ")
#' }
#'
#' @export
genCQ10 <- function(n, p, perc, dist=c("gamma", "norm"),
                    dependence=c("two", "full")) {
  dependence <- match.arg(dependence)
  dist <- match.arg(dist)

  if (dependence == "two") {
    rho <- c(2.883, 2.794, 2.849)
    nAux <- 3 - 1
  } else {
    rho <- runif(p, min=2, max=3)
    nAux <- p - 1
  }

  # true Sigma, considering moving average model with independent variables
  #  : to be band diagonal matrix
  i <- NULL
  Mdim <- length(rho)
  Sigma <- foreach (i=1:(Mdim-1), .combine="+",
           .final=function(res) {
             t(res) + res + diag2(sum(rho*rho), p, p)
             } ) %do% {
    bandVal <- t(rho) %*% diag2(1, Mdim, Mdim, i) %*% rho
    diag2(bandVal, p, p, i)
  }

  trSigma2 <- sum(diag(Sigma %*% Sigma))

  # n of false null
  n_F_null <- p * (100-perc)/100

  # assume mu1 = 0, individual mu2i considering n_F_null
  mu2i <- sqrt(sqrt(trSigma2) / n_F_null)
  mu2pre <- vector("numeric", length = p)

  # ignore 'allocation' options
  # simply sample the n of false null
  if (perc != 100) mu2pre[sample(p, n_F_null)] <- mu2i
  mu2 <- rep(1, times=n) %*% t(mu2pre)

  gen_Z <- function(n, nAux) {
    if (dist == "gamma") {
      matrix(rgamma(n*(p + nAux), shape=4, scale=1), ncol=p + nAux)
    } else {
      matrix(rnorm(n*(p + nAux), mean=0, sd=1), ncol=p + nAux)
    }
  }

  Z1 <- gen_Z(n, nAux)
  Z2 <- gen_Z(n, nAux)

  genX <- function(rhoV, Z) {
    i <- NULL
    iV <- 1:length(rhoV)
    foreach (i=iV, rho=rhoV, .combine="+") %do% {
      rho * Z[, i:(p+i-1)]
    }
  }

  X1 <- genX(rho, Z1)
  X2 <- genX(rho, Z2) + mu2

  invisible(list(X1=X1, X2=X2))
}
