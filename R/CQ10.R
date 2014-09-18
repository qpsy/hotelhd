#' Simulation data for CQ10
#'
#' TEST
#'
#' @param n1 nobs of X1
#' @param n2 nobs of X2
#' @param p number of variables
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
#' res <- genBS96(n1=25, n2=20, p=40, dist="norm")
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
genCQ10 <- function(n, p, perc, dist=c("gamma", "norm"),
                    dependence=c("two", "full"),
                    alloc=c("eq", "inc", "dec")) {
  dependence <- match.arg(dependence)
  dist <- match.arg(dist)
  alloc <- match.arg(alloc)

  gen_Z <- function(n, nAux) {
    if (dist == "gamma") {
      matrix(rgamma(n*(p + nAux-1), shape=4, scale=1), ncol=p + nAuxVar-1)
    } else {
      matrix(rnorm(n*(p + nAux-1), mean=0, sd=1), ncol=p + nAux-1)
    }
  }

  ## TODO: set up Sigma

  ## if (dist == "gamma") {
  ##   Sigma <- diag(rep(4*(1 + rho*rho), p))
  ##   Sigma[(row(Sigma) + 1) == col(Sigma)] <- 4*rho
  ##   Sigma[(row(Sigma) - 1) == col(Sigma)] <- 4*rho
  ## } else {
  ##   Sigma <- (1 - rho) * diag(nrow=p, ncol=p) +
  ##       rho * matrix(1, nrow=p, ncol=p)
  ## }

  genX <- function(rhoV, Z) {
    iV <- 1:length(rhoV)
    foreach (i=iV, rho=rhoV, .combine="+") %do% {
      rho * Z[, i:(p+i-1)]
    }
  }

  if (dependence == "two") {
    rho <- c(2.883, 2.794, 2.849)
    nAux <- 3
  } else {
    rho <- runif(p, min=2, max=3)
    nAux <- p
  }

  Z1 <- gen_Z(n, nAux)
  Z2 <- gen_Z(n, nAux)

  X1 <- genX(rho, Z1)
  X2 <- genX(rho, Z2)

  # assume mu1 = 0
  trSigma2 <- sum(diag(Sigma %*% Sigma))

  if (alloc == "eq") {
    ## TODO: consider percentage
    # assume the p elements of mu2 have the same value mu1i
    # mu1i <- sqrt(eta * sqrt(trSigma2) / p)

    X2 <- X2 + mu2
  } else if (alloc == "inc") {

    X2 <- X2 + mu2
  } else {

    X2 <- X2 + mu2
  }

  invisible(list(d1=d1, d2=d2, phi=phi, lambda=lambda))
}
