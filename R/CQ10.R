#' Simulation data for CQ10
#'
#' TEST
#'
#' @param n nobs of X1 and X2
#' @param p number of variables
#' @param dist gamma or norm
#' @param dependence two or full
#' @param alloc "eq", "inc", "dec"
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
                    dependence=c("two", "full"),
                    alloc=c("eq", "inc", "dec")) {
  dependence <- match.arg(dependence)
  dist <- match.arg(dist)
  alloc <- match.arg(alloc)

  if (dependence == "two") {
    rho <- c(2.883, 2.794, 2.849)
    nAux <- 3 - 1
  } else {
    rho <- runif(p, min=2, max=3)
    nAux <- p - 1
  }

  # no need to calculate Sigma, all are diagonal with sum(rho * rho)
  # Sigma <- diag(var_x, nrow = p)
  # trSigma2 <- sum(diag(Sigma %*% Sigma))
  trSigma2 <- p * sum(rho * rho)

  # n of false null
  n_F_null <- p * (100-perc)/100
  # n of sample in the divided 4 groups
  n_sample <- n_F_null %/% 4

  if (n_F_null %% 4 == 0) sampleV <- rep(n_sample, 4)
  else if (alloc == "inc") sampleV <-  c(rep(n_sample, 3), n_sample+1)
  else sampleV <-  c(n_sample+1, rep(n_sample, 3))

  # assume mu1 = 0, individual mu2i considering n_F_null
  mu2i <- sqrt(sqrt(trSigma2) / n_F_null)
  mu2pre <- vector("numeric", length = p)

  assignMu <- function(cutp) {
    intv <- findInterval(1:p, cutp)
    foreach (i=0:3, sam=sampleV, .combine = c) %do% {
      L <- mu2pre[intv == i]
      L[sample(length(L), sam)] <- mu2i
      L
    }
  }

  if (alloc == "eq") {
    # find 3 cutpoints to divide 4 group
    cutpoints <- floor(p * c(0.25, 0.5, 0.75)) + 1
    mu2 <- assignMu(cutpoints)

  } else if (alloc == "inc") {
    cutpoints <- floor(p * c(0.4, 0.7, 0.9)) + 1
    mu2 <- assignMu(cutpoints)

  } else {
    cutpoints <- floor(p * c(0.1, 0.3, 0.6)) + 1
    mu2 <- assignMu(cutpoints)
  }

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
    iV <- 1:length(rhoV)
    foreach (i=iV, rho=rhoV, .combine="+") %do% {
      rho * Z[, i:(p+i-1)]
    }
  }

  X1 <- genX(rho, Z1)
  X2 <- genX(rho, Z2) + mu2

  invisible(list(X1=X1, X2=X2))
}

