#' Generate simulation data for CLX14
#'
#' TEST
#'
#' @param n nobs of X1 and X2
#' @param p number of variables
#' @param m a character value. Two values of the number of non-zero mu1.
#'   For 'm1', floor(0.05p) of mu1 will be set to non-zeros.
#'   floor(sqrt(p)) for 'm2'. 0 for 'power'.
#' @param muval a character value. The magnitude of the elements of mu1.
#'   If 'muval1', pm sqrt(log(p)/n). 'muval2' for
#'   rnuif(-sqrt(8*log(p)/n), sqrt(8*log(p)/n)).
#' @param model a character value. A model generating covariance matrix.
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
#' res <- genCLX14(n=50, p=50)
#' X1 <- res[["X1"]]
#' X2 <- res[["X2"]]
#' hotelhd(X1, X2, method="H")
#' hotelhd(X1, X2, method="D")
#' hotelhd(X1, X2, method="BS")
#' hotelhd(X1, X2, method="CQ")
#' }
#'
#' @importFrom magic adiag
#' @importFrom MASS mvrnorm
#'
#' @export
genCLX14 <- function(n, p, m=c("m1", "m2", "power"),
                     muval=c("muval1", "muval2"),
                     model=c("model1", "model2", "model3")) {
  m <- match.arg(m)
  if (m == "power") {
    mu1 <- rep(0, p)

  } else {
    ## mu2==0 assumed
    muval <- match.arg(muval)
    logpn <- log(p) / n
    if (muval == "muval1") {
      mu1 <- (2 * rbinom(p, 1, .5) - 1) * # random p of plus, minus signs
          sqrt(logpn)
    } else {
      mu1 <- runif(p, min=-sqrt(8*logpn), max=sqrt(8*logpn))
    }

    if (m == "m1") m <- floor(0.05 * p)
    else m <- floor(sqrt(p))

    mu1[sample(p, p-m)] <- 0
  }

  model <- match.arg(model)
  if (model == "model1") {
    k <- (p %/% 2)
    elem <- matrix(c(1, 0.8, 0.8, 1), ncol=2)
    Sigma <- do.call(adiag, lapply(1:k, function(dummy) elem))

  } else if (model == "model2") {
    Sigma <- matrix(0, nrow=p, ncol=p)
    for (off in 1:(p-1)) {
      bandVal <- 0.6^off
      Sigma <- Sigma + diag2(bandVal, p, p, off)
    }
    Sigma <- Sigma + t(Sigma) + diag(nrow=p, ncol=p)

  } else if (model == "model3") {
    Omega <- matrix(0, nrow=p, ncol=p)
    for (off in 1:4) {
      bandVal <- c(.8, .4, .4, .2)[off]
      Omega <- Omega + diag2(bandVal, p, p, off)
    }
    Omega <- Omega + t(Omega) + 2*diag(nrow=p, ncol=p)
    Sigma <- chol2inv(chol(Omega))
  }

  X1 <- mvrnorm(n=n, mu=mu1, Sigma=Sigma)
  X2 <- mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)

  list(X1=X1, X2=X2)
}
