#' Hotelling's Test with High Dimensional Data
#'
#' TEST
#'
#' Users should take care to ensure that the two dataset have
#' the same variables and the same order of the variables
#'
#'
#' @param d1 a dataset. \code{as.matrix} will be applied.
#' @param d2 a dataset. \code{as.matrix} will be applied.
#'
#' @author Dongjun You
#'
#' @return Test results
#'
#' @references
#' will be added
#'
#'
#' @examples
#' \dontrun{
#' }
#'
#' @importFrom nleqslv nleqslv
#' @import doParallel
#' @import foreach
#'
#' @export
hotelhd <- function(d1, d2, type=c("T2", "nonexact", "BS"), na.rm=TRUE) {
  d1 <- as.matrix(d1)
  d2 <- as.matrix(d2)

  n1 <- NROW(d1)
  n2 <- NROW(d2)
  n <- n1 + n2

  p1 <- NCOL(d1)
  p2 <- NCOL(d2)
  if (p1 != p2)
    stop("Both samples must have the same number of variables (columns)")
  p <- p1

  d1bar <- colMeans(d1, na.rm=na.rm)
  d2bar <- colMeans(d2, na.rm=na.rm)
  meanDiff <- d1bar - d2bar
  S <- ((n1 - 1) * cov(d1) + (n2 - 1) * cov(d2))/(n - 2)

  type <- match.arg(type)

  if (type=="T2") { # Hotelling's T2
    Sinv <- chol2inv(chol(S))
    T2 <- as.vector(n1*n2/n * t(meanDiff) %*% Sinv %*% (meanDiff))
    F_T2 <- (n - p - 1)/(p * (n - 2)) * T2
    pF_T2 <- pf(F_T2, p, n-p-1, lower.tail=FALSE)

    list(statistic=F_T2, pval=pF_T2,
         df=c(p, n - p - 1), n1 = n1, n2 = n2, nvar=p, type=type)

  } else if (type=="nonexact") {
    X <- rbind(d1, d2)
    Htmp <- cbind(1/sqrt(n) * rep(1, n), # 1st row
                    c(rep(sqrt(n2/(n*n1)), n1),
                      -rep(sqrt(n1/(n*n2)), n2))) # 2nd row
    ## n * n orthogonal matrix with given the first two rows
    H <- t(qr.X(qr(Htmp), complete=TRUE))
    Y <- H %*% X

    ## Q vector (Q1,...,Qn)
    Fpre <- Y * Y
    Qs <- rowSums(Fpre)

    ## calculate r1
    t <- (n-2) * log(1/(n-2) * sum(Qs[3:n])) - sum(log(Qs[3:n]))
    r1eq <- function(r1) {
      (1/r1 + (1+1/(n-1))/(3 * r1 * r1)) * (n-3) - t
    }
    r1 <- nleqslv(10, r1eq)$x

    ## calculate r2
    idx <- combn(3:n, 2)
    # if(isTRUE(try(require(doParallel)))) {
    registerDoParallel()
    theta <- foreach (i=1:NCOL(idx)) %dopar% {# angle of every two vectors
      thTmp <- Y[idx[, i], ]
      y1 <- thTmp[1, ]
      y2 <- thTmp[2, ]
      acos(sum(y1 * y2) / (sqrt(sum(y1 * y1) * sum(y2 * y2))))
    }

    sinTheta <- sin(unlist(theta))
    w <- -sum(log(sinTheta * sinTheta))

    r2eq <- function(r2) {
      (1/r2 + (1 + 1/(n-2))/(3*r2*r2)) * (n - 3) +
          (1/r2 + 3/(2*r2*r2)) * choose(n-2, 2) - t - w
    }
    r2 <- nleqslv(10, r2eq)$x

    ## test statistics F
    F_Demp <- (n-2) * sum(Fpre[2, ]) / sum(Fpre[3:n, ])
    pF1_Demp <- pf(F_Demp, r1, r1*(n-2), lower.tail=FALSE)
    pF2_Demp <- pf(F_Demp, r2, r2*(n-2), lower.tail=FALSE)

    list(statistic=F_Demp, r=c(r1=r1, r2=r2),
         pval=c(pval1=pF1_Demp, pval2=pF2_Demp),
         n1=n1, n2=n2, nvar=p, type=type)

  } else if (type=="BS") {
    trS <- sum(diag(S))
    trS2 <- sum(diag(S %*% S))
    Z_BS <- (n1*n2/n * sum(meanDiff * meanDiff) - trS) /
        sqrt(2*(n-1)*(n-2)/(n*(n-3)) * (trS2 - 1/(n-2)*trS*trS))
    pZ_BS <- pnorm(Z_BS, lower.tail=FALSE)

    list(statistics=Z_BS, pval=pZ_BS,
         n1 = n1, n2 = n2, nvar=p, type=type)
  }
}
