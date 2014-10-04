#' Hotelling's Test with High Dimensional Data
#'
#' Users should take care to ensure that the two dataset have
#' the same variables and the same order of the variables
#'
#'
#' @param X1 a n1 by p data matrix. Only a matrix can be used.
#' @param X2 a n2 by p data matrix.
#' @param na.rm logical value of TRUE or FALSE.
#'   option for calculating means of p varuables
#' @param method method to be used. See 'Details'.
#' @param C a constant to decide lambda_n for 'CLX', See Cai et al. (2011).
#' @param omega a character value. the estimation of Omega matrix for 'CLX'.
#'   'clime' for clime, 'ada' for Adaptive Thresholding.
#'   The definition of 'ada' follows CLX14's approach for two-sample test.
#'
#' @author Dongjun You
#'
#' @return Test results
#'
#' @details
#' Method "H" is ...
#'
#' Method "D" is
#'
#' Method "BS" is
#'
#' Method "CQ" is
#'
#' Method "CLX" is using CLIME to estimate Omega. The current version
#'   uses the default options of 'fastclime' function of
#'   'fastclime' package.
#'   Whether the null is rejected or not is only included in the result,
#'   because of the form of asymptotic distribution of the test statistics.
#'   A function 'rejected' for hypothesis test will be returned, which
#'   needs an argument 'alpha' to test Whether the null is rejected or not.
#'
#' Method "Z" is
#'
#' @references
#' will be added
#'
#' @importFrom nleqslv nleqslv
#' @import clime
#' @import Rcpp
#' @useDynLib hotelhd
#'
#' @export
hotelhd <- function(X1, X2, na.rm=TRUE,
                    method=c("H", "D", "BS", "CQ", "CLX", "Z"),
                    C=10, omega=c("clime", "ada")) {
  stopifnot(is.matrix(X1), is.matrix(X1))

  n1 <- NROW(X1)
  n2 <- NROW(X2)
  n <- n1 + n2

  p1 <- NCOL(X1)
  p2 <- NCOL(X2)
  if (p1 != p2)
    stop("Both samples must have the same number of variables (columns)")
  p <- p1

  X1bar <- colMeans(X1, na.rm=na.rm)
  X2bar <- colMeans(X2, na.rm=na.rm)
  meanDiff <- X1bar - X2bar
  S <- ((n1 - 1)*var(X1) + (n2 - 1)*var(X2)) / (n - 2)

  method <- match.arg(method)

  if (method=="H") { # Hotelling's T2
    Sinv <- chol2inv(chol(S))
    H <- as.vector(n1*n2/n * t(meanDiff) %*% Sinv %*% (meanDiff))
    F_H <- (n - p - 1)/(p * (n - 2)) * H
    pF_H <- pf(F_H, p, n-p-1, lower.tail=FALSE)

    list(statistic=F_H, pval=pF_H, df=c(df1=p, df2=n-p-1),
         nobs=c(n1=n1, n2=n2), nvar=p, method=method)

  } else if (method=="D") {
    X <- rbind(X1, X2)
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
    idx_n <- NCOL(idx)

    sinTheta <- vector("numeric", length=idx_n)
    sinTheta_i <- 0
    for (i in 3:(n-1)) {
      y1 <- Y[i, ]
      for (j in (i+1):n) {
        sinTheta_i <- sinTheta_i + 1
        y2 <- Y[j, ]
        sinTheta[sinTheta_i] <- sin(# angle of every two vectors
            acos(sum(y1 * y2) / sqrt(Qs[i] * Qs[j])))
      }
    }
    w <- -sum(log(sinTheta * sinTheta))

    r2eq <- function(r2) {
      (1/r2 + (1 + 1/(n-2))/(3*r2*r2)) * (n - 3) +
          (1/r2 + 3/(2*r2*r2)) * choose(n-2, 2) - t - w
    }
    r2 <- nleqslv(10, r2eq)$x

    ## test statistics F
    F_D <- (n-2) * sum(Fpre[2, ]) / sum(Fpre[3:n, ])
    pF1_D <- pf(F_D, r1, r1*(n-2), lower.tail=FALSE)
    pF2_D <- pf(F_D, r2, r2*(n-2), lower.tail=FALSE)

    list(statistic=F_D, r=c(r1=r1, r2=r2),
         pval=c(pval1=pF1_D, pval2=pF2_D),
         nobs=c(n1=n1, n2=n2), nvar=p, method=method)

  } else if (method=="BS") {
    trS <- sum(diag(S))
    trS2 <- sum(diag(S %*% S))
    Z_Mn <- (n1*n2/n * sum(meanDiff * meanDiff) - trS) /
        sqrt(2*(n-1)*(n-2)/(n*(n-3)) * (trS2 - 1/(n-2)*trS*trS))
    pZ_Mn <- pnorm(Z_Mn, lower.tail=FALSE)

    list(statistics=Z_Mn, pval=pZ_Mn,
         nobs=c(n1=n1, n2=n2), nvar=p, method=method)

  } else if (method=="CQ") {
    prod11 <- tcrossprod(X1)
    prod22 <- tcrossprod(X2)

    Tn <- (sum(prod11) - sum(diag(prod11))) / (n1*(n1-1)) +
        (sum(prod22) - sum(diag(prod22))) / (n2*(n2-1)) -
            2 * sum(tcrossprod(X1, X2)) / (n1*n2)

    ## caltrSigmaSq <- function(X, nn) {
    ##   Sigma1Sq <- matrix(0, nrow=p, ncol=p)
    ##   for (i in 1:(nn-1)) {
    ##     Xj <- X[i, , drop=FALSE]
    ##     for (j in (i+1):nn) {
    ##       Xk <- X[j, , drop=FALSE]
    ##       Xbar_jk <- colMeans(X[-c(i,j), ])
    ##       Sigma1Sq <- Sigma1Sq +
    ##           t(Xj - Xbar_jk) %*% Xj %*% t(Xk - Xbar_jk) %*% Xk
    ##     }
    ##   }
    ##   1/(nn*(nn-1)) * 2*sum(diag(Sigma1Sq))
    ## }

    ## trSigma1Sq <- caltrSigmaSq(X1, n1)
    ## trSigma2Sq <- caltrSigmaSq(X2, n2)

    ## Sigma12 <- matrix(0, nrow=p, ncol=p)
    ## for (j in 1:n1) {
    ##   X1j <- X1[j, , drop=FALSE]
    ##   X1bar_j <- colMeans(X1[-j, ])
    ##   for (k in 1:n2) {
    ##     X2k <- X2[k, , drop=FALSE]
    ##     X2bar_k <- colMeans(X2[-k, ])
    ##     Sigma12 <- Sigma12 +
    ##         t(X1j - X1bar_j) %*% X1j %*% t(X2k - X2bar_k) %*% X2k
    ##   }
    ## }

    ## trSigma12 <- 1/(n1*n2) * sum(diag(Sigma12))

    trSigma1Sq <- caltrSigmaSqC(X1, n1, p)
    trSigma2Sq <- caltrSigmaSqC(X2, n2, p)
    trSigma12 <- caltrSigma12C(X1, X2, n1, n2, p)


    Var_Tn <- 2/(n1*(n1-1)) * trSigma1Sq + 2/(n2*(n2-1)) * trSigma2Sq +
        4/(n1*n2) * trSigma12

    Z_Tn <- Tn / sqrt(Var_Tn)
    pZ_Tn <- pnorm(Z_Tn, lower.tail=FALSE)

    list(statistics=Z_Tn, pval=pZ_Tn,
         nobs=c(n1=n1, n2=n2), nvar=p, method=method)

  } else if (method=="CLX") {
    lambda <- C * (log(p)/n)

    if (omega == "clime") {
      ## clime package
      Omega <- clime(S, sigma=TRUE, lambda=lambda)$Omegalist[[1]]

      ## fastclime package
      #fc <- suppressMessages(fastclime(S))
      #Omega <- fastclime.lambda(fc$lambdamtx, fc$icovlist, lambda)$icov

    } else if (omega == "ada") {
      delta <- 2
      ss1 <- (sweep(X1, 2, X1bar))^2
      ss2 <- (sweep(X2, 2, X2bar))^2
      theta <- (t(ss1) %*% ss1 + (2-n1)*(var(X1))^2 +
                t(ss2) %*% ss2 + (2-n2)*(var(X2))^2) / n

      lambda <- delta * sqrt(log(p)*theta / n)
      Sstar <- S * (abs(S) >= lambda)
      Omega <- solve(Sstar)
    }

    Z <- Omega %*% (X1bar - X2bar)
    omega1 <- var(X1 %*% t(Omega))
    omega2 <- var(X2 %*% t(Omega))
    omega0 <- diag(n1/n * omega1 + n2/n * omega2)

    M <- n1*n2/n * max((Z*Z) / omega0)
    part1 <- 2*log(p) - log(log(p)) - log(pi)
    #rejected <- M >= 2*log(p) - log(log(p)) - log(pi) - 2*log(log(1/(1-alpha)))
    rejected <- function(alpha) {
      M >= part1 - 2*log(log(1/(1-alpha)))
    }

    list(statistic=M, rejected=rejected)

  } else if (method == "Z") {
    lambda <- C * (log(p)/n)

    ## clime package
    Omega <- clime(S, sigma=TRUE, lambda=lambda)$Omegalist[[1]]

    ## fastclime package
    #fc <- suppressMessages(fastclime(S))
    #Omega <- fastclime.lambda(fc$lambdamtx, fc$icovlist, lambda)$icov

    Z <- Omega %*% (X1bar - X2bar)
    omega1 <- var(X1 %*% t(Omega))
    omega2 <- var(X2 %*% t(Omega))
    omega0 <- diag(n1/n * omega1 + n2/n * omega2)

  }
}
