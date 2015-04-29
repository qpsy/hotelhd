#' Hotelling's Test with High Dimensional Data
#'
#' Users should take care to ensure that the two dataset have
#' the same variables and the same order of the variables
#'
#'
#' @param X1,X2 n1(,n2) by p data matrix. See Details.
#' @param na.rm logical value of TRUE or FALSE.
#'   option for calculating means of p varuables
#' @param method method to be used. See 'Details'.
#' @param C a constant to decide lambda_n for 'CLX', See Cai et al. (2011).
#' @param omegaHat A candidate estimate of Omega matrix for the two
#'   maximum
#'   type metods, 'CLX' and 'Z'. When the default 'omega' is used,
#'   Omega will be estimated, See 'omegaHat' argument. Optionally,
#'   'identity' can be choosen to set p by p Identity matrix as
#'   the estimate of Omega matrix.
#' @param omegaEst The estimation method for Omega matrix for 'CLX'.
#'   'clime' for clime, 'ada' for Adaptive Thresholding.
#'   The definition of 'ada' follows CLX14's approach for two-sample test.
#' @param omegaGiven If an omega was given, use it. This is an optional but
#'   is included to improve speed of simulation study.
#'   All the methods using 'clime' package will use the given omega,
#'   so do not need to calculate individually.
#' @param subForM for the method 'M'. 'sub' extract sub-matrix from Omega.
#'   'diag' use diagonals of Omega (to improve efficiency)
#' @param ndim The number of dimensions for generalized maximum type statistic.
#'   Default value is 4.
#' @param R The number of bootstrap statistics for 'Z'.
#' @param block A block size for blockwize multiplier
#'   bootstrap in the method 'Z'. The size should be smaller
#'   than min(n1, n2). The default value is 1.
#' @param alpha a numeric value. The type 1 error.
#'
#' @return Test results
#'
#' @details
#' Data 'X1' and 'X2' only accept matrix as arguments. All the methods
#' implemented in this package assume continuous variables as inputs.
#'
#' Method "H" is for Hotelling's T2 test.
#'
#' Method "D" is for Dempster's (1958) non-exact test.
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
#' Method "Z" is blockwise bootstrap
#'
#' Method "M" is generalized bootstrap
#'
#'
#' @references
#' Dempster, A. P. (1958). A high dimensional two sample significance test. \emph{The Annals of Mathematical Statistics}, 995-1010.
#'
#' @importFrom nleqslv nleqslv
#' @import clime
#' @import Rcpp
#' @useDynLib hotelhd
#'
#' @export
hotelhd <- function(X1, X2, na.rm=TRUE,
                    method=c("H", "D", "BS", "CQ", "CLX", "Z", "M"),
                    C=10, omegaHat=c("omega", "identity"),
                    omegaEst=c("clime", "ada"), omegaGiven=NULL,
                    subForM=c("sub", "diag"), ndim=4,
                    R=500, block=1, alpha=0.05) {
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

  omegaHat <- match.arg(omegaHat)
  omegaEst <- match.arg(omegaEst)
  calcOmegaEst <- function() {
    lambda <- C * (log(p)/n)

    if (omegaEst == "clime") {
      ## clime package
      clime(S, sigma=TRUE, lambda=lambda)$Omegalist[[1]]

      ## fastclime package
      #fc <- suppressMessages(fastclime(S))
      #OmegaHat <- fastclime.lambda(fc$lambdamtx, fc$icovlist, lambda)$icov

    } else if (omegaEst == "ada") {
      delta <- 2
      ss1 <- (sweep(X1, 2, X1bar))^2
      ss2 <- (sweep(X2, 2, X2bar))^2
      theta <- (t(ss1) %*% ss1 + (2-n1)*(var(X1))^2 +
                    t(ss2) %*% ss2 + (2-n2)*(var(X2))^2) / n

      lambda <- delta * sqrt(log(p)*theta / n)
      Sstar <- S * (abs(S) >= lambda)
      solve(Sstar)
    }
  }

  method <- match.arg(method)

  if (method=="H") { # Hotelling's T2
    Sinv <- solve(S)
    H <- as.vector(n1*n2/n * t(meanDiff) %*% Sinv %*% (meanDiff))
    F_H <- (n - p - 1)/(p * (n - 2)) * H
    pF_H <- pf(F_H, p, n-p-1, lower.tail=FALSE)

    list(statistic=F_H, pval=pF_H, df=c(df1=p, df2=n-p-1),
         nobs=c(n1=n1, n2=n2), nvar=p,
         rejected= pF_H < alpha, method=method)

  ##---------------------------------------------------------
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
         nobs=c(n1=n1, n2=n2), nvar=p,
         rejected=c(r1=pF1_D < alpha, r2=pF2_D < alpha), method=method)

  ##-------------------------------------------------------------
  } else if (method=="BS") {
    trS <- sum(diag(S))
    trS2 <- sum(diag(S %*% S))
    Z_Mn <- (n1*n2/n * sum(meanDiff * meanDiff) - trS) /
        sqrt(2*(n-1)*(n-2)/(n*(n-3)) * (trS2 - 1/(n-2)*trS*trS))
    pZ_Mn <- pnorm(Z_Mn, lower.tail=FALSE)

    list(statistics=Z_Mn, pval=pZ_Mn,
         nobs=c(n1=n1, n2=n2), nvar=p,
         rejected=pZ_Mn < alpha, method=method)

  ##--------------------------------------------------
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
         nobs=c(n1=n1, n2=n2), nvar=p,
         rejected=pZ_Tn < alpha, method=method)

  ##--------------------------------------------------------
  } else if (method=="CLX") {
    if (is.null(omegaGiven)) {
      if (omegaHat == "omega") Omega <- calcOmegaEst()
      else Omega <- diag(nrow=p, ncol=p)

    } else {
      Omega <- omegaGiven
    }

    Z <- Omega %*% (X1bar - X2bar)
    omega1 <- (n1-1)/n1 * var(X1 %*% Omega)
    omega2 <- (n2-1)/n2 * var(X2 %*% Omega)
    omega0 <- diag(n1/n * omega1 + n2/n * omega2)

    M <- n1*n2/n * max((Z*Z) / omega0)
    #part1 <- 2*log(p) - log(log(p)) - log(pi)
    rejected <- M >= 2*log(p) - log(log(p)) - log(pi) -
        2*log(log(1/(1-alpha)))
    #rejected <- function(alpha) {
    #  M >= part1 - 2*log(log(1/(1-alpha)))
    #}

    list(statistic=M, nobs=c(n1=n1, n2=n2), nvar=p,
         rejected=rejected, method=method)

  ##--------------------------------------------------
  } else if (method == "Z") {
    if (block > min(n1, n2)) stop("The block size is greater than nobs.")

    if (is.null(omegaGiven)) {
      if (omegaHat == "omega") Omega <- calcOmegaEst()
      else Omega <- diag(nrow=p, ncol=p)

    } else {
      Omega <- omegaGiven
    }

    Z <- Omega %*% (X1bar - X2bar)
    XI1 <- X1 %*% Omega
    XI2 <- X2 %*% Omega
    XI1diff <- sweep(XI1, 2, colMeans(XI1), check.margin=FALSE)
    XI2diff <- sweep(XI2, 2, colMeans(XI2), check.margin=FALSE)
    omega1 <- t(XI1diff) %*% XI1diff / n1
    omega2 <- t(XI2diff) %*% XI2diff / n2
    omega0sqrt <- sqrt(diag(n1/n * omega1 + n2/n * omega2))

    ## test statistics
    sqnn <- sqrt(n1*n2/n)
    T <- max(abs(sqnn * Z))
    Tt <- max(abs(sqnn * Z / omega0sqrt))

    # number of blocks
    if (n1 %% block) l1 <- (n1 %/% block) + 1
    else l1 <- n1 / block

    if (n2 %% block) l2 <- (n2 %/% block) + 1
    else l2 <- n2 / block

    T_boot <- sqnn * quantile(
        vapply(1:R, function(i) {
          max(abs(
              (colMeans(sweep(XI1diff, 1, rep(rnorm(l1), each=block)[1:n1],
                              FUN="*", check.margin=FALSE)) -
               colMeans(sweep(XI2diff, 1, rep(rnorm(l2), each=block)[1:n2],
                              FUN="*", check.margin=FALSE)))))},
               FUN.VALUE=vector("numeric", 1), USE.NAMES=FALSE),
        1 - alpha, names=FALSE)

    Tt_boot <- sqnn * quantile(
        vapply(1:R, function(i) {
          max(abs(
              (colMeans(sweep(XI1diff, 1, rep(rnorm(l1), each=block)[1:n1],
                              FUN="*", check.margin=FALSE)) -
               colMeans(sweep(XI2diff, 1, rep(rnorm(l2), each=block)[1:n2],
                              FUN="*", check.margin=FALSE)))/omega0sqrt))},
               FUN.VALUE=vector("numeric", 1), USE.NAMES=FALSE),
        1 - alpha, names=FALSE)

    list(T=c(T=T, T_boot=T_boot), Tt=c(Tt=Tt, Tt_boot=Tt_boot),
         nobs=c(n1=n1, n2=n2), nvar=p,
         rejected=c(T=T > T_boot, Tt=Tt > Tt_boot), method=method)

  ##--------------------------------------------------------
  } else if (method == "M") {
    if (is.null(omegaGiven)) {
      if (omegaHat == "omega") Omega <- calcOmegaEst()
      else Omega <- diag(nrow=p, ncol=p)

    } else {
      Omega <- omegaGiven
    }

    ## \hat{Z}_{(b)}
    X1b <- sweep(X1, 2, X1bar, check.margin=FALSE) %*% Omega
    X2b <- sweep(X2, 2, X2bar, check.margin=FALSE) %*% Omega
    # omegaInv <- solve(Omega)
    jk <- combn(p, ndim) # column index
    jkc <- NCOL(jk)

    ## M statistic (omitted constant term)
    calcM1 <- function(Z) {
      max(# M(Omega)
          vapply(1:jkc, function(i) {
            jk_i <- jk[, i]
            Z_jk <- Z[jk_i]
            omegaInv_jk <- solve(Omega[jk_i, jk_i])
            t(Z_jk) %*% omegaInv_jk %*% Z_jk
          },
          FUN.VALUE=vector("numeric", 1), USE.NAMES=FALSE))
    }

    ## M statistic (omitted constant term)
    ##  - efficient way of using off-diagonal of omegaInv_jk
    calcM2 <- function(Z) {
      sum(sort(Z / diag(Omega), decreasing=TRUE)[c(1:ndim)])
    }

    subForM <- match.arg(subForM)
    if (subForM == "sub") calcM <- calcM1
    else calcM <- calcM2 # if (subForM == "diag")

    ## test statistic
    Zo <- Omega %*% (X1bar - X2bar)
    M <- calcM(Zo)

    ## constant related to n1, n2 are safely ommitted in calculation and comparison
    M_boot <- quantile(
        vapply(1:R, function(i) {# R of boot statistics
          Zb <-
            colMeans(
                sweep(X1b, 1, rnorm(n1),
                      FUN="*", check.margin=FALSE)) -
            colMeans(
                sweep(X2b, 1, rnorm(n2),
                      FUN="*", check.margin=FALSE))
          calcM(Zb)
        },
        FUN.VALUE=vector("numeric", 1), USE.NAMES=FALSE),
        1 - alpha, names=FALSE)

    list(M=(n1*n2/(n1+n2)) * M,
         nobs=c(n1=n1, n2=n2), nvar=p,
         rejected= M > M_boot, method=method)
  }
}
