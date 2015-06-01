#' Control Parameters for \code{hotelhd}
#'
#' Various parameters that control setting of tests.
#'
#' @param method character. The method to be used.
#'   Each test needs different paremeters. See 'Details'.
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
#'
#' @export
hotelhd_control <- function(method = c("H", "D", "BS", "CQ", "CLX", "Z", "M"),
                            C = 10, omegaHat = c("omega", "identity"),
                            omegaEst = c("clime", "ada"), omegaGiven = NULL,
                            subForM = c("sub", "diag"), ndim=4,
                            R = 500, block = 1)
{
  method <- match.arg(method)

  if (method %in% c("CLX", "Z", "M")) {
    omegaHat <- match.arg(omegaHat)
    omegaEst <- match.arg(omegaEst)
    params <- list(C = C, omegaHat = omegaHat, omegaEst = omegaEst,
                   omegaGiven=omegaGiven)

    if (method == "Z") {
      params <- c(params, list(R = R, block = block))

    } else if (method == "M") {
      subForM <- match.arg(subForM)
      CLXomega <- match.arg(CLXomega)
      params <- c(params, list(R = R, subForM = subForM, ndim = ndim))
    }

  } else if (method %in% c("H", "D", "BS", "CQ")) {# in case for future modification
    NULL
  }

  return(params)
}


## calculation of Omega hat
## if it is useful, export it later.
calcOmega <- function(...)
{
  extraArgs <- list(...)

  if (length(extraArgs)) {
    namesExtraArgs <- names(extraArgs)
    sugmArgs <- names(formals(sugm)) #legal arg names
    indx <- match(namesExtraArgs, sugmArgs, nomatch=0L)

    if (any(indx==0L)) {
      stop(gettextf("Argument %s not matched.", namesExtraArgs[indx==0L]),
           domain = NA)
    } else if (("method" %in% namesExtraArgs) &&
               extraArgs[["method"]] != "clime") {
      stop("Currently, only clime can be used.", domain = NA)
    } else if ("data" %in% namesExtraArgs) {
      stop("Currently, internally caculated data is used.", domain = NA)
    }

    if ("nlambda" %in% namesExtraArgs) {
      omegaList <- sugm(S, method = "clime", ...)$icov
    } else {
      omegaList <- sugm(S, nlambda=50, method = "clime", ...)$icov
    }

    # default: where no additional arguments for 'sugm'
  } else {
    omegaList <- sugm(S, nlambda=50, method="clime")$icov
  }

  dif <- vapply(omegaList,
                function(O) sum(diag(S%*%O))-log(det(O)),
                FUN.VALUE=vector("numeric", 1L))

  omegaList[[which.min(dif)]]
}
