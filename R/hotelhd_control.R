#' Control Parameters for \code{hotelhd}
#'
#' Various parameters that control setting of tests. Currently,
#' the maximum type tests (i.e., 'CLX', 'Z', 'M') only need control parameters.
#'
#' @param C a constant to decide lambda_n for 'CLX', See Cai et al. (2011).
#' @param omegaHat character, "omega" or"identity".
#'   A candidate estimate of Omega matrix for the two
#'   maximum
#'   type metods, 'CLX' and 'Z'. When the default 'omega' is used,
#'   Omega will be estimated, See 'omegaHat' argument. Optionally,
#'   'identity' can be choosen to set p by p Identity matrix as
#'   the estimate of Omega matrix.
#' @param omegaEst character (c("clime")).
#'   The estimation method for \eqn{\Omega} matrix,
#'   which is used for maximum type tests.
#'   Currently, only 'clime' is implemented; see \code{flare::sugm}.
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
hotelhd_control <- function(C = 10, omegaHat = c("omega", "identity"),
                            omegaEst = c("clime"), omegaGiven = NULL,
                            subForM = c("sub", "diag"), ndim=4,
                            R = 500, block = 1)
{
  #method <- match.arg(method)

  #if (method %in% c("CLX", "Z", "M")) {
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

  #} else if (method %in% c("H", "D", "BS", "CQ")) {# in case for future modification
  #  NULL
  #}

  return(params)
}
