#' Control Parameters for \code{hotelhd}
#'
#' Various parameters that control setting of tests. Currently,
#' the maximum type tests (i.e., 'CLX', 'Z', 'M') only need control parameters.
#'
#' @param omegaGiven \code{NULL} or \eqn{p} by \eqn{p} matrix.
#'   When the default \code{NULL} is chosen, \eqn{\Omega} matrix is estimated;
#'   see \code{omegaHat} and \code{omegaEst} parameters.
#'   Optionally, users can provide an estimate of \eqn{\Omega}
#'   from an estimation method that this package does not offer.
#'
#' @param omegaHat character; 'omega' or 'identity'.
#'   A candidate estimate of \eqn{\Omega} matrix for the maximum type methods.
#'   When the default of 'omega' is used, \eqn{\Omega} is estimated.
#'   Users can choose estimation method by \code{omegaEst} parameter.
#'   Optionally,
#'   'identity' can be choosen to set \eqn{p} by \eqn{p} Identity matrix as
#'   an estimate of \eqn{\Omega} matrix.
#'   When the data is not high-dimensional, inverse of sample covariance is offered.
#'
#' @param omegaEst character; 'clime'.
#'   The estimation method for \eqn{\Omega} matrix.
#'   Currently, only 'clime' is implemented; see \code{\link[flare]{sugm}}.
#'
#' @param ndim integer. The number of dimensions for generalized maximum
#'   type statistic, i.e., the method 'M'. The default is 4.
#'
#' @param sub character; 'diag' or 'sub'.
#'   Method about how to form \code{ndim} by \code{ndim} sub-matrix for
#'   the method 'M'. The default 'diag' uses diagonals of \eqn{\Omega}
#'   to improve efficiency, and calculate results of 1 to \code{ndim}.
#'   The 'sub' extracts sub-matrix from \eqn{\Omega}.

#' @param B integer. The number of bootstrap statistics for 'Z'.
#'
#' @param block integer. A block size for blockwize multiplier
#'   bootstrap in the method 'Z'. The size should be smaller
#'   than min(n1, n2). The default value is 1.
#'
#' @export
hotelhd_control <- function(omegaGiven=NULL,
                            omegaHat=c("omega", "identity"),
                            omegaEst=c("clime"),
                            ndim=4L, sub=c("diag", "sub"),
                            B=500L, block=1L)
{
  omegaHat <- match.arg(omegaHat)
  omegaEst <- match.arg(omegaEst)
  sub <- match.arg(sub)

  list(C=C, omegaGiven=omegaGiven, omegaHat=omegaHat,
       omegaEst=omegaEst, ndim=ndim, sub=sub, B=B, block=block)
}
