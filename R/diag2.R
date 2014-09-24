#' A convenient version of diag
#'
#' text
#'
#' @param x
#' @param nrow
#' @param ncol
#' @param off an integer. location of off-diagonal
#'
#' @author Dongjun You
#'
#' @return Test results
#'
#' @details
#' ...
#'
#' @examples
#' \dontrun{
#' }
#'
#' @export

diag2 <- function(x=1, nrow, ncol, off=0) {
  M <- matrix(0, nrow=nrow, ncol=ncol)
  M[col(M) == (row(M) + off)] <- x
  M
}

