#' hamming distance
#'
#' calculates the hamming distance between two matrices or between columns
#' of a matrix using matrix multiplication. Copied from [Johan DeJong's blog](https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/),
#' copyright is ascribed to them.
#'
#' @param X the first matrix
#' @param Y the second matrix
#'
#' @export
#' @return double
hamming <- function(X, Y) {
  if ( missing(Y) ) {
    uniqs <- unique(as.vector(X))
    U <- X == uniqs[1]
    H <- t(U) %*% U
    for ( uniq in uniqs[-1] ) {
      U <- X == uniq
      H <- H + t(U) %*% U
    }
  } else {
    uniqs <- union(X, Y)
    H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
    for ( uniq in uniqs[-1] ) {
      H <- H + t(X == uniq) %*% (Y == uniq)
    }
  }
  nrow(X) - H
}