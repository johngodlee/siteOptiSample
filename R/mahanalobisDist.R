#' Calculate mahalanobis distance between two matrices
#' 
#' @param x numeric matrix
#' @param y numeric matrix
#' 
#' @noRd
#' 

mahanalobisDist <- function(x, y) {
  cov <- cov(rbind(x, y))
  S_inv <- solve(cov)
  out <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      diff <- x[i, ] - y[j, ]
      out[i, j] <- sqrt(t(diff) %*% S_inv %*% diff)
    }
  }
  return(out)
}

