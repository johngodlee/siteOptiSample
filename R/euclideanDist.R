#' Calculate euclidean distance between two matrices
#' 
#' @param x numeric matrix
#' @param y numeric matrix
#' 
#' @noRd
#' 

euclideanDist <- function(x, y) {
  out <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      diff <- x[i, ] - y[j, ]
      out[i, j] <- sqrt(sum(diff^2))
    }
  }
  return(out)
}

