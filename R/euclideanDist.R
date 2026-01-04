#' Calculate euclidean distance between two matrices
#' 
#' @param x numeric matrix of matrix 1
#' @param y numeric matrix of matrix 2
#' @param w optional, numeric vector of weights (one per column). Larger values
#'     increase a variable’s contribution; if NULL, all variables are equally
#'     weighted.
#' 
#' @noRd
#' 

euclideanDist <- function(x, y, w = NULL) {
  # Check matrices have the same number of columns
  if (ncol(x) != ncol(y)) { 
    stop("'x' and 'y' must have the same number of columns")
  }

  # Check weights equal the number of columns as matrices
  if (!is.null(w) && ncol(x) != length(w)) {
    stop("The length of 'w' must equal the number of columns in 'x'")
  }

  # If no weights supplied, use equal weights
  if (is.null(w)) {
    w <- rep(1, ncol(x))
  }

  if (any(w < 0)) {
    stop("'w' must be non-negative")
  }

  out <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      diff <- x[i, ] - y[j, ]
      out[i, j] <- sqrt(sum(w * diff^2))
    }
  }
  return(out)
}

