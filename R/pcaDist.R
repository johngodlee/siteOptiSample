#' Get nearest neighbour distances between two PCAs
#' 
#' @param x PCA scores 
#' @param y PCA scores
#' @param w optional, numeric vector of weights (one per column). Larger values
#'     increase a variable’s contribution; if NULL, all variables are equally
#'     weighted.
#' @param n_pca number of PCA axes used in analysis
#' @param k number of nearest neighbours in y to return for each row in x
#' @param method distance method, either "euclidean" or "mahalanobis"
#'
#' @return if `k = 1` a vector of distances to the nearest neighbour in `y` for
#'     each row in `x`. If `k = 2` a matrix with k columns with ordered nearest
#'     neighbour distances
#'
#' @export
#' 
pcaDist <- function(x, y, w = NULL, n_pca = 3, k = 1, method = "euclidean") {
  # Calculate nearest neighbor distances
  # For each landscape pixel, find distance to nearest plot
  if (method == "mahalanobis") {
    dists_mat <- mahalanobisDist(x[,1:n_pca], y[,1:n_pca], w)
  } else if (method == "euclidean") {
    dists_mat <- euclideanDist(x[,1:n_pca], y[,1:n_pca], w)
  } else {
    stop("method must be either 'euclidean' or 'mahalanobis'")
  }

  min_dists <- apply(dists_mat, 1, function(i) {
    sort(i)[1:k]
  }, simplify = FALSE)
  out <- do.call(rbind, min_dists)

  if (ncol(out) == 1) {
    out <- c(out)
  }

  # Return
  return(out)
}

