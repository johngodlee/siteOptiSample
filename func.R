###########
# UTILITIES
###########

# Resample a raster to any resolution
#
# @param x raster object
# @param res desired resolution, in units of the projection of `x`
# @param ... additional arguments passed to `terra::resample()`
#
# @import terra
#
# @return raster object
# 
rastResample <- function(x, res, ...) { 
  r_rs <- x
  terra::res(r_rs) <- res
  terra::resample(x, r_rs, ...)
}

###################################
# EXTRACT METRICS AT PLOT LOCATIONS
###################################

# Check if object is sf type
#
# @param x object to test
# @param type optional character vector of acceptable sf geometry types
#
# @return logical
# 
# @import sf
# 
isSFType <- function(x, type = NULL) {
  inherits(x, c("sf", "sfc")) && (is.null(type) |
    all(st_geometry_type(x, by_geometry = FALSE) %in% type))
}

# Extract values from raster across polygons
#
# @param r SpatRaster object
# @param p either an sf object containing plot polygons or points, or a
#     dataframe with two columns `x` and `y` describing coordinates in the same
#     coordinate system as `r`
# @param fun function to summarize the extracted data by polygon geometry, see
#     `terra::extract()`
# @param ... additional parameters passed to `terra::extract()`
#
# @return dataframe containing extracted metrics for each polygon
# 
# @import terra
# @import sf
# 
extract_plot_metrics <- function(r, p, fun = mean, ...) {
  # Check type of r
  if (!inherits(r, "SpatRaster")) {
    stop("r must be of class SpatRaster")
  }

  # Convert p to sf if not already
  if (isSFType(p, c("POINT", "POLYGON"))) {
    p_vect <- vect(p)
  } else if (inherits(p, "data.frame") && 
      all(c("x", "y") %in% colnames(p))) {
    p_vect <- st_as_sf(p, coords = c("x", "y"), crs = crs(r))
  } else {
    stop("p must either be an sf POINT/POLYGON object or a data.frame with columns 'x' and 'y'")
  }

  # Extract metrics for each plot polygon
  out <- extract(r, p_vect, fun = fun, 
    na.rm = TRUE, ID = FALSE, ...)
  
  # Return
  return(out)
}

#############################
# REPRESENTATIVENESS ANALYSIS
#############################

# Run a PCA on landscape structural metrics, and optionally place plots in
# that space
# 
# @param r SpatRaster with structural metrics
# @param p optional dataframe containing structural metrics for each plot, e.g.
#     as returned by `extract_plot_metrics()`
# @param ... additional arguments passed to `prcomp()`  
#
# @return list containing: `r_pca`: PCA object from pixel values, and
#     optionally `p_pca`: plot values in PCA space
# 
pca_landscape <- function(r, p = NULL, ...) {
  # Check type of r
  if (!inherits(r, "SpatRaster")) {
    stop("r must be of class SpatRaster")
  }

  # Convert raster to data frame for landscape
  r_df <- as.data.frame(r, xy = FALSE, na.rm = TRUE)
  
  # Select only metric columns 
  metric_cols <- names(r_df)
  
  # PCA to reduce dimensionality
  r_pca <- prcomp(r_df[,metric_cols], ...)
  
  # Project plots into PCA space
  if (!is.null(p)) { 
    p_pca <- scale(p[,metric_cols], 
      center = r_pca$center, scale = r_pca$scale) %*% r_pca$rotation
  } else {
    p_pca <- NULL
  }

  # Return
  return(list(
    "r_pca" = r_pca,
    "p_pca" = p_pca))
}

euclidean_dist <- function(x, y) {
  out <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      diff <- x[i, ] - y[j, ]
      out[i, j] <- sqrt(sum(diff^2))
    }
  }
  return(out)
}

mahalanobis_dist <- function(x, y) {
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

# Get nearest neighbour distances between two PCAs
# 
# @param x PCA scores 
# @param y PCA scores
# @param n_pca number of PCA axes used in analysis
# @param k number of nearest neighbours in y to return for each row in x
# @param method distance method, either "euclidean" or "mahalanobis"
#
# @return if `k = 1` a vector of distances to the nearest neighbour in `y` for
#     each row in `x`. If `k = 2` a matrix with k columns with ordered nearest
#     neighbour distances
#
# @import proxy
# 
pca_dist <- function(x, y, n_pca = 3, k = 1, method = "euclidean") {
  # Calculate nearest neighbor distances
  # For each landscape pixel, find distance to nearest plot
  if (method == "mahalanobis") {
    dists_mat <- mahalanobis_dist(x[,1:n_pca], y[,1:n_pca])
  } else if (method == "euclidean") {
    dists_mat <- euclidean_dist(x[,1:n_pca], y[,1:n_pca])
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


##############################################
# IDENTIFY GAPS AND OPTIMALLY LOCATE NEW PLOTS
##############################################

# Iteratively add candidate plots using the maximin algorithm
# 
# @param r_pca PCA scores of pixels in structural space
# @param p_new_pca optional, PCA scores of candidate plots in structural space.
#     If supplied, selection is restricted to these candidates; otherwise the
#     full set of landscape pixels in `r_pca` is used.
# @param p_pca optional, PCA scores of existing plots in structural space
# @param n_plots maximum number of new plots to add
# @param n_pca number of PCA axes used in analysis
#
# @details The maximin algorithm aims to maximise the minimum distance between
#     proposed plots and existing plots, iteratively placing new plots in
#     pixels with structural attributes most dissimilar to existing plots. For
#     each candidate pixel, compute the distance to the nearest existing plot,
#     then choose the pixel with the highest distance value. As a result, plots
#     occupy structural extremes.

# @return vector of pixel index values for proposed new plots, ordered by
#     largest increase in coverage 
#
# @import proxy
# 
maximin_select <- function(r_pca, p_new_pca = NULL, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[,1:n_pca, drop = FALSE])

  if (!is.null(p_new_pca)) {
    p_new_pca <- as.matrix(p_new_pca[,1:n_pca, drop = FALSE])
  } else {
    p_new_pca <- r_pca
  }

  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca[,1:n_pca, drop = FALSE])
  }
  
  # Initialize minimum distances to existing plots
  if (!is.null(p_pca)) {
    current_plots <- p_pca 
  } else {
    current_plots <- NULL
  }

  # Initialize selected indices
  cand_idx <- seq_len(nrow(p_new_pca))
  sel_idx <- c()

  for (i in seq_len(n_plots)) {
    # Compute minimum distance to nearest current plot for each candidate
    d_cands <- sapply(seq_len(nrow(p_new_pca)), function(j) {
      if (is.null(current_plots)) {
        # If no existing plots, distance to whole raster space
        min(proxy::dist(r_pca, p_new_pca[j, , drop = FALSE]))
      } else {
        # If existing, distance from candidate to existing + selected plots
        min(proxy::dist(p_new_pca[j, , drop = FALSE], current_plots))
      }
    })
    
    # Select candidate with maximum minimum distance
    best_rel_idx <- which.max(d_cands)
    best_idx <- cand_idx[best_rel_idx]   # original index

    # Store original index
    sel_idx <- c(sel_idx, best_idx)
    
    # Update current plots and remove chosen candidate
    current_plots <- rbind(current_plots, p_new_pca[best_rel_idx, , drop = FALSE])
    p_new_pca <- p_new_pca[-best_rel_idx, , drop = FALSE]
    cand_idx <- cand_idx[-best_rel_idx]

    # Stop if no candidates remain
    if (nrow(p_new_pca) == 0) {
      warning("No suitable locations for plot ", i)
      break
    }
  }
  
  # Return
  return(sel_idx)
}
  

# Iteratively add candidate plots using the minimax algorithm
# 
# @param r_pca PCA score of pixels in structural space
# @param p_new_pca optional, PCA scores of candidate plots in structural space.
#     If supplied, selection is restricted to these candidates; otherwise the
#     full set of landscape pixels in `r_pca` is used.
# @param p_pca optional, PCA scores of existing plots in structural space
# @param n_plots maximum number of new plots to add
# @param n_pca number of PCA axes used in analysis
#
# @details
# The minimax algorithm aims to minimize the maximum distance of any pixel to
#     its nearest plot in PCA space, iteratively placing new plots in pixels
#     with structural attributes most dissimilar to their nearest existing
#     plot. For each candidate pixel, compute how adding a plot to that pixel
#     would decrease the maximum pixel-to-nearest-plot distance, then choose
#     the pixel which most effectively minimises the distance value. As a
#     result, plots fill gaps in structural space.
#
# @return vector of pixel index values for proposed new plots, ordered by
# largest increase in coverage 
# 
# @import proxy
# 
minimax_select <- function(r_pca, p_new_pca = NULL, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[,1:n_pca, drop = FALSE])

  if (!is.null(p_new_pca)) {
    p_new_pca <- as.matrix(p_new_pca[,1:n_pca, drop = FALSE])
  } else {
    p_new_pca <- r_pca
  }

  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca[,1:n_pca, drop = FALSE])
  }
  
  # Initialize candidate indices
  cand_idx <- seq_len(nrow(p_new_pca)) 
  sel_idx <- c()
  
  # Initialize minimum distances
  if (!is.null(p_pca)) {
    min_dists <- apply(as.matrix(proxy::dist(r_pca, p_pca)), 1, min)
  } else {
    min_dists <- rep(Inf, nrow(r_pca))
  }
  
  # Precompute distance matrix: pixels x candidate plots
  d_all <- as.matrix(proxy::dist(r_pca, p_new_pca))
  
  for (i in seq_len(n_plots)) {
    # Distances for remaining candidates
    d_cands <- d_all[,cand_idx, drop = FALSE]
    
    # Compute minimax gain
    gains <- apply(d_cands, 2, function(col) {
       max(pmin(min_dists, col))
    })
    
    # Select candidate with minimum gain (tie-break randomly)
    best_rel <- sample(which(gains == min(gains)), 1)
    best_idx <- cand_idx[best_rel]
    sel_idx <- c(sel_idx, best_idx)
    
    # Update minimum distances and candidate pool
    min_dists <- pmin(min_dists, d_all[, best_idx])
    cand_idx <- cand_idx[-best_rel]
    
    if (length(cand_idx) == 0) {
      warning("No suitable locations for plot ", i)
      break
    }
  }
  
  # Return
  return(sel_idx)
}

# Iteratively add candidate plots using the minimax algorithm
# 
# @param r_pca PCA score of pixels in structural space
# @param p_new_pca optional, PCA scores of candidate plots in structural space.
#     If supplied, selection is restricted to these candidates; otherwise the
#     full set of landscape pixels in `r_pca` is used.
# @param p_pca optional, PCA scores of existing plots in structural space
# @param n_plots maximum number of new plots to add
# @param n_pca number of PCA axes used in analysis
#
# @details
# The meanmin algorithm aims to minimize the mean distance of all pixels to
#     their nearest plot in PCA space, iteratively placing new plots in pixels
#     that most reduce the average pixel-to-nearest-plot distance. For each
#     candidate pixel, compute how adding a plot to that pixel would decrease
#     the mean pixel-to-nearest-plot distance, then choose the pixel that
#     yields the largest decrease in mean distance. As a result, plots
#     concentrate in common areas structural space.
#
# @return vector of pixel index values for proposed new plots
# 
# @import proxy
# 
meanmin_select <- function(r_pca, p_new_pca = NULL, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[,1:n_pca, drop = FALSE])

  if (!is.null(p_new_pca)) {
    p_new_pca <- as.matrix(p_new_pca[,1:n_pca, drop = FALSE])
  } else {
    p_new_pca <- r_pca
  }

  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca[, 1:n_pca, drop = FALSE])
  }
  
  # Initialize minimum distances to existing plots
  if (!is.null(p_pca)) {
    d_exist <- as.matrix(proxy::dist(r_pca, p_pca))
    min_dists <- apply(d_exist, 1, min)
  } else {
    min_dists <- rep(Inf, nrow(r_pca))
  }

  # Prepare candidate index mapping to original input rows
  cand_idx <- seq_len(nrow(p_new_pca))
  sel_idx <- c()

  # Initialize selected indices
  for (i in seq_len(n_plots)) {
    n_candidates <- nrow(p_new_pca)
    if (n_candidates == 0) {
      warning("No potential candidates for plot ", i)
      break
    }

    # Compute mean distance if we added each candidate
    mean_after <- sapply(seq_len(n_candidates), function(rel_j) {
      cand_vec <- p_new_pca[rel_j, , drop = TRUE]
      # vectorised Euclidean distance from all pixels to candidate
      cand_mat <- matrix(cand_vec, nrow = nrow(r_pca), ncol = ncol(r_pca), byrow = TRUE)
      new_min <- pmin(min_dists, sqrt(rowSums((r_pca - cand_mat)^2)))
      mean(new_min)
    })

    # choose best candidate (minimizes mean distance); break ties randomly
    best_rel <- sample(which(mean_after == min(mean_after)), 1)
    best_orig_idx <- cand_idx[best_rel]  # original index in input p_new_pca
    sel_idx <- c(sel_idx, best_orig_idx)

    # Update min_dists using the chosen candidate
    chosen_vec <- p_new_pca[best_rel, , drop = TRUE]
    chosen_mat <- matrix(chosen_vec, nrow = nrow(r_pca), ncol = ncol(r_pca), byrow = TRUE)
    d_chosen <- sqrt(rowSums((r_pca - chosen_mat)^2))
    min_dists <- pmin(min_dists, d_chosen)

    # Remove chosen candidate from pool (and update cand_idx mapping)
    p_new_pca <- if (n_candidates > 1) p_new_pca[-best_rel, , drop = FALSE] else matrix(nrow = 0, ncol = n_pca)
    cand_idx <- cand_idx[-best_rel]

    if (nrow(p_new_pca) == 0) {
      if (i < n_plots) warning("No suitable locations for plot ", i + 1)
      break
    }
  }

  return(sel_idx)
}


# Select candidate plots using Latin Hypercube Sampling (LHS)
#
# @param r_pca PCA scores of landscape pixels in structural space
# @param p_new_pca Optional. PCA scores of candidate plots in structural space.
#     If supplied, selection is restricted to these candidates; otherwise the
#     full set of landscape pixels (`r_pca`) is used.
# @param p_pca Optional. PCA scores of existing plots in structural space; used
#     only when enforcing a minimum spatial separation (`min_dist`) so existing
#     plots are treated as already-selected for spatial exclusion.
# @param n_plots Number of plots to select.
# @param n_pca Number of PCA axes used in the sampling space.
# @param min_dist Minimum spatial distance (in map units) required between any
#     newly selected plot and existing/previously selected plots. If `NULL`
#     no spatial exclusion is applied.
#
# @details
# The Latin Hypercube Sampling (LHS) algorithm divides the n_pca-dimensional
#     PCA space into equal-probability intervals along each axis and generates
#     a design that samples one value per interval for each axis, producing
#     target points that ensure uniform coverage of multivariate space. For
#     each LHS target point the algorithm selects the landscape pixel (or a
#     candidate pixel from `p_new_pca` if provided) with PCA values closest to
#     the target. If `min_dist` is supplied and coordinate attributes are
#     attached to the PCA objects, selected plots are constrained to be at
#     least `min_dist` apart (and away from existing plots provided in
#     `p_pca`). The function proceeds greedily through the LHS targets,
#     removing chosen pixels from the candidate pool so that each target is
#     matched to a unique location.
#
# @return Integer vector of pixel index values (rows of `r_pca`) for the
#     proposed new plots, ordered by the sequence in which they were selected.
#
# @import lhs
# @import FNN
# 
lhs_select <- function(r_pca, p_new_pca = NULL, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[, 1:n_pca, drop = FALSE])
  
  if (!is.null(p_new_pca)) {
    p_new_pca <- as.matrix(p_new_pca[, 1:n_pca, drop = FALSE])
  } else {
    p_new_pca <- r_pca
  }
  
  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca[, 1:n_pca, drop = FALSE])
  }
  
  # Determine range for each PCA axis across all available pixels
  pca_ranges <- apply(r_pca, 2, range)
  
  # Create strata for each dimension
  strata_breaks <- lapply(1:n_pca, function(dim) {
    seq(pca_ranges[1, dim], pca_ranges[2, dim], length.out = n_plots + 1)
  })
  
  # Assign each candidate to strata
  candidate_strata <- matrix(0, nrow = nrow(p_new_pca), ncol = n_pca)
  for (dim in 1:n_pca) {
    candidate_strata[, dim] <- findInterval(
      p_new_pca[, dim], strata_breaks[[dim]], rightmost.closed = TRUE)

    # Ensure values are in range 1:n_plots
    candidate_strata[, dim] <- pmax(1, pmin(n_plots, candidate_strata[, dim]))
  }
  
  # Track which strata have been used in each dimension
  used_strata <- vector("list", n_pca)
  for (dim in 1:n_pca) {
    used_strata[[dim]] <- integer(0)
  }
  
  # If existing plots provided, mark their strata as used
  if (!is.null(p_pca)) {
    for (dim in 1:n_pca) {
      existing_strata <- findInterval(
        p_pca[,dim], strata_breaks[[dim]], rightmost.closed = TRUE)
      existing_strata <- pmax(1, pmin(n_plots, existing_strata))
      used_strata[[dim]] <- unique(c(used_strata[[dim]], existing_strata))
    }
  }
  
  # Initialize
  cand_idx <- seq_len(nrow(p_new_pca))
  sel_idx <- c()
  for (i in seq_len(n_plots)) {
    if (length(cand_idx) == 0) {
      warning("No suitable locations for plot ", i)
      break
    }
    
    # Score each candidate based on how many unused strata it occupies
    scores <- numeric(length(cand_idx))
    for (j in seq_along(cand_idx)) {
      orig_idx <- cand_idx[j]
      # Count how many dimensions have unused strata for this candidate
      unused_count <- 0
      for (dim in 1:n_pca) {
        if (!(candidate_strata[orig_idx, dim] %in% used_strata[[dim]])) {
          unused_count <- unused_count + 1
        }
      }
      scores[j] <- unused_count
    }
    
    # Select candidate(s) with highest score (most unused strata)
    best_candidates <- which(scores == max(scores))
    
    # If tie, choose candidate closest to stratum centers
    if (length(best_candidates) > 1) {
      center_dists <- sapply(best_candidates, function(rel_idx) {
        orig_idx <- cand_idx[rel_idx]
        dist_sum <- 0
        for (dim in 1:n_pca) {
          stratum <- candidate_strata[orig_idx, dim]
          center <- mean(strata_breaks[[dim]][c(stratum, stratum + 1)])
          dist_sum <- dist_sum + (p_new_pca[orig_idx, dim] - center)^2
        }
        sqrt(dist_sum)
      })
      best_rel <- best_candidates[which.min(center_dists)]
    } else {
      best_rel <- best_candidates[1]
    }
    
    best_idx <- cand_idx[best_rel]
    sel_idx <- c(sel_idx, best_idx)
    
    # Update used strata
    for (dim in 1:n_pca) {
      used_strata[[dim]] <- unique(c(used_strata[[dim]], 
        candidate_strata[best_idx, dim]))
    }
    
    # Remove selected candidate
    cand_idx <- cand_idx[-best_rel]
  }
  
  return(sel_idx)
}

# Select candidate plots using the k-means centroids method
#
# @param r_pca PCA scores of landscape pixels in structural space.
# @param p_new_pca Optional. PCA scores of candidate plots in structural space.
#     If supplied, selection is restricted to these candidates; otherwise the
#     full set of landscape pixels (`r_pca`) is used.
# @param p_pca Optional. PCA scores of existing plots in structural space. If
#     supplied, these plots are treated as fixed and excluded from clustering
#     to prevent duplication.
# @param n_plots Number of plots (clusters) to select.
# @param n_pca Number of PCA axes used in the clustering space.
#
# @details
# The k-means centroids method partitions the PCA space of all available pixels
#     into `n_plots` clusters and selects the pixel nearest to each cluster
#     centroid as a representative location. This approach samples regions of
#     high data density proportionally, producing a design that represents
#     typical conditions rather than extremes. When `p_pca` is provided, those
#     plots are excluded from clustering so new selections cover the remaining
#     space.
#
# @return Integer vector of pixel index values (rows of `r_pca`) for the
#     proposed new plots, ordered arbitrarily by cluster.
#
# @import stats
# @import FNN
# 
kmeans_select <- function(r_pca, p_new_pca = NULL, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[, 1:n_pca, drop = FALSE])
  
  # Determine candidate pool
  if (!is.null(p_new_pca)) {
    p_new_pca <- as.matrix(p_new_pca[, 1:n_pca, drop = FALSE])
  } else {
    p_new_pca <- r_pca
  }
  
  # Optionally exclude existing plots from candidates
  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca[, 1:n_pca, drop = FALSE])
    # Remove duplicates: any candidate exactly matching existing plot
    keep <- !apply(p_new_pca, 1, function(row)
      any(apply(p_pca, 1, function(p) all(abs(p - row) < .Machine$double.eps^0.5))))
    p_new_pca <- p_new_pca[keep, , drop = FALSE]
  }
  
  # Handle small candidate pool
  if (nrow(p_new_pca) < n_plots) {
    warning("Fewer candidates than requested plots; returning all candidates.")
    return(seq_len(nrow(p_new_pca)))
  }
  
  # Run k-means clustering on candidate PCA scores
  km <- kmeans(p_new_pca, centers = n_plots, nstart = 10)
  
  # Find candidate closest to each centroid
  centers <- km$centers
  nn <- FNN::get.knnx(p_new_pca, centers, k = 1)
  
  # Convert to original indices
  sel_idx <- as.integer(nn$nn.index)
  
  return(sel_idx)
}


###################
# VISUALISE RESULTS
###################

# Create a map of landscape representativeness by plots
#
# @param r raster
# @param r_dist values of structural distance, e.g. as returned by `pca_dist()`
# @param p optional existing plot locations, either an sf object containing plot
#     polygons or points, or a dataframe with two columns `x` and `y`
#     describing coordinates in the same coordinate system as `r`
# @param p_new optional proposed plot locations, either an sf object containing
#     plot polygons or points, or a dataframe with two columns `x` and `y`
#     describing coordinates in the same coordinate system as `r`
#
# @import ggplot2 
# @import scico 
# @import sf 
# 
map_plot_vis <- function(r, r_dist, p = NULL, p_new = NULL) {

  # Convert plot locations to dataframe if not already
  if (!is.null(p)) {
    if (isSFType(p, c("POINT", "POLYGON"))) {
      p_old <- as.data.frame(st_coordinates(st_centroid(p)))
      names(p_old) <- c("x", "y")
    } else if (inherits(p, "data.frame") && all(c("x", "y") %in% colnames(p))) {
      p_old <- p
    } else {
      stop("p must either be an sf POINT/POLYGON object or a data.frame with columns 'x' and 'y'")
    }
    p_old$type <- "Existing plot"
  } else {
    p_old <- NULL
  }

  if (!is.null(p_new)) {
    if (isSFType(p_new, c("POINT", "POLYGON"))) {
      p_new <- as.data.frame(st_coordinates(st_centroid(p_new)))
      names(p_new) <- c("x", "y")
    } else if (inherits(p_new, "data.frame") && all(c("x", "y") %in% colnames(p_new))) {
      p_new <- p_new
    } else {
      stop("p_new must either be an sf POINT/POLYGON object or a data.frame with columns 'x' and 'y'")
    }
    p_new$type <- "Proposed plot"
  } else {
    p_new <- NULL
  }

  p_all <- rbind(p_old, p_new)

  # Prepare raster
  r_df <- as.data.frame(crds(r))
  names(r_df) <- c("x", "y")
  r_df$dist <- r_dist

  # Create map
  out <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = dist)) +
    scale_fill_scico(name = "Relative distance to nearest plot", palette = "bamako") +
    coord_equal() +
    theme_bw() +
    ggtitle("Landscape representativeness") +
    labs(
      x = NULL, 
      y = NULL) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12))

  # Optionally add plots
  if (!is.null(p_all)) { 
    # Define legend items
    colour_map <- c("Existing plot" = "#E74B5E", "Proposed plot" = "#5499DE")
    shape_map <- c("Existing plot" = "circle", "Proposed plot" = "circle")

    out <- out + geom_point(data = p_all, 
      aes(x = x, y = y, shape = type, colour = type)) +
    scale_colour_manual(name = NULL, 
      values = colour_map) + 
    scale_shape_manual(name = NULL, 
      values = shape_map)
  }

  # Return
  return(out)
}

# Create a PCA biplot of landscape representativeness
# 
# @param r_pca PCA scores of pixels in structural space, as returned by
#     `pca_landscape()$r_pca`
# @param p_pca PCA scores of existing plots in structural space, as returned
#     by `pca_landscape()$p_pca`
# @param p_new_pca PCA scores of proposed plots in structural space, as returned
#     by `pca_landscape()$p_pca`
# 
# @import ggplot2
# 
pca_plot_vis <- function(r_pca, p_pca = NULL, p_new_pca = NULL) {

  # Extract PCA values from r_pca
  pca_pt_r <- as.data.frame(r_pca[,paste0("PC", 1:2)])
  pca_pt_r$type <- "Landscape value" 

  # Optionally extract PCA values from p_pca
  if (!is.null(p_pca)) { 
    pca_pt_p <- as.data.frame(p_pca[,paste0("PC", 1:2)])
    pca_pt_p$type <- "Existing plot" 
  } else {
    pca_pt_p <- NULL
  }

  # Optionally extract PCA values from p_new_pca
  if (!is.null(p_new_pca)) { 
    pca_pt_n <- as.data.frame(p_new_pca[,paste0("PC", 1:2)])
    pca_pt_n$type <- "Proposed plot" 
  } else {
    pca_pt_n <- NULL
  }

  # Combine PCA values
  pca_pt_all <- rbind(pca_pt_r, pca_pt_p, pca_pt_n)
  pca_pt_all$type <- factor(pca_pt_all$type, 
    levels = c("Landscape value", "Existing plot", "Proposed plot"))

  # Define legend items 
  size_map <- c("Existing plot" = 3, "Proposed plot" = 3, "Landscape value" = 1)
  opacity_map <- c("Existing plot" = 1, "Proposed plot" = 1, "Landscape value" = 0.5)
  colour_map <- c("Existing plot" = "#E74B5E", "Proposed plot" = "#5499DE", "Landscape value" = "grey")
  shape_map <- c("Existing plot" = "circle", "Proposed plot" = "circle", "Landscape value" = "circle")

  # Create plot
  out <- ggplot() +
    geom_point(data = pca_pt_all, 
      aes(x = PC1, y = PC2, 
        size = type, shape = type, colour = type, alpha = type)) + 
    scale_size_manual(name = NULL, values = size_map) + 
    scale_colour_manual(name = NULL, values = colour_map) + 
    scale_shape_manual(name = NULL, values = shape_map) + 
    scale_alpha_manual(name = NULL, values = opacity_map) + 
    labs(x = "PC1", y = "PC2") + 
    theme_bw() +
    ggtitle("Structural space coverage") +
    theme(plot.title = element_text(hjust = 0.5))

  # Return
  return(out)
}

# Compare histograms of relative distance in structural space before and after
# adding new plots
#
# @param dist values of structural distance to nearest existing plot, as
#     returned by `pca_dist()`
# @param dist_new optional values of structural distance to nearest proposed
#     plot, as returned by `pca_dist()`
# 
# @import ggplot2
# @import dplyr
# @import tidyr
# 
hist_plot_vis <- function(dist, dist_new = NULL) {

  # Create dataframe
  df <- data.frame(
      orig = dist,
      with_new = dist_new) %>%
    pivot_longer(
      everything(), 
      names_to = "scenario", 
      values_to = "distance") %>% 
    mutate(scenario = factor(scenario, 
      levels = c("orig", "with_new"),
      labels = c("Existing plots", "Existing and proposed plots")))
  
  # Create histogram
  out <- ggplot(df, aes(x = distance, fill = scenario)) +
    geom_histogram(alpha = 0.7, position = "identity", bins = 50) +
    scale_fill_manual(name = NULL, 
      values = c("Existing plots" =  "#E74B5E", "Existing and proposed plots" = "#5499DE")) +
    theme_bw() +
    ggtitle("Landscape coverage by plots") +
    labs(
      x = "Relative distance in structural space to nearest plot",
      y = "Number of pixels") + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom")

  # Return
  return(out)
}
  
#######################################
# CLASSIFY PIXELS BY REPRESENTATIVENESS
#######################################

# Classify pixels by whether they are well-represented by existing (and
#     proposed) plots
#
# @param r_pca PCA scores of pixels in structural space, as returned by
#     `pca_landscape()$r_pca`
# @param p PCA scores of candidate plots in structural space, as returned
#     by `pca_landscape()$p_pca`
# @param p_new optional PCA scores of candidate plots in structural space, as
#     returned by `pca_landscape()$p_pca`
# @param n_pca number of PCA axes used in analysis
# @param ci quantile threshold for distances, e.g. 0.95 = 95th percentile of
#     distances among plots 
# 
landscape_classif <- function(r_pca, p, p_new = NULL, 
  n_pca = 3, ci = 0.95) {

  # Extract PCA scores from chosen PCs
  r_df <- as.data.frame(r_pca)[,paste0("PC", 1:n_pca)]
  p <- p[,paste0("PC", 1:n_pca)]

  if (!is.null(p_new)) { 
    p_new <- p_new[,paste0("PC", 1:n_pca)]
  }

  # Compute centroid (mean vector) and covariance matrix of the plots
  p_cent <- colMeans(p)
  p_cov <- cov(p)

  if (!is.null(p_new)) { 
    p_all <- rbind(p, p_new)
    p_all_cent <- colMeans(p_all)
    p_all_cov <- cov(p_all)
  }

  # Compute distances
  # Compute Euclidean distances for plots and pixels using existing and suggested plots
  p_dist <- sqrt(rowSums((p - matrix(p_cent, nrow(p), ncol(p), byrow = TRUE))^2))
  r_dist <- sqrt(rowSums((r_df - matrix(p_cent, nrow(r_df), ncol(r_df), byrow = TRUE))^2))

  if (!is.null(p_new)) { 
    p_all_dist <- sqrt(rowSums((p_all - matrix(p_all_cent, nrow(p_all), ncol(p_all), byrow = TRUE))^2))
    r_all_dist <- sqrt(rowSums((r_df - matrix(p_all_cent, nrow(r_df), ncol(r_df), byrow = TRUE))^2))
  }

  # Choose cutoff — e.g. 95th percentile of plot distances
  p_cutoff <- quantile(p_dist, ci)

  if (!is.null(p_new)) { 
    p_all_cutoff <- quantile(p_all_dist, ci)
  }

  # Flag pixels that are too far from plots
  px_out <- r_dist > p_cutoff

  if (!is.null(p_new)) { 
    px_all_out <- r_all_dist > p_all_cutoff
  }

  # Create output dataframe
  r_df$px_out <- px_out
  r_df$dist <- r_dist

  if (!is.null(p_new)) {
    r_df$px_all_out <- px_all_out
    r_df$group <- case_when(
      r_df$px_out == TRUE & r_df$px_all_out == TRUE ~ "Poorly represented",
      r_df$px_out == FALSE & r_df$px_all_out == FALSE ~ "Well-represented by existing plots",
      r_df$px_out == TRUE & r_df$px_all_out == FALSE ~ "Well-represented by existing and proposed plots",
      TRUE ~ NA_character_)
  } else {
    r_df$group <- case_when(
      r_df$px_out == TRUE ~ "Poorly represented",
      r_df$px_out == FALSE ~ "Well-represented by existing plots")
  }
  r_df$group <- factor(r_df$group, 
    levels = c(
      "Poorly represented", 
      "Well-represented by existing plots", 
      "Well-represented by existing and proposed plots"))

  out <- r_df[,c("dist", "group")]

  # Return
  return(out)
}

# Create a map of pixel classifications 
# 
# @param r raster
# @param classif dataframe with classifications of pixels in column `group`, as
#     returned by `landscape_classif()`
# 
# @import ggplot2 
# @import terra
# 
map_classif_vis <- function(r, classif) { 

  # Check type of r
  if (!inherits(r, "SpatRaster")) {
    stop("r must be of class SpatRaster")
  }

  # Prepare raster
  r_df <- as.data.frame(crds(r))
  names(r_df) <- c("x", "y")
  r_df <- cbind(r_df, classif)

  # Create map
  out <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = group)) +
    scale_fill_discrete(name = NULL, guide = guide_legend(nrow = 3)) +
    coord_equal() +
    theme_bw() +
    ggtitle("Landscape representativeness") +
    labs(
      x = NULL, 
      y = NULL) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12))

  # Return
  return(out)
}

# Create a PCA biplot of pixel classifications 
# 
# @param r_pca PCA scores of pixels in structural space, as returned by
#     `pca_landscape()$r_pca`
# @param classif dataframe with classifications of pixels in column `group`, as
#     returned by `landscape_classif()`
# 
# @import ggplot2
# 
pca_classif_vis <- function(r_pca, classif) {
  r_df <- cbind(as.data.frame(r_pca), classif)

  out <- ggplot(r_df, aes(PC1, PC2, colour = group)) +
    geom_point(alpha = 0.8) +
    scale_colour_discrete(name = NULL) +
    theme_bw() + 
    labs(x = "PC1", y = "PC2") + 
    theme(legend.position = "bottom")

  # Return
  return(out)
}

