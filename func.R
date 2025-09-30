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
#
# @return list containing: `r_pca`: PCA object from pixel values, and
#     optionally `p_pca`: plot values in PCA space
# 
pca_landscape <- function(r, p = NULL) {
  # Check type of r
  if (!inherits(r, "SpatRaster")) {
    stop("r must be of class SpatRaster")
  }

  # Convert raster to data frame for landscape
  r_df <- as.data.frame(r, xy = FALSE, na.rm = TRUE)
  
  # Select only metric columns 
  metric_cols <- names(r_df)
  
  # PCA to reduce dimensionality
  r_pca <- prcomp(r_df[,metric_cols], scale. = TRUE)
  
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

# Get nearest neighbour distances between two PCAs
# 
# @param x PCA scores 
# @param y PCA scores
# @param n_pca number of PCA axes used in analysis
# @param k number of nearest neighbours in y to return for each row in x
# @param method distance method passed to `proxy::dist()`
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
  dists_mat <- proxy::dist(x[,1:n_pca], y[,1:n_pca], method = method)
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
# @param p_new_pca PCA scores of candidate plots in structural space
# @param p_pca optional, PCA scores of existing plots in structural space
# @param n_plots maximum number of new plots to add
# @param n_pca number of PCA axes used in analysis
#
# @details The maximin algorithm iteratively proposes new plots in pixels with
#     structural attributes furthest away from the existing (and previously
#     proposed) plots.
#
# @return vector of pixel index values for proposed new plots, ordered by
#     largest increase in coverage 
#
# @import proxy
# 
maximin_select <- function(r_pca, p_new_pca, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[,1:n_pca, drop = FALSE])
  p_new_pca <- as.matrix(p_new_pca[,1:n_pca, drop = FALSE])
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
# @param p_new_pca PCA scores of candidate plots in structural space
# @param p_pca optional, PCA scores of existing plots in structural space
# @param n_plots maximum number of new plots to add
# @param n_pca number of PCA axes used in analysis
#
# @details
# The minimax algorithm aims to minimize the maximum distance of any pixel to
#     its nearest plot in PCA space.
#
# @return vector of pixel index values for proposed new plots, ordered by
# largest increase in coverage 
# 
# @import proxy
# 
minimax_select <- function(r_pca, p_new_pca, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[,1:n_pca, drop = FALSE])
  p_new_pca <- as.matrix(p_new_pca[,1:n_pca, drop = FALSE])
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
# @param p_new_pca PCA scores of candidate plots in structural space
# @param p_pca optional, PCA scores of existing plots in structural space
# @param n_plots maximum number of new plots to add
# @param n_pca number of PCA axes used in analysis
#
# @details
# The meanmin algorithm aims to minimize the average distance of all pixels to
#     their nearest plot in PCA space, selecting new plots that most reduce the
#     mean distance iteratively.
#
# @return vector of pixel index values for proposed new plots
# 
# @import proxy
# 
meanmin_select <- function(r_pca, p_new_pca, p_pca = NULL, n_plots = 10, n_pca = 3) {
  # Restrict to selected PCA axes
  r_pca <- as.matrix(r_pca[,1:n_pca, drop = FALSE])
  p_new_pca <- as.matrix(p_new_pca[,1:n_pca, drop = FALSE])
  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca[, 1:n_pca, drop = FALSE])
  }
  
  # Initialize minimum distances to existing plots
  if (!is.null(p_pca)) {
    min_dists <- apply(as.matrix(proxy::dist(r_pca, p_pca)), 1, min)
  } else {
    min_dists <- rep(Inf, nrow(r_pca))
  }
  
  # Initialize selected indices
  sel_idx <- c()
  for (i in seq_len(n_plots)) {
    n_candidates <- nrow(p_new_pca)
    if (n_candidates == 0) {
      warning("No potential candidates for plot ", i)
      break
    }
    
    # Compute mean distance if we added each candidate
    mean_after <- sapply(seq_len(n_candidates), function(j) {
      candidate <- matrix(p_new_pca[j,], nrow = nrow(r_pca), ncol = n_pca, byrow = TRUE)
      new_min <- pmin(min_dists, sqrt(rowSums((r_pca - candidate)^2)))
      mean(new_min)
    })
    
    # Select candidate that minimizes mean distance
    best_idx <- which.min(mean_after)
    sel_idx <- c(sel_idx, best_idx)
    
    # Update minimum distances
    candidate <- matrix(p_new_pca[best_idx,], nrow = nrow(r_pca), ncol = n_pca, byrow = TRUE)
    min_dists <- pmin(min_dists, sqrt(rowSums((r_pca - candidate)^2)))
    
    # Remove chosen candidate
    p_new_pca <- p_new_pca[-best_idx, , drop = FALSE]

    # Stop if no candidates remain
    if (nrow(p_new_pca) == 0) {
      warning("No suitable locations for plot ", i)
      break
    }
  }
  
  # Return
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
    colour_map <- c("Existing plot" = "red", "Proposed plot" = "blue")
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
  colour_map <- c("Existing plot" = "red", "Proposed plot" = "blue", "Landscape value" = "grey")
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
      values = c("Existing plots" = "red", "Existing and proposed plots" = "blue")) +
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
# @param method either "mahalanobis" or "euclidean" to specify the distance
#     method used
# 
landscape_classif <- function(r_pca, p, p_new = NULL, 
  n_pca = 3, ci = 0.95, method = "mahalanobis") {

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
  if (method == "mahalanobis") {
    # Compute Mahalanobis distances for plots and pixels using existing and suggested plots
    p_dist <- mahalanobis(p, p_cent, p_cov)
    r_dist <- mahalanobis(r_df, p_cent, p_cov)

    if (!is.null(p_new)) { 
      p_all_dist <- mahalanobis(p_all, p_all_cent, p_all_cov)
      r_all_dist <- mahalanobis(r_df, p_all_cent, p_all_cov)
    }
  } else if (method == "euclidean") { 
    # Compute Euclidean distances for plots and pixels using existing and suggested plots
    p_dist <- sqrt(rowSums((p - matrix(p_cent, nrow(p), ncol(p), byrow = TRUE))^2))
    r_dist <- sqrt(rowSums((r_df - matrix(p_cent, nrow(r_df), ncol(r_df), byrow = TRUE))^2))

    if (!is.null(p_new)) { 
      p_all_dist <- sqrt(rowSums((p_all - matrix(p_all_cent, nrow(p_all), ncol(p_all), byrow = TRUE))^2))
      r_all_dist <- sqrt(rowSums((r_df - matrix(p_all_cent, nrow(r_df), ncol(r_df), byrow = TRUE))^2))
    }
  } else {
    stop("Method must be either 'mahalanobis' or 'euclidean'")
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

