###################################
# EXTRACT METRICS AT PLOT LOCATIONS
###################################

#' Check if object is sf type
#'
#' @param x object to test
#' @param type optional character vector of acceptable sf geometry types
#'
#' @return logical
#' 
#' @import sf
#' 
#' @NoRd
#' 
isSFType <- function(x, type = NULL) {
  inherits(x, c("sf", "sfc")) && (is.null(type) |
    all(st_geometry_type(x, by_geometry = FALSE) %in% type))
}

#' Extract values from raster across polygons
#'
#' @param r SpatRaster object
#' @param p either an sf object containing plot polygons or points, or a
#'     dataframe with two columns `x` and `y` describing coordinates in the same
#'     coordinate system as `r`
#' @param fun function to summarize the extracted data by polygon geometry, see
#'     `terra::extract()`
#' @param ... additional parameters passed to `terra::extract()`
#'
#' @return dataframe containing extracted metrics for each polygon
#' 
#' @import terra
#' @import sf
#' 
#' @export
#' 
extractPlotMetrics <- function(r, p, fun = mean, ...) {
  # Check type of r
  if (!inherits(r, "SpatRaster")) {
    stop("r must be of class SpatRaster")
  }

  # Convert p to sf if not already
  if (isSFType(p, c("POINT", "POLYGON"))) {
    p_vect <- terra::vect(p)
  } else if (inherits(p, c("data.frame", "matrix")) && 
      all(c("x", "y") %in% colnames(p))) {
    p_vect <- sf::st_as_sf(as.data.frame(p), coords = c("x", "y"), crs = sf::crs(r))
  } else {
    stop("p must either be an sf POINT/POLYGON object or a data.frame with columns 'x' and 'y'")
  }

  # Extract metrics for each plot polygon
  out <- as.matrix(terra::extract(r, p_vect, fun = fun, 
    na.rm = TRUE, ID = FALSE, ...))
  
  # Return
  return(out)
}

#############################
# REPRESENTATIVENESS ANALYSIS
#############################

#' Run a PCA on landscape structural metrics, and optionally place plots in
#' that space
#' 
#' @param r SpatRaster with structural metrics
#' @param p optional dataframe containing structural metrics for each plot, e.g.
#'     as returned by `extractPlotMetrics()`
#' @param ... additional arguments passed to `prcomp()`  
#'
#' @return list containing: `r_pca`: PCA object from pixel values, and
#'     optionally `p_pca`: plot values in PCA space
#' 
#' @export
#' 
PCALandscape <- function(r, p = NULL, ...) {
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

#' Calculate euclidean distance between two matrices
#' 
#' @param x numeric matrix
#' @param y numeric matrix
#' 
#' @NoRd
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

#' Calculate mahalanobis distance between two matrices
#' 
#' @param x numeric matrix
#' @param y numeric matrix
#' 
#' @NoRd
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

#' Get nearest neighbour distances between two PCAs
#' 
#' @param x PCA scores 
#' @param y PCA scores
#' @param n_pca number of PCA axes used in analysis
#' @param k number of nearest neighbours in y to return for each row in x
#' @param method distance method, either "euclidean" or "mahalanobis"
#'
#' @return if `k = 1` a vector of distances to the nearest neighbour in `y` for
#'     each row in `x`. If `k = 2` a matrix with k columns with ordered nearest
#'     neighbour distances
#'
#' @import proxy
#' 
#' @export
#' 
pcaDist <- function(x, y, n_pca = 3, k = 1, method = "euclidean") {
  # Calculate nearest neighbor distances
  # For each landscape pixel, find distance to nearest plot
  if (method == "mahalanobis") {
    dists_mat <- mahanalobisDist(x[,1:n_pca], y[,1:n_pca])
  } else if (method == "euclidean") {
    dists_mat <- euclideanDist(x[,1:n_pca], y[,1:n_pca])
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

#######################################
# CLASSIFY PIXELS BY REPRESENTATIVENESS
#######################################

#' Classify pixels by whether they are well-represented by existing (and
#'     proposed) plots
#'
#' @param r_pca PCA scores of pixels in structural space, as returned by
#'     `PCALandscape()$r_pca`
#' @param p PCA scores of candidate plots in structural space, as returned
#'     by `PCALandscape()$p_pca`
#' @param p_new optional PCA scores of candidate plots in structural space, as
#'     returned by `PCALandscape()$p_pca`
#' @param n_pca number of PCA axes used in analysis
#' @param ci quantile threshold for distances, e.g. 0.95 = 95th percentile of
#'     distances among plots 
#' 
#' @export
#' 
classifLandscape <- function(r_pca, p, p_new = NULL, 
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

##############################################
# IDENTIFY GAPS AND OPTIMALLY LOCATE NEW PLOTS
##############################################

#' Recommend locations for additional plots
#' 
#' @param r SpatRaster with structural metrics
#' @param p optional, either an sf object containing plot polygons or points, or
#'     a dataframe with two columns `x` and `y` describing coordinates in the
#'     same coordinate system as `r`
#' @param n_plots maximum number of new plots to add
#' @param p_new_dim optional, dimensions of new plots in the same coordinate 
#'     system as `r`. Either a single value for square plots, or a vector of two
#'     values for rectangular plots. Should be perfectly divisible by the
#'     resolution of `r`
#' @param r_mask optional, a raster which defines a mask of potential plot
#'     locations, e.g. based on accessibility. 
#' @param pca logical, if TRUE (default) the variables in `r` are run through a
#'     PCA to reduce issues of collinearity among variables. If TRUE, provide
#'     `n_pca`
#' @param n_pca optional, the number of PCA axes used to calculate the distance
#'     between pixels and plots. If not provided and `pca` is TRUE, all PCA axes
#'     are used
#' @param method a function to optimally place plots, e.g. `meanminSelect()` or
#'     `minimaxSelect()`
#' 
#' @return sf dataframe with polygons of proposed new plots
#'
#' @import terra
#' @import sf
#' 
#' @export
#'
newPlotSelect <- function(r, p, n_plots, p_new_dim = NULL, r_mask = NULL, 
  pca = TRUE, n_pca = NULL, method = meanminSelect) {

  # If new plot dimensions not provided, set to raster resolution
  if (is.null(p_new_dim)) {
    p_new_dim <- res(r)
  }

  # If new plot dimensions are smaller than the raster resolution
  if (any(p_new_dim < res(r))) {
    warning("'p_new_dim' is smaller than the resolution of 'r'.\nSetting 'p_new_dim' to the resolution of 'r'\n",
      immediate. = TRUE)
    p_new_dim <- res(r)
  }

  if (any(p_new_dim %% res(r) != 0)) {
    p_new_dim <- floor(p_new_dim / res(r)) * res(r)
    warning("'p_new_dim' is not perfectly divisible by the resolution of 'r'.\nRounding 'p_new_dim' down a value which is evenly divisible: ", paste(p_new_dim, collapse = ", "))
  }

  # If new plot dimensions are just a single value, multiply it up
  if (length(p_new_dim) == 1) {
    p_new_dim <- rep(p_new_dim, 2)
  }

  # If mask for new plots is provided, mask the raster
  if (!is.null(r_mask)) { 
    r_mask <- terra::mask(r, r_mask)
  } else {
    r_mask <- r
  }

  # Convert p to spatial vector object 
  # Extract metrics for each plot polygon
  if (!is.null(p)) {
    old_ext <- extractPlotMetrics(r, p)
  } else {
    old_ext <- NULL
  }
  
  # Optional PCA to reduce dimensionality
  if (pca && is.null(n_pca)) {
   n_pca <- ncol(r)
  }

  if (pca) { 
    old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
    r_pca <- rep(r[[1]], n_pca)
    terra::values(r_pca)[complete.cases(terra::values(r)),] <- old_pca$r_pca$x[,1:n_pca]
    names(r_pca) <- colnames(old_pca$r_pca$x)[1:n_pca]
    if (!is.null(p)) {
      p_pca <- old_pca$p_pca[,1:n_pca]
    } else {
      p_pca <- NULL
    }
  } else {
    r_pca <- scale(r)
    if (!is.null(p)) {
      p_pca <- old_ext[,names(r_pca)]
    } else {
      p_pca <- NULL
    }
  }

  # Get index values in PCA for existing plots and potential new plots
  if (!is.null(p)) {
    old_ind <- which(complete.cases(terra::values(terra::mask(r, terra::vect(p)))))
  } else {
    old_ind <- NULL
  }

  # Initialise PCA values for candidate plots 
  new_ind <- which(complete.cases(terra::values(r_mask)))

  # Implement plot selection algorithm
  p_list <- method(
    r_pca = r_pca,
    p_pca = p_pca,
    old_ind = old_ind,
    new_ind = new_ind,
    n_plots = n_plots,
    p_new_dim = p_new_dim)

  # Prepare output sf object
  out <- sf::st_sf(geometry = do.call(c, p_list), crs = terra::crs(r))

  # Return
  return(out)
}


#' Iteratively add candidate plots using the meanmin algorithm
#' 
#' @param r_pca raster of PCA scores of pixels in structural space
#' @param p_pca PCA values of existing plots in structural space
#' @param old_ind optional, pixel IDs in `r_pca` specifying existing plot
#'     locations. 
#' @param new_ind optional, pixel IDs in `r_pca` specifying candidate plot
#'     locations. If not supplied, the full set of landscape pixels in `r_pca`
#'     is used.
#' @param n_plots maximum number of new plots to add
#' @param p_new_dim dimensions of new plots in the same coordinate system as
#'     `r`. A vector of two values. Should be perfectly divisible by the
#'     resolution of `r`.
#'
#' @details
#' The meanmin algorithm aims to minimize the mean distance of all pixels to
#'     their nearest plot in PCA space, iteratively placing new plots in pixels
#'     that most reduce the average pixel-to-nearest-plot distance. For each
#'     candidate pixel, compute how adding a plot to that pixel would decrease
#'     the mean pixel-to-nearest-plot distance, then choose the pixel that
#'     yields the largest decrease in mean distance. As a result, plots
#'     concentrate in common areas structural space.
#'
#' @return vector of pixel index values for proposed new plots
#' 
#' @import proxy
#' @import terra
#' @import sf
#' 
#' @export
#' 
meanminSelect <- function(r_pca, p_pca, old_ind, new_ind, n_plots, p_new_dim) { 
  # Initialise empty list to store new plot polygons
  p_list <- list()

  # For each plot 
  for (i in seq_len(n_plots)) { 
    message(i, "/", n_plots)
    # Get index values in PCA for candidate plots
    # Excluding pixels occupied by existing plots 
    new_ind <- setdiff(new_ind, old_ind)

    # Stop if no potential locations available
    if (length(new_ind) == 0) { 
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
      break
    }

    # Initialise minimum distances between each pixel and existing plots
    min_dists <- r_pca[[1]]
    if (!is.null(p_pca)) {
      terra::values(min_dists) <- apply(as.matrix(
        proxy::dist(terra::values(r_pca), p_pca)), 1, min)
    } else {
      terra::values(min_dists)[complete.cases(terra::values(min_dists))] <- Inf
    }

    # Calculate reduction of mean distance of pixels to nearest plot if each
    #     pixel was added
    # For each candidate pixel:
    mean_after <- sapply(new_ind, function(j) {
      cand_vec <- unlist(terra::values(r_pca)[j,])
      # Euclidean distance from all pixels to candidate
      cand_mat <- matrix(cand_vec, nrow = nrow(terra::values(r_pca)), ncol = ncol(terra::values(r_pca)), byrow = TRUE)
      new_min <- pmin(terra::values(min_dists), sqrt(rowSums((terra::values(r_pca) - cand_mat)^2)))
      mean(new_min, na.rm = TRUE) 
    })
    r_cand <- r[[1]]
    terra::values(r_cand) <- NA_real_
    terra::values(r_cand)[new_ind,] <- mean_after

    # Create a moving window to compute mean distance of other subplots
    p_width_x <- p_new_dim[1]
    p_height_y <- p_new_dim[2]
    res_x <- res(r)[1]
    res_y <- res(r)[2]

    # Number of cells needed in half-window
    n_x <- ceiling((p_width_x / res_x) / 2)
    n_y <- ceiling((p_height_y / res_y) / 2)

    # Full window dimensions
    full_cols <- (2 * n_x) + 1 
    full_rows <- (2 * n_y) + 1

    # Size of the new plot in raster cells
    p_x <- round(p_new_dim[1] / res_x)
    p_y <- round(p_new_dim[2] / res_y)
    P <- c(p_x, p_y)

    # Compute symmetric padding (handles odd differences)
    pad_left_x  <- floor((full_cols - P[1]) / 2)
    pad_right_x <- ceiling((full_cols - P[1]) / 2)
    pad_left_y  <- floor((full_rows - P[2]) / 2)
    pad_right_y <- ceiling((full_rows - P[2]) / 2)

    # Create centered window
    w_matrix <- matrix(NA, nrow = full_rows, ncol = full_cols)

    # Placement indices
    start_row <- pad_left_y + 1
    end_row   <- full_rows - pad_right_y
    start_col <- pad_left_x + 1
    end_col   <- full_cols - pad_right_x

    # Insert centered block
    w_matrix[start_row:end_row, start_col:end_col] <- 1

    # Run focal operation
    r_rs <- terra::focal(r_cand, w = w_matrix, fun = "mean")

    # Stop if there's no potential locations left
    if (all(!is.finite(values(r_rs)))) {
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
      break
    }

    # Extract coordinates of selected plot
    # In the case of ties, pick first 
    min_val <- min(terra::values(r_rs), na.rm = TRUE)
    which_min_val <- which(terra::values(r_rs) == min_val)
    rc <- terra::rowColFromCell(r_rs, which_min_val)
    row0 <- rc[1]
    col0 <- rc[2]

    # Dimensions of window
    h <- nrow(w_matrix)
    w <- ncol(w_matrix)
    pos <- which(w_matrix == 1, arr.ind = TRUE)

    # The cell in the window that corresponds to the focal raster cell:
    # It is ALWAYS the middle cell of the full window matrix
    row_center <- ceiling(h / 2)
    col_center <- ceiling(w / 2)

    # Compute row/col offsets FROM THE CENTER
    row_offsets <- pos[,1] - row_center
    col_offsets <- pos[,2] - col_center

    # Apply offsets to the focal cell
    sel_rows <- row0 + row_offsets
    sel_cols <- col0 + col_offsets

    # Convert to raster cell numbers
    sel_id <- terra::cellFromRowCol(r_rs, sel_rows, sel_cols)
    sel_id <- sel_id[!is.na(sel_id)]

    # Extract polygon of selected plot
    r_sub <- r_pca[[1]]
    terra::values(r_sub) <- NA
    r_sub[sel_id] <- 1
    p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
    p_sf <- sf::st_geometry(sf::st_as_sf(p_sub))

    # Add polygon of selected plot to list
    p_list[[i]] <- p_sf

    # Update cell IDs of selected plots
    old_ind <- c(old_ind, sel_id)
  }

  # Return
  return(p_list)
}


#' Iteratively add candidate plots using the minimax algorithm
#' 
#' @param r_pca raster of PCA scores of pixels in structural space
#' @param p_pca PCA values of existing plots in structural space
#' @param old_ind optional, pixel IDs in `r_pca` specifying existing plot
#'     locations. 
#' @param new_ind optional, pixel IDs in `r_pca` specifying candidate plot
#'     locations. If not supplied, the full set of landscape pixels in `r_pca`
#'     is used.
#' @param n_plots maximum number of new plots to add
#' @param p_new_dim dimensions of new plots in the same coordinate system as
#'     `r`. A vector of two values. Should be perfectly divisible by the
#'     resolution of `r`.
#'
#' @details 
#' The minimax algorithm aims to minimise the maximum distance between
#'     proposed plots and existing plots, iteratively placing new plots in
#'     pixels with structural attributes most dissimilar to existing plots. For
#'     each candidate pixel, compute the distance to the nearest existing plot,
#'     then choose the pixel with the highest distance value. As a result, plots
#'     occupy structural extremes.
#'
#' @return vector of pixel index values for proposed new plots
#' 
#' @import proxy
#' @import terra
#' @import sf
#' 
#' @export
#' 
minimaxSelect <- function(r_pca, p_pca, old_ind, new_ind, n_plots, p_new_dim) { 
  # Initialise empty list to store new plot polygons
  p_list <- list()

  # For each plot 
  for (i in seq_len(n_plots)) { 
    message(i, "/", n_plots)
    # Get index values in PCA for candidate plots
    # Excluding pixels occupied by existing plots 
    new_ind <- setdiff(new_ind, old_ind)

    # Stop if no potential locations available
    if (length(new_ind) == 0) { 
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
      break
    }

    # Add PCA values for new plots to p_pca
    p_pca <- values(r_pca)[old_ind,]

    # Initialise minimum distances between each pixel and existing plots
    r_cand <- r_pca[[1]]
    if (!is.null(p_pca)) {
      terra::values(r_cand) <- apply(as.matrix(
        proxy::dist(terra::values(r_pca), p_pca)), 1, min)
    } else {
      terra::values(r_cand)[complete.cases(terra::values(r_cand))] <- Inf
    }

    # Create a moving window to compute mean distance of other subplots
    p_width_x <- p_new_dim[1]
    p_height_y <- p_new_dim[2]
    res_x <- res(r)[1]
    res_y <- res(r)[2]

    # Number of cells needed in half-window
    n_x <- ceiling((p_width_x / res_x) / 2)
    n_y <- ceiling((p_height_y / res_y) / 2)

    # Full window dimensions
    full_cols <- (2 * n_x) + 1 
    full_rows <- (2 * n_y) + 1

    # Size of the new plot in raster cells
    p_x <- round(p_new_dim[1] / res_x)
    p_y <- round(p_new_dim[2] / res_y)
    P <- c(p_x, p_y)

    # Compute symmetric padding (handles odd differences)
    pad_left_x  <- floor((full_cols - P[1]) / 2)
    pad_right_x <- ceiling((full_cols - P[1]) / 2)
    pad_left_y  <- floor((full_rows - P[2]) / 2)
    pad_right_y <- ceiling((full_rows - P[2]) / 2)

    # Create centered window
    w_matrix <- matrix(NA, nrow = full_rows, ncol = full_cols)

    # Placement indices
    start_row <- pad_left_y + 1
    end_row   <- full_rows - pad_right_y
    start_col <- pad_left_x + 1
    end_col   <- full_cols - pad_right_x

    # Insert centered block
    w_matrix[start_row:end_row, start_col:end_col] <- 1

    # Run focal operation
    r_rs <- terra::focal(r_cand, w = w_matrix, fun = "mean")

    # Find pixel with largest value
    max_val <- max(terra::values(r_rs, na.rm = TRUE))
    which_max_val <- which(terra::values(r_rs) == max_val)

    # Extract coordinates of selected plot
    rc <- terra::rowColFromCell(r_rs, which_max_val)
    row0 <- rc[1]
    col0 <- rc[2]

    # Dimensions of window
    h <- nrow(w_matrix)
    w <- ncol(w_matrix)
    pos <- which(w_matrix == 1, arr.ind = TRUE)

    # The cell in the window that corresponds to the focal raster cell:
    # It is ALWAYS the middle cell of the full window matrix
    row_center <- ceiling(h / 2)
    col_center <- ceiling(w / 2)

    # Compute row/col offsets FROM THE CENTER
    row_offsets <- pos[,1] - row_center
    col_offsets <- pos[,2] - col_center

    # Apply offsets to the focal cell
    sel_rows <- row0 + row_offsets
    sel_cols <- col0 + col_offsets

    # Convert to raster cell numbers
    sel_id <- terra::cellFromRowCol(r_rs, sel_rows, sel_cols)
    sel_id <- sel_id[!is.na(sel_id)]

    # Extract polygon of selected plot
    r_sub <- r_pca[[1]]
    terra::values(r_sub) <- NA
    r_sub[sel_id] <- 1
    p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
    p_sf <- sf::st_geometry(sf::st_as_sf(p_sub))

    # Add polygon of selected plot to list
    p_list[[i]] <- p_sf

    # Update cell IDs of selected plots
    old_ind <- c(old_ind, sel_id)

  }

  # Return
  return(p_list)
}

