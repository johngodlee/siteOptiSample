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
#' @importFrom proxy dist
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


