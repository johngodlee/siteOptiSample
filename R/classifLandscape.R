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

