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


