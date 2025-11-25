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

