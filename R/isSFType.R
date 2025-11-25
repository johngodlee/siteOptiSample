#' Check if object is sf type
#'
#' @param x object to test
#' @param type optional character vector of acceptable sf geometry types
#'
#' @return logical
#' 
#' @import sf
#' 
#' @noRd
#' 

isSFType <- function(x, type = NULL) {
  inherits(x, c("sf", "sfc")) && (is.null(type) |
    all(st_geometry_type(x, by_geometry = FALSE) %in% type))
}

