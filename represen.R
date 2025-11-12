# How representative of forest structure are plots in a site?
# John L. Godlee (johngodlee@gmail.com)
# Last updated: 2025-09-19

# Purpose:
# To what extent are existing plots representative of the broader landscape forest structure? 
# Identify optimal locations for new plots to fill gaps in multidimensional forest structure space.

# Packages 
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(terra)
library(sf)
library(scico)  # Color palettes
library(proxy)  # Fast distances
library(GGally)  # Correlation plot

source("./func.R")

#############
# IMPORT DATA
#############

# Import ALS metrics
als <- rast("./dat/als/san_lorenzo_predictors_50m.tiff")

# Import 50x50 m plot polygons
subplots <- st_read("./dat/plots/sub_summ.gpkg")

# Filter subplots to ALS data
als_poly <- st_as_sf(as.polygons(als, na.rm = TRUE))
subplots_als <- st_join(subplots, als_poly, join = st_intersects) %>% 
  filter(!is.na(pulse_density)) 

# Only include "bottomleft" subplots
p <- subplots_als %>% 
  filter(grepl("bottomleft", subplot_id))

# Select relevant metrics from ALS data
# names(als)
r <- als[[c(
  "zmax.canopy.vox",
  "zmean.canopy.vox",
  "zcv.canopy.vox",
  "zskew.canopy.vox",
  "zkurt.canopy.vox",
  "zq10.canopy.vox",
  "zq30.canopy.vox",
  "zq50.canopy.vox",
  "zq70.canopy.vox",
  "zq90.canopy.vox",
  "lad_cv.aboveground.vox",
  "lad_sum.aboveground.vox")]]

# Check for collinearity
als_vals <- as.data.frame(r)

# Compute correlation between layers
# Visualize correlation
ggcorr(als_vals, label = TRUE)
# if very collinear maybe just use biplot of metrics rather than a PCA

########################### 
# TEST INDIVIDUAL FUNCTIONS
########################### 

# Extract metrics from raster for each existing plot
old_ext <- extractPlotMetrics(r, p, fun = mean)
# Returns a dataframe of values in the same row order as p
# Optional, could also use the computed values for each subplot, 
# rather than extracting point estimates from a grid.

# Put pixels and existing plots in the same PCA space
old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
# Returns two dataframes of PCA scores, for pixels (r_pca) and plots (p_pca)
# r_pca is a full PCA object
# p_pca is only PCA scores
# PCA is run on variables scaled to unit variance shifted to center on zero

# Calculate distance from each pixel to nearest plot in PCA space
old_dist_euclidean <- pcaDist(old_pca$r_pca$x, old_pca$p_pca, 
  n_pca = 3, k = 1, method = "euclidean")

old_dist_mahalanobis <- pcaDist(old_pca$r_pca$x, old_pca$p_pca, 
  n_pca = 3, k = 1, method = "mahalanobis")
# Returns vector of distances for each pixel in r_pca

old_dist_list <- list(
  "euclidean" = old_dist_euclidean,
  "mahalanobis" = old_dist_mahalanobis)

# Compare Euclidean and Mahalanobis distances
data.frame(
  euclidean = c(old_dist_euclidean), 
  mahalanobis = c(old_dist_mahalanobis)) %>% 
  ggplot(., aes(x = euclidean, y = mahalanobis)) + 
  geom_abline(linetype = 2, colour = "red") + 
  geom_point() + 
  geom_smooth(colour = "blue") + 
  geom_smooth(method = "lm", colour = "green") + 
  theme_classic()

data.frame(
  euclidean = c(old_dist_euclidean), 
  mahalanobis = c(old_dist_mahalanobis)) %>% 
  pivot_longer(everything()) %>% 
  ggplot(., aes(x = value, fill = name)) + 
    geom_histogram(position = "identity", alpha = 0.5) + 
    theme_classic()

# Create map of the dissimilarity of pixels from existing plots
visMap_out_list <- lapply(names(old_dist_list), function(x) {  
  visMap(r, old_dist_list[[x]], p = p) + 
    scale_fill_scico(
      name = paste0("Relative ", x, " distance to nearest plot"), 
      palette = "bamako",
      limits = c(0, 8),
      oob = scales::squish
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        barwidth = 20,
        barheight = 1
      )
    )
})
names(visMap_out_list) <- names(old_dist_list)

lapply(names(visMap_out_list), function(x) { 
  ggsave(visMap_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/visMap_out_", x, ".png"))
})

# Create plot of existing plots within landscape PCA space 
visPCA_out <- visPCA(old_pca$r_pca$x, old_pca$p_pca)

ggsave(visPCA_out, width = 6, height = 4, 
  file = "./img/visPCA_out.png")

# Classify pixels by how well-represented they are by the plots
old_pca_classif <- classifLandscape(
  r_pca = old_pca$r_pca$x, 
  p = old_pca$p_pca, 
  n_pca = 3, 
  ci = 0.95
)
# Returns dataframe with the distance of each pixel and its classification

# Create map of pixel classification 
visMapClassif_out <- visMapClassif(r, old_pca_classif)

ggsave(visMapClassif_out, width = 8, height = 6, file = 
  "./img/visMapClassif_out.png")

# Create PCA biplot of pixel classification 
visPCAClassif_out <- visPCAClassif(old_pca$r_pca$x, old_pca_classif)

ggsave(visPCAClassif_out, width = 6, height = 4, file = 
  "./img/visPCAClassif_out.png")

# Extract candidate plot locations (50% sample of pixels in the landscape)
# Excluding existing plots
# XY coordinates of each candidate location
subplots_cent <- p %>% 
  st_centroid() %>% 
  vect() %>% 
  geom(.) %>% 
  as.data.frame() %>%
  dplyr::select(x, y) %>% 
  cellFromXY(r, .)
r_noplot <- r
r_noplot[subplots_cent] <- NA
cand_cds <- as.data.frame(crds(r_noplot))
cand <- cand_cds[sample(seq_len(nrow(cand_cds)), floor(nrow(cand_cds) / 2)),]
# This could also be, for example, locations based on accessibility within the site
# or existing plots that could be re-measured

# Map the candidate plots (pixels)
visMap_out_list$euclidean + 
  geom_point(data = cand, aes(x = x, y = y), colour = "blue")

# Extract raster values for candidate plot locations
cand_ext <- extractPlotMetrics(r, cand)

# Put candidate plots in PCA space
cand_pca <- PCALandscape(r, cand_ext, center = TRUE, scale. = TRUE)
# Returns PCA objects for pixels and candidate plots

# Check original landscape PCA should be identical to new landscape PCA
stopifnot(all(old_pca$r_pca$x == cand_pca$r_pca$x))

# Create a PCA biplot to check data
ggplot() + 
  geom_point(data = cand_pca$r_pca$x, aes(x = PC1, y = PC2), 
    size = 2) + 
  geom_point(data = cand_pca$p_pca, aes(x = PC1, y = PC2), 
    colour = "red", size = 0.5) + 
  geom_point(data = old_pca$p_pca, aes(x = PC1, y = PC2), 
    colour = "green", size = 0.5) + 
  theme_bw()
# Candidate plots are a subset of the original pixels

# Optimally locate new plots within landscape based on PCA scores using
# different methods
opt_maximin <- maximinSelect(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_minimax <- minimaxSelect(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_meanmin <- meanminSelect(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_maxilhs <- lhsSelect(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_kmeans <- kmeansSelect(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

# Create list of results
opt_list <- list(
  maximin = opt_maximin,
  minimax = opt_minimax,
  meanmin = opt_meanmin,
  maxilhs = opt_maxilhs,
  kmeans = opt_kmeans)

# Check no selections are duplicated
lapply(opt_list, function(x) { 
  stopifnot(all(!duplicated(x)))
})

# Create PCA plots to check results are valid for each method
wrap_plots(lapply(names(opt_list), function(x) { 

  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  # Combine PCA scores of existing and proposed plots
  old_cand_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel)

  # Create a PCA biplot to check data
  ggplot() + 
    geom_point(data = cand_pca$r_pca$x, aes(x = PC1, y = PC2), size = 2) + 
    geom_point(data = old_pca$p_pca, aes(x = PC1, y = PC2), colour = "green", size = 0.6) + 
    geom_point(data = cand_pca_sel, aes(x = PC1, y = PC2), colour = "red", size = 0.5) + 
    theme_bw() + 
    ggtitle(x)
}))
# Candidate plots are a subset of the original pixels 
# for all methods

# Create maps with results for each method
visMap_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  # Combine PCA scores of existing and proposed plots
  old_cand_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel)

  # Analyse representativeness of existing and proposed plots
  old_cand_dist <- pcaDist(cand_pca$r_pca$x, old_cand_pca_sel, n_pca = 3, k = 1)

  # Create map of the dissimilarity of pixels from proposed and existing plots 
  visMap(
    r = r, 
    r_dist = old_cand_dist, 
    p = p, 
    p_new = cand[opt_list[[x]],]
  ) +
    scale_fill_scico(
      name = "Relative distance to nearest plot", 
      palette = "bamako",
      limits = c(0, 8),
      oob = scales::squish
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        barwidth = 20, 
        barheight = 1
      )
    )
})
names(visMap_new_out_list) <- names(opt_list)

lapply(names(visMap_new_out_list), function(x) { 
  ggsave(visMap_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/visMap_new_out_", x, ".png"))
})

# Create biplots of proposed and existing plots within landscape PCA space 
visPCA_new_out_list <- lapply(names(opt_list), function(x) { 

  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  visPCA(
    r_pca = cand_pca$r_pca$x, 
    p_pca = as.data.frame(old_pca$p_pca),
    p_new_pca = as.data.frame(cand_pca_sel))
})
names(visPCA_new_out_list) <- names(opt_list)

lapply(names(visPCA_new_out_list), function(x) { 
  ggsave(visPCA_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/visPCA_new_out_", x, ".png"))
})

# Create histograms of pixel-plot dissimilarity with and without proposed plots
visHist_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  # Combine PCA scores of existing and proposed plots
  old_cand_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel)

  # Analyse representativeness of existing and proposed plots
  old_cand_dist <- pcaDist(cand_pca$r_pca$x, old_cand_pca_sel, n_pca = 3, k = 1)

  # Create plot
  visHist(old_dist_list$euclidean, old_cand_dist)
})
names(visHist_new_out_list) <- names(opt_list)

lapply(names(visHist_new_out_list), function(x) { 
  ggsave(visHist_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/visHist_new_out_", x, ".png"))
})

# Create maps which classify pixels by how well-represented they are by the
# existing and new plots
visMapClassif_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  cand_classif <- classifLandscape(
    r_pca = old_pca$r_pca$x, 
    p = old_pca$p_pca, 
    p_new = cand_pca_sel,
    n_pca = 3, 
    ci = 0.95
  )

  # Create map of pixel classification 
  visMapClassif_new_out <- visMapClassif(r, cand_classif)
})
names(visMapClassif_new_out_list) <- names(opt_list)

lapply(names(visMapClassif_new_out_list), function(x) { 
  ggsave(visMapClassif_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/visMapClassif_new_out_", x, ".png"))
})

# Create PCA biplots of pixel classification 
visPCAClassif_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  cand_classif <- classifLandscape(
    r_pca = old_pca$r_pca$x, 
    p = old_pca$p_pca, 
    p_new = cand_pca_sel,
    n_pca = 3, 
    ci = 0.95
  )

  visPCAClassif_new_out <- visPCAClassif(old_pca$r_pca$x, cand_classif) + 
    guides(colour = guide_legend(nrow = 2))
})
names(visPCAClassif_new_out_list) <- names(opt_list)

lapply(names(visPCAClassif_new_out_list), function(x) { 
  ggsave(visPCAClassif_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/visPCAClassif_new_out_", x, ".png"))
})


###############################################
# TRY COMBINED WORKFLOW WITH TOP-LEVEL FUNCTION
###############################################


# Run function
wout <- newPlotSelect(
  r = r,
  p = p,
  p_new_dim = 100,
  n_plots = 3,
  n_pca = 3,
  pcaDist_method = "euclidean",
  resample_fun = mean,
  select_fun = meanminSelect)

# Check visualisation functions can run with outputs from wrapper function

# Create PCA biplot of proposed and existing plots
visPCA(
  r_pca = wout$r_pca$x,
  p_pca = wout$p_pca,
  p_new_pca = wout$p_new_pca)

# Create map of the dissimilarity of pixels from proposed and existing plots 
wrap_plots(
  visMap(
    r = r,
    r_dist = wout$r_new_dist, 
    p = p, 
    p_new = wout$p_new_coords) +
    scale_fill_scico(
      name = "Relative distance to nearest plot", 
      palette = "bamako") +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        barwidth = 20,
        barheight = 1)),
  visMap(
    r = r,
    r_dist = wout$r_dist, 
    p = p, 
    p_new = NULL) +
    scale_fill_scico(
      name = "Relative distance to nearest plot", 
      palette = "bamako") +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        barwidth = 20,
        barheight = 1)))

# Create histogram of pixel-plot dissimilarity with and without proposed plots
visHist(wout$r_dist, wout$r_new_dist)

# Compare distance with old plots only vs. including new plots
plot(wout$r_dist, wout$r_new_dist)
