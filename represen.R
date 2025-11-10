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
subplots_fil <- subplots_als %>% 
  filter(grepl("bottomleft", subplot_id))

# Select relevant metrics from ALS data
# names(als)
als_sel <- als[[c(
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
als_vals <- as.data.frame(als_sel)

# Compute correlation between layers
# Visualize correlation
ggcorr(als_vals, label = TRUE)
# if very collinear maybe just use biplot of metrics rather than a PCA

############## 
# RUN ANALYSIS
############## 

# Extract metrics from raster for each existing plot
old_ext <- extract_plot_metrics(als_sel, subplots_fil, fun = mean)
# Returns a dataframe of values in the same row order as subplots_fil
# Optional, could also use the computed values for each subplot, 
# rather than extracting point estimates from a grid.

# Put pixels and existing plots in the same PCA space
old_pca <- pca_landscape(als_sel, old_ext, center = TRUE, scale. = TRUE)
# Returns two dataframes of PCA scores, for pixels (r_pca) and plots (p_pca)
# r_pca is a full PCA object
# p_pca is only PCA scores
# PCA is run on variables scaled to unit variance shifted to center on zero

# Calculate distance from each pixel to nearest plot in PCA space
old_dist_euclidean <- pca_dist(old_pca$r_pca$x, old_pca$p_pca, 
  n_pca = 3, k = 1, method = "euclidean")

old_dist_mahalanobis <- pca_dist(old_pca$r_pca$x, old_pca$p_pca, 
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
map_plot_vis_out_list <- lapply(names(old_dist_list), function(x) {  
  map_plot_vis(als_sel, old_dist_list[[x]], p = subplots_fil) + 
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
names(map_plot_vis_out_list) <- names(old_dist_list)

lapply(names(map_plot_vis_out_list), function(x) { 
  ggsave(map_plot_vis_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/map_plot_vis_out_", x, ".png"))
})

# Create plot of existing plots within landscape PCA space 
pca_plot_vis_out <- pca_plot_vis(old_pca$r_pca$x, old_pca$p_pca)

ggsave(pca_plot_vis_out, width = 6, height = 4, 
  file = "./img/pca_plot_vis_out.png")

# Classify pixels by how well-represented they are by the plots
old_pca_classif <- landscape_classif(
  r_pca = old_pca$r_pca$x, 
  p = old_pca$p_pca, 
  n_pca = 3, 
  ci = 0.95
)
# Returns dataframe with the distance of each pixel and its classification

# Create map of pixel classification 
map_classif_vis_out <- map_classif_vis(als_sel, old_pca_classif)

ggsave(map_classif_vis_out, width = 8, height = 6, file = 
  "./img/map_classif_vis_out.png")

# Create PCA biplot of pixel classification 
pca_classif_vis_out <- pca_classif_vis(old_pca$r_pca$x, old_pca_classif)

ggsave(pca_classif_vis_out, width = 6, height = 4, file = 
  "./img/pca_classif_vis_out.png")

# Extract candidate plot locations (50% sample of pixels in the landscape)
# Excluding existing plots
# XY coordinates of each candidate location
subplots_cent <- subplots_fil %>% 
  st_centroid() %>% 
  vect() %>% 
  geom(.) %>% 
  as.data.frame() %>%
  dplyr::select(x, y) %>% 
  cellFromXY(als_sel, .)
als_sel_noplot <- als_sel
als_sel_noplot[subplots_cent] <- NA
cand_cds <- as.data.frame(crds(als_sel_noplot))
cand <- cand_cds[sample(seq_len(nrow(cand_cds)), floor(nrow(cand_cds) / 2)),]
# This could also be, for example, locations based on accessibility within the site
# or existing plots that could be re-measured

# Map the candidate plots (pixels)
map_plot_vis_out_list$euclidean + 
  geom_point(data = cand, aes(x = x, y = y), colour = "blue")

# Extract raster values for candidate plot locations
cand_ext <- extract_plot_metrics(als_sel, cand)

# Put candidate plots in PCA space
cand_pca <- pca_landscape(als_sel, cand_ext, center = TRUE, scale. = TRUE)
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
opt_maximin <- maximin_select(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_minimax <- minimax_select(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_meanmin <- meanmin_select(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_maxilhs <- lhs_select(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

opt_kmeans <- kmeans_select(
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
map_plot_vis_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  # Combine PCA scores of existing and proposed plots
  old_cand_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel)

  # Analyse representativeness of existing and proposed plots
  old_cand_dist <- pca_dist(cand_pca$r_pca$x, old_cand_pca_sel, n_pca = 3, k = 1)

  # Create map of the dissimilarity of pixels from proposed and existing plots 
  map_plot_vis(
    r = als_sel, 
    r_dist = old_cand_dist, 
    p = subplots_fil, 
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
names(map_plot_vis_new_out_list) <- names(opt_list)

lapply(names(map_plot_vis_new_out_list), function(x) { 
  ggsave(map_plot_vis_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/map_plot_vis_new_out_", x, ".png"))
})

# Create biplots of proposed and existing plots within landscape PCA space 
pca_plot_vis_new_out_list <- lapply(names(opt_list), function(x) { 

  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  pca_plot_vis(
    r_pca = cand_pca$r_pca$x, 
    p_pca = as.data.frame(old_pca$p_pca),
    p_new_pca = as.data.frame(cand_pca_sel))
})
names(pca_plot_vis_new_out_list) <- names(opt_list)

lapply(names(pca_plot_vis_new_out_list), function(x) { 
  ggsave(pca_plot_vis_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/pca_plot_vis_new_out_", x, ".png"))
})

# Create histograms of pixel-plot dissimilarity with and without proposed plots
hist_plot_vis_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  # Combine PCA scores of existing and proposed plots
  old_cand_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel)

  # Analyse representativeness of existing and proposed plots
  old_cand_dist <- pca_dist(cand_pca$r_pca$x, old_cand_pca_sel, n_pca = 3, k = 1)

  # Create plot
  hist_plot_vis(old_dist_list$euclidean, old_cand_dist)
})
names(hist_plot_vis_new_out_list) <- names(opt_list)

lapply(names(hist_plot_vis_new_out_list), function(x) { 
  ggsave(hist_plot_vis_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/hist_plot_vis_new_out_", x, ".png"))
})

# Create maps which classify pixels by how well-represented they are by the
# existing and new plots
map_classif_vis_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  cand_classif <- landscape_classif(
    r_pca = old_pca$r_pca$x, 
    p = old_pca$p_pca, 
    p_new = cand_pca_sel,
    n_pca = 3, 
    ci = 0.95
  )

  # Create map of pixel classification 
  map_classif_vis_new_out <- map_classif_vis(als_sel, cand_classif)
})
names(map_classif_vis_new_out_list) <- names(opt_list)

lapply(names(map_classif_vis_new_out_list), function(x) { 
  ggsave(map_classif_vis_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/map_classif_vis_new_out_", x, ".png"))
})

# Create PCA biplots of pixel classification 
pca_classif_vis_new_out_list <- lapply(names(opt_list), function(x) { 
  # Extract PCA scores for proposed plots
  cand_pca_sel <- cand_pca$p_pca[opt_list[[x]],]

  cand_classif <- landscape_classif(
    r_pca = old_pca$r_pca$x, 
    p = old_pca$p_pca, 
    p_new = cand_pca_sel,
    n_pca = 3, 
    ci = 0.95
  )

  pca_classif_vis_new_out <- pca_classif_vis(old_pca$r_pca$x, cand_classif) + 
    guides(colour = guide_legend(nrow = 2))
})
names(pca_classif_vis_new_out_list) <- names(opt_list)

lapply(names(pca_classif_vis_new_out_list), function(x) { 
  ggsave(pca_classif_vis_new_out_list[[x]], width = 8, height = 6, 
    file = paste0("./img/pca_classif_vis_new_out_", x, ".png"))
})


###################################################
# TRY WORKFLOW WITH PLOTS LARGER THAN RASTER PIXELS
###################################################

# Create a raster of landscape PCA values 
r_pca_rast <- rast(replicate(ncol(old_pca$r_pca$x), als_sel[[1]]))
values(r_pca_rast)[complete.cases(values(r_pca_rast)),] <- old_pca$r_pca$x
names(r_pca_rast) <- colnames(old_pca$r_pca$x)

# Resample landscaep PCA raster to 4 ha (200x200 m)
r_pca_rast_big <- rastResample(r_pca_rast, 200)

# Extract values from landscape PCA raster 
r_pca_rast_big_val <- values(r_pca_rast_big)[complete.cases(values(r_pca_rast_big)),]

# Optionally mask pixels containing existing plots
r_pca_rast_big_fil <- mask(r_pca_rast_big, vect(subplots_fil), inverse = TRUE)

# Extract values from masked landscape PCA
r_pca_rast_big_fil_val <- values(r_pca_rast_big_fil)[
  complete.cases(values(r_pca_rast_big_fil)),]

# Identify optimal locations for additional 4 ha (200x200 m) plots 
# Using maximin algorithm. Could be any of the algorithms.
opt_maximin_big <- maximin_select(
  r_pca = r_pca_rast_big_val,  # PCA values resampled to 200x200 m
  p_new_pca = r_pca_rast_big_fil_val,  # Candidate plot PCA values
  p_pca = old_pca$p_pca,  # PCA values of existing plots
  n_plots = 5, 
  n_pca = 3
)

# Extract PCA scores of proposed plots
cand_pca_sel_big <- r_pca_rast_big_val[opt_maximin_big,]

# Combine PCA scores of existing and proposed plots
old_cand_big_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel_big)

# Analyse representativeness of existing and proposed plots for each pixel
cand_big_dist <- pca_dist(r_pca_rast_big_val, old_cand_big_pca_sel, 
  n_pca = 3, k = 1)

# Extract coordinates for candidate plots
cand_pca_sel_big_cds <- as.data.frame(crds(r_pca_rast_big))[opt_maximin_big,]

# Create PCA biplot of proposed and existing plots within 200x200 m 
# landscape PCA space
pca_plot_vis(
  r_pca = r_pca_rast_big_val,
  p_pca = as.data.frame(old_pca$p_pca),
  p_new_pca = as.data.frame(cand_pca_sel_big))

# Create map of the dissimilarity of pixels from proposed and existing plots 
map_plot_vis_new_big <- map_plot_vis(
  r = r_pca_rast_big, 
  r_dist = cand_big_dist, 
  p = subplots_fil, 
  p_new = cand_pca_sel_big_cds 
  ) +
  scale_fill_scico(
    name = "Relative distance to nearest plot", 
    palette = "bamako"
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = 20,
      barheight = 1 
    )
  )

ggsave(map_plot_vis_new_big, width = 8, height = 6, 
  file = "./img/map_plot_vis_new_big.png")
