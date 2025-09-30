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

# Only include bottomleft subplots
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
old_pca <- pca_landscape(als_sel, old_ext)
# Returns two dataframes of PCA scores, for pixels (r_pca) and plots (p_pca)
# r_pca is a full PCA object
# p_pca is only PCA scores

# Calculate distance from each pixel to nearest plot in PCA space
old_dist <- pca_dist(old_pca$r_pca$x, old_pca$p_pca, n_pca = 3, k = 1)
# Returns vector of distances for each pixel in r_pca

# Create map of the dissimilarity of pixels from existing plots
map_plot_vis_out <- map_plot_vis(als_sel, old_dist, p = subplots_fil)

ggsave(map_plot_vis_out, width = 6, height = 6, 
  file = "./img/map_plot_vis_out.png")

# Create plot of existing plots within landscape PCA space 
pca_plot_vis_out <- pca_plot_vis(old_pca$r_pca$x, old_pca$p_pca)

ggsave(pca_plot_vis_out, width = 6, height = 4, 
  file = "./img/pca_plot_vis_out.png")

# Classify pixels by how well-represented they are by the plots
old_pca_classif <- landscape_classif(
  r_pca = old_pca$r_pca$x, 
  p = old_pca$p_pca, 
  n_pca = 3, 
  ci = 0.95, 
  method = "mahalanobis"
)
# Returns dataframe with the distance of each pixel and its classification

# Create map of pixel classification 
map_classif_vis_out <- map_classif_vis(als_sel, old_pca_classif)

ggsave(map_classif_vis_out, width = 6, height = 6, file = 
  "./img/map_classif_vis_out.png")

# Create PCA biplot of pixel classification 
pca_classif_vis_out <- pca_classif_vis(old_pca$r_pca$x, old_pca_classif)

ggsave(pca_classif_vis_out, width = 6, height = 4, file = 
  "./img/pca_classif_vis_out.png")

# Extract candidate plot locations (50% sample of pixels in the landscape)
# XY coordinates of each candidate location
cand_cds <- as.data.frame(crds(als_sel))
cand <- cand_cds[sample(seq_len(nrow(cand_cds)), floor(nrow(cand_cds) / 2)),]
# This could also be, for example, locations based on accessibility within the site
# or existing plots that could be re-measured

# Extract raster values for candidate plot locations
cand_ext <- extract_plot_metrics(als_sel, cand)

# Put candidate plots in PCA space
cand_pca <- pca_landscape(als_sel, cand_ext)
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

# Optimally locate new plots within landscape based on PCA scores
opt_maximin <- maximin_select(
  r_pca = cand_pca$r_pca$x, 
  p_new_pca = cand_pca$p_pca,
  p_pca = old_pca$p_pca,
  n_plots = 12, 
  n_pca = 3
)

# Do the same but using the minimax algorithm
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

# Check no selections are duplicated
stopifnot(all(!duplicated(opt_maximin)))
stopifnot(all(!duplicated(opt_minimax)))
stopifnot(all(!duplicated(opt_meanmin)))

# Select one method
opt_sel <- opt_meanmin

# Extract PCA scores for proposed plots
cand_pca_sel <- cand_pca$p_pca[opt_sel,]

# Combine PCA scores of existing and proposed plots
old_cand_pca_sel <- rbind(old_pca$p_pca, cand_pca_sel)

# Create a PCA biplot to check data
ggplot() + 
  geom_point(data = cand_pca$r_pca$x, aes(x = PC1, y = PC2), size = 2) + 
  geom_point(data = old_pca$p_pca, aes(x = PC1, y = PC2), colour = "green", size = 0.6) + 
  geom_point(data = cand_pca_sel, aes(x = PC1, y = PC2), colour = "red", size = 0.5) + 
  theme_bw()
# Candidate plots are a subset of the original pixels 

# Analyse representativeness of existing and proposed plots
old_cand_dist <- pca_dist(cand_pca$r_pca$x, old_cand_pca_sel, n_pca = 3, k = 1)

# Create map of the dissimilarity of pixels from proposed and existing plots 
map_plot_vis_new_out <- map_plot_vis(
  r = als_sel, 
  r_dist = old_cand_dist, 
  p = subplots_fil, 
  p_new = cand[opt_sel,]
)

ggsave(map_plot_vis_new_out, width = 6, height = 6, 
  file = "./img/map_plot_vis_new_out.png")

# Create plot of proposed and existing plots within landscape PCA space 
pca_plot_vis_new_out <- pca_plot_vis(
  r_pca = cand_pca$r_pca$x, 
  p_pca = as.data.frame(old_pca$p_pca),
  p_new_pca = as.data.frame(cand_pca_sel))

ggsave(pca_plot_vis_new_out, width = 6, height = 4, 
  file = "./img/pca_plot_vis_new_out.png")

# Create histograms of pixel-plot dissimilarity with and without proposed plots
hist_plot_vis_new_out <- hist_plot_vis(old_dist, old_cand_dist)

ggsave(hist_plot_vis_new_out, width = 6, height = 4, 
  file = "./img/hist_plot_vis_new_out.png")

# Classify pixels by how well-represented they are by the existing and new plots
cand_classif <- landscape_classif(
  r_pca = old_pca$r_pca$x, 
  p = old_pca$p_pca, 
  p_new = cand_pca_sel,
  n_pca = 3, 
  ci = 0.95, 
  method = "mahalanobis"
)

# Create map of pixel classification 
map_classif_vis_new_out <- map_classif_vis(als_sel, cand_classif)

ggsave(map_classif_vis_new_out, width = 6, height = 6, 
  file = "./img/map_classif_vis_new_out.png")

# Create PCA biplot of pixel classification 
pca_classif_vis_new_out <- pca_classif_vis(old_pca$r_pca$x, cand_classif)

ggsave(pca_classif_vis_new_out, width = 6, height = 4, 
  file = "./img/pca_classif_vis_new_out.png")

