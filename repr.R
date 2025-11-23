# How representative of forest structure are plots in a site?
# John L. Godlee (johngodlee@gmail.com)
# Last updated: 2025-11-18

# Purpose:
# To what extent are existing plots representative of the broader landscape forest structure? 
# Identify optimal locations for new plots to fill gaps in multidimensional forest structure space.

# Packages
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(scico)
library(GGally)
library(patchwork)

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

###########################################################
# TEST INDIVIDUAL FUNCTIONS FOR REPRESENTATIVENESS ANALYSIS
###########################################################

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
euclid_mahal_dist_comp <- data.frame(
  euclidean = c(old_dist_euclidean), 
  mahalanobis = c(old_dist_mahalanobis)) %>% 
  ggplot(., aes(x = euclidean, y = mahalanobis)) + 
  geom_abline(linetype = 2, colour = "red") + 
  geom_point(alpha = 0.5) + 
  geom_smooth(colour = "blue") + 
  geom_smooth(method = "lm", colour = "green") + 
  theme_classic() +
  labs(
    x = "Euclidean distance",
    y = "Mahalanobis distance") + 
  ggtitle("Distance in structural space from each pixel to nearest plot")

ggsave(euclid_mahal_dist_comp, width = 8, height = 5, 
  file = "./img/euclid_mahal_dist_comp.png")

euclid_mahal_dist_hist <- data.frame(
  euclidean = c(old_dist_euclidean), 
  mahalanobis = c(old_dist_mahalanobis)) %>% 
  pivot_longer(everything()) %>% 
  ggplot(., aes(x = value, fill = name)) + 
    geom_histogram(position = "identity", alpha = 0.5) + 
    scale_fill_discrete(name = "Distance") + 
    theme_classic() + 
    labs(
      x = "Distance",
      y = "Number of pixels")

ggsave(euclid_mahal_dist_hist, width = 8, height = 5, 
  file = "./img/euclid_mahal_dist_hist.png")

# Prepare rasters with distances
old_dist_r_list <- lapply(old_dist_list, function(x) { 
  out <- r[[1]]
  values(out)[!is.na(values(out))] <- x
  return(out)
})
names(old_dist_r_list) <- names(old_dist_list)

# Define colour and shape palettes for plot points
size_map <- c(
  "Existing plot" = 3, 
  "Proposed plot" = 3, 
  "Landscape value" = 1,
  "meanmin" = 3, 
  "minimax" = 3) 

alpha_map <- c(
  "Existing plot" = 0.5, 
  "Proposed plot" = 0.5, 
  "Landscape value" = 1,
  "meanmin" = 0.5, 
  "minimax" = 0.5)

colour_map <- c(
  "Existing plot" = "#E74B5E", 
  "Proposed plot" = "#5499DE", 
  "Existing plot" = "#6bab62", 
  "meanmin" = "#9149a3", 
  "minimax" = "#bf9f3c", 
  "Landscape value" = "darkgrey")

# Create map of the dissimilarity of pixels from existing plots
vis_map_list <- lapply(names(old_dist_r_list), function(x) { 
  ggplot() +
    geom_spatraster(data = old_dist_r_list[[x]]) + 
    scale_fill_scico(name = paste0("Relative ", x, " distance to nearest plot"), 
      palette = "bamako",
      limits = c(0, 8),
      oob = scales::squish) +
    geom_sf(data = p, 
      colour = colour_map["Existing plot"], fill = NA) + 
    guides(
      fill = guide_colourbar(
        title.position = "top",
        barwidth = 20,
        barheight = 1)) + 
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
})
names(vis_map_list) <- names(old_dist_r_list)

lapply(names(vis_map_list), function(x) { 
  ggsave(vis_map_list[[x]], width = 8, height = 6, 
    file = paste0("./img/vis_map_out_", x, ".png"))
})

# Extract PCA values from r_pca
pca_pt_r <- as.data.frame(old_pca$r_pca$x)
pca_pt_r$type <- "Landscape value" 

# Extract PCA values from p_pca
pca_pt_p <- as.data.frame(old_pca$p_pca)
pca_pt_p$type <- "Existing plot" 

# Combine PCA values
pca_pt_all <- rbind(pca_pt_r, pca_pt_p)
pca_pt_all$type <- factor(pca_pt_all$type, 
  levels = c("Landscape value", "Existing plot"))

# Create a PCA biplot of landscape representativeness
vis_pca <- ggplot() +
  geom_point(data = pca_pt_all, 
    aes(x = PC1, y = PC2, 
      size = type, colour = type, alpha = type),
      shape = "circle") + 
  scale_size_manual(name = NULL, values = size_map) + 
  scale_colour_manual(name = NULL, values = colour_map) + 
  scale_alpha_manual(name = NULL, values = alpha_map) + 
  labs(x = "PC1", y = "PC2") + 
  theme_bw() +
  ggtitle("Structural space coverage") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(vis_pca, width = 6, height = 4, 
  file = "./img/vis_pca.png")

# Classify pixels by how well-represented they are by the plots
old_pca_classif <- classifLandscape(
  r_pca = old_pca$r_pca$x, 
  p = old_pca$p_pca, 
  n_pca = 3, 
  ci = 0.95
)
# Returns dataframe with the distance of each pixel and its classification

# Prepare raster with pixel classifications
r_classif <- r[[1]]
values(r_classif)[!is.na(values(r_classif))] <- old_pca_classif$group
r_classif_cls <- data.frame(
  id = c(1, 2), 
  classif = c("Poorly represented", "Well-represented by existing plots"))
levels(r_classif) <- r_classif_cls

# Create map of pixel classification 
vis_map_classif <- ggplot() +
  geom_spatraster(data = r_classif) + 
  scale_fill_discrete(name = NULL, guide = guide_legend(nrow = 3),
    na.value = NA, na.translate = FALSE) +
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

ggsave(vis_map_classif, width = 8, height = 6, 
  file = "./img/vis_map_classif.png")

# Create PCA biplot of pixel classification 
old_pca_classif_df <- as.data.frame(old_pca$r_pca$x)
old_pca_classif_df$group <- old_pca_classif$group

vis_pca_classif <- ggplot(old_pca_classif_df, aes(PC1, PC2, colour = group)) +
  geom_point(alpha = 0.8) +
  scale_colour_discrete(name = NULL) +
  theme_bw() + 
  labs(x = "PC1", y = "PC2") + 
  theme(legend.position = "bottom")

ggsave(vis_pca_classif, width = 6, height = 4, 
  file = "./img/vis_pca_classif.png")

##############################################
# IDENTIFY GAPS AND OPTIMALLY LOCATE NEW PLOTS
##############################################

# Test with existing plots
test1 <- newPlotSelect(
  r = r,
  p = p,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 3,
  method = meanminSelect)

# Test without existing plots
test2 <- newPlotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 3,
  method = meanminSelect)

# Test without PCA 
test3 <- newPlotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = FALSE, 
  n_pca = NULL,
  method = meanminSelect)

# Test with PCA using all columns
test4 <- newPlotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 12,
  method = meanminSelect)

# Output should be identical
ggplot() + 
  geom_sf(data = test3, colour = "green", linewidth = 2) + 
  geom_sf(data = test4, colour = "red", fill = NA) 

# Test with mask layer
r_mask <- r[[1]]
r_mask[values(r_mask) < 35] <- NA_real_
test5 <- newPlotSelect(
  r = r,
  p = NULL,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = r_mask,
  pca = TRUE, 
  n_pca = 3,
  method = meanminSelect)

# test5 (red) should always overlap mask layer   
ggplot() + 
  geom_spatraster(data = r_mask) + 
  geom_sf(data = test4, colour = "green", fill = NA) + 
  geom_sf(data = test5, colour = "red", fill = NA) 

# Test with minimax algorithm
test6 <- newPlotSelect(
  r = r,
  p = p,
  n_plots = 10,
  p_new_dim = c(100, 100),
  r_mask = NULL,
  pca = TRUE, 
  n_pca = 3,
  method = minimaxSelect)

# Compare plot placement uing minimax and meanmin algorithms
ggplot() + 
  geom_spatraster(data = r[[1]]) + 
  geom_sf(data = test1, colour = "green", fill = NA) + 
  geom_sf(data = test6, colour = "red", fill = NA) 

# Compare PCA positions using minimax and meanmin algorithms

# Extract metrics from raster for each existing plot
test1_ext <- extractPlotMetrics(r, test1, fun = NULL)
test6_ext <- extractPlotMetrics(r, test6, fun = NULL)
# Using FUN = NULL returns the values of all values within each pixel

# Put pixels and existing plots in the same PCA space
test1_pca <- as.data.frame(PCALandscape(r, test1_ext, center = TRUE, scale. = TRUE)$p_pca)
test1_pca$type <- "meanmin"
test6_pca <- as.data.frame(PCALandscape(r, test6_ext, center = TRUE, scale. = TRUE)$p_pca)
test6_pca$type <- "minimax"

pca_pt_test_1_6 <- rbind(pca_pt_r, test1_pca, test6_pca, pca_pt_p)
pca_pt_test_1_6$type <- factor(pca_pt_test_1_6$type,
  levels = c("Landscape value", "meanmin", "minimax", "Existing plot"))

vis_pca_meanmin_minimax_comp <- ggplot() +
  geom_point(data = pca_pt_test_1_6, 
    aes(x = PC1, y = PC2, 
      size = type, colour = type, alpha = type),
      shape = "circle") + 
  scale_size_manual(name = NULL, values = size_map) + 
  scale_colour_manual(name = NULL, values = colour_map) + 
  scale_alpha_manual(name = NULL, values = alpha_map) + 
  labs(x = "PC1", y = "PC2") + 
  theme_bw() +
  ggtitle("Structural space coverage") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(vis_pca_meanmin_minimax_comp, width = 6, height = 4, 
  file = "./img/vis_pca_meanmin_minimax_comp.png")

# Calculate distance from each pixel to nearest plot in PCA space
# Using meanmin and minimax algorithms
test1_dist_euclidean <- pcaDist(old_pca$r_pca$x, 
  as.matrix(test1_pca[,grepl("PC", colnames(test1_pca))]), 
  n_pca = 3, k = 1, method = "euclidean")

test6_dist_euclidean <- pcaDist(old_pca$r_pca$x, 
  as.matrix(test6_pca[,grepl("PC", colnames(test6_pca))]), 
  n_pca = 3, k = 1, method = "euclidean")
 
# Create histograms to compare resulting distances
test1_6_dist_hist <- data.frame(
  existing = c(old_dist_euclidean), 
  meanmin = c(test1_dist_euclidean),
  minimax = c(test6_dist_euclidean)) %>% 
  pivot_longer(everything()) %>% 
  mutate(name = factor(name, 
    levels = c("existing", "meanmin", "minimax"),
    labels = c("Existing plot", "meanmin", "minimax"))) %>% 
  ggplot(., aes(x = value, fill = name)) + 
    geom_histogram(position = "identity", colour = "black") + 
    scale_fill_manual(name = NULL, values = colour_map) + 
    facet_wrap(~name, scales = "fixed", ncol = 1) + 
    theme_classic() + 
    theme(legend.position = "none") + 
    labs(
      x = "Distance",
      y = "Number of pixels")

ggsave(test1_6_dist_hist, width = 6, height = 5, 
  file = "./img/test1_6_dist_hist.png")
