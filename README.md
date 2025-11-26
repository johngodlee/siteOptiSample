# siteOptiSample 

`siteOptiSample` is an R package developed for the GEO-TREES initiative. It provides tools to design optimized tree-inventory plot sampling strategies for forest biomass reference sites. By analysing raster datasets such as structural metrics from airborne LiDAR data, this package helps select plot locations that best represent the heterogeneity of a given site.

## Overview

Designing a sampling strategy for ground-truth data collection is critical for accurate biomass estimation. Random sampling often leads to poor representation of the variation of forest structure present in the site.

## Installation

You can install the development version of `siteOptiSample` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("GEO-TREES/siteOptiSample")
```

## Dependencies

This package relies heavily on the modern R spatial ecosystem. Ensure you have the following installed:

* `terra` - for spatial raster data
* `sf` - for spatial vector data
* `ggplot2` - for visualisation

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.
