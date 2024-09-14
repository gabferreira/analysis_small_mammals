#' ---
#' title: recorting PD by PAs and remnants
#' author: gabriela alves ferreira
#' date: 2023-12-15
#' ---
#' 

# setwd("/media/gabriela/Gabi_HD/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises/enm_mammals")
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")

library(terra)
# install.packages("phyloraster")
library(phyloraster)

# calculate the delta PD
# load phylogenetic diversity rasters
pd <- rast("pd_mammals.tif")
pd
plot(pd)

pas <- rast("E:/shapes_rasters/rasters/atlantic_forest/497_atlantic_spatial_protected_areas_binary_cropped.tif")
plot(pas)
# rem <- rast("/media/gabriela/Gabi_HD/shapes_rasters/rasters/atlantic_forest/003_atlantic_spatial_forest_vegetation_binary.tif")
rem <- rast("E:/shapes_rasters/rasters/atlantic_forest/003_atlantic_spatial_reprojected_remnants.tif")
plot(rem)

# rem_agg <- terra::aggregate(rem, res(pd)[1]*3600, fun = "modal", cores = 10)
# rem_agg
# plot(rem_agg)
# 
# rem_agg_proj <- terra::project(rem_agg, pd, method = "near")
# rem_agg_proj
# plot(rem_agg_proj)
# 
# rem_agg_proj_na <- rem_agg_proj
# rem_agg_proj_na[rem_agg_proj_na == 0] <- NA
# plot(rem_agg_proj_na)
# 
# writeRaster(rem_agg_proj_na, "/media/gabriela/Gabi_HD/shapes_rasters/rasters/atlantic_forest/003_atlantic_spatial_reprojected_remnants.tif")

##############################
# delta
# delta_pd2040_370 <- delta.grid(pd$pd_pres, pd$pd_2040_SSP370)
# delta_pd2040_585 <- delta.grid(pd$pd_pres, pd$pd_2040_SSP585)
# delta_pd2070_370 <- delta.grid(pd$pd_pres, pd$pd_2070_SSP370)
delta_pd2070_585 <- delta.grid(pd$pd_pres, pd$pd_2070_SSP585)

# delta_pd <- c(delta_pd2040_370, delta_pd2040_585, delta_pd2070_370, delta_pd2070_585)
# plot(delta_pd)

#############################
# mask the delta_pd by the PAs and remnants
# using only the worse scenario (delta_pd2070_585)

# pas
X11()
# align rast 'pas' to raster 'delta_pd2070_585'
pas_aligned <- resample(pas, delta_pd2070_585, method = "near")

# crop
delta_pa <- crop(delta_pd2070_585, pas_aligned, mask = TRUE)
plot(delta_pa)

# remnants
# align rast 'rem' to raster 'delta_pd2070_585'
rem_aligned <- resample(rem, delta_pd2070_585, method = "near")

delta_rem <- crop(delta_pd2070_585, rem_aligned, mask = TRUE)
plot(delta_rem)

# save the rasters
writeRaster(delta_pa, "delta_pd_2070_585_PAs.tif", overwrite=TRUE)
writeRaster(delta_rem, "delta_pd_2070_585_REMNANTS.tif", overwrite=TRUE)
