#' ---
#' title: phylogenetic diversity calculation
#' author: gabriela alves ferreira
#' date: 2023-12-11
#' ---

setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/04_presence_absence_rasters/")

# packages
library(terra)
library(phyloraster)
library(ape)

###########################
ma <- vect("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/01_data/02_variables/00_limit/merge_limites_MA_buffer20km_WGS84.shp")
plot(ma)

# load binary rasters
bi_pres <- rast("rast_pres.tif")
# bi_2050_SSP370 <- rast("rast_2011_2040_SSP370.tif")
# bi_2050_SSP585 <- rast("rast_2011_2040_SSP585.tif")
# bi_2070_SSP370 <- rast("rast_2041_2070_SSP370.tif")
# bi_2070_SSP585 <- rast("rast_2041_2070_SSP585.tif")
bi_2050_SSP370 <- rast("rast_2041_2070_SSP370.tif")
bi_2050_SSP585 <- rast("rast_2041_2070_SSP585.tif")
bi_2070_SSP370 <- rast("rast_2071_2100_SSP370.tif")
bi_2070_SSP585 <- rast("rast_2071_2100_SSP585.tif")

## substitui a primeira letra do nome pra virar maiuscula e bater com arvore
names(bi_pres) <- gsub("(^[a-z])", "\\U\\1", names(bi_pres), perl = TRUE)
names(bi_2050_SSP370) <- gsub("(^[a-z])", "\\U\\1", names(bi_2050_SSP370), perl = TRUE)
names(bi_2050_SSP585) <- gsub("(^[a-z])", "\\U\\1", names(bi_2050_SSP585), perl = TRUE)
names(bi_2070_SSP370) <- gsub("(^[a-z])", "\\U\\1", names(bi_2070_SSP370), perl = TRUE)
names(bi_2070_SSP585) <- gsub("(^[a-z])", "\\U\\1", names(bi_2070_SSP585), perl = TRUE)

###########################
# replace NA by 0 and mask by AF
bi_pres <- subst(bi_pres, NA, 0)
bi_pres <- mask(bi_pres, ma, 0)

bi_2050_SSP370 <- subst(bi_2050_SSP370, NA, 0)
bi_2050_SSP370 <- mask(bi_2050_SSP370, ma, 0)

bi_2050_SSP585 <- subst(bi_2050_SSP585, NA, 0)
bi_2050_SSP585 <- mask(bi_2050_SSP585, ma, 0)

bi_2070_SSP370 <- subst(bi_2070_SSP370, NA, 0)
bi_2070_SSP370 <- mask(bi_2070_SSP370, ma, 0)

bi_2070_SSP585 <- subst(bi_2070_SSP585, NA, 0)
bi_2070_SSP585 <- mask(bi_2070_SSP585, ma, 0)

###########################
# calculate PD for each scenario and year
# load phylogenetic tree
tree <- read.tree("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises/enm_mammals/01_data/03_phylo_tree/Tree_bootstrap_500.txt")
setdiff(names(bi_pres), tree$tip.label)

# removing species that are not in the tree
data_pres <- phylo.pres(bi_pres, tree)
data_50_370 <- phylo.pres(bi_2050_SSP370, tree)
data_50_585 <- phylo.pres(bi_2050_SSP585, tree)
data_70_370 <- phylo.pres(bi_2070_SSP370, tree)
data_70_585 <- phylo.pres(bi_2070_SSP585, tree)


###########################
# calculate richness
sr_pres <- rast.sr(data_pres$x)
sr_2050_SSP370 <- rast.sr(data_50_370$x)
sr_2050_SSP585 <- rast.sr(data_50_585$x)
sr_2070_SSP370 <- rast.sr(data_70_370$x)
sr_2070_SSP585 <- rast.sr(data_70_370$x)

# calculate pd
pd_pres <- rast.pd(data_pres$x, data_pres$tree)
pd_2050_SSP370 <- rast.pd(data_50_370$x, data_50_370$tree)
pd_2050_SSP585 <- rast.pd(data_50_585$x, data_50_585$tree)
pd_2070_SSP370 <- rast.pd(data_70_370$x, data_70_370$tree)
pd_2070_SSP585 <- rast.pd(data_70_370$x, data_70_370$tree)


# visualize the results
# pdf("pd_mammals.pdf", height = 6, width = 6)
par(mfrow = c(3, 2))

plot(pd_pres, main = "PD Baseline")
plot(pd_2050_SSP370, main = "PD 2050 SSP370")
plot(pd_2050_SSP585, main = "PD 2050 SSP585")
plot(pd_2070_SSP370, main = "PD 2070 SSP370")
plot(pd_2070_SSP585, main = "PD 2070 SSP585")
dev.off()

###########################
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")
# save the rasters for PD
writeRaster(pd_pres, "phylo_div_Baseline.tif",
            overwrite=TRUE)
writeRaster(pd_2050_SSP370, "phylo_div_2011_2040_SSP370.tif",
            overwrite=TRUE)
writeRaster(pd_2050_SSP585, "phylo_div_2011_2040_SSP585.tif",
            overwrite=TRUE)
writeRaster(pd_2070_SSP370, "phylo_div_2041_2070_SSP370.tif",
            overwrite=TRUE)
writeRaster(pd_2070_SSP585, "phylo_div_2041_2070_SSP585.tif", 
            overwrite=TRUE)

# or
pd <- c(pd_pres, pd_2050_SSP370, pd_2050_SSP585, pd_2070_SSP370, pd_2070_SSP585)
names(pd) <- c("pd_pres", "pd_2050_SSP370", "pd_2050_SSP585", "pd_2070_SSP370",
               "pd_2070_SSP585")
x11()
plot(pd)
writeRaster(pd, "pd_mammals.tif")

###########################
# save the rasters for SR
writeRaster(sr_pres, "richness_Baseline.tif",
            overwrite=TRUE)
writeRaster(sr_2050_SSP370, "richness_2011_2040_SSP370.tif",
            overwrite=TRUE)
writeRaster(sr_2050_SSP585, "richness_2011_2040_SSP585.tif",
            overwrite=TRUE)
writeRaster(sr_2070_SSP370, "richness_2041_2070_SSP370.tif",
            overwrite=TRUE)
writeRaster(sr_2070_SSP585, "richness_2041_2070_SSP585.tif", 
            overwrite=TRUE)

# or
sr <- c(sr_pres, sr_2050_SSP370, sr_2050_SSP585, sr_2070_SSP370, sr_2070_SSP585)
names(sr) <- c("sr_pres", "sr_2050_SSP370", "sr_2050_SSP585", "sr_2070_SSP370",
               "sr_2070_SSP585")
x11()
plot(sr)
writeRaster(sr, "sr_mammals.tif")
