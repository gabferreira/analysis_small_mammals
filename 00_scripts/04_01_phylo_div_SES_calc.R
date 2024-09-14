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
bi_2050_SSP370 <- rast("rast_2011_2040_SSP370.tif")
bi_2050_SSP585 <- rast("rast_2011_2040_SSP585.tif")
bi_2070_SSP370 <- rast("rast_2041_2070_SSP370.tif")
bi_2070_SSP585 <- rast("rast_2041_2070_SSP585.tif")

## first letter in caps lock
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

# calculate ses pd
library(SESraster)
pd_pres_ses <- rast.pd.ses(x = data_pres$x, tree = data_pres$tree,
                           random = "spat", aleats = 10)
pd_2050_SSP370_ses <- rast.pd.ses(x = data_pres$x, tree = data_pres$tree,
                                  random = "spat", aleats = 10)
pd_2050_SSP585_ses <- rast.pd.ses(x = data_pres$x, tree = data_pres$tree,
                                  random = "spat", aleats = 10)
pd_2070_SSP370_ses <- rast.pd.ses(x = data_pres$x, tree = data_pres$tree,
                                  random = "spat", aleats = 10)
pd_2070_SSP585_ses <- rast.pd.ses(x = data_pres$x, tree = data_pres$tree,
                                  random = "spat", aleats = 10)

##########################################################################
## calculating sesPD with canaper
library(canaper)

####################################
## Present
# transform in dataframe to use in the canaper package
x.m <- as.data.frame(data_pres$x, xy = TRUE)

# Filtrar células que têm pelo menos uma espécie presente
x.m <- x.m[rowSums(x.m[, -c(1:2)]) > 0, ]  # Exclui as colunas de coordenadas XY para o filtro

# Coordenadas das células filtradas
xy <- x.m[, c("x", "y")]

# getting coordinates for each cell
names(xy) <- c("lon", "lat")

# removing the xy columns from the community matrix 
x.m <- x.m[,-c(1,2)]
head(x.m)

# run randomization test
set.seed(071421)
rand_test_results <- cpr_rand_test(
    x.m, tree,
    null_model = "curveball", metrics = c("pd"), n_reps = 499)
#7:09 comecou
#7:42 terminou

rand_test_results[, 1:9]
View(rand_test_results)

# plot the results
# (add lat/long columns)
ses_resu <- data.frame(rand_test_results, xy)

# Plot the results
# rasterizing to save the results for bufonidae
r <- rast(ses_resu[,c("lon", "lat", "pd_obs")], type="xyz")
plot(r)

pd_pres_ses <- rast(ses_resu[,c("lon", "lat", "pd_obs_z")], type="xyz")
plot(pd_pres_ses)

## saving the resul for SES PD for present
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")
writeRaster(pd_pres_ses, "pd_pres_SES.tif")

##########################################################################
## 2050 SSP370
library(canaper)

# transform in dataframe to use in the canaper package
x.m <- as.data.frame(data_50_370$x, xy = TRUE)

# Filtrar células que têm pelo menos uma espécie presente
x.m <- x.m[rowSums(x.m[, -c(1:2)]) > 0, ]  # Exclui as colunas de coordenadas XY para o filtro

# Coordenadas das células filtradas
xy <- x.m[, c("x", "y")]

# getting coordinates for each cell
names(xy) <- c("lon", "lat")

# removing the xy columns from the community matrix 
x.m <- x.m[,-c(1,2)]
head(x.m)

# run randomization test
set.seed(071421)
rand_test_results <- cpr_rand_test(
    x.m, tree,
    null_model = "curveball", metrics = c("pd"), n_reps = 999)
#7:09 comecou
#7:42 terminou

rand_test_results[, 1:9]
View(rand_test_results)

# plot the results
# (add lat/long columns)
ses_resu <- data.frame(rand_test_results, xy)

# Plot the results
# rasterizing to save the results for bufonidae
r <- rast(ses_resu[,c("lon", "lat", "pd_obs")], type="xyz")
plot(r)

pd_50_370_ses <- rast(ses_resu[,c("lon", "lat", "pd_obs_z")], type="xyz")
plot(pd_50_370_ses)

## saving the resul for SES PD for present
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")
writeRaster(pd_50_370_ses, "pd_2011_2040_370_SES.tif")

####################################
## 2050 SSP 585
# transform in dataframe to use in the canaper package
x.m <- as.data.frame(data_50_585$x, xy = TRUE)

# Filtrar células que têm pelo menos uma espécie presente
x.m <- x.m[rowSums(x.m[, -c(1:2)]) > 0, ]  # Exclui as colunas de coordenadas XY para o filtro

# Coordenadas das células filtradas
xy <- x.m[, c("x", "y")]

# getting coordinates for each cell
names(xy) <- c("lon", "lat")

# removing the xy columns from the community matrix 
x.m <- x.m[,-c(1,2)]
head(x.m)

# run randomization test
set.seed(071421)
rand_test_results <- cpr_rand_test(
    x.m, tree,
    null_model = "curveball", metrics = c("pd"), n_reps = 999)
#7:09 comecou
#7:42 terminou

# rand_test_results[, 1:9]
# View(rand_test_results)

# plot the results
# (add lat/long columns)
ses_resu <- data.frame(rand_test_results, xy)

# Plot the results
# rasterizing to save the results for bufonidae

pd_2050_585_ses <- rast(ses_resu[,c("lon", "lat", "pd_obs_z")], type="xyz")
# plot(pd_2050_585_ses)

## saving the resul for SES PD for present
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")
writeRaster(pd_2050_585_ses, "pd_2011_2040_585_SES.tif", overwrite = T)

##########################################################################
## 2070 SSP370
library(canaper)

# transform in dataframe to use in the canaper package
x.m <- as.data.frame(data_70_370$x, xy = TRUE)

# Filtrar células que têm pelo menos uma espécie presente
x.m <- x.m[rowSums(x.m[, -c(1:2)]) > 0, ]  # Exclui as colunas de coordenadas XY para o filtro

# Coordenadas das células filtradas
xy <- x.m[, c("x", "y")]

# getting coordinates for each cell
names(xy) <- c("lon", "lat")

# removing the xy columns from the community matrix 
x.m <- x.m[,-c(1,2)]
# head(x.m)

# run randomization test
set.seed(071421)
rand_test_results <- cpr_rand_test(
    x.m, tree,
    null_model = "curveball", metrics = c("pd"), n_reps = 999)
#7:09 comecou
#7:42 terminou

# rand_test_results[, 1:9]
# View(rand_test_results)

# plot the results
# (add lat/long columns)
ses_resu <- data.frame(rand_test_results, xy)

# Plot the results
# rasterizing to save the results for bufonidae
r <- rast(ses_resu[,c("lon", "lat", "pd_obs")], type="xyz")
# plot(r)

pd_70_370_ses <- rast(ses_resu[,c("lon", "lat", "pd_obs_z")], type="xyz")
# plot(pd_70_370_ses)

## saving the resul for SES PD for present
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")
writeRaster(pd_70_370_ses, "pd_2041_2070_370_SES.tif")

#########################################################################
## 2070 SSP 585
# transform in dataframe to use in the canaper package
x.m <- as.data.frame(data_70_585$x, xy = TRUE)

# Filtrar células que têm pelo menos uma espécie presente
x.m <- x.m[rowSums(x.m[, -c(1:2)]) > 0, ]  # Exclui as colunas de coordenadas XY para o filtro

# Coordenadas das células filtradas
xy <- x.m[, c("x", "y")]

# getting coordinates for each cell
names(xy) <- c("lon", "lat")

# removing the xy columns from the community matrix 
x.m <- x.m[,-c(1,2)]
# head(x.m)

# run randomization test
set.seed(071421)
rand_test_results <- cpr_rand_test(
    x.m, tree,
    null_model = "curveball", metrics = c("pd"), n_reps = 999)
#7:09 comecou
#7:42 terminou

# rand_test_results[, 1:9]
# View(rand_test_results)

# plot the results
# (add lat/long columns)
ses_resu <- data.frame(rand_test_results, xy)

# Plot the results
# rasterizing to save the results for bufonidae

pd_2070_585_ses <- rast(ses_resu[,c("lon", "lat", "pd_obs_z")], type="xyz")
# plot(pd_2070_585_ses)

## saving the resul for SES PD for present
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")
writeRaster(pd_2070_585_ses, "pd_2041_2070_585_SES.tif")

