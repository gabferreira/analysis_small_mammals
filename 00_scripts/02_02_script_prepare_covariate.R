#' ---
#' title: var - download
#' author: mauricio vancine
#' date: 2022-04-21
#' ---

# prepare r ---------------------------------------------------------------

# packages
library(tidyverse)
library(rnaturalearth)
library(sf)
library(terra)
library(doParallel)
library(parallelly)
library(foreach)

# options
options(timeout = 1e6)

# prepare -----------------------------------------------------------------

# directory
dir.create("01_data/02_variables/01_climate_present/01_cropped")
dir.create("01_data/02_variables/02_climate_future/01_cropped")

### import vector ----
sa <- sf::st_union(
    rnaturalearth::ne_countries(scale = 10, continent = "South America"),
    rnaturalearth::ne_countries(scale = 10, country = "France")) %>% 
    sf::st_crop(rnaturalearth::ne_countries(scale = 10, continent = "South America")) %>% 
    sf::st_union() %>% 
    terra::vect()
sa

### wc present ----
files <- dir(path = "/media/mude/afe69132-ffdb-4892-b809-a0f7d2b8f423/geospatial_data_base/02_raster/chelsa/bioclimatic/", 
             pattern = "1981-2010", full.names = TRUE)
files

chelsa_p <- terra::rast(files)
chelsa_p

names(chelsa_p) <- c(paste0("bio0", 1:9), paste0("bio", 10:19))
chelsa_p

chelsa_p_sa <- chelsa_p %>% 
    terra::crop(sa) %>% 
    terra::mask(sa)
chelsa_p_sa

chelsa_p_sa_5km <- terra::aggregate(chelsa_p_sa, fact = 5, cores = 6)
chelsa_p_sa_5km

plot(chelsa_p_sa_5km[[1]])
lines(sa)

terra::writeRaster(chelsa_p_sa_5km, 
                   paste0("01_data/02_variables/01_climate_present/01_cropped/", names(chelsa_p_sa_5km), "_sa.tif"), 
                   overwrite = TRUE)

### wc future ----
pattern <- paste0(c("ssp370", "ssp585"), collapse = "|")

files <- dir(path = "/media/mude/afe69132-ffdb-4892-b809-a0f7d2b8f423/geospatial_data_base/02_raster/chelsa/bioclimatic", 
             pattern = pattern, full.names = TRUE) %>% 
    grep(paste0("bio17|bio18|bio19"), ., value = TRUE) %>% 
    grep(paste0(".tif$"), ., value = TRUE)
files

length(files)

doParallel::registerDoParallel(parallelly::availableCores(omit = 6))
foreach::foreach(i=files) %do% {
    
    print(i)
    
    r <- terra::rast(i)
    
    r %>% 
        terra::crop(sa) %>%
        terra::mask(sa) %>%
        terra::aggregate(fact = 5) %>% 
        terra::writeRaster(paste0("01_data/02_variables/02_climate_future/01_cropped/", sub(".tif", "_sa.tif", basename(i))), overwrite = TRUE)
    
}
doParallel::stopImplicitCluster()

# mean
fut_scenarios <- expand.grid(
    year = c("2011-2040", "2041-2070", "2071-2100"),
    ssp = c("ssp370", "ssp585")) %>% 
    apply(MARGIN = 1, FUN = paste0, collapse = "_")
fut_scenarios

var_f <- dir(path = "01_data/02_variables/02_climate_future/01_cropped", 
             pattern = ".tif$", full.names = TRUE) %>% 
    terra::rast()
var_f

doParallel::registerDoParallel(parallelly::availableCores(omit = 6))
for(i in fut_scenarios){
    
    print(i)
    
    fut_scenarios_i <- stringr::str_split(i, "_", simplify = TRUE)
    
    var_f_year <- terra::subset(var_f, grep(fut_scenarios_i[, 1], names(var_f), value = TRUE))   
    var_f_year_ssp <- terra::subset(var_f_year, grep(fut_scenarios_i[, 2], names(var_f_year), value = TRUE))   
    
    foreach::foreach(j = 1:19) %dopar% {

        print(j)
        bio <- ifelse(j < 10, paste0("bio0", j), paste0("bio", j))
        var_f_year_spp_bio <- terra::subset(var_f_year_ssp, grep(bio, names(var_f_year_ssp), value = TRUE))   
        var_f_year_spp_bio_mean <- terra::mean(var_f_year_spp_bio)      
        names(var_f_year_spp_bio_mean) <- paste0("CHELSA_", bio, "_", i, "_V.2.1")
        terra::writeRaster(var_f_year_spp_bio_mean, paste0("01_data/02_variables/02_climate_future/02_cropped_mean/CHELSA_", bio, "_", i, "_V.2.1_sa.tif"), overwrite = TRUE)
        
    }
    
}
doParallel::stopImplicitCluster()

# end ---------------------------------------------------------------------