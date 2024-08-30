#' ---
#' title: occurrences - cleaning
#' author: mauricio vancine
#' date: 2024-04-10
#' ---

# prepare r ---------------------------------------------------------------

# packages
library(tidyverse)
library(lubridate)
library(rnaturalearth)
library(sf)
library(terra)
library(tmap)
library(CoordinateCleaner)

# options
tmap_options(check.and.fix = TRUE)
sf::sf_use_s2(FALSE)



# import ----------------------------------------------------------------

# import
occ <- readr::read_csv("01_data/01_occurrences/01_raw/00_occ_raw.csv") %>% 
    tibble::rowid_to_column(var = "id")
occ

# species list
species_list <- readr::read_csv2("01_data/01_occurrences/00_species_list/species_list.csv") %>% 
    dplyr::arrange(species) %>% 
    dplyr::pull(species)
species_list

# information
nrow(occ)
length(unique(occ$species))

# south america
sa <- rnaturalearth::ne_countries(scale = 10, continent = "South America") %>% 
    sf::st_union(rnaturalearth::ne_countries(scale = 10, country = "France")) %>% 
    sf::st_crop(rnaturalearth::ne_countries(scale = 10, continent = "South America"))
sa

tm_shape(sa) +
    tm_polygons()

# date filter -------------------------------------------------------------

# temporal

# fauna
occ_filter_date <- occ %>% 
    dplyr::mutate(date_filter = ifelse(year >= 1970 & year <= 2024 | is.na(year), TRUE, FALSE))
occ_filter_date

# precision ---------------------------------------------------------------

# fauna
occ_filter_date_precision <- occ_filter_date %>% 
    dplyr::mutate(precision_filter = ifelse(longitude %>% as.character() %>% stringr::str_split_fixed(., pattern = "[.]", n = 2) %>% .[, 2] %>% stringr::str_length() >= 3 &
                                                latitude %>% as.character() %>% stringr::str_split_fixed(., pattern = "[.]", n = 2) %>% .[, 2] %>% stringr::str_length() >= 3, 
                                            TRUE, FALSE)) %>% 
    dplyr::mutate(precision_filter_date = ifelse(date_filter, TRUE, ifelse(is.na(year) & precision_filter, TRUE, FALSE)))
occ_filter_date_precision

# south america filter ---------------------------------------------------

# vector
sa <- sf::st_union(
    rnaturalearth::ne_countries(scale = 10, continent = "South America"),
    rnaturalearth::ne_countries(scale = 10, country = "France")) %>% 
    sf::st_crop(rnaturalearth::ne_countries(scale = 10, continent = "South America")) %>% 
    sf::st_union() %>% 
    sf::st_as_sf() %>% 
    dplyr::mutate(fid = 1)
sa
plot(sa, col = "gray")

# fauna
occ_filter_date_precision_sa <- occ_filter_date_precision %>% 
    dplyr::mutate(lon = longitude, lat = latitude) %>% 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
    sf::st_join(sa) %>% 
    dplyr::mutate(america_filter = ifelse(is.na(fid), FALSE, TRUE)) %>% 
    dplyr::select(-fid) %>% 
    sf::st_drop_geometry()
occ_filter_date_precision_sa

# spatial filter -----------------------------------------------------

# spatial distance filter
occ_filter_date_precision_sa_spatial <- tibble::tibble()
for(i in species_list){
    
    # information
    print(i)
    
    # filter
    occ_i <- dplyr::filter(occ_filter_date_precision_sa, species == i)
    
    # thin
    thinned_data_i <- spThin::thin(loc.data = occ_i, 
                                   lat.col = "latitude", 
                                   long.col = "longitude", 
                                   spec.col = "species", 
                                   thin.par = 5, 
                                   reps = 1, 
                                   locs.thinned.list.return = TRUE, 
                                   write.files = FALSE, 
                                   write.log.file = FALSE) %>% 
        .[[1]] %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(spatial_filter = TRUE) %>% 
        dplyr::rename(longitude = Longitude,
                      latitude = Latitude)
    thinned_data_i
    
    # join
    occ_filter_date_precision_sa_spatial_i <- occ_i %>% 
        dplyr::left_join(thinned_data_i) %>% 
        dplyr::mutate(spatial_filter = ifelse(is.na(spatial_filter), FALSE, spatial_filter))
    
    # bind
    occ_filter_date_precision_sa_spatial <- dplyr::bind_rows(occ_filter_date_precision_sa_spatial, 
                                                             occ_filter_date_precision_sa_spatial_i)
    
}
occ_filter_date_precision_sa_spatial

# bias filter -------------------------------------------------------------

# bias
occ_filter_date_precision_sa_spatial_bias <- CoordinateCleaner::clean_coordinates(
    x = occ_filter_date_precision_sa_spatial, 
    species = "species",
    lon = "longitude", 
    lat = "latitude",
    tests = c("capitals", # radius around capitals
              "centroids", # radius around country and province centroids
              "duplicates", # records from one species with identical coordinates
              "equal", # equal coordinates
              "gbif", # radius around GBIF headquarters
              "institutions", # radius around biodiversity institutions
              # "outliers", # remove outliers
              "seas", # in the sea
              "urban", # within urban area
              "validity", # outside reference coordinate system
              "zeros" # plain zeros and lat = lon
    ),
    capitals_rad = 2000,
    centroids_rad = 2000,
    centroids_detail = "both",
    inst_rad = 100,
    outliers_method = "quantile",
    outliers_mtp = 5,
    outliers_td = 1000,
    outliers_size = 10,
    range_rad = 0,
    zeros_rad = 0.5,
    capitals_ref = NULL,
    centroids_ref = NULL,
    country_ref = NULL,
    country_refcol = "countryCode",
    inst_ref = NULL,
    range_ref = NULL,
    # seas_ref = continent_border,
    seas_scale = 110,
    urban_ref = NULL,
    value = "spatialvalid") %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(.cen = case_when(longitude == -52.8731 & latitude == -10.8339 ~ FALSE, .default = .cen),
                  .summary = case_when(longitude == -52.8731 & latitude == -10.8339 ~ FALSE, .default = .summary)) %>%     dplyr::mutate(lon = longitude, lat = latitude, bias_filter = .summary) %>% 
    dplyr::mutate(bias_filter = .summary)
occ_filter_date_precision_sa_spatial_bias

# filter ------------------------------------------------------------------

# filter
occ_cleaned <- occ_filter_date_precision_sa_spatial_bias %>% 
    dplyr::filter(precision_filter_date == TRUE,
                  america_filter == TRUE,
                  spatial_filter == TRUE,
                  bias_filter == TRUE) %>% 
    dplyr::select(-id) %>% 
    tibble::rowid_to_column(var = "id")
occ_cleaned

# export ------------------------------------------------------------------

# export
readr::write_csv(occ_cleaned, "01_data/01_occurrences/02_clean/00_occ_cleaned.csv")

occ_cleaned <- readr::read_csv("01_data/01_occurrences/02_clean/00_occ_cleaned.csv")
occ_cleaned_n <- occ_cleaned %>% 
    count(species, sort = TRUE)

# end ---------------------------------------------------------------------