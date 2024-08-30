#' ---
#' title: occurrences - download - portal da biodiversidade
#' author: mauricio vancine
#' date: 2024-01-30
#' ---

# prepare r -------------------------------------------------------------

# packages
library(tidyverse)
library(janitor)
library(parallelly)
library(future)
library(furrr)
library(tmap)

# options
options(timeout = 3e5)

# import data -------------------------------------------------------------

# import south america limite
li <- spData::world
li

tm_shape(li) +
    tm_polygons()

# species
fauna_list <- readr::read_csv2("01_data/01_occurrences/00_species_list/species_list.csv") %>% 
    dplyr::pull(species)
fauna_list

# import
list_files_portalbio <- dir(path = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/05_plataforma_biodiversidade", pattern = ".csv", recursive = TRUE, full.names = TRUE) %>% 
    stringr::str_subset(pattern = "portalbio_export")
list_files_portalbio

occ_portalbio <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/05_plataforma_biodiversidade/portalbio_one_species.csv")
for(i in list_files_portalbio){
    
    occ_portalbio_i <- vroom::vroom(i, delim = ";", col_types = cols()) %>% 
        janitor::clean_names() %>%
        dplyr::mutate(species = especie,
                      longitude = as.numeric(longitude),
                      latitude = as.numeric(latitude),
                      year = lubridate::year(lubridate::dmy(data_do_registro)),
                      source = "portalbio") %>% 
        dplyr::select(species, longitude, latitude, year, source)
    occ_portalbio <- dplyr::bind_rows(occ_portalbio, occ_portalbio_i)
    
}
occ_portalbio

# filter
occ_portalbio_fauna <- occ_portalbio %>% 
    dplyr::filter(species %in% fauna_list)
occ_portalbio_fauna

# vector
occ_portalbio_fauna_v <- occ_portalbio_fauna %>% 
    dplyr::mutate(lon = longitude,
                  lat = latitude) %>% 
    dplyr::filter(lon > -180 & lon < 180) %>% 
    dplyr::filter(lat > -90 & lat < 90) %>% 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
occ_portalbio_fauna_v

# maps
tm_shape(li[li$iso_a2 == "BR", ]) +
    tm_polygons() +
    tm_shape(occ_portalbio_fauna_v) +
    tm_dots()

# export
readr::write_csv(occ_portalbio_fauna, "01_data/01_occurrences/01_raw/03_occ_raw_portalbio_fauna.csv")

# end ---------------------------------------------------------------------