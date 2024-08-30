#' ---
#' title: occurrences - download - sibbr
#' author: mauricio vancine
#' date: 2024-01-30
#' ---

# prepare r -------------------------------------------------------------

# packages
library(tidyverse)

# options
options(timeout = 3e5)

# import data -------------------------------------------------------------

# species list
fauna_list <- readr::read_csv2("01_data/01_occurrences/01_raw/01_fauna/00_species_list/fauna_species_list_frugivore.csv") %>% 
    dplyr::pull(species)
fauna_list

# download ----------------------------------------------------------------

# download
# download.file(url = "https://ipt.sibbr.gov.br/sibbr/resource?r=sibbr_mamiferos_01", 
#               destfile = "01_data/00_occurrences/raw/03_sibbr/sibbr_mamiferos.zip", 
#               mode = "wb")

# unzip
# directories <- sub("sibbr_", "", sub(".zip", "", dir(path = "01_data/00_occurrences/raw/03_sibbr/", pattern = ".zip$")))
# directories
# 
# for(i in directories){
#     
#     print(i)
#     unzip(zipfile = paste0("01_data/00_occurrences/raw/03_sibbr/sibbr_", i, ".zip"), 
#           exdir = paste0("01_data/00_occurrences/raw/03_sibbr/", i))
# }


# filter ------------------------------------------------------------------

# mamiferos
occ_data_sibbr_mamiferos <- readr::read_delim("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/03_sibbr/mamiferos/occurrence.txt", delim = "\t") %>% 
    dplyr::mutate(species_searched = scientificName,
                  name = scientificName,
                  longitude = as.numeric(decimalLongitude),
                  latitude = as.numeric(decimalLatitude),
                  prov = "sibbr",
                  date = as.numeric(year),
                  key = as.character(catalogNumber)) %>% 
    tidyr::drop_na(longitude, latitude) %>% 
    dplyr::select(species_searched, name, longitude, latitude, prov, date, key) %>% 
    dplyr::filter(species_searched %in% fauna_list) %>% 
    dplyr::filter(longitude > -180, 
                  longitude < 180,
                  latitude > -90,
                  latitude < 90)
occ_data_sibbr_mamiferos

# export
readr::write_csv(occ_data_sibbr_mamiferos, "01_data/01_occurrences/01_raw/02_occ_raw_sibbr_fauna.csv")

# end ---------------------------------------------------------------------