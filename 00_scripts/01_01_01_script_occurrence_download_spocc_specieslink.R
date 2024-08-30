#' ---
#' title: occurrences - download - spocc and specieslink 
#' author: mauricio vancine
#' date: 2024-01-30
#' ---

# prepare r -------------------------------------------------------------

# packages
library(tidyverse)
library(janitor)
library(jsonlite)
library(lubridate)
library(parallelly)
library(doParallel)
library(foreach)
library(furrr)
library(future)
library(spocc)
library(rnaturalearth)
library(spData)
library(sf)
library(tmap)

# options
options(timeout = 3e5)

# import data -------------------------------------------------------------

# import south america limite
li <- spData::world
li

tm_shape(li) +
    tm_polygons()

# species list
species_list <- readr::read_csv2("01_data/01_occurrences/00_species_list/species_list.csv") %>% 
    dplyr::pull(species)
species_list

# download fauna --------------------------------------------------------

# occ
# doParallel::registerDoParallel(parallelly::availableCores(omit = 2))
foreach::foreach(i=species_list) %do% {
    
    # species
    print(i)
    
    ## specieslink ----
    # information
    print("splink")
    
    # download
    occ_splink <- jsonlite::fromJSON(
        paste0("https://specieslink.net/ws/1.0/search?scientificname=", 
               tolower(gsub(" ", "+", i)), 
               "&apikey=aXGEJtnQW12sPuSyKMX7&offset=0&limit=50000"))$features$properties
    
    # conditional without data
    if (length(occ_splink) == 0){
        
        occ_splink_data <- tibble::tibble(species_searched = i,
                                          name = NA,
                                          longitude = NA,
                                          latitude = NA,
                                          prov = "specieslink",
                                          date = NA,
                                          key = NA)
        
        # conditional with data and year
    } else{
        if ("yearcollected" %in% colnames(occ_splink)) {
            occ_splink_data <- occ_splink %>% 
                tidyr::drop_na(decimallongitude, decimallatitude) %>% 
                dplyr::mutate(species_searched = i,
                              name = scientificname,
                              longitude = as.numeric(decimallongitude),
                              latitude = as.numeric(decimallatitude),
                              prov = "specieslink",
                              date = as.numeric(yearcollected),
                              key = as.character(catalognumber)) %>% 
                dplyr::select(species_searched, name, longitude, latitude, prov, date, key)
            
        } else {
            # Se 'yearcollected' não estiver presente, faça algo apropriado (crie uma coluna de data vazia, por exemplo)
            occ_splink_data <- occ_splink %>% 
                tidyr::drop_na(decimallongitude, decimallatitude) %>% 
                dplyr::mutate(species_searched = i,
                              name = scientificname,
                              longitude = as.numeric(decimallongitude),
                              latitude = as.numeric(decimallatitude),
                              prov = "specieslink",
                              date = NA,
                              key = as.character(catalognumber)) %>% 
                dplyr::select(species_searched, name, longitude, latitude, prov, date, key)
        }
    }
    
    
    # spocc ----
    # information
    print("spocc")
    
    # download
    occ_spocc <- spocc::occ(query = i, 
                            from = c("gbif", "vertnet", "idigbio", "ecoengine"),
                            has_coords = TRUE,
                            limit = 1e6,
                            throw_warnings = FALSE)
    
    # data
    occ_spocc_data <- spocc::occ2df(occ_spocc)
    
    # conditional without data
    if(nrow(occ_spocc_data) == 0){
        
        occ_spocc_data <- tibble::tibble(species_searched = i,
                                         name = NA,
                                         longitude = NA,
                                         latitude = NA,
                                         prov = "spocc",
                                         date = NA,
                                         key = NA)
        
        # conditional without year  
    } else if(!"date" %in% colnames(occ_spocc_data)){
        
        occ_spocc_data <- occ_spocc_data %>% 
            dplyr::mutate(species_searched = i, .before = 1) %>% 
            dplyr::mutate(longitude = as.numeric(longitude),
                          latitude = as.numeric(latitude),
                          date = NA,
                          key = as.character(key)) %>% 
            dplyr::select(species_searched, name, longitude, latitude, prov, date, key)
        
        # conditional with data and year
    } else{
        
        occ_spocc_data <- occ_spocc_data %>% 
            dplyr::mutate(species_searched = i, .before = 1) %>% 
            dplyr::mutate(longitude = as.numeric(longitude),
                          latitude = as.numeric(latitude),
                          date = lubridate::year(occ_spocc_data$date),
                          key = as.character(key)) %>% 
            dplyr::select(species_searched, name, longitude, latitude, prov, date, key)
        
    }
    
    # combine data ----
    occ_data <- dplyr::bind_rows(occ_splink_data, occ_spocc_data)
    
    # export 
    readr::write_csv(occ_data, 
                     paste0("01_data/01_occurrences/01_raw/spocc_specieslink/01_spocc_specieslink_", sub(" ", "_", tolower(i)), ".csv"))
    
}
# doParallel::stopImplicitCluster()

# integrated --------------------------------------------------------------

# import
plan(multisession, workers = parallelly::availableCores(omit = 2))
occ_data <- dir(path = "01_data/01_occurrences/01_raw/01_spocc_specieslink/", 
                pattern = ".csv", full.names = TRUE) %>% 
    furrr::future_map_dfr(readr::read_csv, col_types = "ccddcdc")
occ_data

# vector
occ_data_v <- occ_data %>% 
    tidyr::drop_na(longitude, latitude) %>% 
    dplyr::mutate(lon = longitude,
                  lat = latitude) %>% 
    dplyr::filter(lon > -180 & lon < 180) %>% 
    dplyr::filter(lat > -90 & lat < 90) %>% 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
occ_data_v

# map
tm_shape(li) +
    tm_polygons() +
    tm_shape(occ_data_v) +
    tm_bubbles(size = .2, 
               col = "species_searched", 
               col.legend = tm_legend_hide())

# export
readr::write_csv(occ_data, "01_data/01_occurrences/01_raw/01_occ_raw_splink_spocc.csv")

# end ---------------------------------------------------------------------