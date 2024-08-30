#' ---
#' title: occurrences - download - data papers
#' author: mauricio vancine
#' date: 2024-01-30
#' ---

# prepare r -------------------------------------------------------------

# packages
library(tidyverse)
library(janitor)
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
fauna_list <- readr::read_csv2("01_data/01_occurrences/00_species_list/species_list.csv") %>% 
    dplyr::pull(species)
fauna_list

# download ----------------------------------------------------------------

## atlantic camtrap ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.1998&file=ecy1998-sup-0001-DataS1.zip"
 
# # unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy1998-sup-0001-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")
 
# import sites
ca_si <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_CAMTRAPS_1-0_LOCATION.csv") %>%
    dplyr::select(location_id, X, Y)
ca_si

# import sites names
ca_si_na <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_CAMTRAPS_1-0_SURVEY.csv") %>% 
    dplyr::left_join(., ca_si, by = "location_id") %>% 
    dplyr::select(location_id, survey_id, yearfinish,  X, Y)
ca_si_na

# import species names
ca_sp_na <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_CAMTRAPS_1-0_SPECIES.csv") %>% 
    dplyr::select(order:genus, species_name, species_code)
ca_sp_na

# import species
ca_sp <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_CAMTRAPS_1-0_RECORDS.csv") %>% 
    dplyr::select(survey_id, species_code, presence_absence) %>% 
    dplyr::filter(presence_absence == 1) %>% 
    dplyr::mutate(species_code = ifelse(species_code == "Potu_flav", "Poto_flav", species_code)) %>% 
    dplyr::left_join(., ca_sp_na, by = "species_code") %>% 
    dplyr::select(-c(presence_absence,  species_code))
ca_sp

# join data
ca_occ <- ca_sp %>% 
    dplyr::left_join(ca_si_na, by = "survey_id") %>% 
    dplyr::mutate(species = species_name,
                  longitude = as.numeric(X), 
                  latitude = as.numeric(Y), 
                  year = as.numeric(yearfinish),
                  source = "atlantic_camtrap") %>%
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
ca_occ

## atlantic mammals ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2785&file=ecy2785-sup-0001-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy2785-sup-0001-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import data
lm_occ <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_MAMMAL_MID_LARGE _assemblages_and_sites.csv") %>% 
    dplyr::mutate(species = Actual_species_Name,
                  longitude = as.numeric(Longitude), 
                  latitude = as.numeric(Latitude), 
                  year = as.numeric(Year_finish),
                  source = "atlantic_large_mammals") %>%
    dplyr::select(species, longitude, latitude, year, source, Precision) %>% 
    tidyr::drop_na(longitude, latitude) %>% 
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::filter(species %in% fauna_list,
                  Precision == "precise") %>% 
    tidyr::drop_na(longitude, latitude) %>% 
    dplyr::select(-Precision)
lm_occ

## atlantic small mammals ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.1893&file=ecy1893-sup-0002-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy1893-sup-0002-datas1.zip",
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import sites
sm_occ_si <- data.table::fread("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_SM_Study_Site.csv") %>% 
    dplyr::rename(ID = `ID\xa0`) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(Precision == "Precise") %>% 
    dplyr::select(ID, Reference_number, Latitude, Longitude) %>% 
    dplyr::mutate(id = paste(ID, sub(" / ", "_", Reference_number), sep = "_")) %>% 
    dplyr::select(id, Latitude, Longitude)
sm_occ_si

# import species
sm_occ_sp <- data.table::fread("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_SM_Capture.csv") %>%
    dplyr::rename(ID = `ID\xa0`) %>% 
    tibble::as_tibble() %>%
    dplyr::select(ID, Reference_number, Order, Genus, Actual_species_name, Year_finish) %>% 
    dplyr::mutate(id = paste(ID, sub(" / ", "_", Reference_number), sep = "_")) %>% 
    dplyr::select(id, Order, Genus, Actual_species_name, Year_finish)
sm_occ_sp

# join data
sm_occ <- sm_occ_sp %>% 
    dplyr::left_join(sm_occ_si, by = "id") %>%
    dplyr::mutate(order = Order, 
                  genus = str_split(Actual_species_name, " ", simplify = TRUE)[, 1],
                  species = Actual_species_name,
                  longitude = as.numeric(Longitude), 
                  latitude = as.numeric(Latitude), 
                  year = as.numeric(Year_finish),
                  source = "atlantic_small_mammals") %>% 
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
sm_occ

## atlantic small mammals abundance ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2005&file=ecy2005-sup-0001-datas1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy2005-sup-0001-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import sites
sm_ab_si <- data.table::fread("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/Localities.csv") %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(Study_year = as.numeric(ifelse(str_length(Study_year) > 4, 
                                                 stringr::str_sub(Study_year, -4, -1),
                                                 Study_year))) %>% 
    dplyr::select(SampleID, Latitude, Longitude, Study_year)
sm_ab_si

# import species
sm_ab_sp <- data.table::fread("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/Mammal_Communities.csv") %>% 
    tibble::as_tibble() %>% 
    dplyr::select(SampleID, Valid_Species) %>% 
    dplyr::filter(Valid_Species != "")
sm_ab_sp

# join data
sm_ab_occ <- dplyr::left_join(sm_ab_sp, sm_ab_si, by = "SampleID") %>% 
    dplyr::mutate(species = Valid_Species, 
                  genus = str_split(Valid_Species, " ", simplify = TRUE)[, 1],
                  name = Valid_Species, 
                  longitude = Longitude, 
                  latitude = Latitude, 
                  year = as.numeric(Study_year),
                  source = "atlantic_small_mammals_abu") %>% 
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
sm_ab_occ

## atlantic mammals uprb ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2107&file=ecy2107-sup-0002-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy2107-sup-0002-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import sites
muprb_si <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/Mammals_UPRB_study_sites.csv")
muprb_si

muprb_sp <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/Mammals_UPRB_species.csv") %>% 
    dplyr::select(-c(order:genus, frequency_of_occurrence:Brazilian_status)) %>% 
    t() %>% 
    as.data.frame()
colnames(muprb_sp) <- muprb_sp[1, ]
muprb_sp <- muprb_sp[-1, ]
muprb_sp$site <- rownames(muprb_sp)
muprb_sp <- tidyr::pivot_longer(muprb_sp, cols = -site) %>% 
    dplyr::filter(value > 0) %>% 
    dplyr::select(-value) %>% 
    dplyr::mutate(site = as.numeric(stringr::str_replace_all(site, "site_", "")))
muprb_sp

# join data
muprb_occ <- dplyr::left_join(muprb_sp, muprb_si) %>% 
    dplyr::mutate(species = name, 
                  year = as.numeric(year_finish),
                  source = "atlantic_mammal_uprb") %>% 
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
muprb_occ

## atlantic mammal traits ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2106&file=ecy2106-sup-0002-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy2106-sup-0002-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import
mt_occ <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_TR_all_data.csv") %>% 
    janitor::clean_names() %>% 
    dplyr::mutate(species = binomial,
                  source = "atlantic_mammal_traits") %>%
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
mt_occ

## atlantic frugivory ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.1818&file=ecy1818-sup-0002-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy1818-sup-0002-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import
fru_fauna_occ <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ATLANTIC_frugivory.csv") %>% 
    janitor::clean_names() %>% 
    dplyr::filter(precision == "preciso") %>% 
    dplyr::mutate(species = frugivore_species,
                  year = NA,
                  source = "atlantic_frugivory") %>%
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
fru_fauna_occ

## brazil roadkill ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2464&file=ecy2464-sup-0001-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy2464-sup-0001-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

br_occ <- readr::read_csv("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/Brazil_Roadkill_20180527.csv") %>% 
    dplyr::mutate(species = Scientific_name, 
                  longitude = Long, 
                  latitude = Lat,
                  year = Year,
                  source = "brasil_roadkill") %>%
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
br_occ

## neotropical alien mammals ----

# download data
# "https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3115&file=ecy3115-sup-0001-DataS1.zip"

# unzip
# unzip(zipfile = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/ecy3115-sup-0001-datas1.zip", 
#       exdir = "../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers")

# import
neo_alien_occ <- readr::read_delim("../LEDDiv/ms-atlantic-forest-networks/01_data/01_occurrences/01_raw/02_fauna/04_data_papers/DataS1/NEOTROPICAL_ALIEN_MAMMALS_OCCURENCE.csv", delim = ";") %>% 
    dplyr::filter(PRECISION %in% c(as.character(0:100))) %>%
    dplyr::mutate(species = SPECIES,
                  longitude = LONG_X, 
                  latitude = LAT_Y,
                  year = as.numeric(RECORD_YEAR),
                  source = "neotropical_alien_mammals") %>%
    dplyr::select(species, longitude, latitude, year, source) %>% 
    dplyr::mutate(species = str_trim(species)) %>% 
    dplyr::filter(species %in% fauna_list) %>% 
    tidyr::drop_na(longitude, latitude)
neo_alien_occ

## combine ----
occ_data_papers_fauna <- dplyr::bind_rows(ca_occ, lm_occ,
                                          sm_occ, sm_ab_occ, muprb_occ, 
                                          mt_occ, fru_fauna_occ, 
                                          br_occ, neo_alien_occ)
occ_data_papers_fauna

occ_data_papers_fauna_v <- sf::st_as_sf(occ_data_papers_fauna, coords = c("longitude", "latitude"), crs = 4326)
occ_data_papers_fauna_v

tm_shape(li[li$iso_a2 == "BR", ]) +
    tm_polygons() +
    tm_shape(occ_data_papers_fauna_v) +
    tm_dots()

## export ----
readr::write_csv(occ_data_papers_fauna, "01_data/01_occurrences/01_raw/04_occ_raw_data_paper_fauna.csv")

# end ---------------------------------------------------------------------