#' ---
#' title: occurrences - download - integrate
#' author: mauricio vancine
#' date: 2024-08-22
#' aim: integrated occurrences
#' ---

# prepare r ---------------------------------------------------------------

# packages
library(tidyverse)
library(furrr)
library(future)
library(tmap)
library(viridis)

# import data -------------------------------------------------------------

## specieslink and gbif ----
occ_specieslink_spocc <- readr::read_csv("01_data/01_occurrences/01_raw/01_occ_raw_splink_spocc.csv") %>% 
    dplyr::rename(species = species_searched,
                  year = date,
                  source = prov) %>% 
    dplyr::select(species, longitude, latitude, year, source)
occ_specieslink_spocc

## sibbr ----
occ_sibbr <- readr::read_csv("01_data/01_occurrences/01_raw/02_occ_raw_sibbr.csv") %>% 
    dplyr::rename(species = name,
                  year = date) %>% 
    dplyr::select(species, longitude, latitude, year) %>% 
    dplyr::mutate(source = "sibbr")
occ_sibbr

## portal biodiversidade ----
occ_portalbio <- readr::read_csv("01_data/01_occurrences/01_raw/03_occ_raw_portalbio.csv") %>% 
    dplyr::select(species, longitude, latitude, year) %>% 
    dplyr::mutate(source = "portalbio")
occ_portalbio

## data papers ----
occ_data_papers <- readr::read_csv("01_data/01_occurrences/01_raw/04_occ_raw_data_paper.csv")
occ_data_papers

## paper ----
occ_paper <- readr::read_csv("01_data/01_occurrences/01_raw/05_occs_paper.csv") %>% 
    dplyr::mutate(year = NA, source = "paper")
occ_paper

## combine ----
occ <- dplyr::bind_rows(occ_specieslink_spocc, occ_data_papers,
                        occ_paper, occ_sibbr, occ_portalbio) %>% 
    tidyr::drop_na(longitude, latitude) %>% 
    dplyr::filter(longitude > -180 & longitude < 180) %>% 
    dplyr::filter(latitude > -90 & latitude < 90)
occ

## export ----
readr::write_csv(occ, "01_data/01_occurrences/01_raw/00_occ_raw.csv")

# end ---------------------------------------------------------------------
