#' ----
#' title: flexsdm
#' author: mauricio vancine
#' date: 30/08/2024
#' operational system: gnu/linux - ubuntu - pop_os
#' link: https://sjevelazco.github.io/flexsdm/articles/v02_modeling.html
#' ----

# prepare r -------------------------------------------------------------

# packages
library(tidyverse)
library(terra)
library(flexsdm)
library(tmap) # remotes::install_github("r-tmap/tmap@v4")

# import data -------------------------------------------------------------

## af
af <- terra::vect("01_data/02_variables/00_limit/merge_limites_MA_buffer20km_WGS84.shp")
af

## occurrences ----
occ <- readr::read_csv("01_data/01_occurrences/02_clean/00_occ_cleaned.csv") %>% 
    dplyr::rename(x = longitude, y = latitude) %>% 
    dplyr::select(species, x, y) %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(n = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(n >= 10)
occ

tm_shape(terra::vect(occ, geom = c("x", "y"), crs = "EPSG:4326")) +
    tm_dots() +
    tm_shape(af) +
    tm_borders(col = "red")

## variables ----

### current ----
var_c <- dir(path = "01_data/02_variables/01_climate_present/01_cropped", 
             pattern = "_sa.tif", full.names = TRUE) %>% 
    terra::rast()
var_c

tm_shape(var_c[[1]]) +
    tm_raster(col.legend = tm_legend(, show = FALSE),
              col.scale = tm_scale_continuous(values = "Spectral")) +
    tm_shape(af) +
    tm_borders(col = "red")

### future ----
var_f <- dir(path = "01_data/02_variables/02_climate_future/02_cropped_mean/", 
             pattern = ".tif$", full.names = TRUE) %>% 
    terra::rast()
var_f

tm_shape(var_f[[1]]) +
    tm_raster(col.legend = tm_legend(, show = FALSE),
              col.scale = tm_scale_continuous(values = "Spectral")) +
    tm_shape(af) +
    tm_borders(col = "red")

# sdm ---------------------------------------------------------------------

# parameters
ca_buffer <- 300000
n_cores <- 6
fut <- expand.grid(
    year = c("2011-2040", "2041-2070", "2071-2100"),
    ssp = c("ssp370", "ssp585")) %>% 
    apply(MARGIN = 1, FUN = paste0, collapse = "_")

# loop for each species
for(i in unique(occ$species)){
    
    # information
    print(i)
    sp_name <- sub(" ", "_", tolower(i))
    sp_dir <- paste0("02_results/", sp_name, "/")
    
    # directory
    dir.create(path = sp_dir)
    
    ## species filter ----
    occ_i <- occ %>% 
        dplyr::filter(species == i) %>% 
        dplyr::mutate(pr_ab = 1) %>% 
        dplyr::select(x, y, pr_ab) %>% 
        tibble::rowid_to_column()
    readr::write_csv(occ_i, paste0(sp_dir, "01_occ_pr", sp_name, ".csv"))
    
    ## calibration area ----
    ca <- flexsdm::calib_area(data = occ_i, 
                              x = "x", 
                              y = "y", 
                              method = c("buffer", width = ca_buffer),
                              crs = crs(var_c))
    terra::writeVector(ca, paste0(sp_dir, "02_var_calibration_area_", sp_name, ".gpkg"), overwrite = TRUE)
    
    map_occ <- tm_shape(var_f[[1]]) +
        tm_raster(col.legend = tm_legend(show = FALSE),
                  col.scale = tm_scale_categorical(values = "gray80")) +
        tm_shape(af) +
        tm_borders() +
        tm_shape(ca) +
        tm_borders(col = "red") +
        tm_shape(terra::vect(occ_i, geom = c("x", "y"), crs = "EPSG:4326")) +
        tm_dots() +
        tm_title(text = i, fontface = "italic")
    tmap::tmap_save(map_occ, paste0(sp_dir, "01_occ_map_", sp_name, ".png"))
    
    ## adjust variables ----
    var_c_ca <- terra::crop(var_c, ca, mask = TRUE)
    
    ## collinearity variables ----
    var_c_ca_colin <- flexsdm::correct_colinvar(env_layer = var_c_ca, method = c("vif", th = "2"))
    var_c_ca_sel_names <- var_c_ca_colin$vif_table$Variables
    var_c_ca_sel <- var_c_ca_colin$env_layer
    readr::write_csv(var_c_ca_colin$vif_table, paste0(sp_dir, "02_var_vif_", sp_name, ".csv"))
    
    var_f_sel <- terra::subset(var_f, grep(paste0(names(var_c_ca_sel), collapse = "|"), names(var_f), value = TRUE))
    var_f_ca_sel <- terra::crop(var_f_sel, ca, mask = TRUE)
    
    ## occurrences environmental filtering ----
    occ_i_filt <- flexsdm::occfilt_env(data = occ_i, 
                                       x = "x", 
                                       y = "y", 
                                       id = "rowid", 
                                       env_layer = var_c_ca_sel, 
                                       nbins = 12) %>% 
        dplyr::mutate(pr_ab = 1)
    readr::write_csv(occ_i_filt, paste0(sp_dir, "01_occ_", sp_name, "_filt.csv"))    
    
    ## more than 20 occs ----
    if(nrow(occ_i_filt) >= 20){
        
        ### data partitioning ----
        n_part <- 4
        occ_i_filt_sblock <- NA
        while(any(is.na(occ_i_filt_sblock))){
            occ_i_filt_sblock <- flexsdm::part_sblock(
                data = occ_i_filt,
                env_layer = var_c_ca_sel,
                pr_ab = "pr_ab",
                x = "x",
                y = "y",
                n_part = n_part)
            n_part <- n_part - 1
        }
        
        occ_i_filt_part <- occ_i_filt_sblock$part
        occ_i_filt_block_layer <- flexsdm::get_block(env_layer = var_c_ca_sel, best_grid = occ_i_filt_sblock$grid)
        
        readr::write_csv(occ_i_filt_part, paste0(sp_dir, "01_occ_", sp_name, "_filt_part.csv"))    
        terra::writeRaster(occ_i_filt_block_layer, paste0(sp_dir, "01_occ_", sp_name, "_block_layer.tif"), overwrite = TRUE)
        
        cl <- c("#64146D", "#9E2962", "#F47C15", "#FCFFA4")
        map_occ_block <- tm_shape(var_f[[1]]) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = "gray80")) +
            tm_shape(occ_i_filt_block_layer) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = cl)) +
            tm_shape(af) +
            tm_borders() +
            tm_shape(ca) +
            tm_borders(col = "red") +
            tm_shape(terra::vect(occ_i_filt, geom = c("x", "y"), crs = "EPSG:4326")) +
            tm_dots() +
            tm_title(text = i)
        tmap::tmap_save(map_occ_block, paste0(sp_dir, "01_occ_map_", sp_name, "_block.png"))
        
        ### pseudo-absence and background sampling ----
        psa <- lapply(1:4, function(x){
            flexsdm::sample_pseudoabs(
                data = occ_i_filt_part,
                x = "x",
                y = "y",
                n = sum(occ_i_filt_part$.part == x),
                method = "random",
                rlayer = occ_i_filt_block_layer,
                maskval = x,
                calibarea = ca)}) %>%
            dplyr::bind_rows()
        psa <- flexsdm::sdm_extract(data = psa, x = "x", y = "y", env_layer = occ_i_filt_block_layer)
        readr::write_csv(psa, paste0(sp_dir, "01_occ_psa_", sp_name, ".csv"))    
        
        bg <- lapply(1:4, function(x){
            flexsdm::sample_background(
                data = occ_i_filt_part,
                x = "x",
                y = "y",
                n = sum(occ_i_filt_part$.part == x) * 10,
                method = "random",
                rlayer = occ_i_filt_block_layer,
                maskval = x,
                calibarea = ca)}) %>%
            dplyr::bind_rows()
        bg_part <- flexsdm::sdm_extract(data = bg, x = "x", y = "y", env_layer = occ_i_filt_block_layer)
        readr::write_csv(bg, paste0(sp_dir, "01_occ_bg_", sp_name, ".csv"))   
        
        map_psa_block <- tm_shape(var_c[[1]]) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = "gray80")) +
            tm_shape(occ_i_filt_block_layer) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = cl)) +
            tm_shape(af) +
            tm_borders() +
            tm_shape(ca) +
            tm_borders(col = "red") +
            tm_shape(terra::vect(psa, geom = c("x", "y"), crs = "EPSG:4326")) +
            tm_dots(fill = "gray") +
            tm_title(text = i)
        tmap::tmap_save(map_psa_block, paste0(sp_dir, "01_occ_map_", sp_name, "_psa_block.png"))
        
        map_bg_block <- tm_shape(var_c[[1]]) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = "gray80")) +
            tm_shape(occ_i_filt_block_layer) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = cl)) +
            tm_shape(af) +
            tm_borders() +
            tm_shape(ca) +
            tm_borders(col = "red") +
            tm_shape(terra::vect(bg, geom = c("x", "y"), crs = "EPSG:4326")) +
            tm_dots(fill = "gray") +
            tm_title(text = i)
        tmap::tmap_save(map_bg_block, paste0(sp_dir, "01_occ_map_", sp_name, "_bg_block.png"))
        
        ### extracting environmental values ----
        occ_i_filt_part_pa <- dplyr::bind_rows(occ_i_filt_part, psa)
        
        occ_i_filt_part_pa_data <- flexsdm::sdm_extract(
            data = occ_i_filt_part_pa,
            x = "x",
            y = "y",
            env_layer = var_c_ca_sel,
            filter_na = TRUE)
        readr::write_csv(occ_i_filt_part_pa_data, paste0(sp_dir, "02_var_data_occ_psa_", sp_name, ".csv"))    
        
        bg_part_data <- flexsdm::sdm_extract(
            data = bg_part,
            x = "x",
            y = "y",
            env_layer = var_c_ca_sel,
            filter_na = TRUE)
        readr::write_csv(bg_part_data, paste0(sp_dir, "02_var_data_bg_", sp_name, ".csv"))    
        
        ### models ----
        
        #### fit models ----
        glm_f <- flexsdm::fit_glm(
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            thr = "max_sens_spec")
        glm_f_plot <- flexsdm::p_pdp(model = glm_f$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_glm_", sp_name, ".png"), 
               plot = glm_f_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        k <- -1
        gam_f <- NULL
        while(is.null(gam_f)){
            gam_f <- flexsdm::fit_gam(
                data = occ_i_filt_part_pa_data,
                response = "pr_ab",
                predictors = var_c_ca_sel_names,
                partition = ".part",
                thr = "max_sens_spec",
                k = k)
            k <- k + 1
        }
        gam_f_plot <- flexsdm::p_pdp(model = gam_f$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_gam_", sp_name, ".png"), 
               plot = gam_f_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        gau_f <- flexsdm::fit_gau(
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            background = bg_part_data,
            thr = "max_sens_spec")
        gau_f_plot <- flexsdm::p_pdp(model = gau_f$model, training_data = occ_i_filt_part_pa, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_gau_", sp_name, ".png"), 
               plot = gau_f_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        #### tune models ----
        gbm_tune <- expand.grid(
            n.trees = c(20, 50, 100),
            shrinkage = c(0.1, 0.5, 1),
            n.minobsinnode = c(1, 3, 5, 7, 9))
        gbm_t <- flexsdm::tune_gbm(            
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            grid = gbm_tune,
            thr = "max_sens_spec",
            metric = "TSS",
            n_cores = n_cores)
        gbm_t_plot <- flexsdm::p_pdp(model = gbm_t$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_gbm_", sp_name, ".png"), 
               plot = gbm_t_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        max_tune <- expand.grid(
            regmult = seq(.5, 4, .5),
            classes = c("l", "lq", "h", "lqh", "lqhp", "lqhpt"))
        max_t <- flexsdm::tune_max(            
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            background = bg_part_data,
            grid = max_tune,
            thr = "max_sens_spec",
            metric = "TSS",
            n_cores = n_cores)
        max_t_plot <- flexsdm::p_pdp(model = max_t$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_max_", sp_name, ".png"), 
               plot = max_t_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        net_tune <- expand.grid(
            size = c(2, 4, 6, 8, 10),
            decay = c(0.001, 0.05, 0.1, 1, 3, 4, 5, 10))
        net_t <- flexsdm::tune_net(           
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            grid =  net_tune,
            thr = "max_sens_spec",
            metric = "TSS",
            n_cores = n_cores)
        net_t_plot <- flexsdm::p_pdp(model = net_t$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_net_", sp_name, ".png"), 
               plot = net_t_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        raf_tune <- expand.grid(mtry = seq(1, 7, 1))
        raf_t <- flexsdm::tune_raf(           
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            grid = raf_tune,
            thr = "max_sens_spec",
            metric = "TSS",
            n_cores = n_cores)
        raf_t_plot <- flexsdm::p_pdp(model = raf_t$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_raf_", sp_name, ".png"), 
               plot = raf_t_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        svm_tune <- expand.grid(
            C = c(2, 4, 8, 16, 20),
            sigma = c(0.01, 0.1, 0.2, 0.3, 0.4))
        svm_t <-  flexsdm::tune_svm(            
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            grid = svm_tune,
            thr = "max_sens_spec",
            metric = "TSS",
            n_cores = n_cores)
        svm_t_plot <- flexsdm::p_pdp(model = svm_t$model, training_data = occ_i_filt_part_pa_data, 
                                     colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
        ggsave(filename = paste0(sp_dir, "03_mod_response_svm_", sp_name, ".png"), 
               plot = svm_t_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
        ### ensemble ----
        model_list <- Filter(function(x) x$performance$TSS_mean > 0, 
                             list(glm_f, gam_f, gau_f, 
                                  gbm_t, max_t, net_t, raf_t, svm_t))
        
        ens_m <- flexsdm::fit_ensemble(
            models = model_list,
            ens_method = "meanw",
            thr_model = "max_sens_spec",
            metric = "TSS")
        readr::write_csv(ens_m$performance, paste0(sp_dir, "04_ens_performance_", sp_name, ".csv"))
        
        model_performance <- flexsdm::sdm_summarize(list(glm_f, gam_f, gau_f, gbm_t, max_t, net_t, raf_t, svm_t))
        readr::write_csv(model_performance, paste0(sp_dir, "03_fit_performance_", sp_name, ".csv"))
        
        ### predict ----
        
        #### current ----
        predict_c <- flexsdm::sdm_predict(
            models = ens_m,
            pred = var_c_ca_sel,
            nchunk = 10,
            thr = "max_sens_spec",
            con_thr = TRUE)
        
        names(predict_c$meanw) <- c(paste0(sp_name, "_meanw_con_c"), paste0(sp_name, "_meanw_con_thr_c"))
        
        terra::writeRaster(predict_c$meanw[[1]], paste0(sp_dir, "05_ens_con_c_", sp_name, ".tif"), overwrite = TRUE)
        terra::writeRaster(predict_c$meanw[[2]], paste0(sp_dir, "05_ens_thr_c_", sp_name, ".tif"), overwrite = TRUE)
        
        #### future ----
        for(j in fut){
            
            var_f_ca_sel_j <- terra::subset(var_f_ca_sel, grep(j, names(var_f_ca_sel)))
            names(var_f_ca_sel_j) <- names(var_c_ca_sel)
            
            predict_f <- flexsdm::sdm_predict(
                models = ens_m,
                pred = var_f_ca_sel_j,
                thr = "max_sens_spec",
                con_thr = TRUE,
                predict_area = NULL)
            
            names(predict_f$meanw) <- c(paste0(sp_name, "_meanw_con_f_", j), paste0(sp_name, "_meanw_con_thr_f_", j))
            
            terra::writeRaster(predict_f$meanw[[1]], paste0(sp_dir, "05_ens_con_f_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            terra::writeRaster(predict_f$meanw[[2]], paste0(sp_dir, "05_ens_thr_f_", sp_name, "_", j, ".tif"), overwrite = TRUE)
        }
        
        ### msdm ----
        
        #### current ----
        predict_c_m <- flexsdm::msdm_posteriori(
            records = occ_i_filt_part_pa,
            x = "x",
            y = "y",
            pr_ab = "pr_ab",
            method = "pres",
            cont_suit = predict_c$meanw[[1]],
            thr = "max_sens_spec",
            buffer = NULL)
        
        names(predict_c_m[[1]]) <- paste0(names(predict_c_m)[1], "_msdm")
        names(predict_c_m[[2]]) <- paste0(names(predict_c_m)[1], "_msdm_thr")
        
        terra::writeRaster(predict_c_m[[1]], paste0(sp_dir, "05_ens_con_c_msdm_", sp_name, ".tif"), overwrite = TRUE)
        terra::writeRaster(predict_c_m[[2]], paste0(sp_dir, "05_ens_con_c_msdm_thr_", sp_name, ".tif"), overwrite = TRUE)
        
        #### future ----
        pred_f <- dir(path = sp_dir, pattern = "ssp", full.names = TRUE) %>% 
            grep("_thr", ., value = TRUE, invert = TRUE) %>% 
            terra::rast()
        
        for(j in fut){
            
            pred_f_j <- terra::subset(pred_f, grep(j, names(pred_f)))
            
            predict_f_m <- flexsdm::msdm_posteriori(
                records = occ_i_filt_part_pa,
                x = "x",
                y = "y",
                pr_ab = "pr_ab",
                method = "pres",
                cont_suit = predict_f$meanw[[1]],
                thr = "max_sens_spec",
                buffer = NULL)
            
            names(predict_f_m[[1]]) <- paste0(names(predict_f_m)[1], "_msdm")
            names(predict_f_m[[2]]) <- paste0(names(predict_f_m)[1], "_msdm_thr")
            
            terra::writeRaster(predict_f_m[[1]], paste0(sp_dir, "05_ens_con_f_msdm_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            terra::writeRaster(predict_f_m[[2]], paste0(sp_dir, "05_ens_con_f_msdm_thr_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            
        }
        
    } else{
        
        ## less than 20 and more than 10 occs ----
        
        ### pseudo-absence and background sampling ----
        psa <- flexsdm::sample_pseudoabs(
            data = occ_i_filt,
            x = "x",
            y = "y",
            n = sum(hespero$pr_ab), 
            method = "random",
            rlayer = var_c_ca_sel,
            calibarea = ca)
        
        bg <- flexsdm::sample_background(
            data = occ_i_filt,
            x = "x",
            y = "y",
            n = sum(hespero$pr_ab) * 10, 
            method = "random",
            rlayer = var_c_ca_sel,
            calibarea = ca)
        
        map_psa <- tm_shape(var_c[[1]]) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = "gray80")) +
            tm_shape(af) +
            tm_borders() +
            tm_shape(ca) +
            tm_borders(col = "red") +
            tm_shape(terra::vect(psa, geom = c("x", "y"), crs = "EPSG:4326")) +
            tm_dots(fill = "gray50") +
            tm_title(text = i)
        tmap::tmap_save(map_psa, paste0(sp_dir, "01_occ_map_", sp_name, "_psa.png"))
        
        map_bg <- tm_shape(var_c[[1]]) +
            tm_raster(col.legend = tm_legend(show = FALSE),
                      col.scale = tm_scale_categorical(values = "gray80")) +
            tm_shape(af) +
            tm_borders() +
            tm_shape(ca) +
            tm_borders(col = "red") +
            tm_shape(terra::vect(bg, geom = c("x", "y"), crs = "EPSG:4326")) +
            tm_dots(fill = "gray50") +
            tm_title(text = i)
        tmap::tmap_save(map_bg, paste0(sp_dir, "01_occ_map_", sp_name, "_bg.png"))
        
        ### data partitioning ----
        occ_i_filt_part <- flexsdm::part_random(
            data = occ_i_filt,
            pr_ab = "pr_ab",
            method = c(method = "rep_kfold", folds = 3, replicates = 5))
        occ_i_filt_part
        
        psa_part <- flexsdm::part_random(
            data = psa,
            pr_ab = "pr_ab",
            method = c(method = "rep_kfold", folds = 3, replicates = 5))
        psa_part
        
        bg_part <- flexsdm::part_random(
            data = bg,
            pr_ab = "pr_ab",
            method = c(method = "rep_kfold", folds = 3, replicates = 5))
        bg_part
        
        readr::write_csv(occ_i_filt_part, paste0(sp_dir, "01_occ_", sp_name, "_filt_part.csv"))    
        readr::write_csv(psa_part, paste0(sp_dir, "01_occ_", sp_name, "_psa_part.csv"))    
        readr::write_csv(bg_part, paste0(sp_dir, "01_occ_", sp_name, "_bg_part.csv"))    
        
        ### extracting environmental values ----
        occ_i_filt_part_pa <- dplyr::bind_rows(occ_i_filt_part, psa_part)
        
        occ_i_filt_part_pa_data <- flexsdm::sdm_extract(
            data = occ_i_filt_part_pa,
            x = "x",
            y = "y",
            env_layer = var_c_ca_sel,
            filter_na = TRUE)
        readr::write_csv(occ_i_filt_part_pa_data, paste0(sp_dir, "02_var_data_occ_psa_", sp_name, ".csv"))    
        
        bg_part_data <- flexsdm::sdm_extract(
            data = bg_part,
            x = "x",
            y = "y",
            env_layer = var_c_ca_sel,
            filter_na = TRUE)
        readr::write_csv(bg_part_data, paste0(sp_dir, "02_var_data_bg_", sp_name, ".csv"))    
        
        ### models ----
        
        #### esm models ----
        glm_esm <- flexsdm::esm_glm(
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            thr = "max_sens_spec")
        if(any(is.na(glm_esm))){
        } else{
            for(k in 1:length(glm_esm$esm_model)){
                glm_esm_plot <- flexsdm::p_pdp(model = glm_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_glm_", sp_name, "_", k, ".png"), 
                       plot = glm_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        gam_esm <- flexsdm::esm_gam(
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            thr = "max_sens_spec")
        if(any(is.na(gam_esm))){
        } else{
            for(k in 1:length(gam_esm$esm_model)){
                gam_esm_plot <- flexsdm::p_pdp(model = gam_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_gam_", sp_name, "_", k, ".png"), 
                       plot = gam_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        gau_esm <- flexsdm::esm_gau(
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            background = bg_part_data,
            thr = "max_sens_spec")
        if(any(is.na(gau_esm))){
        } else{
            for(k in 1:length(gau_esm$esm_model)){
                gau_esm_plot <- flexsdm::p_pdp(model = gau_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_gau_", sp_name, "_", k, ".png"), 
                       plot = gau_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        gbm_esm <- flexsdm::esm_gbm(            
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            thr = "max_sens_spec")
        if(any(is.na(gbm_esm))){
        } else{
            for(k in 1:length(gbm_esm$esm_model)){
                gbm_esm_plot <- flexsdm::p_pdp(model = gbm_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_gbm_", sp_name, "_", k, ".png"), 
                       plot = gbm_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        max_esm <- flexsdm::esm_max(            
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            background = bg_part_data,
            thr = "max_sens_spec")
        if(any(is.na(max_esm))){
        } else{
            for(k in 1:length(max_esm$esm_model)){
                max_esm_plot <- flexsdm::p_pdp(model = max_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_max_", sp_name, "_", k, ".png"), 
                       plot = max_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        net_esm <- flexsdm::esm_net(           
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            thr = "max_sens_spec")
        if(any(is.na(net_esm))){
        } else{
            for(k in 1:length(net_esm$esm_model)){
                net_esm_plot <- flexsdm::p_pdp(model = net_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_net_", sp_name, "_", k, ".png"), 
                       plot = net_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        svm_esm <- flexsdm::esm_svm(            
            data = occ_i_filt_part_pa_data,
            response = "pr_ab",
            predictors = var_c_ca_sel_names,
            partition = ".part",
            thr = "max_sens_spec")
        if(any(is.na(svm_esm))){
        } else{
            for(k in 1:length(svm_esm$esm_model)){
                svm_esm_plot <- flexsdm::p_pdp(model = svm_esm$esm_model[[k]], training_data = occ_i_filt_part_pa_data, 
                                               colorl = "red", resid = TRUE, colorp = "gray", alpha = 1)
                ggsave(filename = paste0(sp_dir, "03_mod_response_svm_", sp_name, "_", k, ".png"), 
                       plot = svm_esm_plot, width = 20, height = 15, units = "cm", dpi = 300)
            }
        }
        
        ### ensemble ----
        model_list <- Filter(function(x) !any(is.na(x)), 
                             list(glm_esm, gam_esm, gau_esm, gbm_esm, max_esm, 
                                  net_esm, svm_esm))
        
        ### performance ----
        model_performance <- flexsdm::sdm_summarize(model_list)
        readr::write_csv(model_performance, paste0(sp_dir, "03_fit_performance_", sp_name, ".csv"))
        
        ### predict ----
        
        #### current ----
        for(e in model_list){
            
            predict_c_e <- flexsdm::sdm_predict(
                models = e,
                pred = var_c_ca_sel,
                nchunk = 10,
                thr = "max_sens_spec",
                con_thr = TRUE)
            
            terra::writeRaster(predict_c_e[[1]][[1]], paste0(sp_dir, "05_ens_con_c_", names(predict_c_e), "_", sp_name, ".tif"), overwrite = TRUE)
            terra::writeRaster(predict_c_e[[1]][[2]], paste0(sp_dir, "05_ens_con_thr_c_", names(predict_c_e), "_", sp_name, ".tif"), overwrite = TRUE)
            
        }
        
        predict_c_con_mean <- dir(sp_dir, pattern = "05_ens_con_c_esm_", full.names = TRUE) %>% 
            terra::rast() %>% 
            terra::mean()
        names(predict_c_con_mean) <- paste0(sp_name, "_mean_con_c")    
        terra::writeRaster(predict_c_con_mean, paste0(sp_dir, "05_ens_con_c_", sp_name, ".tif"), overwrite = TRUE)
        
        predict_c_con_thr_mean <- dir(sp_dir, pattern = "05_ens_con_thr_c_esm_", full.names = TRUE) %>% 
            terra::rast() %>% 
            terra::mean()
        names(predict_c_con_thr_mean) <- paste0(sp_name, "_mean_con_thr_c")
        terra::writeRaster(predict_c_con_thr_mean, paste0(sp_dir, "05_ens_con_thr_c_", sp_name, ".tif"), overwrite = TRUE)
        
        #### future ----
        for(j in fut){
            
            var_f_ca_sel_j <- terra::subset(var_f_ca_sel, grep(j, names(var_f_ca_sel)))
            names(var_f_ca_sel_j) <- names(var_c_ca_sel)
            
            for(e in model_list){
                
                predict_f_e <- flexsdm::sdm_predict(
                    models = e,
                    pred = var_f_ca_sel_j,
                    nchunk = 10,
                    thr = "max_sens_spec",
                    con_thr = TRUE)
                
                terra::writeRaster(predict_f_e[[1]][[1]], paste0(sp_dir, "05_ens_con_f_", names(predict_f_e), "_", sp_name, "_", j, ".tif"), overwrite = TRUE)
                terra::writeRaster(predict_f_e[[1]][[2]], paste0(sp_dir, "05_ens_con_thr_f_", names(predict_f_e), "_", sp_name, "_", j, ".tif"), overwrite = TRUE)
                
            }
            
            predict_f_con_mean <- dir(sp_dir, pattern = "05_ens_con_f_esm_", full.names = TRUE) %>% 
                grep(j, ., value = TRUE) %>% 
                terra::rast() %>% 
                terra::mean()
            names(predict_f_con_mean) <- paste0(sp_name, "_mean_con_f_", j)    
            terra::writeRaster(predict_f_con_mean, paste0(sp_dir, "05_ens_con_f_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            
            predict_f_con_thr_mean <- dir(sp_dir, pattern = "05_ens_con_thr_f_esm_", full.names = TRUE) %>% 
                grep(j, ., value = TRUE) %>% 
                terra::rast() %>% 
                terra::mean()
            names(predict_f_con_thr_mean) <- paste0(sp_name, "_mean_con_thr_f_", j)
            terra::writeRaster(predict_f_con_thr_mean, paste0(sp_dir, "05_ens_con_thr_f_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            
        }
        
        ### msdm ----
        
        #### current ----
        predict_c_con_mean_m <- flexsdm::msdm_posteriori(
            records = occ_i_filt_part_pa,
            x = "x",
            y = "y",
            pr_ab = "pr_ab",
            method = "pres",
            cont_suit = predict_c_con_mean,
            thr = "max_sens_spec",
            buffer = NULL)
        
        names(predict_c_con_mean_m[[1]]) <- paste0(sp_name, "_mean_con_c_msdm")
        names(predict_c_con_mean_m[[2]]) <- paste0(sp_name, "_mean_con_c_msdm_thr")
        
        terra::writeRaster(predict_c_con_mean_m[[1]], paste0(sp_dir, "05_ens_con_c_msdm_", sp_name, ".tif"), overwrite = TRUE)
        terra::writeRaster(predict_c_con_mean_m[[2]], paste0(sp_dir, "05_ens_con_c_msdm_thr_", sp_name, ".tif"), overwrite = TRUE)
        
        #### future ----
        predict_f_con_mean <- dir(path = sp_dir, pattern = "ssp", full.names = TRUE) %>% 
            grep("esm", ., value = TRUE, invert = TRUE) %>% 
            grep("_thr", ., value = TRUE, invert = TRUE) %>% 
            terra::rast()
        
        for(j in fut){
            
            predict_f_con_mean_j <- terra::subset(predict_f_con_mean, grep(j, names(predict_f_con_mean)))
            
            predict_f_con_mean_j_m <- flexsdm::msdm_posteriori(
                records = occ_i_filt_part_pa,
                x = "x",
                y = "y",
                pr_ab = "pr_ab",
                method = "pres",
                cont_suit = predict_f_con_mean_j,
                thr = "max_sens_spec",
                buffer = NULL)
            
            names(predict_f_con_mean_j_m[[1]]) <- paste0(sp_name, "_mean_con_f_msdm_", j)
            names(predict_f_con_mean_j_m[[2]]) <- paste0(sp_name, "_mean_con_f_msdm_thr_", j)
            
            terra::writeRaster(predict_f_con_mean_j_m[[1]], paste0(sp_dir, "05_ens_con_f_msdm_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            terra::writeRaster(predict_f_con_mean_j_m[[2]], paste0(sp_dir, "05_ens_con_f_msdm_thr_", sp_name, "_", j, ".tif"), overwrite = TRUE)
            
        }
        
    }
    
}

# end ---------------------------------------------------------------------