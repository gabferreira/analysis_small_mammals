#' ---
#' title: figures manuscript mammals
#' author: gabriela alves ferreira
#' date: 2024-09-06
#' ---

# setwd("/media/gabriela/Gabi_HD/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises/enm_mammals")
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_diversity_metrics")

# packages
library(terra)
library(dplyr)
library(tmap) 

# load richness rasters
sr_pres <- rast("richness_Baseline.tif")
sr_2040_SSP370 <- rast("richness_2011_2040_SSP370.tif")
sr_2040_SSP585 <- rast("richness_2011_2040_SSP585.tif")
sr_2070_SSP370 <- rast("richness_2041_2070_SSP370.tif")
sr_2070_SSP585 <- rast("richness_2041_2070_SSP585.tif")

# load phylogenetic diversity rasters
pd_pres <- rast("phylo_div_Baseline.tif")
pd_2040_SSP370 <- rast("phylo_div_2011_2040_SSP370.tif")
pd_2040_SSP585 <- rast("phylo_div_2011_2040_SSP585.tif")
pd_2070_SSP370 <- rast("phylo_div_2041_2070_SSP370.tif")
pd_2070_SSP585 <- rast("phylo_div_2041_2070_SSP585.tif")

# load SES phylogenetic diversity rasters
pd_pres_ses <- rast("pd_pres_SES.tif")
pd_2040_SSP370_ses <- rast("pd_2011_2040_370_SES.tif")
pd_2040_SSP585_ses <- rast("pd_2011_2040_585_SES.tif")
pd_2070_SSP370_ses <- rast("pd_2041_2070_370_SES.tif")
pd_2070_SSP585_ses <- rast("pd_2041_2070_585_SES.tif")

# plot(pd_pres/sr_pres)
# plot(pd_pres_ses)

# load phylogenetic delta croped by PAs and remnants
delta_pd_rem <- rast("delta_pd_2070_585_REMNANTS.tif")
delta_pd_pas <- rast("delta_pd_2070_585_PAs.tif")

# shapefile south america
amer <- vect("E:/shapes_rasters/shapes/americadoSul/South_America_Countries/South_America.shp")
# ext(amer)
bras <- vect("E:/shapes_rasters/shapes/UFs_Brasil/estados_brasileiros_2010.shp")
af <- terra::vect("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/01_data/02_variables/00_limit/merge_limites_MA_buffer20km_WGS84.shp")


##############################
# richness maps

# breaks for richness maps
breaks_rq <- seq(0, 57, by = 0.01)
library(viridis)
colors_rq <- viridis::magma(7, direction = -1) # accessible pallete

colors_rq <- colorRampPalette(c("#FCFDBFFF", "#FEAF77FF" ,"#F1605DFF", "#B63679FF",
                                "#721F81FF", "#2D1160FF","#000004FF"))(length(breaks_rq)-1)

# ploting and saving
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Figures")
# jpeg("richness_Present.jpg", width = 560, height = 640, quality = 100)
tiff("richness_Present.tiff", width = 10, height = 10, unit = "in", res = 300)

plot(sr_pres, col = c(colors_rq), range = c(0, 57),
     axes = T,
     main = "(a) TD - Present",
     plg = list( # parameters for drawing legend
       # title = "b)",
       # title.cex = 2, # Legend title size
       cex = 1.8), # Legend text size
     pax=list( # parameters for drawing axes
       cex.axis = 1.5), # Axis text size
     cex.main = 1.75) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

# ploting
# jpeg("richness_future.jpeg", width = 820, height = 620, quality = 100)
tiff("richness_future.tiff", width = 10, height = 10, unit = "in", res = 300)
par(mfrow = c(2,2))
# jpeg("richness_SSP585.jpg", width = 1020, height = 620, quality = 100)
# par(mfrow = c(1,2))

plot(sr_2040_SSP370, col = c(colors_rq), range = c(0, 57),
     axes = T,
     main = "(b) TD - 2050 SSP370",
     plg = list( # parameters for drawing legend
       # title = "b)",
       # title.cex = 2, # Legend title size
       cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
       cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(sr_2040_SSP585, col = c(colors_rq), range = c(0, 57),
     axes = T,
     main = "(c) TD - 2050 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
                cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

# dev.off()

# jpeg("richness_SSP370.jpg", width = 1020, height = 620, quality = 100)
# par(mfrow = c(1,2))


plot(sr_2070_SSP370, col = c(colors_rq), range = c(0, 57),
     axes = T,
     main = "(d) TD - 2070 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
                cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(sr_2070_SSP585, col = c(colors_rq), range = c(0, 57),
     axes = T,
     main = "(e) TD - 2070 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

##############################
# delta
library(phyloraster)
delta2040_370 <- delta.grid(sr_pres, sr_2040_SSP370)
delta2040_585 <- delta.grid(sr_pres, sr_2040_SSP585)

delta2070_370 <- delta.grid(sr_pres, sr_2070_SSP370)
delta2070_585 <- delta.grid(sr_pres, sr_2070_SSP585)

## breaks para os mapas de delta riqueza
minmax(c(delta2040_370, delta2040_585,
         delta2070_370, delta2070_585))

breaks_dt <- seq(-21, 21, by = 0.01)
colors_dt <- colorRampPalette(c("#c44601", "#f57600","#faaf90","#C0C6CB",
                                "#8babf1", "#0073e6", "#054fb9"))(length(breaks_dt))
# colors_dt <- colorRampPalette(c("#A20A19","#BB061C","#ED0024","#f07470",
#                                 "#C0C6CB","#00A6D7","#0058B3","darkblue","black"))(length(breaks_dt))

# plotando
# tiff("delta_richness2.tiff", width = 920, height = 920, res = 600)
tiff("delta_richness.tiff", width = 10, height = 10, res = 300, units = "in")
# jpeg("delta_richness.jpeg", width = 820, height = 620, quality = 100)
# x11()
par(mfrow = c(2,2))
plot(delta2040_370, col = c(colors_dt), range = c(-21, 21),
     axes = T,
     main = "(b) Delta TD - 2050 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(delta2040_585, col = c(colors_dt), range = c(-21, 21),
     axes = T,
     main = "(c) Delta TD - 2050 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

# breaks_dt <- seq(-17, 7, by = 0.1)
# colors_dt <- colorRampPalette(c("#450B0C","#660000","#B30D02","#FF0800","#C0C6CB","#00A6D7","#0058B3","darkblue"))(length(breaks_dt))


plot(delta2070_370, col = c(colors_dt), range = c(-21, 21),
     axes = T,
     main = "(d) Delta TD - 2070 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(delta2070_585, col = c(colors_dt), range = c(-21, 21),
     axes = T,
     main = "(e) Delta TD - 2070 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

##############################
# phylogenetic diversity maps

# breaks for phylodiv maps
breaks_pd <- seq(0, 12.27122 , by = 0.01)
library(viridis)
colors_pd <- viridis::viridis(7, direction = -1) # accessible pallete
# colors_rq <- colorRampPalette(c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF",
#                                 "#35B779FF", "#8FD744FF", "#FDE725FF"))(length(breaks_rq)-1)
colors_pd <- colorRampPalette(c("#FDE725FF", "#8FD744FF", "#35B779FF", "#21908CFF", "#31688EFF",
                                "#443A83FF", "#440154FF"))(length(breaks_pd)-1)

# ploting and saving
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Figures")
# jpeg("phylodiv_Present.jpg", width = 560, height = 640, quality = 100)
tiff("phylodiv_present.tiff", width = 10, height = 10, res = 300, units = "in")

plot(pd_pres, col = c(colors_pd), range = c(0, 12.27122),
     axes = T,
     main = "(a) PD - Present",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 1.8), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.5), # Axis text size
     cex.main = 1.75) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

# ploting
# jpeg("phylodiv_future.jpeg", width = 820, height = 620, quality = 100)
tiff("phylodiv_future.tiff", width = 10, height = 10, res = 300, units = "in")
par(mfrow = c(2,2))
# jpeg("richness_SSP585.jpg", width = 1020, height = 620, quality = 100)
# par(mfrow = c(1,2))
plot(pd_2040_SSP370, col = c(colors_pd), range = c(0, 12.27122),
     axes = T,
     main = "(b) PD - 2050 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(pd_2040_SSP585, col = c(colors_pd), range = c(0, 12.27122),
     axes = T,
     main = "(c) PD - 2050 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(pd_2070_SSP370, col = c(colors_pd), range = c(0, 12.27122),
     axes = T,
     main = "(d) PD - 2070 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(pd_2070_SSP585, col = c(colors_pd), range = c(0, 12.27122),
     axes = T,
     main = "(e) PD - 2070 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

##############################
# delta
library(phyloraster)
delta_pd2040_370 <- delta.grid(pd_pres, pd_2040_SSP370)
delta_pd2040_585 <- delta.grid(pd_pres, pd_2040_SSP585)

delta_pd2070_370 <- delta.grid(pd_pres, pd_2070_SSP370)
delta_pd2070_585 <- delta.grid(pd_pres, pd_2070_SSP585)

## breaks para os mapas de delta riqueza
minmax(c(delta_pd2040_370, delta_pd2040_585,
         delta_pd2070_370, delta_pd2070_585))

breaks_dt <- seq(-4.779682, 3.930551, by = 0.01)
# colors_dt <- colorRampPalette(c("#993636","#AC5341","#BF6F4D","#D18C58","#F7C56F","#C0C6CB","#00A6D7","#0058B3","darkblue"))(length(breaks_dt))
# colors_dt <- colorRampPalette(c("#A20A19","#BB061C","#D40320","#ED0024","#f07470","#C0C6CB","#00A6D7","#0058B3","darkblue"))(length(breaks_dt))
colors_dt <- colorRampPalette(c("#c44601", "#f57600","#faaf90","#C0C6CB",
                                "#8babf1", "#0073e6", "#054fb9"))(length(breaks_dt))

# plotando
# tiff("delta_richness2.tiff", width = 920, height = 920, res = 600)
tiff("delta_phylodiv.tiff", width = 10, height = 10, res = 300, units = "in")
# jpeg("delta_phylodiv.jpeg", width = 820, height = 620, quality = 100)

par(mfrow = c(2,2))
plot(delta_pd2040_370, col = c(colors_dt), range = c(-4.421559, 3.930551),
     axes = T,
     main = "(b) Delta PD - 2050 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2.2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(delta_pd2040_585, col = c(colors_dt), range = c(-4.421559, 3.930551),
     axes = T,
     main = "(c) Delta PD - 2050 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2.2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

# breaks_dt <- seq(-17, 7, by = 0.1)
# colors_dt <- colorRampPalette(c("#450B0C","#660000","#B30D02","#FF0800","#C0C6CB","#00A6D7","#0058B3","darkblue"))(length(breaks_dt))

plot(delta_pd2070_370, col = c(colors_dt), range = c(-4.421559, 3.930551),
     axes = T,
     main = "(d) Delta PD - 2070 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2.2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(delta_pd2070_585, col = c(colors_dt), range = c(-4.421559, 3.930551),
     axes = T,
     main = "(e) Delta PD - 2070 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2.2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

##############################
# SES
# breaks for phylodiv SES maps
minmax(c(pd_pres_ses, pd_2040_SSP370_ses, pd_2040_SSP585_ses,
         pd_2070_SSP370_ses, pd_2070_SSP585_ses))
breaks_pd_ses <- seq(-0.4303128, 0.4092114 , by = 0.001)
library(viridis)
colors_pd_ses <- viridis::rocket(7, direction = -1) # accessible pallete
# colors_rq <- colorRampPalette(c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF",
#                                 "#35B779FF", "#8FD744FF", "#FDE725FF"))(length(breaks_rq)-1)
colors_pd_ses <- colorRampPalette(c("#FAEBDDFF", "#F6AA82FF", "#F06043FF",
                                "#CB1B4FFF", "#841E5AFF", "#3F1B44FF",
                                 "#03051AFF"))(length(breaks_pd_ses)-1)

# ploting and saving
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Figures")
# jpeg("phylodiv_Present.jpg", width = 560, height = 640, quality = 100)
tiff("phylodiv_present_SES.tiff", width = 10, height = 10, res = 300, units = "in")

plot(pd_pres_ses, col = c(colors_pd_ses), range = c(-0.4303128, 0.4092114),
     axes = T,
     main = "(a) PD SES - Present",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 1.8), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.5), # Axis text size
     cex.main = 1.75) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

# ploting
# jpeg("phylodiv_future.jpeg", width = 820, height = 620, quality = 100)
tiff("phylodiv_future_SES.tiff", width = 10, height = 10, res = 300, units = "in")
par(mfrow = c(2,2))
# jpeg("richness_SSP585.jpg", width = 1020, height = 620, quality = 100)
# par(mfrow = c(1,2))
plot(pd_2040_SSP370_ses, col = c(colors_pd_ses), range = c(-0.4303128, 0.4092114),
     axes = T,
     main = "(b) PD SES - 2050 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(pd_2040_SSP585_ses, col = c(colors_pd_ses), range = c(-0.4303128, 0.4092114),
     axes = T,
     main = "(c) PD SES - 2050 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(pd_2070_SSP370_ses, col = c(colors_pd_ses), range = c(-0.4303128, 0.4092114),
     axes = T,
     main = "(d) PD SES - 2070 SSP370",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(pd_2070_SSP585_ses, col = c(colors_pd_ses), range = c(-0.4303128, 0.4092114),
     axes = T,
     main = "(e) PD SES - 2070 SSP585",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 2), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1.8), # Axis text size
     cex.main = 2) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

##############################
# delta PAs and REMNANTS
delta_pd_pas
delta_pd_rem

# breaks for delta PA and delta remnant maps
breaks_delta_pa_rm <- seq(-4.336499, 3.151682 , by = 0.01)
colors_delta_pa_rm <- colorRampPalette(c("#c44601", "#f57600","#faaf90","#C0C6CB",
                                          "#0073e6", "#054fb9"))(length(
                                           breaks_delta_pa_rm)-1)
# colors_delta_pa_rm <- colorRampPalette(c("brown","#836953",
#                                          "#0b9abd"))(length(breaks_delta_pa_rm)-1)

# ploting and saving
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Figures")

# jpeg("delta_PD_PAs_remnants_2070_585.jpg", width = 1020, height = 520, 
#      quality = 100)
tiff("delta_PD_PAs_remnants_2070_585.tiff", width = 10, height = 6, 
     res = 300, unit = "in")
par(mfrow = c(1,2))

plot(delta_pd_pas, col = c(colors_delta_pa_rm), range = c(-4.336499, 3.151682),
     axes = T,
     main = "(a) Delta PD - Protected areas",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 1.1), # Legend text size
     pax=list( # parameters for drawing axes
                  cex.axis = 1), # Axis text size
     cex.main = 1.1) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)

plot(delta_pd_rem, col = c(colors_delta_pa_rm), range = c(-4.336499, 3.151682 ),
     axes = T,
     main = "(b) Delta PD - Atlantic Forest Remnants",
     plg = list( # parameters for drawing legend
         # title = "b)",
         # title.cex = 2, # Legend title size
         cex = 1.1), # Legend text size
     pax=list( # parameters for drawing axes
         cex.axis = 1), # Axis text size
     cex.main = 1.1) # Title text size
terra::lines(bras, lwd = 0.3, alpha = 0.9)
terra::lines(amer, lwd = 0.3, alpha = 0.9)
dev.off()

##############################
# supplementary figure 1 with margin limits of atlantic forest and occurrence records
# ploting and saving

## occurrences ----
occ <- readr::read_csv("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/01_data/01_occurrences/02_clean/00_occ_cleaned.csv") %>% 
    dplyr::rename(x = longitude, y = latitude) %>% 
    dplyr::select(species, x, y) %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(n = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(n >= 10)
occ

setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Figures")
tiff("supplementary_AF_occs.tiff", width = 10, height = 10, unit = "in", res = 300)

ext(amer) <- c(-82, -26.2413902282715, -58.4986114501953, 12.5902767181396)

plot(amer, col = "lightgray", 
     xlim = c(-88, -30.24),  
     ylim = c(-58.5, 12.59), 
     axes = TRUE, 
     # main = "(a) South America with Atlantic Forest Limits",
     plg = list(
         cex = 1.8), 
     pax = list(
         cex.axis = 1.5), 
     cex.main = 1.75) 

# South america
terra::lines(bras, lwd = 0.3, alpha = 0.9)  
terra::lines(amer, lwd = 0.3, alpha = 0.9)  

# AF
plot(af, add = TRUE, col = "darkgreen", border = "black") 

# occs
occ_vect <- terra::vect(occ, geom = c("x", "y"), crs = "EPSG:4326")
plot(occ_vect, add = TRUE, col = "blue", pch = 19, cex = 0.1) 

dev.off()