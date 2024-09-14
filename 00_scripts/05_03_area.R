library(terra)

binPath <- "E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/04_presence_absence_rasters/"
bin.files <- list.files(path = binPath, pattern='.tif', all.files=TRUE)

# loading rasters
fut_2011_2040_SSP370 <- rast(paste0(binPath, bin.files[[4]]))
fut_2011_2040_SSP585 <- rast(paste0(binPath, bin.files[[5]]))
fut_2041_2070_SSP370 <- rast(paste0(binPath, bin.files[[11]]))
fut_2041_2070_SSP585 <- rast(paste0(binPath, bin.files[[12]]))
Baseline <- rast(paste0(binPath, bin.files[[13]]))

# function to calculate area
calc_area_km2 <- function(raster) {
    # Calcula a área de cada célula em km²
    cell_area_km2 <- prod(res(raster)) * 111.32^2  # Convertendo graus para km²
    # Conta o número de células com presença (valor 1)
    num_cells_pres <- global(raster, fun = "sum", na.rm = TRUE)[[1]]
    # Calcula a área total de presença
    area_km2 <- num_cells_pres * cell_area_km2
    return(area_km2)
}

# raster list and species name
raster_list <- list(Baseline, fut_2011_2040_SSP370, fut_2011_2040_SSP585, fut_2041_2070_SSP370, fut_2041_2070_SSP585)
scenario_names <- c("Baseline", "2050 SSP370", "2050 SSP585", "2070 SSP370", "2070 SSP585")

# Calculating area
species_names <- names(Baseline)  # species names
area_data <- data.frame(Species = species_names)

for (i in seq_along(scenario_names)) {
    areas <- sapply(1:nlyr(raster_list[[i]]), function(j) calc_area_km2(raster_list[[i]][[j]]))
    
    # area as colummn in each scenario
    area_data[[scenario_names[i]]] <- areas
}

print(area_data)

# proportion of loss in area
area_data_prop <- area_data
area_data_prop$`2050 SSP370` <- area_data_prop$`2050 SSP370`/area_data_prop$Baseline
area_data_prop$`2050 SSP585` <- area_data_prop$`2050 SSP585`/area_data_prop$Baseline
area_data_prop$`2070 SSP370` <- area_data_prop$`2070 SSP370`/area_data_prop$Baseline
area_data_prop$`2070 SSP585` <- area_data_prop$`2070 SSP585`/area_data_prop$Baseline

# Reshape the data for plotting
melted_data <- melt(area_data_prop[,-2], id.vars = "Species", variable.name = "Scenario", value.name = "Area_km2")

# Load required libraries
library(ggplot2)
library(viridis)  # For colorblind-friendly palettes

# Plotting the violin plot for each scenario with colorblind-friendly palette and no legend
ggplot(melted_data, aes(x = Scenario, y = Area_km2, fill = Scenario)) +
    geom_violin(trim = FALSE) +  # Create the violin plot
    geom_jitter(width = 0.2, size = 0.8) +  # Add jitter to show individual data points
    labs(title = "Change in Area per Scenario", x = "Scenario", y = "Proportion of change") +
    scale_fill_viridis(discrete = TRUE, alpha = 0.9, begin = 0.2) +  # Use viridis color palette for colorblind-friendliness
    theme_classic() +  # Minimal theme for clean visuals
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
          legend.position = "none")  # Remove the legend

setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Figures")

png("prop_change.png", width = 8, height = 6, 
     res = 300, unit = "in")
tiff("prop_change.tiff", width = 8, height = 6, 
     res = 300, unit = "in")
# par(mfrow = c(1,2))

ggplot(melted_data, aes(x = Scenario, y = Area_km2, fill = Scenario)) +
    geom_violin(trim = FALSE) +  # Create the violin plot
    geom_jitter(width = 0.2, size = 0.8) +  # Add jitter to show individual data points
    labs(title = NULL, x = "Year and Scenario", y = "Proportion of change") +
    scale_fill_viridis(discrete = TRUE, alpha = 0.9, begin = 0.2) +  # Use viridis color palette for colorblind-friendliness
    theme_classic() +  # Minimal theme for clean visuals
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
          legend.position = "none")  # Remove the legend

dev.off()

# saving the csv
write.csv(area_data, "E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/Tables/area_scenarios.csv", row.names = FALSE)
