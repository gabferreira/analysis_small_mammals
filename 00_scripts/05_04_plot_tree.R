## Phylogenetic tree 
## Gabriela Alves-Ferreira (gabriela-alves77@hotmail.com)
## 09.09.2024
## R version 4.3.2

# packages
library(terra)
library(phyloraster)
library(ape)
library(ggtree)

###########################
# calculate PD for each scenario and year
# load phylogenetic tree
tree <- read.tree("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises/enm_mammals/01_data/03_phylo_tree/Tree_bootstrap_500.txt")
setdiff(names(bi_pres), tree$tip.label)

# load binary rasters
setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/04_presence_absence_rasters/")
bi_pres <- rast("rast_pres.tif")

# removing species that are not in the tree
setdiff(names(bi_pres), tree$tip.label)
names(bi_pres) <- gsub("(^[a-z])", "\\U\\1", names(bi_pres), perl = TRUE)
data_pres <- phylo.pres(bi_pres, tree)

data_pres$tree$tip.label <- gsub("_", " ", data_pres$tree$tip.label)

# library(ggplot2)
# 
# tr <- ggtree(data_pres$tree, layout = "fan") + 
#     geom_tiplab(size = 3, fontface = "italic") + # Rotular as espécies alinhadas
#     theme_tree2() +  # Usar um tema mais elaborado
#     theme(plot.margin = margin(1, 1, 1, 1, "cm")) # Ajustar as margens

setwd("E:/Manuscritos_em_producao/Manuscrito_Ricardo_Div_Filogenetica/novas_analises_round_2/01_data_small_mammals/05_01_final_maps")
tiff("phylo_tree_mammals.tiff", width = 10, height = 8, 
     res = 300, unit = "in")
# tr
plot.phylo(data_pres$tree, type = "phylogram", cex = 0.5, font = 4, 
           tip.color = "black", edge.color = "grey20", edge.width = 1,
           label.offset = 0.01)  # Ajuste do espaço entre os rótulos e a árvore

dev.off()
