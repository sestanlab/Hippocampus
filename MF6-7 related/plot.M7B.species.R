## Plot the expression of DCX across all cluster in mouse, rhesus, Human DG and human EC cells.
source("../scripts/hip.fun.R") 
library(ggpubr)
library(schex)

rdg <- readRDS(paste0("../MF2_final/load_files/", "Rhesus_DG_seu.rds"))
hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds"))
pdg <- readRDS(paste0("../MF2_final/load_files/", "Pig_DG_seu.rds"))
mdg <- readRDS(paste0("../MF2_final/load_files/", "Mouse_DG_seu.rds"))


seulist <- list(Human = hipec, Rhesus = rdg, Pig = pdg, Mouse = mdg) %>%
			lapply(., function(x) make_hexbin(x, nbins = 100, dimension_reduction = "UMAP"))



p <- plot_hexbin_feature_byregion(sce = seulist, assay = "RNA", slot = "data", features = "METTL7B", action = "mean", nrow = 1, ncol = 4, colors = c("#f0f0f0", "#d7301f", "#b30000", "#7f0000"), file_name = "MF5.METTL7B.hex", output_dir = outputdir, pdf_size = c(16,4), region_order = names(seulist), legend = "bottom", return_rawp = TRUE, exp_ceiling = 1.5)[["METTL7B"]]
pdf(paste0(outputdir, "MF6.M7B.Exp.species.hex.pdf"), width = 20, height = 6)
plot_grid(plotlist = p, nrow = 1, ncol = 4) %>% print()
dev.off() 

jpeg(paste0(outputdir, "MF6.M7B.Exp.species.hex.jpeg"), width = 20, height = 6, res = 300, unit = "in")
plot_grid(plotlist = p, nrow = 1, ncol = 4) %>% print()
dev.off() 



for (ii in names(seulist)){
	DimFig(seulist[[ii]], group.by = "cluster", file_name = paste0("UMAP_cluster_", ii))
}






