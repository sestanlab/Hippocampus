## Organize the dataset 
source("../scripts/hip.fun.R")

## Further remove low-quality cells identified during the integration
##filterFile <- paste0(inputdir, "Inte.InN.4reg.harmony.slim.filtered.rds")
##if (!file.exists(filterFile)){
##	seu <- readRDS(paste0(inputdir, "Inte.InN.4reg.harmony.slim.rds")) %>%
##				subset(seurat_clusters != "15") %>%
##				RunUMAP(., dims = 1:30, reduction = "harmony")
##	##DimFig(seu, group.by = "subtypes", file_name = "Inte.InN.4reg.harmony.filter", split.by = "merge_region", split.order = c("HIP", "EC", "MTG", "PFC"), plot.scale = 0.7)
##  saveRDS(seu, file = filterFile)
##}
for (mm in c("seurat", "MNN", "harmony")){ ## 2 methods
	seu <- readRDS(paste0(inputdir, "Inte.InN.4reg.", mm, ".slim.rds")) %>%
					subset(subtypes != "CR RELN NDNF")


	colors <- c("#FF420E","#89DA59","#4CB5F5","#FFBB00", "grey90") %>% 
					setNames(., c("HIP", "EC", "PFC", "MTG", "bg"))

	plot_data <- data.frame(seurat_clusters = as.character(seu@meta.data$seurat_clusters),
						region = seu@meta.data$merge_region,
						xaxis = seu$umap@cell.embeddings[, 1], 
						yaxis = seu$umap@cell.embeddings[, 2], 
						stringsAsFactors = FALSE)


	plist <- lapply(c("HIP", "EC", "MTG", "PFC"), function(reg) {
		plot_data <- plot_data %>%
						mutate(region_new = ifelse(region == reg, region, "bg"))
		new_data <- plot_data[lapply(c("bg", reg), function(x) which(plot_data$region_new == x)) %>% unlist(), ]

		p <- ggplot(new_data, aes(x = xaxis, y = yaxis, color = region_new)) +
	      	geom_point(size = 0.02) +
	      	scale_color_manual(values = colors[c("bg", reg)]) +
	        theme_classic() + 
	        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")
	    return(p)
		})


	outPrefix <- switch(mm, harmony = "MF4.harmony", seurat = "SF4.seurat", MNN = "SF4.MNN")
	jpeg(paste0(outputdir, outPrefix, ".inte.InN.4reg.region.UMAP.jpeg"), width = 16, height = 4, units = "in", res = 300)
	plot_grid(plotlist = plist, nrow = 1, ncol = 4) %>% print()
	dev.off()
}







