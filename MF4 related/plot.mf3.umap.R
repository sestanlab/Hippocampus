## Plot the MainFigure and supplementary figure
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")

color_vec <- readRDS(file = paste0(dataDir, "cluster.color.rds"))

seulist <- list(ExN = readRDS(file = paste0(inputdir, "HIPEC.", "ExN", ".final.rename.seu.rds")), 
				InN = readRDS(file = paste0(inputdir, "HIPEC.", "InN", ".final.rename.seu.rds")), 
				NNC = readRDS(file = paste0(inputdir, "HIPEC.", "NNC", ".final.rename.seu.rds")))

## remove CR from the InN groups
seulist$InN <- seulist$InN[, seulist$InN@meta.data$fig1cluster != "CR RELN NDNF"]
group.by <- "fig1cluster"

##------------------------------------------------------------------
## plot the MF umap plots
theme_empty <- function(x) {
	theme_classic() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "none")
}
plist <- lapply(names(seulist), function(ctp) {
	plot_data <- seulist[[ctp]]@meta.data[, c("region_new", group.by)] %>%
				cbind(., seulist[[ctp]]$umap@cell.embeddings) %>%
				mutate(region = ifelse(region_new == "EC", "EC", "HIP"))

	pp <- ggplot(plot_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = group.by)) +
		geom_point(size = 0.01) + 
		scale_color_manual(values = color_vec[levels(as.factor(plot_data[, group.by]))])+
		theme_empty() 
	qq <- ggplot(plot_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = "region")) +
		geom_point(size = 0.01) + 
		scale_color_manual(values = c(HIP = "#FF420E", EC = "#89DA59"))+
		theme_empty() 
	pq <- plot_grid(pp, qq, nrow = 1, ncol = 2, align = "h")
	}) %>% 
		setNames(., names(seulist))

plot.scale <- 0.9
for (ctp in  names(plist)){
	jpeg(paste0(outputdir, "MF3.UMAP.", ctp, ".jpeg"), width = 20 * plot.scale, height = 10 * plot.scale, res = 400, units = "in")
	print(plist[[ctp]])
	dev.off()
}



 
