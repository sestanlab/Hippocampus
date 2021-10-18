## Plot the UMAP of the dataset
source("../scripts/hip.fun.R")


## Read the dataset
hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds"))


set.seed <- 42
meta_data <- hipec@meta.data[sample(colnames(hipec)), ]
plot_data <- meta_data[, c("group", "region_new", "samplename")] %>%
				cbind(., hipec$umap@cell.embeddings[rownames(meta_data), ])


reg_colors <- paste0(c("#b2182b", "#b8e186", "#4d9221", "#80ffff", "#fee090"), "") %>%
				setNames(., c("DG", "CA24", "CA1", "SUB", "EC"))

theme_empty <- function(x) {
	theme_classic() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "none")
}


## Plot the subregion
p1 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = group)) +
		geom_point(size = 0.01) + 
		scale_color_manual(values = c(ExN = "#C51B7D", InN = "#4D9221", NNC = "#999999"))+
		theme_empty()
p2 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = region_new)) +
		geom_point(size = 0.01) + 
		scale_color_manual(values = reg_colors)+
		theme_empty() 
p3 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = samplename)) +
		geom_point(size = 0.01) + 
		scale_color_manual(values = gg_color_hue(6) %>% setNames(., levels(as.factor(plot_data$samplename))))+
		theme_empty()

plot_data <- plot_data %>%
			mutate(region = ifelse(region_new == 'DG', "DG", "bg"))
dg_data <-  plot_data[lapply(c("bg", "DG"), function(x) which(plot_data$region == x)) %>% unlist(), ]		
p4 <- ggplot(dg_data, aes(x = UMAP_1, y = UMAP_2, color = region)) +
		geom_point(size = 0.01) + 
		scale_color_manual(values = c(DG = "#b2182b", bg = "#e5e5e5"))+
		theme_empty()

plot.scale <- 1.35
jpeg(paste0(outputdir, "MF1.Cluster.UMAP.jpeg"), width = plot.scale * 4* 8, height = plot.scale * 8, units = "in", res = 300)
plot_grid(plotlist = list(p1, p2, p3, p4), nrow = 1, ncol = 4) %>% print()
dev.off()



## plot the legend
lp_data <- data.frame(cluster = c(levels(as.factor(plot_data$samplename)), names(reg_colors), c("ExN", "InN", "NNC")), 
					color = c(gg_color_hue(6), reg_colors, c("#C51B7D", "#4D9221", "#999999")),
					stringsAsFactors = FALSE) %>%
				mutate(value = 1)
lp <- ggplot(lp_data, aes(y = cluster,  x = 1, color = cluster))+
		geom_point(size = 3) + 
		##scale_color_identity() +
		scale_color_manual(values = setNames(lp_data$color, lp_data$cluster)) +
		theme_classic() +
		theme(legend.position = "right")
pdf(paste0(outputdir, "MF1.Cluster.UMAP.legend.pdf"), width = 4, height = 6)
print(lp)
dev.off()






