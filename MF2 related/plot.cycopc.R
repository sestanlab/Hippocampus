## Plot the expression of cycOPC
args <- commandArgs(trailingOnly = TRUE) 
source("../scripts/hip.fun.R")
library(ggpubr)


hrpm <- readRDS(file = paste0(inputdir, "HRPM.all.seurat.expHVG.1500.slim.rds"))

FeatureFig(hrpm, features = c("PDGFRA", "SOX6", "OLIG2"), file_name = "cycOPC.exp", split.by = "species", split.order = c("Mouse", "Pig", "Rhesus", "Human"))


sq <- c(7, 9, 2.8, 4.5)
p <- FeatureFig(hrpm, features = "PDGFRA", file_name = "cycOPC.exp.UMAP", output.ggplot = TRUE, pt.size = 0.1)[[1]] +
					geom_vline(xintercept = sq[1:2]) +
					geom_hline(yintercept = sq[3:4])
jpeg(paste0(outputdir, "cycOPC.exp.UMAP.jpeg"), width = 8, height = 8, units = "in", res = 300)
print(p)
dev.off()

submeta <- colnames(hrpm)[hrpm$umap@cell.embeddings[, 1] >= sq[1] & 
							hrpm$umap@cell.embeddings[, 1] <= sq[2] &
							hrpm$umap@cell.embeddings[, 2] >= sq[3] & 
							hrpm$umap@cell.embeddings[, 2] <= sq[4]] %>%
				hrpm@meta.data[., ]
opc_meta <- submeta[submeta$fig2cluster == "OPC", ]


cls_use <- c("nIPC", "NB", "OPC", "Oligo", "GC")
subseu <- hrpm[, hrpm@meta.data$fig2cluster %in% cls_use]
subseu@meta.data[rownames(opc_meta), "fig2cluster"] <- "cycOPC"
seulist <- SplitObject(subseu, split.by = "species")
seulist$Pig <- subset(seulist$Pig, fig2cluster != "nIPC")

all_genes <- c("PDGFRA", "OLIG2", "SLC1A3", "MKI67", "TOP2A", "CENPF", "PCNA", "PROX1", "DCX")


sp_cols <- setNames(c("#031a7f", "#ffa600", "#89DA59", "#FF420E"), c("Mouse", "Pig", "Rhesus", "Human"))
plot.margin <- unit(c(-0, 0, -0, 0), "inch")
dotlist <- lapply(names(sp_cols), function(sp) {
	new_order <- c("OPC", "cycOPC", "nIPC", "NB", "GC")
	cls_order <- new_order[new_order %in% unique(seulist[[sp]]@meta.data$fig2cluster)]
	Idents(seulist[[sp]]) <- "fig2cluster"
	p <- DotPlot(seulist[[sp]], features = rev(all_genes), cols = c("lightgrey", sp_cols[sp]), dot.scale = 4, dot.min = 0.025, scale.by = "size")+
	        RotatedAxis() +
			scale_y_discrete(limits = rev(cls_order)) +
			scale_size(range = c(0, 4), limits = c(2.5, 100)) +
	        theme(axis.text.x=element_text(size = 7), axis.text.y=element_text(size = 6),legend.position = "right", axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid.major = element_line(size = 0.25, color = "grey"), panel.grid.minor = element_line(size = 0.25, color = "grey"), plot.margin = plot.margin)
	if (sp %in% c("Mouse", "Pig", "Rhesus")){
		p <- p + 
			theme(axis.text.x=element_blank(), legend.position = "right") +
			scale_color_gradient(low = "lightgrey", high = sp_cols[sp], guide = "none") +
			guides(color = "none")
	}
	return(p)
	})

layout1 <- "
A
B
C
D
"
pdf(paste0(outputdir, "SF2.cycOPC.exp.pdf"), useDingbats = FALSE, width = 5, height = 3.5)
patchwork::wrap_plots(dotlist, guides = "collect", design = layout1)
dev.off()



###------------------------------------------------------------------------------------
## FeaturePlot showing the PDGFRA expression
FeatureFig(hrpm, features = "PDGFRA", file_name = "SF2.cycOPC.exp.UMAP.PDGFRA", output.ggplot = FALSE, pt.size = 0.01, cols = c("lightgrey", "darkred"), plot.scale = 0.8, ngradient = 40)
FeatureFig(hrpm, features = "PDGFRA", file_name = "SF2.cycOPC.exp.UMAP.PDGFRA.bysp", output.ggplot = FALSE, pt.size = 0.01, cols = c("lightgrey", "darkred"), plot.scale = 0.6, ngradient = 40, split.by = "species", split.order = c("Mouse", "Pig", "Rhesus", "Human"))















