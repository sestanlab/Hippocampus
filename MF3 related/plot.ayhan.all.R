source("../scripts/hip.fun.R")



## Prepare
type <- "All"
seu <- switch(type, All = readRDS(file = paste0(inputdir, "Konopka.Inte.Seurat.2000.v2.slim.rds")), 
				Sub = readRDS(file = paste0(inputdir, "Konopka.Inte.Seurat.slim.2500.slim.rds")))
p.scale <- switch(type, All = 1,
				Sub = 0.8)
file_name <- switch(type, All = "SF.all",
				Sub = "SF.Sub")



## ADD the doublet scores into the meta.data
meta <- readRDS(file = paste0(inputdir, "Konopka_all_db_res.rds"))
seu@meta.data <- cbind(seu@meta.data, meta[colnames(seu), ,drop = FALSE])
	
## ADD the Oligo AUC into the meta.data
auc_res <- readRDS(file = paste0(inputdir, "Konopka_Oligo_auc_org.rds"))
seu@meta.data <- cbind(seu@meta.data, auc_res[colnames(seu), ,drop = FALSE])


## Do the Hex
library(schex)
seu <- make_hexbin(seu, nbins = 150, dimension_reduction = "UMAP")

source("./vis.fun.R")
genes <- c("LPAR1", "MOBP", "Oligo", "doublet_scores")
plist <- lapply(genes, function(x) {
	p <- plot_hex_cbn_v2(sce = seu, feature = x, action = "mean", cols = c("#f5f5dc", "#31a354","#253494"))
	return(p)
	})


plot.scale <- 1
nrow <- 1
ncol <- 4

library(ggpubr)
ggarrange(plotlist = plist, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = "right") %>%
	ggexport(filename = paste0(outputdir, "SF2.Exp_scores.pdf"), width = plot.scale * ncol * 5 + 2, height = plot.scale * nrow * 5)












plot_hex_cbn_v2(sce, feature, action = "mean", cols = viridis(3), assay = "RNA")



