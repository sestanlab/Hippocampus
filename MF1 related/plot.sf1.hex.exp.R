## Visualize gene expression on HEX plots
source("../scripts/hip.fun.R")
library(schex)
library(ggpubr)

hipec_hex <- readRDS(file = paste0(dataDir, "HIPEC_hex.cbn.rds"))


genes <- c("ZBTB20", "SLC17A7", "GAD1", "AQP4", "PDGFRA", "MOBP", "PTPRC", "CLDN5", "CEMIP")
plist <- lapply(genes, function(x) {
	p <- plot_hex_cbn(sce = hipec_hex, feature = x, action = "mean", cols = c("grey90",  "red"))
	return(p)
	})


plot.scale <- 1
nrow <- 3
ncol <- 3

ggarrange(plotlist = plist, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = "right") %>%
	ggexport(filename = paste0(outputdir, "SF1.marker.hex.exp.withtitle.pdf"), width = plot.scale * ncol * 5 + 2, height = plot.scale * nrow * 5)


