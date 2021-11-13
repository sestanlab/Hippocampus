## Plot the DCX expression in pig in 0h, 1h ,7h
source("../scripts/hip.fun.R")


## Load the pig integration dataset
pig <- readRDS(paste0("~/project/OrganEx/Hippocampus/load_files/OE.HIP.all.seurat.filter.final.highReso.rds"))
pig0 <- readRDS(paste0(inputdir, "zf.HRPM.all.seurat.20210407.rds")) %>%
			subset(species == "Pig")

rm_cells <- setdiff(colnames(pig)[pig@meta.data$condition == "h0"], colnames(pig0))
pig <- pig[, setdiff(colnames(pig), rm_cells)]
pig <- subset(pig, condition == c("h0", "h1", "h7"))

share_cells <- intersect(colnames(pig), colnames(pig0))
pig@meta.data[share_cells, "newcls"] <- pig0@meta.data[share_cells, "fig2cluster"]

DimFig(pig, group.by = c("newcls", "hres"), file_name = "Pig_all", plot.scale = 0.8)


## Get the expression plots using schex
library(schex)
pig_list <- SplitObject(pig, split.by = "condition") %>%
			lapply(., function(x) make_hexbin(x, nbins = 100, dimension_reduction = "UMAP"))


gene <- "METTL7B"
p <- plot_hexbin_feature_byregion(sce = pig_list, assay = "RNA", slot = "data", features = gene, action = "mean", nrow = 1, ncol = 3, colors = c("#f0f0f0", "#d9d9d9", "#fdbb84", "#b30000","#7f0000"), file_name = paste0("SF3.Pig.PMI.", gene, ".hex"), output_dir = outputdir, pdf_size = c(12,4), region_order = c("h0", "h1", "h7"), legend = "bottom", return_rawp = TRUE, exp_ceiling = 2.5)[[gene]]
pdf(paste0(outputdir, "SF3.Pig.PMI.", gene, ".Exp.hex.pdf"), width = 15, height = 6)
plot_grid(plotlist = p, nrow = 1, ncol = 3) %>% print()
dev.off() 


##viridis(4)




