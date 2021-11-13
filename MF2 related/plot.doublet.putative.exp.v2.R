## Plot the integration UMAP figures
args <- commandArgs(trailingOnly = TRUE) 
source("../scripts/hip.fun.R")
library(ggpubr)



##------------------------------------------------------------------------
## 3. putative human cells [Heatmap exp]  
hm <- readRDS(file = paste0("../SF2_doublets/load_files/", "HM.Inte.GCDB.v2.slim.rds"))
pt_cells <- colnames(hm)[hm$umap@cell.embeddings[, 1] <= 1 & hm$umap@cell.embeddings[, 1] >= -0.8 & hm$umap@cell.embeddings[, 2] >= 4.6  & hm$umap@cell.embeddings[, 2] <= 8 & hm@meta.data$species == "Human"]
pt_cells <- c(pt_cells, c("HSB179_6_eDG_AATTTCCGTCCTACAA"))##"HSB628_1_DG_CCGATGGGTGTTAAAG", 


hdg <- readRDS(paste0(inputdir, "Human_DG_seu.rds"))
hdg@meta.data$mres <- gsub("InN CGE", "InN", hdg@meta.data$mres) %>%
						gsub("InN MGE", "InN", .)
hdg_withdb <- readRDS(paste0("../SF2_doublets/load_files/", "HDG.inte.GCdoublets.seurat.slim.rds"))



## Get the bg cells
set.seed(42)
bg_cols <- c("GC", "InN", "OPC", "Astro") ## immnue
bg_cells <- lapply(bg_cols, function(x) {
	sample(setdiff(colnames(hdg)[hdg@meta.data$mres == x], pt_cells), length(pt_cells))
	}) %>%
			setNames(., bg_cols)
all_cells <- c(pt_cells, unlist(bg_cells, use.names = FALSE))

cell_cate <- lapply(names(bg_cells), function(x) {
	setNames(rep(x, length(bg_cells[[x]])), bg_cells[[x]])
	}) %>%
			unlist() %>%
			c(., setNames(rep("putative", length(pt_cells)), pt_cells))


all_genes <- c("SLC17A7", "GAD1", "GAD2", "LHX6", "SATB1", "SST", "PVALB", "NR2F2", "VIP", "RELN", "OLIG2", "PDGFRA", "MOBP", "SLC1A3", "AQP4", "PTPRC", "CLDN5", "CEMIP", "PROX1", "DCX", "CALB2", "DPYSL3")
exp_data <- hdg_withdb$RNA@counts[all_genes, pt_cells] %>%
				cbind(., hdg$RNA@counts[all_genes, setdiff(all_cells, pt_cells)]) %>%
				.[, all_cells]
plot_data <- exp_data %>%
					as.matrix() %>%
					reshape2::melt() %>%
					setNames(., c("gene", "cell", "exp")) %>%
					mutate(gene = as.character(gene), cell = as.character(cell)) %>%
					##mutate(ctp = ifelse(cell %in% pt_cells, "yes", "no")) %>%
					##mutate(ctp = factor(ctp, levels = c("yes",  "no"))) %>%
					mutate(ctp = cell_cate[cell]) %>%
					mutate(ctp = factor(as.character(ctp), levels = c("putative", bg_cols))) %>%
					mutate(exp = MinMax(exp, min = 0, max = 3)) %>%
					mutate(exp = as.character(exp)) 

colors <- viridis(4)
putap <- ggplot(plot_data) + 
				geom_tile(aes_string(y = "gene", x = "cell", fill = "exp"), color = "black", height = 1, width = 1) + 
				theme_classic() +
				##scale_x_discrete(limits = all_cells) +
				scale_y_discrete(limits = rev(all_genes)) +
				scale_fill_manual(values = colors) +
				scale_color_identity() + 
				facet_wrap(vars(ctp), nrow = 1, ncol = length(bg_cols) + 1, scales = "free_x") +
	            theme(panel.spacing = unit(0.02, "in"), strip.background = element_blank(), strip.text = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 10), axis.line = element_blank(), legend.position = "bottom", axis.title = element_blank())
pdf(paste0(outputdir, "SF2.doublet.putative.exp.v2.pdf"), width = 6, height = 8)
print(putap)
dev.off()

