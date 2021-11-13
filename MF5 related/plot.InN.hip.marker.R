## Find the HIP-specific markers in the cluster "InN SST ADAMTS12"
source("../scripts/hip.fun.R")
library(ggpubr)


if (FALSE){
seulist <- readRDS(file = paste0(inputdir, "Inte.InN.4reg.harmony.slim.rds")) %>%
		SplitObject(., split.by = "merge_region")

slim_list <- lapply(c("HIP", "EC", "MTG", "PFC"), function(reg) {
	subseu <- seulist[[reg]]
	all_cls <- table(subseu@meta.data$subtypes) %>%
					.[. >= 30] %>%
					names()
	if (reg %in% c("MTG", "PFC")){
		all_cls <- grep("SST", all_cls, value = TRUE)
	}
	seuuse <- subseu[, subseu@meta.data$subtypes %in% all_cls]
	return(seuuse)
	}) %>%
		setNames(., c("HIP", "EC", "MTG", "PFC"))


sst_ins <- lapply(c("EC", "MTG", "PFC"), function(reg) {
	subseu <- seulist[[reg]]
	all_cls <- table(subseu@meta.data$subtypes) %>%
					.[. >= 30] %>%
					names()
	sst_cls <- grep("SST", all_cls, value = TRUE)
	seuuse <- subseu[, subseu@meta.data$subtypes %in% sst_cls]
	return(seuuse)
	}) %>%
		setNames(., c("EC", "MTG", "PFC"))


mar_res <- lapply(names(sst_ins), function(reg) {
	subseu <- sst_ins[[reg]]
	Idents(subseu) <- "subtypes"
	mars <- FindAllMarkers(subseu, max.cells.per.ident = 2500, min.pct = 0.15, only.pos = TRUE, logfc.threshold = 0.15) %>%
						mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01))
	return(mars)
	}) %>%
		setNames(., names(sst_ins))

sst_markers <- FindMarkers(seulist$HIP, group.by = "subtypes", ident.1= "InN SST ADAMTS12", max.cells.per.ident = 2500, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.1)

ratio_list <- lapply(names(slim_list), function(reg) {
	subseu <- slim_list[[reg]]
	all_cls <- levels(as.factor(subseu@meta.data$subtypes))
	ratio <- lapply(all_cls, function(cls) 
		Matrix::rowMeans(subseu$RNA@data[, subseu@meta.data$subtypes == cls] != 0)
		) %>% 
			setNames(., all_cls) %>%
			as.data.frame(., check.names = FALSE)
	return(ratio)
	}) %>%
		setNames(., names(slim_list))
save(slim_list, mar_res, sst_markers, ratio_list, file = paste0(inputdir, "HIP.specific.InN.marker.res.Rdata"))
}



load(file = paste0(inputdir, "HIP.specific.InN.marker.res.Rdata"))

markers <- lapply(c("MTG", "PFC"), function(x){
	genes <- mar_res[[x]] %>%
				subset(ratio_fc >= 1.1) %>%
				subset(p_val_adj <= 0.05 & avg_logFC >= 0.1) %>%
				.$gene %>%
				unique()
	genes
	}) %>% unlist() %>% unique()


idx <- which(colnames(ratio_list$HIP) == "InN SST ADAMTS12")
max_ratio <- apply(ratio_list$HIP, 1, function(x) x[idx]/sort(x, decreasing = TRUE)[2])
max_ratio[is.na(max_ratio)] <- 0
max_genes <- names(max_ratio)[which(max_ratio >= 1.1)]

sst_fmar <- sst_markers %>%
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				subset(ratio_fc >= 2.5 & pct.2 <= 0.1)
##top_genes
hip_genes <- setdiff(sst_fmar$gene[sst_fmar$gene %in% max_genes], markers %>% unique())
max_ratio <- sapply(ratio_list, function(x) x[hip_genes, ] %>% apply(., 1, max) %>% max())
sum(c("EVC2", "EVC") %in% hip_genes)
sum(c("EVC2", "EVC") %in% max_genes)
length(hip_genes)


## Check the UMI numbers for each dataset
colors <- c("#FF420E","#89DA59","#4CB5F5","#FFBB00", "grey90") %>% 
				setNames(., c("HIP", "EC", "PFC", "MTG", "bg"))
plist <- lapply(names(slim_list), function(reg) {
	p <- DotPlot(slim_list[[reg]], features = hip_genes, group.by = "subtypes", cols = c("lightgrey", colors[[reg]]), dot.min = 0.05, dot.scale = 3) +
			coord_flip() + 
			theme_cowplot() +
			scale_radius(limits = c(5, 60), range = c(0, 2.5)) +
			RotatedAxis() +
			theme(legend.position = "none", axis.line = element_line(size = 0.2), axis.ticks.y = element_line(size = 0.2), axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6))
	if (reg != "HIP"){
		p <- p + 
			theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
	}
	return(p)
	})
ggarrange(plotlist = plist, nrow = 1, ncol = 4, align = "h", common.legend = TRUE, legend = "top", widths = c(1.3, 1, 1, 1)) %>%
		ggexport(filename = paste0(outputdir, "MF4.InN.hip.makers.pdf"), width = 8.5, height = 2)












