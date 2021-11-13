## Plot the heatmap of 4 region markers
source("../scripts/hip.fun.R")
library(gridExtra)

cls.order <- list(HIP = c("DG GC PROX1 PDLIM5", "DG GC PROX1 SGCZ", "CA2 CFAP299 HGF", "SUB proximal ROBO1 COL5A2", "CA3 CFAP299 SYN3", "DG MC ARHGAP24 DLC1", "CA1 dorsal GRIK1 GRM3", "CA1 ventral ACVR1C SYT13", "SUB distal FN1 NTNG1", "SUB proximal ROBO1 SEMA3E"), 
					EC = c("EC L2 CUX2 PDGFD","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L2 CUX2 LAMA3","EC L2 RELN BMPR1B","EC L2 RELN BCL11B","EC L3 PCP4 ADARB2","EC L5 RORB TPBG","EC L5 RORB TLL1","EC L6 THEMIS CDH13","EC L6 THEMIS RGS12","EC L5 BCL11B ADRA1A","EC L6 TLE4 SULF1","EC L56 TLE4 NXPH2","EC L6b TLE4 CCN2"), 
					MTG = c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1"), 
					PFC = readRDS(file = "~/project/PFC/data/PFC.cluster.order.rds")[c(1:14, 16:22, 15, 23:40)][c(1:20, 22:26, 21, 31:34, 27:30, 35:40)])


## Get the average expression of all shared genes across 4 regions
avgFile <- paste0(inputdir, "AVG.ExN.4reg.rds")
if (!file.exists(avgFile)){
	seulist <- readRDS(file = paste0(dataDir, "FourReg.ExN.final.01172021.rds"))
	exp_genes <- Reduce("intersect", lapply(seulist, function(x) rownames(x$RNA)))


	avglist <- lapply(names(seulist), function(reg) {
		xx <- seulist[[reg]]
		yy <- xx[, xx@meta.data$subtypes %in% cls.order[[reg]]]
		Idents(yy) <- "subtypes"
		avg <- log(AverageExpression(yy, features = exp_genes)$RNA + 1)
		return(avg)
		}) %>%
			setNames(., names(seulist))
	saveRDS(avglist, file = avgFile)
}


## Get the marker list
all.mar <- readRDS(file = paste0("../MF3_cellregion/load_files/", "Markers.4reg.ExN.rds"))
max.marker <- 50
mars <- lapply(names(cls.order), function(reg) {
	reg.mars <- lapply(cls.order[[reg]], function(x) {
		xx <- all.mar %>%
				subset(region == reg & cluster == x) %>%
				subset(pct.1 >= 0.3 & avg_logFC >= log(2) & ratio_fc >= 1.1) %>%
				mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
				subset(p_val_adj <= 0.01)
		yy <- setNames(xx$avg_logFC, xx$gene)
		gg <- names(yy)[order(yy, decreasing = TRUE)[1:min(length(yy), max.marker)]]
		if (is.na(gg)){
			print(paste0("region, ", reg, " cluster ", x))
		}
		return(gg)
		}) %>% setNames(., cls.order[[reg]])
	return(reg.mars)
	}) %>% setNames(., names(cls.order))

plot_glist <- lapply(mars, function(x) unlist(x) %>% unique())


avglist <- readRDS(file = paste0(inputdir, "AVG.ExN.4reg.rds"))
plot_data <- lapply(names(avglist), function(x) {
					yy <- avglist[[x]] %>% t() %>%
								scale() %>% t() %>%
								MinMax(., min = -1.8, max = 1.8) %>%
								.[, cls.order[[x]]]
					colnames(yy) <- paste0(x, "|", colnames(yy))
					return(yy)
					}) %>%
						do.call(cbind, .)



plot_list=list()
col_anno <- data.frame(row.names = colnames(plot_data), 
				region = extract_field(colnames(plot_data), 1, "|"), 
				stringsAsFactors = FALSE)
colvec <- c("#FF420E","#89DA59","#4CB5F5","#FFBB00") %>% 
				setNames(., c("HIP", "EC", "PFC", "MTG"))
for (reg in c("HIP", "EC", "MTG", "PFC")){
	x <- pheatmap::pheatmap(plot_data[plot_glist[[reg]], ], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(c("#fcfbfd", "#dadaeb", "#3f007d"))(30), border_color = NA, show_rownames = FALSE, show_colnames = FALSE, fontsize_col = 8, gaps_col = cumsum(sapply(avglist, ncol)[1:3]), labels_col = extract_field(colnames(plot_data), "rm_start", "|"), annotation_col = col_anno, annotation_colors = list(region = colvec), annotation_legend = FALSE)
	plot_list[[reg]] = x[[4]]
}

gg <- arrangeGrob(grobs = plot_list, nrow = 1, ncol = 4)
ggsave(paste0(outputdir, "SF4.marker.heatmap.pdf"), gg, width = 16, height = 4)



