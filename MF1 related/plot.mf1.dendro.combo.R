## Plot the dendrogram of all clusters & related bar & dotplots
source("../scripts/hip.fun.R")
library(ggdendro)

## Read the dataset
hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds"))
group.by <- "fig1cluster"
   


cls_order <- c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","EC L2 CUX2 PDGFD","EC L2 CUX2 LAMA3","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L6 THEMIS CDH13","EC L5 RORB TPBG","EC L6 THEMIS RGS12","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","SUB distal FN1 NTNG1","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","SUB proximal ROBO1 SEMA3E","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2","CR RELN NDNF","InN MEIS2 SHISAL2B","InN SST NPY","InN SST ADAMTS12","InN SST EPB41L4A","InN SST OTOF","InN PVALB PLEKHH2","InN LHX6 AC008415.1","InN PVALB MEPE","InN PVALB PLCL1","InN VIP NOX4","InN VIP SCTR","InN VIP ABI3BP","InN VIP SCML4","InN VIP CHRNA2","InN VIP PENK","InN NR2F2 PTPRK","InN NR2F2 MIR4300HG","InN NR2F2 SLC17A8","InN NR2F2 ANO2","InN NR2F2 DDR2","InN LAMP5 CHST9","InN LAMP5 KIT","InN LAMP5 NMBR","Astro AQP4 CHRDL1","Astro AQP4 GFAP","OPC PDGFRA EGR1","OPC PDGFRA GRIA4","COP GPR17 ADAM33","Oligo CPXM2 KANK4","Oligo OPALIN LAMA2","Oligo OPALIN LINC01098","Oligo OPALIN SLC5A11","Micro C1QB CD83","Micro C1QB P2RY12","Macro F13A1 COLEC12","Myeloid LSP1 LYZ","T SKAP1 CD247","aEndo DKK2 FBLN5","Endo CLDN5 VWF","PC CLDN5 ABCC9","vSMC ABCC9 P2RY14","aSMC ACTA2 TAGLN","VLMC COL1A1 COL1A2")    



## get the average expression of the dataset
if (FALSE){
	Idents(hipec) <- group.by
	avgs <- log(AverageExpression(hipec)$RNA + 1)
	saveRDS(avgs, file = paste0(inputdir, "Avg.HIPEC.subtypes.rds"))
}


## Use the merge HVG from 3 cell types
if (FALSE){
	inn <- readRDS(file = paste0("../MF3_cellregion/load_files/", "HIPEC.", "InN", ".final.rename.seu.rds"))
	exn <- readRDS(file = paste0("../MF3_cellregion/load_files/", "HIPEC.", "ExN", ".final.rename.seu.rds"))
	nnc <- readRDS(file = paste0("../MF3_cellregion/load_files/", "HIPEC.", "NNC", ".final.rename.seu.rds"))

	hvg <- c(rownames(inn$pca@feature.loadings), rownames(exn$pca@feature.loadings)) %>% 
				c(., rownames(nnc$pca@feature.loadings)) %>%
					unique()
	saveRDS(hvg, file = paste0(inputdir, "HVG.HIPEC.subtypes.rds"))
}


avgs <- readRDS(file = paste0(inputdir, "Avg.HIPEC.subtypes.rds"))
hvg  <- readRDS(file = paste0(inputdir, "HVG.HIPEC.subtypes.rds"))
avg_use <- avgs[hvg, ]


## Build the dendrogram
reg_colors <- paste0(c("#b2182b", "#b8e186", "#4d9221", "#80ffff", "#fee090"), "") %>%
				setNames(., c("DG", "CA24", "CA1", "SUB", "EC"))
meta_data <- hipec@meta.data


theme_empty <- function(x) {
	theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), legend.position = "none")
}


# Rectangular lines 
dendro <- dist(MinMax(scale(t(avg_use)), min = -4, max = 4)) %>% 
			hclust(., method = "ward.D2") %>% 
			as.dendrogram()


dendro <- reorder(dendro, match(colnames(avgs), rev(cls_order)), agglo.FUN=mean)
pdf(paste0(outputdir, "cluster.dendrogram.pdf"), width = 7, height = 12)
ggdendrogram(dendro, rotate = TRUE)
dev.off()


ddata <- dendro_data(dendro, type = "rectangle")
den_p <- ggplot(segment(ddata)) + 
			geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
			coord_flip() + 
			scale_y_reverse(expand = c(0.2, 0)) +
			theme_empty()


## Plot the region contribution
clssize_data <- meta_data %>%
			group_by(!!sym(group.by), group) %>%
			summarize(ncells = n()) %>%
			mutate(ncells = sqrt(ncells)) %>%
			ungroup()
regcontri_data <- meta_data %>%
			group_by(!!sym(group.by), region_new) %>%
			summarize(ncells = n()) %>%
			ungroup() %>%
			group_by(region_new) %>%
			mutate(ratio = ncells/sum(ncells)) %>%
			ungroup() %>%
			group_by(!!sym(group.by)) %>%
			mutate(rratio = ratio * 100/sum(ratio))

indicontri_data <- meta_data %>%
			group_by(!!sym(group.by), samplename) %>%
			summarize(ncells = n()) %>%
			ungroup() %>%
			group_by(samplename) %>%
			mutate(ratio = ncells/sum(ncells)) %>%
			ungroup() %>%
			group_by(!!sym(group.by)) %>%
			mutate(rratio = ratio * 100/sum(ratio))

size_bks <- seq(0, 200, 50)
size_p <- ggplot(clssize_data, aes_string(x = "ncells", y = group.by, fill = "group")) +
		geom_bar(positio = "stack", color = NA, stat = "identity") +
		scale_fill_manual(values = c(ExN = "#C51B7D", InN = "#4D9221", NNC = "#999999")) +
		scale_y_discrete(limits = rev(cls_order))+
		scale_x_continuous(position = "top", breaks = size_bks, labels = as.character(size_bks^2)) +
		theme_empty() +
		theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45,hjust = 0), axis.ticks.x = element_line(size = 0.2), axis.line.x = element_line(size = 0.2))
reg_p <- ggplot(regcontri_data, aes_string(x = "rratio", y = group.by, fill = "region_new")) +
		geom_bar(positio = "stack", color = NA, stat = "identity") +
		scale_fill_manual(values = reg_colors) +
		scale_y_discrete(limits = rev(cls_order))+
		scale_x_continuous(position = "top") +
		theme_empty() +
		theme(axis.text.x = element_text(size = 7, angle = 45,hjust = 0), axis.ticks.x = element_line(size = 0.2), axis.line.x = element_line(size = 0.2))
ind_p <- ggplot(indicontri_data, aes_string(x = "rratio", y = group.by, fill = "samplename")) +
		geom_bar(positio = "stack", color = NA, stat = "identity") +
		scale_fill_manual(values = gg_color_hue(6) %>% setNames(., levels(as.factor(indicontri_data$samplename)))) +
		scale_y_discrete(limits = rev(cls_order))+
		scale_x_continuous(position = "top") +
		theme_empty() +
		theme(axis.text.x = element_text(size = 7, angle = 45,hjust = 0), axis.ticks.x = element_line(size = 0.2), axis.line.x = element_line(size = 0.2))


genes <- c("RBFOX3", "SLC17A7", "PROX1", "CFAP299", "CUX2", "RELN", "FN1", "RORB","TLE4", "GAD1", "GAD2", "MEIS2", "LHX6", "SST", "PVALB", "VIP", "NR2F2", "LAMP5", "AQP4", "PDGFRA", "MOBP", "PTPRC", "C1QB", "CLDN5", "ACTA2", "CEMIP")

Idents(hipec) <- group.by
exp_p <- DotPlot(hipec, features = genes, cols = c("lightgrey", "red"), dot.min = 0.05, dot.scale = 3) +
			scale_x_discrete(limits = genes, position = "top") +
			scale_y_discrete(limits = rev(cls_order)) +
			theme(legend.position = "right", axis.title = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text.x = element_text(size = 7, angle = 45,hjust = 0), axis.text.y = element_blank())
pdf(paste0(outputdir, "MF1.Cluster.combo.pdf"), width = 12, height = 7)
plot_grid(den_p, size_p, reg_p, ind_p, exp_p, nrow = 1, ncol = 5, rel_widths = c(1, 1.2, 0.7, 0.7, 2.5), align = "h") %>% print()
dev.off()

 







