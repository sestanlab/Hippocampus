## Find HIP-specific markers and visualize it.
source("../scripts/hip.fun.R")
library(ggpubr)
seulist <- readRDS(file = paste0(dataDir, "FourReg.ExN.final.01172021.rds"))

exp_genes <- Reduce("intersect", lapply(seulist, function(x) rownames(x$RNA)))

## Get the expression ratio
cls.order <- list(HIP = c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","SUB proximal ROBO1 SEMA3E","SUB distal FN1 NTNG1"),
					EC = c("EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L2 CUX2 LAMA3","EC L2 CUX2 PDGFD","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L5 RORB TPBG","EC L6 THEMIS CDH13","EC L6 THEMIS RGS12","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2"), 
					MTG = c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1"), 
					PFC = readRDS(file = "~/project/PFC/data/PFC.cluster.order.rds")[c(1:14, 16:22, 15, 23:40)])
if (FALSE) {
	ratio_list <- lapply(names(cls.order), function(reg) {
		seu <- seulist[[reg]]
		cls <- cls.order[[reg]]
		ratio <- lapply(cls, function(x) {
			xx <- seu[, seu@meta.data$subtypes == x]
			message(paste0("Finish cluster ",  x))
			return(Matrix::rowMeans(xx$RNA@data != 0))
			}) %>% setNames(., cls) %>%
					as.data.frame(., check.names = FALSE) %>%
					as.matrix()
		return(ratio)
		}) %>% 
				setNames(., names(cls.order))
	saveRDS(ratio_list, file = paste0(inputdir, "ExN.4reg.exp.ratio.rds"))
}


## Region-specific gene expression. 
ratio_list <- readRDS(file = paste0(inputdir, "ExN.4reg.exp.ratio.rds"))

hip_x <- ratio_list$HIP
hip_yy <- apply(hip_x[exp_genes, ], 1, max)
hip_glist <- apply(hip_x[exp_genes, ], 1, function(y) colnames(x)[which(y == max(y))[1]]) 
hip_order <- lapply(colnames(hip_x), function(mm) names(hip_glist)[hip_glist == mm]) %>% unlist()


reg_max <- lapply(ratio_list, function(x) apply(x[exp_genes, ], 1, max)) %>%
			as.data.frame(., check.names = FALSE)
reg_max <- reg_max[apply(reg_max,1,max) >= 0.1, ]
##reg_max$ratio_fc <- reg_max$HIP/apply(reg_max[, 2:4], 1, max)

hip_genes <- rownames(reg_max)[reg_max$HIP >= 0.3 & 
					reg_max$HIP/apply(reg_max[, 2:4], 1, max) >= 2.2 & 
					apply(reg_max[, 2:4], 1, max) <= 0.22]
ec_genes <- rownames(reg_max)[reg_max$EC >= 0.3 & 
					reg_max$EC/apply(reg_max[, c(1,3,4)], 1, max) >= 2.2 & 
					apply(reg_max[, c(1,3,4)], 1, max) <= 0.22]
mtg_genes <- rownames(reg_max)[reg_max$MTG >= 0.3 & 
					reg_max$MTG/apply(reg_max[, c(1,2)], 1, max) >= 2.2 & 
					apply(reg_max[, c(1,2)], 1, max) <= 0.22]
pfc_genes <- rownames(reg_max)[reg_max$PFC >= 0.3 & 
					reg_max$PFC/apply(reg_max[, c(1,2)], 1, max) >= 2.2 & 
					apply(reg_max[, c(1,2)], 1, max) <= 0.22]
neo_genes <- intersect(mtg_genes, pfc_genes)
genes <- c(hip_order[hip_order %in% hip_genes], ec_genes, neo_genes[1:5])


## Plot the expression
color_vec <- c("#FF420E","#89DA59","#4CB5F5","#FFBB00") %>% 
				setNames(., c("HIP", "EC", "PFC", "MTG"))
plot_genes <- c("CUX2", "FGFR2", "SH3RF2", "PDGFD", "GLP2R","TLE4", "BCL11B", "RNF152", "FEZF2", "TRPC6", genes) %>% unique()
max_scale <- reg_max[plot_genes, 1:4] %>% max()
plist <- lapply(c("HIP", "EC", "MTG", "PFC"), function(reg) {
	seu <- seulist[[reg]]
	seu <- seu[, seu@meta.data$subtypes %in% cls.order[[reg]]]
	seu@meta.data$subtypes <- factor(as.character(seu@meta.data$subtypes), levels = cls.order[[reg]])
	p <- DotPlot(seu, features = plot_genes, cols = c("#D3D3D3", color_vec[reg]), dot.min = 0.10, dot.scale = 2.5, group.by = "subtypes") +
			coord_flip() +
			theme_cowplot()  +
			RotatedAxis() +
			scale_size(limits = c(10, 100 *max_scale), range = c(0, 2.5)) + 
			scale_x_discrete(limits = plot_genes %>% rev()) + 
			theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 7), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_line(size = 0.1))
	if (reg != "HIP"){
		p <- p + 
			theme(axis.text.y = element_blank())
	}
	return(p)
	})

##emps <- function(x) {return(paste(rep(" ", x), collapse = ""))}
ggarrange(plotlist = plist, nrow = 1, ncol = 4, align = "h", common.legend = TRUE, legend = "bottom", widths = c(0.8, 0.8, 1, 1.2)) %>%
annotate_figure(., 
		top = text_grob(paste0("Temporal lobe: allo-, meso- and neocortex", paste(rep(" ", 10), collapse = ""), "Prefrontal neocortex"), color = "black", hjust = 0.2, size = 10), 
		left = text_grob(paste0("Region-specific", paste(rep(" ", 5), collapse = ""), "Layer-specific"), hjust = 0.25, size = 10, color = "black", rot = 90)) %>%
ggexport(., filename = paste0(outputdir, "MF4.HIP.specific.genes.pdf"), width = 10, height = 4)
saveRDS(hip_genes, file = paste0(inputdir, "HIP-specific-genes.rds"))





hip_genes <- readRDS(file = paste0(inputdir, "HIP-specific-genes.rds"))
source("./bs.fun.R") ##"~/project/side_analysis/Kahle/SMARCC1/plot.gene.fun.R")
for (gene in hip_genes){
	PlotBSLine(gene = gene, file_name= "BS_Trend", output_dir = "./report/")
}











