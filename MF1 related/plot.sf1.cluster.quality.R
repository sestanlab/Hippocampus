## Plot the number of cells in each sample 
source("../scripts/hip.fun.R")
library(lemon)

hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds")) 
cls_order <- c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","EC L2 CUX2 CALB1","EC L2 CUX2 PDGFD","EC L2 CUX2 LAMA3","EC L2 CUX2 IL1RAPL2","EC L3 PCP4 ADARB2","EC L6 THEMIS CDH13","EC L5 RORB TLL1","EC L5 RORB TPBG","EC L6 THEMIS RGS12","EC L56 TLE4 NXPH2","SUB proximal ROBO1 SEMA3E","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2","EC L5 BCL11B ADRA1A","SUB distal FN1 NTNG1","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","CR RELN NDNF","InN SST NPY","InN SST ADAMTS12","InN SST EPB41L4A","InN SST OTOF","InN PVALB PLEKHH2","InN PVALB MEPE","InN PVALB PLCL1","InN LHX6 AC008415.1","InN VIP NOX4","InN VIP SCTR","InN VIP ABI3BP","InN VIP SCML4","InN VIP CHRNA2","InN VIP PENK","InN NR2F2 PTPRK","InN NR2F2 MIR4300HG","InN LAMP5 CHST9","InN LAMP5 KIT","InN LAMP5 NMBR","InN NR2F2 SLC17A8","InN NR2F2 ANO2","InN NR2F2 DDR2","InN MEIS2 SHISAL2B", "Astro AQP4 CHRDL1","Astro AQP4 GFAP","OPC PDGFRA EGR1","OPC PDGFRA GRIA4","COP GPR17 ADAM33","Oligo CPXM2 KANK4","Oligo OPALIN LAMA2","Oligo OPALIN LINC01098","Oligo OPALIN SLC5A11","Micro C1QB CD83","Micro C1QB P2RY12","Macro F13A1 COLEC12","Myeloid LSP1 LYZ","T SKAP1 CD247","aEndo DKK2 FBLN5","Endo CLDN5 VWF","PC CLDN5 ABCC9","vSMC ABCC9 P2RY14","aSMC ACTA2 TAGLN","VLMC COL1A1 COL1A2")    


meta_use <- hipec@meta.data
rep.by <- "repname"
ind.by <- "samplename"
reg.by <- "region_new"
group.by <- "fig1cluster"



###----------------------------------------------------------------------------------------
## Plot nUMI, nGene, mapped reads by cluster
plot_cols <- c("nFeature_RNA", "nCount_RNA")

plot_data <- meta_use[, c(plot_cols, group.by)] %>%
				reshape2::melt(id.vars = group.by, measure.vars = plot_cols, variable.name = "quality") %>%
				mutate(!!group.by := as.character(!!sym(group.by))) %>%
				mutate(!!group.by := factor(!!sym(group.by), levels = rev(cls_order)))
cluster_colors <- rep(c("#C51B7D", "#4D9221", "#999999"), c(25, 24, 20))

plot_data$quality <- plot_data$quality %>%
						gsub("nCount_RNA", "nUMIs", .) %>%
						gsub("nFeature_RNA", "nGenes", .)
plot_data$quality <- factor(plot_data$quality, levels = c("nUMIs", "nGenes"))


p <- ggplot(plot_data, aes_string(y = group.by, x = "value")) +
		geom_violin(fill = "grey", scale = "width", size = 0.05, adjust = 1.5, trim =TRUE, alpha = 0.8) + 
		geom_boxplot(width=0.3, outlier.shape = NA, lwd= 0.05, alpha = 0.3) + 
		coord_capped_cart(top='both', left='both') +
		theme_classic() + 
		scale_x_continuous(position = "top") + 
		theme(panel.border=element_blank(), axis.line=element_line(size = 0.25), axis.ticks = element_line(size = 0.25),axis.text.x = element_text(size = 10, hjust = 0, angle = 45), axis.text.y = element_text(size = 7, colour = rev(cluster_colors)), axis.title = element_blank()) + 
		facet_wrap(facets = vars(quality), nrow = 1, ncol = length(plot_cols), strip.position = "top", scales = "free_x") + 
		theme(strip.background = element_blank(), strip.text = element_text(size = 11, face = "bold"), panel.spacing = unit(0.05, "in"), legend.title = element_text(size = 10), legend.text = element_text(size = 12), strip.placement = "outside")
pdf(paste0(outputdir, "SF1.Cluster_quality.pdf"), width = 5, height = 7)##, units = "in",  res = 300)
print(p)
dev.off()
