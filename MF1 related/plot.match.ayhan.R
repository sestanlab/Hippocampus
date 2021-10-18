## Plot the UMAP of the dataset
source("../scripts/hip.fun.R")
library(ggdendro)

hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds"))



avgs <- readRDS(file = paste0("../MF1_overview/load_files/", "Avg.HIPEC.subtypes.rds"))
#hvg  <- readRDS(file = paste0(inputdir, "HVG.HIPEC.subtypes.rds"))
cls_order <- c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","EC L2 CUX2 PDGFD","EC L2 CUX2 LAMA3","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L6 THEMIS CDH13","EC L5 RORB TPBG","EC L6 THEMIS RGS12","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","SUB distal FN1 NTNG1","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","SUB proximal ROBO1 SEMA3E","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2","CR RELN NDNF","InN MEIS2 SHISAL2B","InN SST NPY","InN SST ADAMTS12","InN SST EPB41L4A","InN SST OTOF","InN PVALB PLEKHH2","InN LHX6 AC008415.1","InN PVALB MEPE","InN PVALB PLCL1","InN VIP NOX4","InN VIP SCTR","InN VIP ABI3BP","InN VIP SCML4","InN VIP CHRNA2","InN VIP PENK","InN NR2F2 PTPRK","InN NR2F2 MIR4300HG","InN NR2F2 SLC17A8","InN NR2F2 ANO2","InN NR2F2 DDR2","InN LAMP5 CHST9","InN LAMP5 KIT","InN LAMP5 NMBR","Astro AQP4 CHRDL1","Astro AQP4 GFAP","OPC PDGFRA EGR1","OPC PDGFRA GRIA4","COP GPR17 ADAM33","Oligo CPXM2 KANK4","Oligo OPALIN LAMA2","Oligo OPALIN LINC01098","Oligo OPALIN SLC5A11","Micro C1QB CD83","Micro C1QB P2RY12","Macro F13A1 COLEC12","Myeloid LSP1 LYZ","T SKAP1 CD247","aEndo DKK2 FBLN5","Endo CLDN5 VWF","PC CLDN5 ABCC9","vSMC ABCC9 P2RY14","aSMC ACTA2 TAGLN","VLMC COL1A1 COL1A2")    

seu <- readRDS(file = paste0(inputdir, "Konopka.Inte.all.Seurat.2500.slim.rds"))
Idents(seu) <- "Cluster"
pub_avg <- log(AverageExpression(seu)$RNA + 1)
pub_exp <- rownames(seu)[Matrix::rowSums(seu$RNA@counts!=0) >= 20]




pub_order <- c("Den.Gyr1", "Den.Gyr2", "Den.Gyr3", "Pyr1", "Pyr2", "In1", "In2", "In3", "Astro1", "Astro2", "Astro3","OPC1", "OPC2", "OPC3", "OPC4", "Olig1", "Olig2", "Olig3", "Olig4", "Olig5", "Micro1", "Micro2", "Micro3", "Endo")

hvg <- rownames(hipec$pca@feature.loadings) 
gene_use <- intersect(hvg, pub_exp)
pub_avg <- pub_avg[gene_use, pub_order]
colnames(pub_avg) <- gsub("Den.Gyr", "DG GC", colnames(pub_avg)) %>%
						gsub("Pyr", "Pyr ExN", .) %>%
						gsub("Olig", "Oligo", .) %>%
						gsub("^In", "InN", .)


cor_mat <- cor(avgs[gene_use, ], pub_avg, method = "p")

pdf(paste0(outputdir, "SF1.Match.ayhan.heatmap.pdf"), width = 7, height = 7)
pheatmap::pheatmap(cor_mat[cls_order, ], cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(3))(30), border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 8, fontsize_row = 4)
dev.off()


