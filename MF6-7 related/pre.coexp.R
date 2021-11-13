## Prepare the co-expression seurat object list
source("../scripts/hip.fun.R")
source("./pie.fun.R")


order_list <- list(HIP = c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","SUB proximal ROBO1 SEMA3E","SUB distal FN1 NTNG1"), 
					EC = c("EC L2 CUX2 PDGFD","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L2 CUX2 LAMA3","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L5 RORB TPBG","EC L6 THEMIS CDH13","EC L6 THEMIS RGS12","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2"), 
					MTG = c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1"), 
					PFC = readRDS(file = "~/project/PFC/data/PFC.cluster.order.rds")[c(1:14, 16:22, 15, 23:40)])


sel_genes <- c("METTL7B", "APP", "RTN3", "RTN4", "NAE1", "LRP1", "APBB1", "PPP3CA", "UQCRFS1", "CDK5", "CALM1", "EIF2AK3", "ADAM17", "CASP3", "COX6B1", "CYCS", "ATP5F1B", "ITPR1", "ITPR2", "ATF6", "PDIA3")


seulist <- readRDS(file = paste0(dataDir, "FourReg.ExN.final.01172021.rds"))
hdfc <- readRDS(file = "~/project/PFC/side/load_files/Human.noLiftOver.recount.filt.seu.rds")
hdfc@meta.data$subtypes <- hdfc@meta.data$fig1cluster
newlist <- lapply(names(order_list), function(reg) {
	new_genes <- ifelse_check(reg == "MTG", gsub("ATP5F1B", "ATP5B", sel_genes), sel_genes)
	subseu <- seulist[[reg]][new_genes, ]
	return(subseu)
	}) %>% 
			setNames(., names(order_list))



match_genes <- sapply(sel_genes, function(x) grep(paste0("-", x, "$"), rownames(hdfc), value = TRUE))
new_counts <- hdfc$RNA@counts[match_genes, ]; rownames(new_counts) <- extract_field(rownames(new_counts), "rm_start", "-") %>% setNames(., NULL)
new_data <- hdfc$RNA@data[match_genes, ]; rownames(new_data) <- extract_field(rownames(new_data), "rm_start", "-") %>% setNames(., NULL)
new_pfc <- hdfc[match_genes, intersect(colnames(seulist$PFC), colnames(hdfc))]
new_pfc[["RNA"]] <- CreateAssayObject(data = new_data[, colnames(new_pfc)], min.cells = 0, min.features = 0)


newlist$PFC <- new_pfc
saveRDS(newlist, file = paste0(inputdir, "ExN.4reg.mettl7b.coexp.slim.seu.rds"))










