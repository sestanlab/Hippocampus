## Get the co-expression of METTL7B & its interactors in different cell types
source("../scripts/hip.fun.R")
source("./pie.fun.R")


order_list <- list(HIP = c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","SUB proximal ROBO1 SEMA3E","SUB distal FN1 NTNG1"), 
					EC = c("EC L2 CUX2 PDGFD","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L2 CUX2 LAMA3","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L5 RORB TPBG","EC L6 THEMIS CDH13","EC L6 THEMIS RGS12","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2"), 
					MTG = c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1"), 
					PFC = readRDS(file = "~/project/PFC/data/PFC.cluster.order.rds")[c(1:14, 16:22, 15, 23:40)])

cls_order <- unlist(order_list, use.names = FALSE)
sel_genes <- c("METTL7B", "APP", "RTN3", "RTN4", "NAE1", "LRP1", "APBB1", "PPP3CA", "UQCRFS1", "CDK5", "CALM1", "EIF2AK3", "ADAM17", "CASP3", "COX6B1", "CYCS", "ATP5F1B", "ITPR1", "ITPR2", "ATF6", "PDIA3")


##----------------------------------------------------------------------------
## Prepare the pie plot
coexpFile <- paste0(inputdir, "METTL7B_coexp_res.rds")
if (!file.exists(coexpFile)){
	seulist <- readRDS(file = paste0(inputdir, "ExN.4reg.mettl7b.coexp.slim.seu.rds"))

	## Get the co-expression for each gene pair across clusters
	###!!!!!!!!!!!make sure the group.by ==  "subtypes"
	coexp_res <- lapply(names(seulist), function(reg) {
		print(paste0("Working on region ", reg))
		message(paste0("The following genes are not present in region: ", reg))
		print(setdiff(sel_genes, rownames(seulist[[reg]])))
		reg_res <- lapply(order_list[[reg]], function(cls) {
						seu <- subset(seulist[[reg]], subtypes == cls)
						data <- seu$RNA@data
						mett7b_idx <- data["METTL7B", ] != 0
						bg_ratio <- data.frame(gene = "METTL7B", 
											mett7b = mean(mett7b_idx),
											coexp = 0, 
											other = 0, 
											noexp = mean(!mett7b_idx), 
											cluster = cls, 
											stringsAsFactors = FALSE)

						ratios <- lapply(sel_genes, function(gene) {
									if ( (!gene %in% rownames(data)) | gene == "METTL7B"){
										if (gene == "ATP5F1B" & reg %in% c("MTG")){
											gene_idx <- data["ATP5B", ] != 0
											ratio <- data.frame(gene = "ATP5F1B", 
												mett7b = mean(mett7b_idx & (!gene_idx)),
												coexp = mean(mett7b_idx & gene_idx), 
												other = mean((!mett7b_idx) & gene_idx), 
												noexp = mean(!(mett7b_idx | gene_idx)), 
												stringsAsFactors = FALSE)
										} else {
											ratio <- data.frame(gene = gene, 
												mett7b = mean(mett7b_idx),
												coexp = 0, 
												other = 0, 
												noexp = mean(!mett7b_idx), 
												stringsAsFactors = FALSE)
										}
										
									} else {
										gene_idx <- data[gene, ] != 0
										ratio <- data.frame(gene = gene, 
												mett7b = mean(mett7b_idx & (!gene_idx)),
												coexp = mean(mett7b_idx & gene_idx), 
												other = mean((!mett7b_idx) & gene_idx), 
												noexp = mean(!(mett7b_idx | gene_idx)), 
												stringsAsFactors = FALSE)
									}
									return(ratio)
							}) %>% 
								do.call(rbind, .) %>%
								mutate(cluster = cls)
						ratios <- rbind(ratios, bg_ratio)
						print(paste0("Finish cluster ", cls))
						return(ratios)
					}) %>% 
						do.call(rbind,.) %>%
						mutate(region = reg)
		return(reg_res)
		}) %>% 
					do.call(rbind,.)
	saveRDS(coexp_res, file = coexpFile)
}

## Missing genes in MTG: "ATP5F1B", replaced with "ATP5B" based on genecards annotation: https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATP5F1B&keywords=ATP5B
## Missing genes in PFC: "ATP5F1B" "ITPR2" 


##----------------------------------------------------------------------------
## Plot the scatter Pie
coexp_res <- readRDS(file = paste0(inputdir, "METTL7B_coexp_res.rds"))


## Set some parameters defining the intervals among pies
small_gap <- 0.6
big_gap <- 1
small_pie <- 0.2
##big_cls <- c("DG GC PROX1 EGR1", "DG GC PROX1 PDLIM5", "DG GC PROX1 SGCZ", "CA2 CFAP299 HGF", "CA3 CFAP299 SYN3", "SUB distal FN1 NTNG1", "SUB proximal ROBO1 COL5A2")
big_cls <- cls_order[sapply(cls_order, function(x) rowSums(coexp_res[coexp_res$cluster == x, c("mett7b", "coexp")])[1] >= 0.075)] ## Only enlarge the METTL7B-expressed clusters
coexp_res$radius <- ifelse(coexp_res$cluster %in% big_cls, big_gap/2, small_pie)



## Generate the x-axis breaks for the clusters
cls_vec <- rep(small_gap, length(cls_order)-1)
for (i in big_cls){
	idx <- which(cls_order == i)
	if (idx != 1){
		cls_vec[idx-1] <- big_gap
	} 
	cls_vec[idx] <- big_gap
}
x_vec <- round(cumsum(c(0, cls_vec)), digits = 2) %>% setNames(., cls_order)


## X axis text color
text_color <- ifelse(cls_order %in% big_cls, "red", "black")


## Where to put the dashed line to separate different regions
reg_sep <- c()
sep_idx <- sapply(order_list[1:3], function(x) which(cls_order == rev(x)[1]))
for (i in 1:length(sep_idx)){
	loc <- sep_idx[i]
	diff <- setNames((x_vec[loc+1] - x_vec[loc]), NULL)
	if (isTRUE(all.equal(diff, small_gap))){
		reg_sep[i] <- mean(c(x_vec[loc+1], x_vec[loc]))
	} else {
		reg_sep[i] <- (3/4) * big_gap - small_pie/2 + x_vec[loc]
	}
}
##reg_sep <- c(10, 15, 20)


p <- PlotScatterPie(coexp_res, group.col = "cluster", feature.col = "gene", r.col = "radius", pie_scale = 1, split.order = c("mett7b", "coexp", "other", "noexp"), group.order = cls_order, feature.order = sel_genes, rsf = 1, cluster.coord = x_vec) +
		geom_vline(xintercept = reg_sep, linetype = "dashed", size = 0.2) +
		theme(axis.text.x = element_text(size = 3.5, color = text_color), axis.text.y = element_text(size = 6)) +
		theme(legend.position = "bottom", panel.grid.major = element_blank())

pdf(paste0(outputdir, "SF7.METTL7B_interactors_coexp.pdf"), width = 8.5, height = 5)
print(p)
dev.off()





