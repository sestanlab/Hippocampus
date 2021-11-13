## Compare the mouse and Rhesus markers
source("../scripts/hip.fun.R") 
library(ggpubr)


##------------------------------------------------------------------------
## Calculate the markers for each cluster in rhesus & mouse & pig
marFile <- paste0(inputdir, "Markers.RPM.res.rds")
if (!file.exists(marFile)){
	## Find Marker analysis
	allsps <- c("Mouse", "Pig", "Rhesus", "Human")
	allcls <- c("Astro", "RGL", "nIPC", "NB", "GC")
	seulist <- lapply(allsps, function(sp) {
					xx <- readRDS(paste0(inputdir, sp, "_DG_seu.rds")) %>%
							subset(., fig2cluster %in% allcls)
					return(xx)
					}) %>%
					setNames(., allsps)
	share_genes <- lapply(seulist, rownames) %>%
					Reduce("intersect", .)


	group.by <- "fig2cluster"
	allcls <- levels(as.factor(seulist$Mouse@meta.data[, group.by]))
	marres <- lapply(c("Mouse", "Pig", "Rhesus"), function(sp) {
		subseu <- seulist[[sp]][share_genes, ]
		Idents(subseu) <- group.by
		mar <- lapply(allcls, function(cls) {
			min_thre <- switch(sp, Mouse = 0.05, Pig = 0.05, Rhesus = 0.05)
			if (sum(subseu@meta.data[, group.by] == cls) < 10){
				res <- NULL
			} else {
				res <- FindMarkers(subseu, ident.1 = cls, ident.2 = NULL, logfc.threshold = 0.1, min.pct = min_thre, only.pos = TRUE, max.cells.per.ident = 1000) %>%
					rownames_to_column("gene") %>%
					mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
					mutate(species = sp)
			}
			return(res)
			}) %>% setNames(., allcls)
		return(mar)
		}) %>% setNames(., c("Mouse", "Pig", "Rhesus"))
	saveRDS(marres, file = marFile)
}



##------------------------------------------------------------------------
## Get the average expression for downstream visualization
avgFile <- paste0(inputdir, "AVG.HRPM.GC.Rdata")
if (!file.exists(avgFile)){
	group.by <- "fig2cluster"
	allcls <- c("Astro", "RGL", "nIPC", "NB", "GC")
	allsps <- c("Mouse", "Pig", "Rhesus", "Human")
	seulist <- lapply(allsps, function(sp) readRDS(paste0(inputdir, sp, "_DG_seu.rds"))) %>%
					setNames(., allsps)
	share_genes <- lapply(seulist, rownames) %>%
					Reduce("intersect", .)


	res <- lapply(c("Mouse", "Pig", "Rhesus", "Human"), function(sp) {
				subseu <- seulist[[sp]]
				cls_use <- allcls[allcls %in% unique(subseu@meta.data[, group.by])]
				Idents(subseu) <- group.by
				avg <- log(AverageExpression(subseu)$RNA + 1)[, cls_use]


				exp_ratio <- lapply(cls_use, function(cls) 
					Matrix::rowMeans(subseu$RNA@data[, subseu@meta.data[, group.by] == cls] != 0)) %>%
								setNames(., cls_use) %>%
								as.data.frame()


				## Remove the cluster with limited number of cells
				rm_cls <- table(subseu@meta.data[, group.by]) %>%
							.[. < 20] %>%
							names()
				print(rm_cls)
				avg <- avg[share_genes, setdiff(colnames(avg), rm_cls)]
				exp_ratio <- exp_ratio[share_genes, setdiff(colnames(exp_ratio), rm_cls)]

				## Make column names unique
				colnames(avg) <- paste0(sp, "_", colnames(avg))
				colnames(exp_ratio) <- paste0(sp, "_", colnames(exp_ratio))
				return(list(avg = avg, ratio = exp_ratio))
				})
	avgs <- lapply(res, function(x) x$avg) %>% 
				do.call(cbind, .)
	ratios <- lapply(res, function(x) x$ratio) %>% 
				do.call(cbind, .)
	save(avgs, ratios, file = avgFile)
}
load(file = avgFile)



##------------------------------------------------------------------------
## Get the specific markers for DotPlot visualization
marres <- readRDS(file = paste0(inputdir, "Markers.RPM.res.rds"))
load(file = paste0(inputdir, "AVG.HRPM.GC.Rdata"))
sp_order <- c("Mouse", "Pig", "Rhesus")
exclu_markers <- lapply(sp_order, function(sp) {
	res <- marres[[sp]]
	exp_ratio <- ratios[grep(sp, colnames(ratios))]
	colnames(exp_ratio) <- extract_field(colnames(exp_ratio), 2, "_")

	spres <- lapply(c("nIPC", "NB"), function(cls) {
		if (is.null(res[[cls]])){
			gg <- c()
		} else {
			gg <- res[[cls]] %>%
					subset(ratio_fc >= 1.2) %>%
					mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
					subset(p_val_adj <= 0.05) %>%
					arrange(desc(avg_logFC)) %>%
					.$gene
			idx1 <- exp_ratio[gg, cls]/apply(exp_ratio[gg, setdiff(colnames(exp_ratio), cls), drop = FALSE], 1, max) >= 1.2
			idx2 <- apply(exp_ratio[gg, setdiff(colnames(exp_ratio), cls), drop = FALSE], 1, max) <= 0.125
			gg <- gg[idx1 & idx2]

			if (cls == "nIPC"){
				all_cc <- unlist(cc.genes) %>% unique()
				non_cc_gg <- gg[!gg %in% all_cc] %>% .[1:10]
				cc_gg <- gg[gg %in% all_cc] %>% .[1:10]
				gg <- c(non_cc_gg, cc_gg)
			} else {
				gg <- gg[1:min(length(gg), 20)]
			}
		}
		return(gg)
		}) %>% setNames(., c("nIPC", "NB"))
	return(spres)
	}) %>% setNames(., sp_order)
saveRDS(exclu_markers, file = paste0(inputdir, "Markers.RPM.exclu.rds"))



##----------------------------------------------------------------------------------
## Get the conserved markers (>= 2 species)
marres <- readRDS(file = paste0(inputdir, "Markers.RPM.res.rds"))
load(file = paste0(inputdir, "AVG.HRPM.GC.Rdata"))
genelist <- lapply(names(marres), function(sp) {
	res <- marres[[sp]]
	exp_ratio <- ratios[grep(sp, colnames(ratios))]
	colnames(exp_ratio) <- extract_field(colnames(exp_ratio), 2, "_")
	cls_use <- c("nIPC", "NB")
	spres <- lapply(cls_use, function(cls) {
		if (is.null(res[[cls]])){
			new_gg <- c()
		} else {
			gg <- res[[cls]] %>%
					subset(ratio_fc >= 1.2) %>%
					mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
					subset(p_val_adj <= 0.05) %>%
					arrange(desc(avg_logFC)) %>%
					.$gene
			idx1 <- (exp_ratio[gg, cls]/apply(exp_ratio[gg, setdiff(colnames(exp_ratio), cls_use), drop = FALSE], 1, max)) >= 1.25
			##idx1 <- (exp_ratio[gg, cls]/apply(exp_ratio[gg, ], 1, function(mm) sort(mm, decreasing = TRUE)[3])) >= 2
			##idx1[is.na(idx1)] <- FALSE
			##idx2 <- apply(exp_ratio[gg, setdiff(colnames(exp_ratio), cls), drop = FALSE], 1, max) <= 0.4
			gg <- gg[idx1]
			new_gg <- setNames(1:length(gg), gg)
		}
		return(new_gg)
		}) %>% setNames(., cls_use)
	return(spres)
	}) %>% setNames(., names(marres)) %>%
		do.call(c, .)

get_con <- function(x, y, z = NULL, n = 20){
	if (is.null(z)){
		aa <- intersect(names(x), names(y))
		bb <- (x[aa] + y[aa])
		cc <- sort(bb, decreasing = FALSE)
	} else {
		aa <- intersect(names(x), names(y)) %>% intersect(., names(z))
		bb <- sapply(1:length(aa), function(pp) c(x[aa[1]], y[aa[1]], z[aa[1]]) %>% sort(., decreasing = FALSE) %>% .[1:2] %>% sum()) %>% 
				setNames(., aa)
		cc <- sort(bb, decreasing = FALSE)
	}
	dd <- cc[1:min(n, length(cc))]
	return(names(dd))
}


con_list <- list(nIPC_RM = get_con(genelist$Rhesus.nIPC, genelist$Mouse.nIPC), 
				NB_RMP = get_con(genelist$Rhesus.NB, genelist$Mouse.NB, genelist$Pig.NB), 
				NB_RM = get_con(genelist$Rhesus.NB, genelist$Mouse.NB), 
				NB_RP = get_con(genelist$Rhesus.NB, genelist$Pig.NB), 
				NB_PM = get_con(genelist$Mouse.NB, genelist$Pig.NB))
con_list$NB_RM <- setdiff(con_list$NB_RM, con_list$NB_RMP)
con_list$NB_RP <- setdiff(con_list$NB_RP, con_list$NB_RMP)
con_list$NB_PM <- setdiff(con_list$NB_PM, con_list$NB_RMP)
con_list$common <- c("SLC17A7", "GAD1", "OLIG2", "PDGFRA", "MOBP", "SLC1A3", "AQP4", "PTPRC", "CLDN5", "CEMIP", "PROX1")
con_list <- con_list[c("common", "nIPC_RM", "NB_RMP", "NB_PM", "NB_RP", "NB_RM")]
save(genelist, con_list, file = paste0(inputdir, "Markers.RPM.conserved.Rdata"))


























