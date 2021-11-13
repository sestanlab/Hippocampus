## Compare the mouse and Rhesus markers
source("../scripts/hip.fun.R") 
library(ggpubr)



##------------------------------------------------------------------------
## Get rhesus and mouse marker sets
## Get the conserved markers (>= 2 species)
marres <- readRDS(file = paste0(inputdir, "Markers.RPM.res.rds"))
load(file = paste0(inputdir, "AVG.HRPM.GC.Rdata"))
genelist <- lapply(names(marres), function(sp) {
	res <- marres[[sp]]
	exp_ratio <- ratios[grep(sp, colnames(ratios))]
	colnames(exp_ratio) <- extract_field(colnames(exp_ratio), 2, "_")
	spres <- lapply(names(res), function(cls) {
		if (is.null(res[[cls]])){
			new_gg <- c()
		} else {
			gg <- res[[cls]] %>%
					subset(ratio_fc >= 1.2) %>%
					mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
					subset(p_val_adj <= 0.05) %>%
					arrange(desc(avg_logFC)) %>%
					.$gene
			##idx1 <- (exp_ratio[gg, cls]/apply(exp_ratio[gg, setdiff(colnames(exp_ratio), cls), drop = FALSE], 1, max)) >= 1.25
			idx1 <- (exp_ratio[gg, cls]/apply(exp_ratio[gg, ], 1, function(mm) sort(mm, decreasing = TRUE)[3])) >= 2
			##idx1[is.na(idx1)] <- FALSE
			##idx2 <- apply(exp_ratio[gg, setdiff(colnames(exp_ratio), cls), drop = FALSE], 1, max) <= 0.4
			gg <- gg[idx1]
			new_gg <- setNames(1:length(gg), gg)
		}
		return(new_gg)
		}) %>% setNames(., names(res))
	return(spres)
	}) %>% setNames(., names(marres)) %>%
		do.call(c, .)



load(file = paste0(inputdir, "AVG.HRPM.GC.Rdata"))
sp_order <- c("Mouse", "Pig", "Rhesus")
avg_use <- avgs[, lapply(c(sp_order, "Human"), function(x) {
			xx <- grep(x, colnames(avgs), value = TRUE)
			bg <- paste0(x, "_", c("Astro", "RGL", "nIPC", "NB", "GC"))
			yy <- bg[bg %in% xx]
			yy
			}) %>% unlist()]


list_use <- lapply(genelist, function(x) names(x)[1:min(length(x), 75)])
names(list_use) <- gsub("\\.", "_", names(list_use))
exp_list <- lapply(sp_order, function(sp) {
	data <- avg_use[unlist(list_use[grep(sp, colnames(avg_use), value = TRUE)]) %>% unique(), ]
	rownames(data) <- paste0(sp, "_", rownames(data))
	data <- t(data) %>% scale() %>% t() %>%
				MinMax(., min = -2, max = 3)
	return(data)
	})



pdf(paste0(outputdir, "SF2.marker.shared.heatmap.pdf"), width = 4, height = 3)
pheatmap::pheatmap(do.call(rbind, exp_list), cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(viridis(5))(30), border_color = NA, show_rownames = FALSE, show_colnames = TRUE, fontsize_col = 8, gaps_row = cumsum(sapply(exp_list, nrow))[1:2], gaps_col = c(5, 8, 12), labels_col = extract_field(colnames(exp_list[[1]]), -1, "_"), annotation_colors = list())
dev.off()



## Check the marker similarity
## For each cell type in species A, find the top 50 markers and check whether these 50 markers also serve as markers in species B
pairs  <- c("Mouse|Rhesus", "Rhesus|Mouse", "Mouse|Pig", "Pig|Mouse", "Rhesus|Pig", "Pig|Rhesus")
sp_use <- c("Pig", "Mouse", "Rhesus")
marker_simi <- lapply(pairs, function(pair) {
	sp1 <- strsplit(pair, "|", fixed = TRUE)[[1]][1]
	sp2 <- strsplit(pair, "|", fixed = TRUE)[[1]][2]
	spres <- sapply(c("Astro", "nIPC", "NB", "GC"), function(cls) {
		if ((!is.null(marres[[sp1]][[cls]]))  && (!is.null(marres[[sp2]][[cls]]))){
			tops <- marres[[sp1]][[cls]] %>%
						subset(ratio_fc >= 1.2 & avg_logFC >= 0.5) %>%
						mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
						subset(p_val_adj <= 0.05) %>%
						arrange(desc(avg_logFC)) %>%
						top_n(n = 100, wt = avg_logFC) %>%
						.$gene
			othermarkers <- marres[[sp2]][[cls]] %>%
						subset(ratio_fc >= 1.2 & avg_logFC >= 0.5) %>%
						mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
						subset(p_val_adj <= 0.01) %>%
						.$gene
			simi <- sum(tops %in% othermarkers)/length(tops)
		}else {
			simi <- NA
		}
		
		return(simi)
		})
	return(spres)
	}) %>% 
			setNames(., pairs) %>%
			as.data.frame()


simi_data <- data.frame(MR = (marker_simi$Mouse.Rhesus + marker_simi$Rhesus.Mouse)/2, 
					MP = (marker_simi$Mouse.Pig + marker_simi$Pig.Mouse)/2, 
					RP = (marker_simi$Rhesus.Pig + marker_simi$Pig.Rhesus)/2, 
					row.names = rownames(marker_simi),
					cluster = rownames(marker_simi),
					stringsAsFactors = FALSE) %>%
					tidyr::gather(., "type", "similarity", c("MR", "MP", "RP")) %>%
					subset(!is.na(similarity)) %>%
					mutate(type = factor(as.character(type), levels = rev(c("MR", "MP", "RP"))))



library(lemon)
cls_cols <- c("#fcb5c5","#e01134", "#d834b7", "#7dc6e8", "#08519c") %>% 
				setNames(., c("Astro", "RGL", "nIPC", "NB", "GC"))
ss <- ggplot(simi_data, aes(y = similarity, x = cluster, fill = type)) + 
		geom_bar(stat = "identity", position = position_dodge2(preserve = "total"), color = "black") +
		theme_classic() +
		scale_x_discrete(limits = c("Astro", "nIPC", "NB", "GC")) +
		coord_capped_cart(left='both', bottom = "both") +
		theme(legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2))
pdf(paste0(outputdir, "SF2.marker.shared.similarity.pdf"), width = 5, height = 3)
print(ss)
dev.off()







