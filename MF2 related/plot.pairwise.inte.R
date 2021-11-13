## Organize all the pair-wise integration analysis
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")


##------------------------------------------------------------------------------
## Load the pair-wise integration dataset
icells <- c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG") ##The first is NB, second is nIPC
inte.met <- c("seurat", "harmony")[as.numeric(args[1])]
seulist <- list(`Human.Mouse young` = readRDS(paste0(inputdir, "HM.pw.GC.", inte.met,".expHVG.1500.slim.rds")),
			`Human.Mouse middle` = readRDS(paste0(inputdir, "HRPM.pw.GC.", inte.met,".expHVG.1500.Human.Mouse.slim.rds")), 
			`Human.Pig` = readRDS(paste0(inputdir, "HRPM.pw.GC.", inte.met,".expHVG.1500.Human.Pig.slim.rds")), 
			`Human.Rhesus` = readRDS(paste0(inputdir, "HRPM.pw.GC.", inte.met,".expHVG.1500.Human.Rhesus.slim.rds")))


## Load the cell identity
data_list <- lapply(names(seulist), function(pair) {
	sp1 <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
	sp2 <- strsplit(pair, ".", fixed = TRUE)[[1]][2]

	pdata <- seulist[[pair]]@meta.data %>%
				rownames_to_column("cell") %>%
				mutate(cluster = fig2cluster) %>%
				select(cell, vis.batch, cluster) %>%
				subset(vis.batch %in% c(sp1, sp2)) %>%
				mutate(xaxis = seulist[[pair]]$umap@cell.embeddings[, 1]) %>%
				mutate(yaxis = seulist[[pair]]$umap@cell.embeddings[, 2]) %>%
				#mutate(dsize = ifelse(cell %in% icells, 10, 0.01))
				mutate(dsize = ifelse(cell %in% icells, 0.01, 0.01))
				
	return(pdata)
	}) %>%
		setNames(., names(seulist))



## Set the colors
cls_colors <- c("#fcb5c5","#e01134", "#d834b7", "#7dc6e8", "#08519c") %>% 
			setNames(., c("Astro","RGL","nIPC","NB","GC"))
sp_cols <- setNames(c("#ffa600", "#031a7f", "#031a7f", "#89DA59", "#FF420E"), c("Pig", "Mouse young", "Mouse middle", "Rhesus", "Human"))


## Plot the UMAP (colored by species & clusters)
plist <- lapply(names(data_list), function(pair) {
	bb <- 0
	xrange <- c(min(data_list[[pair]]$xaxis)- bb, 
			max(data_list[[pair]]$xaxis) + bb)
	yrange <- c(min(data_list[[pair]]$yaxis)- bb, 
			max(data_list[[pair]]$yaxis) + bb)

	sp1 <- "Human"
	sp2 <- strsplit(pair, ".", fixed = TRUE)[[1]] %>%
			setdiff(., sp1)

	p1 <- ggplot(data_list[[pair]], aes(x = xaxis, y = yaxis, color = vis.batch)) +
      		##geom_point(size = 0.01) +
      		geom_point(aes(size = dsize)) +
      		scale_color_manual(values = sp_cols) +
        	theme_classic() + 
        	scale_size_identity() +
        	theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none", plot.title = element_blank())  +
        	xlim(limits = xrange) +
			ylim(limits = yrange)

	p2 <- ggplot(subset(data_list[[pair]], vis.batch == sp1), aes(x = xaxis, y = yaxis, color = cluster)) +
      		##geom_point(size = 0.01) +
      		geom_point(aes(size = dsize)) +
      		scale_color_manual(values = cls_colors) +
        	theme_classic() + 
        	scale_size_identity() +
        	theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none", plot.title = element_blank())  +
        	xlim(limits = xrange) +
			ylim(limits = yrange)

	p3 <- ggplot(subset(data_list[[pair]], vis.batch == sp2), aes(x = xaxis, y = yaxis, color = cluster)) +
      		##geom_point(size = 0.01) +
      		geom_point(aes(size = dsize)) +
      		scale_color_manual(values = cls_colors) +
        	theme_classic() + 
        	scale_size_identity() +
        	theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none", plot.title = element_blank())  +
        	xlim(limits = xrange) +
			ylim(limits = yrange)
	newp <- list(p1, p2, p3) %>%
				lapply(., function(x) x + theme(plot.margin = unit(c(0, 0, 0, 0), "inch")))
	return(newp)
	})


jpeg(paste0(outputdir, "SF2.pairwise.GC.1.species.", inte.met, ".jpeg"), height = 4 * 3, width = 4 * 4, units = "in", res = 300)
patchwork::wrap_plots(do.call(c, plist), nrow = 3, ncol = 4, byrow = FALSE)
dev.off()


##------------------------------------------------------------------------------
## Identify the matched human cells 
## Only run once
if (FALSE){
	inte.met <- c("seurat", "harmony")[as.numeric(args[1])]
	seulist <- list(`Human.Mouse young` = readRDS(paste0(inputdir, "HM.pw.GC.", inte.met,".expHVG.1500.slim.rds")),
			`Human.Mouse middle` = readRDS(paste0(inputdir, "HRPM.pw.GC.", inte.met,".expHVG.1500.Human.Mouse.slim.rds")), 
			`Human.Pig` = readRDS(paste0(inputdir, "HRPM.pw.GC.", inte.met,".expHVG.1500.Human.Pig.slim.rds")), 
			`Human.Rhesus` = readRDS(paste0(inputdir, "HRPM.pw.GC.", inte.met,".expHVG.1500.Human.Rhesus.slim.rds")))

	pt_vec <- list()

	for (pair in names(seulist)){
		##pair <- names(seulist)[3]
		## The order: x-min, x-max, y-min, y-max
		sq_list1 <- list(`Human.Mouse young` = c(2.5, 15, 5, 11),
						`Human.Mouse middle` = c(0, 10, 6, 15),
						`Human.Pig` = c(-1, 9, 6, 15), 
						`Human.Rhesus` = c(1.75, 7, 8.5, 11.5))
		sq_list2 <- list(`Human.Mouse young` = c(-1.5, 11, 3.8, 15),
						`Human.Mouse middle` = c(0, 10, 6, 15),
						`Human.Pig` = c(-0.8, 9, 4.3, 15), 
						`Human.Rhesus` = c(1, 5, -1, 1.5))
		sq <- switch(inte.met, seurat = sq_list1[[pair]], harmony = sq_list2[[pair]])
		p <- DimFig(seulist[[pair]], group.by = "mres", file_name = paste0("final.HRM.inte.GC.ident.", "human_match"), output.ggplot = TRUE, pt.size = 0.1, cells.highlight = c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG"), sizes.highlight = 4)[[3]] +
					geom_vline(xintercept = sq[1:2]) +
					geom_hline(yintercept = sq[3:4])
		jpeg(paste0(outputdir, "Intemediate.Putative.", inte.met, ".", pair, ".jpeg"), width = 4, height = 4, units = "in", res = 300)
		print(p)
		dev.off()
		pt_vec[[pair]] <- colnames(seulist[[pair]])[seulist[[pair]]$umap@cell.embeddings[, 1] >= sq[1] & 
											seulist[[pair]]$umap@cell.embeddings[, 1] <= sq[2] &
											seulist[[pair]]$umap@cell.embeddings[, 2] >= sq[3] & 
											seulist[[pair]]$umap@cell.embeddings[, 2] <= sq[4] & 
											seulist[[pair]]@meta.data$species == "Human"]
		print(pt_vec[[pair]])
	}
	saveRDS(pt_vec, file = paste0(inputdir, "Putative.", inte.met, ".cells.rds"))
}


## Additionally, there is one cell when doing the four-species integration
## Only run once
if (FALSE){
	hrpm <- readRDS(file = paste0(inputdir, "HRPM.GC.seurat.expHVG.1500.slim.rds"))

	sq <- c(0, 10, 6, 13)
	p <- DimFig(hrpm, group.by = "mres", file_name = paste0("final.HRM.inte.GC.ident.", "human_match"), output.ggplot = TRUE, pt.size = 0.1, cells.highlight = c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG"), sizes.highlight = 4)[[3]] +
				geom_vline(xintercept = sq[1:2]) +
				geom_hline(yintercept = sq[3:4])
	jpeg(paste0(outputdir, "Intemediate.Putative.seurat.HRPM.jpeg"), width = 4, height = 4, units = "in", res = 300)
	print(p)
	dev.off()
	pt <- colnames(hrpm)[hrpm$umap@cell.embeddings[, 1] >= sq[1] & 
										hrpm$umap@cell.embeddings[, 1] <= sq[2] &
										hrpm$umap@cell.embeddings[, 2] >= sq[3] & 
										hrpm$umap@cell.embeddings[, 2] <= sq[4] & 
										hrpm@meta.data$species == "Human"]
	saveRDS(pt, file = paste0(inputdir, "Putative.seurat.HRPM.cells.rds"))
}


## Combine all the putaive cells together
## Only run once
if (FALSE) {
	all_pt <- lapply(c("seurat", "harmony"), function(x) readRDS(file = paste0(inputdir, "Putative.", x, ".cells.rds"))) %>%
				do.call(c, .) %>%
				unlist() %>%
				c(., readRDS(file = paste0(inputdir, "Putative.seurat.HRPM.cells.rds"))) %>%
				unique()
	saveRDS(all_pt, file = paste0(inputdir, "Putative.cbn.cells.rds"))
}












