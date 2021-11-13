## Compare the DCX positive vs negative cells and check for the enrichment of rhesus/mouse NB markers.
source("../scripts/hip.fun.R")
library(ggpubr)



## Calculate the AUROC scores
auc_file <- paste0(inputdir, "DCX.auc.res.rds")
if (!file.exists(auc_file)){
	seulist <- readRDS(file = paste0(inputdir, "Subset.data.Cbn.GC.rds"))
	mars <- readRDS(paste0(inputdir, "Markers.RPM.exclu.rds"))  %>%
			do.call(c, .) %>%
			c(., list(proliferative = c("GFAP", "HOPX", "NES", "TFAP2C", "ASCL1", "EOMES", "NEUROG2"), neuroblast = c("DCX", "NEUROD1", "SOX11", "CALB2", "IGFBPL1", "FXYD7"))) %>%
			lapply(., function(x) setdiff(x, "DCX"))


	all_auc <- lapply(names(seulist)[1:4], function(sp) {
		auc <- GetModuleScore(assay.data = seulist[[sp]]$RNA@data, features = mars, seed = 42, method = "aucell", input_dir = inputdir, file_name = paste0("DCX.AUC.", sp), output_dir = outputdir)$auc[colnames(seulist[[sp]]), ]
		df <- as.data.frame(auc, check.names = FALSE) %>%
				rownames_to_column("cell") %>%
				mutate(DCX = as.numeric(seulist[[sp]]$RNA@data["DCX", ] > 0)) %>%
				mutate(species = sp) %>%
				mutate(cluster = seulist[[sp]]@meta.data$fig2cluster)
		return(df)
		}) %>%
			do.call(rbind, .)
	saveRDS(all_auc, file = auc_file)
}



all_auc <- readRDS(file = auc_file)
## Further scale the values
qt <- 0.99
all_fauc <- all_auc %>%
				mutate(Pig.NB.scaled = MinMax(Pig.NB/quantile(Pig.NB, qt), min = 0, max = 1)) %>%
				mutate(Mouse.NB.scaled = MinMax(Mouse.NB/quantile(Mouse.NB, qt), min = 0, max = 1)) %>%
				mutate(Rhesus.NB.scaled = MinMax(Rhesus.NB/quantile(Rhesus.NB, qt), min = 0, max = 1)) %>%
				mutate(neuroblast.scaled = MinMax(neuroblast/quantile(neuroblast, qt), min = 0, max = 1)) %>%
				mutate(DCX = as.character(DCX))



## Select the cell list
sp_order <- c("Mouse", "Pig", "Rhesus", "Human")
cls_list <- lapply(sp_order, function(sp) {
	ll <- list(c("NB"), c("NB", "GC")) %>%
			setNames(., paste0(sp, ".", c("NB", "NG")))
	ll
	}) %>%
		do.call(c, .)
cls_list["Human.NB"] <- NULL


## Re-organize the AUC results for plot
new_auc <- lapply(names(cls_list), function(x){
	sp <- extract_field(x, 1, ".")
	df <- all_fauc %>%
			subset(species == sp) %>%
			mutate(group = x)
	dcx_neg <- which(df$cluster == "GC" & df$DCX == "0")
	dcx_pos <- which(df$cluster %in% cls_list[[x]] & df$DCX != "0")
	df <- df[c(dcx_pos, dcx_neg), ]
	return(df)
	}) %>%
			do.call(rbind, .) %>%
			select(cell, DCX, species, cluster, group, Pig.NB.scaled, Mouse.NB.scaled, Rhesus.NB.scaled, neuroblast.scaled) %>%
			tidyr::gather(., "type", "AUC", c("Mouse.NB.scaled", "Pig.NB.scaled", "Rhesus.NB.scaled", "neuroblast.scaled")) %>%
			mutate(type = factor(as.character(type), levels = c("Mouse.NB.scaled", "Pig.NB.scaled", "Rhesus.NB.scaled", "neuroblast.scaled"))) %>%
			mutate(DCX = factor(as.character(DCX), levels = c("0", "1"))) %>%
			mutate(newcls = paste0(type, "|", group))


new_order <- paste0(rep(c("Mouse.NB.scaled", "Pig.NB.scaled", "Rhesus.NB.scaled", "neuroblast.scaled"), each = length(cls_list)), "|", rep(names(cls_list), 4))

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


p <- ggplot(new_auc, aes(x = newcls, fill = DCX, y = AUC)) +
		##geom_boxplot(position = "dodge2", color = "black", outlier.shape = NA,size = 0.1) +
		##geom_point(data = function(x) dplyr::filter_(x, ~ outlier), position = 'jitter') +
		##geom_violinhalf(position = "dodge2", trim = TRUE, scale = "width", size = 0.2) +
		geom_split_violin(scale = "width", size = 0.2, trim = TRUE, adjust = 1.5) +
		stat_compare_means(label = "p.signif", method = "wilcox.test", size = 1.5, method.args = list(alternative = "less"), symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", "ns"))) +
		scale_fill_manual(values = setNames(c("#ffffff","#666464"), c("0", "1"))) +
		theme_cowplot() + 
		scale_y_continuous(limits = c(0, 1)) +
		scale_x_discrete(limits = new_order, labels = setNames(rep(c(rep(c("group1", "group2"), 3), "group2"), 4), names(new_order))) +
		geom_vline(xintercept = c(7.5, 14.5, 21.5), linetype = "dashed", size = 0.2) +
		##facet_wrap(vars(type), nrow = 1, ncol = 4) +
		theme(axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1), axis.title = element_blank(), axis.text.y = element_text(size = 3),axis.text.x = element_text(angle = 30, hjust = 1, size = 3, color = rep(c("#031a7f", "#ffa600", "#89DA59", "#FF420E"), c(2, 2, 2, 1))), legend.position = "right")

pdf(paste0(outputdir, "SF2.DCX.pos.neg.comparison.v2.pdf"), width = 4, height = 1.2)
print(p)
dev.off() 




