## Plot the astrocytes laminar distribution across regions
source("../scripts/hip.fun.R")
library(ggpubr)

nnc_4reg <- readRDS(file = paste0(dataDir, "FourReg.NNC.final.01172021.rds"))
cls_list <- list(HIP = c("Astro AQP4 GFAP", "Astro AQP4 CHRDL1"), 
					EC  = c("Astro AQP4 GFAP", "Astro AQP4 CHRDL1"), 
					MTG = c("Astro L1-2 FGFR3 GFAP", "Astro L1-6 FGFR3 SLC14A1"),
					PFC = c("iAstro GFAP FABP7", "pAstro AQP4 SLC1A2"))
slim_list <- lapply(names(nnc_4reg), function(x) 
	nnc_4reg[[x]][, nnc_4reg[[x]]@meta.data$subtypes %in% cls_list[[x]]]
	) %>%
			setNames(., names(nnc_4reg))


newlist <- SplitObject(slim_list$HIP, split.by = "region_new") %>%
				c(., slim_list[c("EC", "MTG", "PFC")])




## Check the UMI numbers for each dataset
colors <- c("#FF420E","#FF420E","#FF420E","#FF420E","#FF420E","#89DA59","#4CB5F5","#FFBB00", "grey90") %>% 
				setNames(., c("HIP", "CA1", "CA24", "DG", "SUB", "EC", "PFC", "MTG", "bg"))
layer_genes <- c("AQP4", "GFAP", "ID3", "WDR49", "CHRDL1", "GRM3")
plist <- lapply(c("DG", "CA24", "CA1", "SUB", "EC", "MTG", "PFC"), function(reg) {
	cls_order <- ifelse_check(reg %in% c("HIP", "CA1", "CA24", "DG", "SUB"), cls_list$HIP, cls_list[[reg]])
	p <- DotPlot(newlist[[reg]], features = layer_genes, group.by = "subtypes", cols = c("lightgrey", colors[[reg]]), dot.min = 0.05, dot.scale = 3) +
			coord_flip() + 
			theme_cowplot() +
			scale_y_discrete(limits = cls_order) +
			scale_radius(limits = c(5, 100), range = c(0, 3)) +
			RotatedAxis() +
			theme(legend.position = "none", axis.line = element_line(size = 0.2), axis.ticks.y = element_line(size = 0.2), axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 6))
	if (reg != "DG"){
		p <- p + 
			theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
	}
	return(p)
	})
ggarrange(plotlist = plist, nrow = 1, ncol = length(plist), align = "h", common.legend = TRUE, legend = "top", widths = c(1.3, rep(1, length(plist)-1))) %>%
		ggexport(filename = paste0(outputdir, "SF4.NNC.layer.makers.pdf"), width = 6, height = 2.5)



##table(slim_list$HIP@meta.data$region_new, slim_list$HIP@meta.data$subtypes)
##Astro AQP4 CHRDL1 Astro AQP4 GFAP
##  CA1               1733            1712
##  CA24               573            2044
##  DG                5421           11787
##  SUB               1960             690
##  EC                2851             1023







