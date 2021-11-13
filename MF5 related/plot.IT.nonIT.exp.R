## Plot the E of 4 region markers
source("../scripts/hip.fun.R")
library(ggrepel)

## First identify the exclusive IT markers & non-IT markers
all.mar <- readRDS(file = paste0("../MF3_cellregion/load_files/", "Markers.4reg.ExN.rds"))
mtg_order <- c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1")
max.marker <- 750
reg <- "MTG"
mars <- lapply(cls.order$MTG, function(x) {
	xx <- all.mar %>%
				subset(region == reg & cluster == x) %>%
				subset(pct.1 >= 0.25 & avg_logFC >= log(2.5) & ratio_fc >= 1.2) %>%
				mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
				subset(p_val_adj <= 0.01)
	yy <- setNames(xx$avg_logFC, xx$gene)
	gg <- names(yy)[order(yy, decreasing = TRUE)[1:min(length(yy), max.marker)]]
	return(gg)
	}) %>% setNames(., cls.order$MTG)

raw_it <- mars[grepl('L2 |L2-|RORB',cls.order$MTG)] %>% unlist() %>% table()
raw_nit <- mars[grepl('FEZF2|SLC17A7',cls.order$MTG)] %>% unlist() %>% table() 


exclu_it <- setdiff(raw_it %>% .[.>=2] %>% names(), names(raw_nit))
exclu_nit <- setdiff(raw_nit %>% .[.>=2] %>% names(), names(raw_it))


## Then plot the scaled AVG of IT & non-IT markers for each region
avglist <- readRDS(file = paste0(inputdir, "AVG.ExN.4reg.rds"))
avglist <- lapply(names(avglist), function(x) {
	xx <- avglist[[x]] %>%
			t() %>% scale() %>% t() 
	xx[is.na(xx)] <- -2.5
	xx <- xx %>%
			MinMax(., min = -2.5, max = 2.5)
	return(xx)
	}) %>%
		setNames(., names(avglist))


plot_data <- lapply(names(avglist), function(reg) {
	data <- data.frame(cluster = colnames(avglist[[reg]]), stringsAsFactors = FALSE) %>%
					mutate(ITexp = colMeans(avglist[[reg]][exclu_it, colnames(avglist[[reg]])])) %>%
					mutate(NonITexp = colMeans(avglist[[reg]][exclu_nit, colnames(avglist[[reg]])])) %>%
					mutate(region = reg)
	return(data)
	})  %>%
		do.call(rbind, .)


## Plot the results
set.seed(42)
maxy <- max(c(max(plot_data$ITexp), max(plot_data$NonITexp)))
xrange <- c(0, max(plot_data$ITexp)+0.0)
yrange <- c(0, max(plot_data$NonITexp)+0.0)
plist <- lapply(names(avglist), function(reg) {
	text_size <- ifelse(reg %in% c("MTG", "PFC"), 2, 3)
	p <- ggplot(subset(plot_data, region == reg), aes_string(x = "ITexp", y = "NonITexp")) +
			geom_point(size = 3, color = "red") +
			geom_text_repel(aes_string(label = "cluster"), size = text_size) +
			theme_bw() + 
			scale_x_continuous(limits = c(0, maxy)) +
			scale_y_continuous(limits = c(0, maxy)) +
			geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.4) +
			labs(x = "AVG EXP of IT markers",  y = "AVG EXP of non-IT markers") +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	if (reg != "HIP"){
		p <- p +
			theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank())
	}
	return(p)
	})

pdf(paste0(outputdir, "SF4.IT.nonIT.AVG.exp.pdf"), width = 16, height = 4)
plot_grid(plotlist = plist, nrow = 1, ncol = 4, rel_widths = c(1.1, 1, 1, 1)) %>% print()
dev.off()




