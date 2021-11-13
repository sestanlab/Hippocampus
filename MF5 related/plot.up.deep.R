## Plot the E of 4 region markers
source("../scripts/hip.fun.R")
library(tidyr)
library(ggpubr)

## First identify the exclusive IT markers & non-IT markers
all.mar <- readRDS(file = paste0("../MF3_cellregion/load_files/", "Markers.4reg.ExN.rds"))
cls.order <- list(HIP = c("DG GC PROX1 PDLIM5", "DG GC PROX1 SGCZ", "DG MC ARHGAP24 DLC1", "CA2 CFAP299 HGF", "CA3 CFAP299 SYN3", "CA1 dorsal GRIK1 GRM3", "CA1 ventral ACVR1C SYT13", "SUB proximal ROBO1 COL5A2", "SUB distal FN1 NTNG1", "SUB proximal ROBO1 SEMA3E"), 
					EC = c("EC L2 CUX2 PDGFD","EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L2 CUX2 LAMA3","EC L2 RELN BMPR1B","EC L2 RELN BCL11B","EC L3 PCP4 ADARB2","EC L5 RORB TPBG","EC L5 RORB TLL1","EC L6 THEMIS CDH13","EC L6 THEMIS RGS12","EC L5 BCL11B ADRA1A","EC L6 TLE4 SULF1","EC L56 TLE4 NXPH2","EC L6b TLE4 CCN2"), 
					MTG = c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1"), 
					PFC = readRDS(file = "~/project/PFC/data/PFC.cluster.order.rds")[c(1:14, 16:22, 15, 23:40)][c(1:20, 22:26, 21, 31:34, 27:30, 35:40)])

max.marker <- 750
reg <- "MTG"
mars <- lapply(cls.order$MTG, function(x) {
	xx <- all.mar %>%
				subset(region == reg & cluster == x) %>%
				subset(pct.1 >= 0.25 & avg_logFC >= log(3) & ratio_fc >= 1.25) %>%
				mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
				subset(p_val_adj <= 0.01)
	yy <- setNames(xx$avg_logFC, xx$gene)
	gg <- names(yy)[order(yy, decreasing = TRUE)[1:min(length(yy), max.marker)]]
	return(gg)
	}) %>% setNames(., cls.order$MTG)


raw_it <- mars[grepl('L2 |L2-|L3-',cls.order$MTG, invert = TRUE)] %>% unlist() %>% table()


exclu_it <- setdiff(raw_it %>% .[.>=2] %>% names(), names(raw_nit))
exclu_nit <- setdiff(raw_nit %>% .[.>=2] %>% names(), names(raw_it))


## Check the percentage of cells expressing 
avglist <- readRDS(file = paste0(inputdir, "AVG.ExN.4reg.rds"))


qt <- 0.4
avgf <- lapply(names(avglist), function(reg) {
	xx <- avglist[[reg]][, cls.order[[reg]]]
	thre <- quantile(as.matrix(xx), probs = qt)
	xx[xx < thre] <- 0
	return(xx)
	}) %>% 
			setNames(., names(avglist))


plot_data <- lapply(names(avgf), function(reg) {
	xx <- avgf[[reg]]
	data <- data.frame(cluster = cls.order[[reg]], stringsAsFactors= FALSE) %>%
				mutate(up = colMeans(avgf[[reg]][exclu_it, ] !=0 )) %>%
				mutate(deep = colMeans(avgf[[reg]][exclu_nit, ] !=0 )) %>%
				mutate(region = reg) %>%
				gather(., "type", "ratio", c("up", "deep"))

	return(data)
	}) %>%
		do.call(rbind, .)

maxy <- max(plot_data$ratio)##ceiling(max(plot_data$ratio)*10)/10
plist <- lapply(names(avgf), function(reg) {
	p <- ggplot(subset(plot_data, region == reg), aes_string(x = "cluster", y = "ratio", fill = "type")) +
			geom_bar(position = "dodge2", stat = "identity") + 
			theme_classic() + 
			scale_x_discrete(limits = cls.order[[reg]]) + 
			RotatedAxis() +
			scale_fill_manual(values = c(up = "#bababa", deep = "#4d4d4d")) +
			scale_y_continuous(limits = c(0, maxy)) +
			theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.line.y = element_line(size =  0.2), axis.ticks.y = element_line(size = 0.2))
	if (reg != "HIP"){
		p <- p +
			theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank())
	}
	return(p)
	})


ggarrange(plotlist = plist, nrow = 1, ncol = 4, align = "h", widths = c(1.5, 1.5, 1.8, 2)) %>%
	ggexport(filename = paste0(outputdir, "SF4.up.deep.exp.ratio.", as.character(qt), ".pdf"), width = 16, height = 4)










