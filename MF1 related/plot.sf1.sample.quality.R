## Plot the number of cells in each sample
source("../scripts/hip.fun.R")
library(lemon)

hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds")) 


meta_use <- hipec@meta.data
rep.by <- "repname"
ind.by <- "samplename"
reg.by <- "region_new"


###----------------------------------------------------------------------------------------
## Plot number of cells by samples
plot_data <- meta_use %>%
				rownames_to_column("cell") %>%
				group_by(!!sym(rep.by), !!sym(ind.by)) %>%
				summarize(size = n()) %>%
				ungroup() %>% 
				mutate(!!rep.by := extract_field(!!sym(rep.by), -1, "_")) %>%
				mutate(!!rep.by := as.numeric(as.factor(!!sym(rep.by)))) %>%
				mutate(!!rep.by := paste0("REP-", !!sym(rep.by))) 


sample_order <- c("HSB179","HSB181","HSB282","HSB231","HSB237","HSB628")
ind.cells <- plot_data %>% 
				group_by(!!sym(ind.by)) %>%
				summarize(size = sum(size)) 


mm <- ggplot(plot_data, aes_string(x = ind.by, y = "size")) +
		geom_bar(aes_string(x = ind.by, y = "size"), fill = "grey", color = "black", position = position_stack(reverse = FALSE), stat = "identity", lwd = 0.3) +
		geom_text(data = ind.cells, mapping = aes_string(x = ind.by, y = "size", label = "size"), nudge_y = 1500, angle = 45, vjust = 0.5, hjust = 0.1, size = 2.8) +
		coord_capped_cart(left='both') +
		scale_x_discrete(limits = sample_order) +  	
		theme_classic() + 
		RotatedAxis() + 
		labs(y = "Sample size", x = "Individual") +
		theme(axis.line=element_line(size = 0.25), axis.ticks=element_line(size = 0.25)) +
		theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),  axis.title.x= element_blank())



plot_data <- meta_use %>%
				rownames_to_column("cell") %>%
				group_by(!!sym(reg.by), !!sym(ind.by)) %>%
				summarize(size = n()) %>%
				ungroup()


sample_order <- c("HSB179","HSB181","HSB282","HSB231","HSB237","HSB628")
ind.cells <- plot_data %>% 
				group_by(!!sym(ind.by)) %>%
				summarize(size = sum(size)) 
reg_colors <- paste0(c("#b2182b", "#b8e186", "#4d9221", "#80ffff", "#fee090"), "") %>%
				setNames(., c("DG", "CA24", "CA1", "SUB", "EC"))

nn <- ggplot(plot_data, aes_string(x = ind.by, y = "size")) +
		geom_bar(aes_string(x = ind.by, y = "size", fill = reg.by), color = "black", position = position_stack(reverse = FALSE), stat = "identity", lwd = 0.3) +
		geom_text(data = ind.cells, mapping = aes_string(x = ind.by, y = "size", label = "size"), nudge_y = 1500, angle = 45, vjust = 0.5, hjust = 0.1, size = 2.8) +
		coord_capped_cart(bottom='both', left='both') +
		scale_x_discrete(limits = sample_order) +  
		scale_fill_manual(values = reg_colors)+
		theme_classic() + 
		RotatedAxis() + 
		labs(y = "Sample size", x = "Individual") +
		theme(axis.line=element_line(size = 0.25), axis.ticks=element_line(size = 0.25), legend.position = "bottom")

pdf(paste0(outputdir, "SF1.Sample_size.pdf"), width = 7, height = 7)
plot_grid(mm,nn,nrow = 2, ncol = 1, align = 'v', rel_heights = c(1, 1.3))
dev.off()



###----------------------------------------------------------------------------------------
## Plot quality by samples
plot_cols <- c("nFeature_RNA", "nCount_RNA")
plot_data <- meta_use[, c(plot_cols, ind.by)] %>%
				reshape2::melt(id.vars = ind.by, measure.vars = plot_cols, variable.name = "quality")


plot_data$quality <- plot_data$quality %>%
						gsub("nCount_RNA", "nUMIs", .) %>%
						gsub("nFeature_RNA", "nGenes", .)
plot_data$quality <- factor(plot_data$quality, levels = c("nUMIs", "nGenes"))
##plot_data <- plot_data[plot_data$individual != "CJB1680", ]


qqq <- ggplot(plot_data, aes_string(x = ind.by, y = "value")) +
		geom_violin(scale = "width", fill = "grey", size = 0.1,adjust = 2,trim =TRUE) + 
		geom_boxplot(fill = "white", color = "black", width=0.15, outlier.shape = NA, lwd= 0.2, alpha = 0.5) + 
		coord_capped_cart(top='both', left='both') +
		theme_classic() + 
		scale_x_discrete(limits = sample_order) +
		##RotatedAxis() + 
		theme(panel.border=element_blank(), axis.line=element_line(size = 0.25), axis.ticks = element_line(size = 0.25),axis.text.x = element_text(size = 10, hjust = 1, angle = 45), axis.text.y = element_text(size = 9), axis.title = element_blank()) + 
		facet_wrap(facets = vars(quality), nrow = length(plot_cols), ncol = 1, strip.position = "left", scales = "free_y") + 
		theme(strip.background = element_blank(), strip.text = element_text(size = 8, face = "bold"), panel.spacing = unit(0.05, "in"), legend.title = element_text(size = 8), legend.text = element_text(size = 12), strip.placement = "outside")
pdf(paste0(outputdir, "SF1.Sample_quality.pdf"), width = 7, height = 7)#, units = "in",  res = 300)
print(qqq)
dev.off()






