## Plot the number of cells in each sample
source("../scripts/hip.fun.R")
library(lemon)

hrp_dg <- readRDS(file = paste0(inputdir, "zf.HRPM.all.seurat.20210407.rds")) %>%
			subset(species %in% c("Human", "Rhesus", "Pig"))


meta_use <- hrp_dg@meta.data
meta_use$samplename <- gsub("RMB([0-9]*)_.*", "RMB\\1", meta_use$samplename) %>%
						gsub("RMB([0-9]*)_.*", "RMB\\1", .)
rep.by <- "repname"
split.by <- "species"
ind.by <- "samplename"
sp_colors <- c(Human = "#FF420E", Rhesus = "#89DA59", Pig = "#ffa600")


###----------------------------------------------------------------------------------------
## Plot number of cells by samples
plot_data <- meta_use %>%
				rownames_to_column("cell") %>%
				group_by(!!sym(rep.by), !!sym(ind.by)) %>%
				summarize(size = n(), !!split.by := unique(!!sym(split.by))) %>%
				ungroup() %>% 
				mutate(!!split.by := factor(!!sym(split.by), levels = names(sp_colors))) %>%
				mutate(!!rep.by := extract_field(!!sym(rep.by), -1, "_")) %>%
				mutate(!!rep.by := as.numeric(as.factor(!!sym(rep.by)))) %>%
				mutate(!!rep.by := paste0("REP-", !!sym(rep.by))) 


all_samples <- plot_data[, ind.by][[1]] %>% unique()
sample_order <- c("HSB179","HSB181","HSB282","HSB231","HSB237","HSB628","RMB1","RMB2","RMB3","p41_h0","p42_h0","p43_h0")


sp_ncells <- plot_data %>% 
				group_by(!!sym(split.by)) %>%
				summarize(size = sum(size)) %>%
				mutate(y_loc = 52000) %>% 
				mutate(!!split.by := factor(!!sym(split.by), levels = names(sp_colors))) %>%
				mutate(cell_origin = c("HSB282", "RMB2", "p42_h0"))


ind.cells <- plot_data %>% 
				group_by(!!sym(ind.by)) %>%
				summarize(size = sum(size), species = unique(species)) 


mm <- ggplot(plot_data, aes_string(x = ind.by, y = "size")) +
		geom_bar(aes_string(x = ind.by, y = "size", fill = split.by), color = "black", position = position_stack(reverse = FALSE), stat = "identity", lwd = 0.3) +
		geom_text(data = ind.cells, mapping = aes_string(x = ind.by, y = "size", label = "size"), nudge_y = 1500, angle = 45, vjust = 0.5, hjust = 0.1, size = 2.8) +
		coord_capped_cart(bottom='both', left='both') +
		scale_fill_manual(values = sp_colors) +
		scale_x_discrete(limits = sample_order) +  
		geom_label(data = sp_ncells, aes_string(x = "cell_origin", y = "y_loc", label = "size", fill = split.by), nudge_x = 0.5) + 		
		theme_classic() + 
		RotatedAxis() + 
		labs(y = "Sample size", x = "Individual") +
		theme(axis.line=element_line(size = 0.25), axis.ticks=element_line(size = 0.25))


pdf(paste0(outputdir, "SF2.Sample_size.pdf"), width = 7, height = 7)
print(mm)
dev.off()













