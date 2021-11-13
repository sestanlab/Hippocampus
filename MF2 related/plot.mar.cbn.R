## Annotate the species markers
source("../scripts/hip.fun.R") 
library(ggpubr)


##------------------------------------------------------------------------
## Prepare the markers to plot
sp_order <- c("Mouse.nIPC", "Rhesus.nIPC", "proliferative", "Mouse.NB", "Pig.NB", "Rhesus.NB", "neuroblast")
mars <- readRDS(paste0(inputdir, "Markers.RPM.exclu.rds"))  %>%
			do.call(c, .) %>%
			c(., list(proliferative = c("GFAP", "HOPX", "NES", "TFAP2C", "ASCL1", "EOMES", "NEUROG2"), neuroblast = c("DCX", "NEUROD1", "SOX11", "CALB2", "IGFBPL1", "FXYD7"))) %>%
			c(list(common = c("SLC17A7", "GAD1", "OLIG2", "PDGFRA", "MOBP", "SLC1A3", "AQP4", "PTPRC", "CLDN5", "CEMIP", "PROX1", "DCX")), .)
mars <- mars[c("common", sp_order)]
all_genes <- unlist(mars) %>% unique()



##------------------------------------------------------------------------
## 1. marker annotation plot
load(file = paste0(inputdir, "Markers.RPM.conserved.Rdata")) ## genelist
genelist <- lapply(genelist, function(x) names(x))
gene2df <- lapply(sp_order, function(x) {
	aa <- as.numeric(all_genes %in% genelist[[x]]) ##genelist doesn't contain "proliferative" or "neuroblast",so it's 0 for these two list. And accordingly, only bb values. will be included.
	bb <- as.numeric(all_genes %in% mars[[x]]) 
	cc <- aa + bb
	if (!x %in% c("proliferative",  "neuroblast")){
		cc[bb == 1] <- 2 ## in the exclusive marker list
	}
	return(setNames(cc, all_genes))
	}) %>% 
			setNames(., sp_order) %>%
			as.data.frame(, check.names = FALSE) %>%
			as.matrix() %>%
			t() %>%
			.[setdiff(sp_order, "common"), ]


plot_data <- gene2df %>%
				as.data.frame(., check.names = FALSE) %>%
				rownames_to_column("type") %>%
				tidyr::gather(., "gene", "value", colnames(gene2df)) %>%
				mutate(gene = factor(as.character(gene), levels = colnames(gene2df))) %>%
				mutate(type = factor(as.character(type), levels = rev(rownames(gene2df)))) %>%
				mutate(value = as.character(value)) %>%
				mutate(isstroke = ifelse(value == 2, "1", "0"))

colors <- colorRampPalette(c("grey95", "black"))(30)

p.anno <- ggplot(plot_data) + 
				geom_tile(aes_string(x = "gene", y = "type", fill = "value"), color = "white", height = 0.9, width = 0.9, size = 0.1) + 
				theme_classic() +
				scale_x_discrete(limits = all_genes) +
				scale_fill_manual(values = setNames(colors[c(1, 10, 30)], c("0", "1", "2"))) +
				RotatedAxis() + 
	            theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 4.25), axis.line = element_blank(), legend.position = "none", axis.title = element_blank(), axis.ticks.x = element_blank())


##------------------------------------------------------------------------
## 2. marker Dotplot
##sublist <- readRDS(file = paste0(inputdir, "Subset.data.Cbn.GC.rds"))
sublist <- readRDS(file = paste0(inputdir, "HRPM.GC.seurat.expHVG.1500.slim.rds")) %>%
			SplitObject(., split.by = "species")

sp_cols <- setNames(c("#031a7f", "#ffa600", "#89DA59", "#FF420E"), c("Mouse", "Pig", "Rhesus", "Human"))
plot.margin <- unit(c(0.05, 0, 0.05, 0), "inch")
dotlist <- lapply(names(sp_cols), function(sp) {
	cls_order <- switch(sp, Human = c("Astro", "GC"), 
					Rhesus = c("Astro", "nIPC", "NB", "GC"),
					Pig = c("Astro", "NB", "GC"),
					Mouse = c("Astro", "RGL", "nIPC", "NB", "GC"))
	Idents(sublist[[sp]]) <- "fig2cluster"
	p <- DotPlot(sublist[[sp]], features = rev(all_genes), cols = c("lightgrey", sp_cols[sp]), dot.scale = 3, dot.min = 0.025, scale.by = "size")+
	        RotatedAxis() +
			scale_y_discrete(limits = rev(cls_order)) +
			##scale_size(range = c(0, 3), limits = c(2.5, 100)) +
			scale_radius(range = c(0, 4), limits = c(2.5, 100)) +
	        theme(axis.text.x=element_text(size = 4.25), axis.text.y=element_text(size = 6),legend.position = "right", axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid.major = element_line(size = 0.25, color = "grey"), panel.grid.minor = element_line(size = 0.25, color = "grey"), plot.margin = plot.margin)
	if (sp %in% c("Mouse", "Pig", "Rhesus")){
		p <- p + 
			theme(axis.text.x=element_blank(), legend.position = "right")
	}
	return(p)
	})



pp <- ggarrange(plotlist = dotlist, nrow = 4, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(0.9, 0.65, 0.8, 1))
p.anno <- p.anno + 
			theme(axis.text.x = element_blank())
pdf(paste0(outputdir, "MF2.mar.anno_dot.pdf"), useDingbats = FALSE, width = 8.5, height = 3.5)
plot_grid(p.anno, pp, nrow = 2, ncol = 1, rel_heights = c(1, 3), align = "v") %>% print()
dev.off()
dev.off()


##------------------------------------------------------------------------
## 3. putative human cells [Heatmap exp]  
pt_cells <- readRDS(file = paste0(inputdir, "Putative.cbn.cells.rds"))
hdg <- readRDS(paste0(inputdir, "Human_DG_seu.rds"))
hgc <- hdg %>%
			subset(mres %in% "GC")


## Get the bg cells
set.seed(42)
bg_cells <- sample(setdiff(colnames(hgc), pt_cells), length(pt_cells))

pt_cells <- setdiff(pt_cells, c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG")) %>%
				##.[c(1:10, 18, 12:17, 11)] %>%
				c(., c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG"))

all_cells <- c(pt_cells, bg_cells)
high_cell <- c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG")
plot_data <- hdg$RNA@counts[all_genes, all_cells] %>%
					as.matrix() %>%
					reshape2::melt() %>%
					setNames(., c("gene", "cell", "exp")) %>%
					mutate(gene = as.character(gene), cell = as.character(cell)) %>%
					mutate(ctp = ifelse(cell %in% pt_cells, "yes", "no")) %>%
					mutate(ctp = factor(ctp, levels = c("yes",  "no"))) %>%
					mutate(exp = MinMax(exp, min = 0, max = 3)) %>%
					mutate(exp = as.character(exp)) %>%
					mutate(color = ifelse(cell %in% high_cell, "black", "black"))##"#440154FF"

##colors <- colorRampPalette(colors = c("grey90", "red"))(4)[1:4]#440154FF
##colors <- colorRampPalette(colors = c("grey95", "red3"))(100)[c(1, 5,35,100)]
colors <- viridis(4)
#colors[1] <- "#FFFFFF"
putap <- ggplot(plot_data) + 
				geom_tile(aes_string(x = "gene", y = "cell", fill = "exp", color = "color"), height = 1, width = 1) + 
				theme_classic() +
				scale_y_discrete(limits = all_cells) +
				scale_x_discrete(limits = all_genes) +
				scale_fill_manual(values = colors) +
				scale_color_identity() + 
				RotatedAxis() + 
				facet_wrap(vars(ctp), nrow = 2, ncol = 1, scales = "free_y") +
	            theme(panel.spacing = unit(0.02, "in"), strip.background = element_blank(), strip.text = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 4.5), axis.line = element_blank(), legend.position = "bottom", axis.title = element_blank(), axis.ticks.x = element_blank())
pdf(paste0(outputdir, "MF2.mar.putative.exp.pdf"), width = 8.5, height = 6)
print(putap)
dev.off()



















