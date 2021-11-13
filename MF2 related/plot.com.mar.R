## Plot the expression of shared markers
source("../scripts/hip.fun.R")
library(ggpubr)


##-----------------------------------------------------------------------------------
## Plot gene annotation
load(paste0(inputdir, "Markers.RPM.conserved.Rdata"))
all_genes <- unlist(con_list[c("common", "nIPC_RM", "NB_RMP", "NB_PM", "NB_RM", "NB_RP")]) %>% unique()


gene2df <- lapply(con_list, function(x)
	as.numeric(all_genes %in% x) %>% setNames(., all_genes)
	) %>% 
			as.data.frame(, check.names = FALSE) %>%
			.[, setdiff(names(con_list), "common")] %>%
			rownames_to_column("gene") %>%
			mutate(nIPC_Rhesus = nIPC_RM, nIPC_Mouse = nIPC_RM) %>%
			mutate(NB_Pig = as.numeric((NB_RMP + NB_PM + NB_RP) != 0)) %>%
			mutate(NB_Mouse = as.numeric((NB_RMP + NB_PM + NB_RM) != 0)) %>%
			mutate(NB_Rhesus = as.numeric((NB_RMP + NB_RP + NB_RM) != 0)) %>%
			select(-c(nIPC_RM, NB_RMP, NB_PM, NB_RP, NB_RM))

sp_cols <- setNames(c("#ffa600", "#031a7f", "#89DA59", "#FF420E"), c("Pig", "Mouse", "Rhesus", "Human"))
plot_data <- gene2df %>%
				tidyr::gather(., "type", "value", c("nIPC_Rhesus", "nIPC_Mouse", "NB_Mouse", "NB_Pig", "NB_Rhesus")) %>%
				mutate(gene = factor(as.character(gene), levels = all_genes)) %>%
				mutate(value = as.character(value)) %>%
				mutate(species = extract_field(type, 2, "_")) %>%
				mutate(color = ifelse(value == "1", sp_cols[species], "#D3D3D3")) %>%
				mutate(type = factor(as.character(type), levels = c("nIPC_Mouse", "nIPC_Rhesus", "NB_Mouse", "NB_Pig", "NB_Rhesus") %>% rev()))


p <- ggplot(plot_data) + 
				geom_tile(aes_string(x = "gene", y = "type", fill = "color"), color = "white", height = 0.9, width = 0.9, size = 0.1) + 
				theme_classic() +
				scale_fill_identity() + 
				RotatedAxis() + 
	            theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 6), axis.line = element_blank(), legend.position = "bottom", axis.title = element_blank(), axis.ticks.x = element_blank())
pdf(paste0(outputdir, "SF3.conserved.markers.anno.pdf"), width = 8, height = 2)
print(p)
dev.off()



##-----------------------------------------------------------------------------------
## Plot gene dot plot
sublist <- readRDS(file = paste0(inputdir, "HRPM.GC.seurat.expHVG.1500.slim.rds")) %>%
			SplitObject(., split.by = "species")

load(paste0(inputdir, "Markers.RPM.conserved.Rdata"))
all_genes <- unlist(con_list[c("common", "nIPC_RM", "NB_RMP", "NB_PM", "NB_RM", "NB_RP")]) %>% unique()


sp_cols <- setNames(c("#031a7f", "#ffa600", "#89DA59", "#FF420E"), c("Mouse", "Pig", "Rhesus", "Human"))
plot.margin <- unit(c(0.05, 0, 0.05, 0), "inch")
plist <- lapply(names(sp_cols), function(sp) {
	cls_order <- switch(sp, Human = c("Astro", "GC"), 
					Rhesus = c("Astro", "nIPC", "NB", "GC"),
					Pig = c("Astro", "NB", "GC"),
					Mouse = c("Astro", "RGL", "nIPC", "NB", "GC"))
	Idents(sublist[[sp]]) <- "fig2cluster"
	p <- DotPlot(sublist[[sp]], features = rev(all_genes), cols = c("lightgrey", sp_cols[sp]), dot.scale = 3, dot.min = 0.025, scale.by = "size")+
	        RotatedAxis() +
			scale_y_discrete(limits = rev(cls_order)) +
			scale_size(range = c(0, 3), limits = c(2.5, 100)) +
	        theme(axis.text.x=element_text(size = 6), axis.text.y=element_text(size = 6),legend.position = "right", axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid.major = element_line(size = 0.25, color = "grey"), panel.grid.minor = element_line(size = 0.25, color = "grey"), plot.margin = plot.margin)
	if (sp %in% c("Mouse", "Pig", "Rhesus")){
		p <- p + 
			theme(axis.text.x=element_blank(), legend.position = "right")
	}
	return(p)
	})

pp <- ggarrange(plotlist = plist, nrow = 4, ncol = 1, common.legend = TRUE, legend = "bottom", heights = c(0.9, 0.65, 0.8, 1.125))
pdf(paste0(outputdir, "SF3.conserved.markers.Dotplot.pdf"), useDingbats = FALSE, width = 8.5, height = 2.8)

##patchwork::wrap_plots(plist, nrow = 4, ncol = 1, guides = "collect")
print(pp)
dev.off()
dev.off()



##-----------------------------------------------------------------------------------
## Plot putative cell expression
load(paste0(inputdir, "Markers.RPM.conserved.Rdata"))
all_genes <- unlist(con_list[c("common", "nIPC_RM", "NB_RMP", "NB_PM", "NB_RM", "NB_RP")]) %>% unique()



hdg <- readRDS(paste0(inputdir, "Human_DG_seu.rds"))
hgc <- hdg %>%
			subset(mres %in% "GC")


all_cells <- c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG")
plot_data <- hdg$RNA@counts[all_genes, all_cells] %>%
					as.matrix() %>%
					reshape2::melt() %>%
					setNames(., c("gene", "cell", "exp")) %>%
					mutate(gene = as.character(gene), cell = as.character(cell)) %>%
					mutate(exp = MinMax(exp, min = 0, max = 3)) %>%
					mutate(exp = as.character(exp)) ##%>%
					##mutate(color = ifelse(cell %in% high_cell, "black", "black"))##"#440154FF"

colors <- viridis(4)
p <- ggplot(plot_data) + 
				geom_tile(aes_string(x = "gene", y = "cell", fill = "exp"), color = "black", height = 1, width = 1) + 
				theme_classic() +
				scale_y_discrete(limits = all_cells) +
				scale_x_discrete(limits = all_genes) +
				scale_fill_manual(values = colors) +
				scale_color_identity() + 
				RotatedAxis() +
	            theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 6), axis.line = element_blank(), legend.position = "bottom", axis.title = element_blank(), axis.ticks.x = element_blank())
pdf(paste0(outputdir, "SF3.conserved.markers.putative.exp.pdf"), width = 8.5, height = 1.5)
print(p)
dev.off()


##-----------------------------------------------------------------------------------
## Do the classification of these conserved markers 
load(paste0(inputdir, "Markers.RPM.conserved.Rdata"))
all_genes <- unlist(con_list[c("common", "nIPC_RM", "NB_RMP", "NB_PM", "NB_RM", "NB_RP")]) %>% unique()



hdg <- readRDS(paste0(inputdir, "Human_DG_seu.rds"))
hgc <- hdg %>%
			subset(mres %in% "GC")


bg_ratio <- Matrix::rowMeans(hgc$RNA@data != 0)
cell_use <- c("HSB179_6_eDG_AATTTCCGTCCTACAA", "HSB628_1_DG_CCGATGGGTGTTAAAG")[c(2, 1)]
plot_data <- data.frame(gene = all_genes, 
					inBG = as.numeric(bg_ratio[all_genes] >= 0.1), 
					inIPC = as.numeric(hdg$RNA@data[all_genes, cell_use[1]] != 0) * 2, 
					inNB = as.numeric(hdg$RNA@data[all_genes, cell_use[2]] != 0) * 2, 
					stringsAsFactors = FALSE) %>%
			mutate(sumIPC = as.character(inBG+inIPC), sumNB = as.character(inBG+inNB))
plot_data$c1 <- plot_data$c2 <- plot_data$c3 <- plot_data$c4 <- 0

ipcidx <- length(unlist(con_list[c("common", "nIPC_RM")]) %>% unique())
nidx <- 12:ipcidx
iidx <- (ipcidx+1):length(all_genes)


plot_data$c1[nidx] <- ifelse(plot_data$sumIPC[nidx] == "2", 1, 0) 
plot_data$c2[nidx] <- ifelse(plot_data$sumIPC[nidx] == "0", 1, 0) 
plot_data$c3[nidx] <- ifelse(plot_data$sumIPC[nidx] == "3", 1, 0) 
plot_data$c4[nidx] <- ifelse(plot_data$sumIPC[nidx] == "1", 1, 0) 
plot_data$c1[iidx] <- ifelse(plot_data$sumNB[iidx] == "2", 1, 0) 
plot_data$c2[iidx] <- ifelse(plot_data$sumNB[iidx] == "0", 1, 0) 
plot_data$c3[iidx] <- ifelse(plot_data$sumNB[iidx] == "3", 1, 0) 
plot_data$c4[iidx] <- ifelse(plot_data$sumNB[iidx] == "1", 1, 0) 
plot_data <- plot_data %>%
				tidyr::gather("group", "value", c("c1", "c2", "c3", "c4")) %>%
				mutate(color = ifelse(value == 0, "grey85", "black")) %>%
				mutate(group = factor(as.character(group), levels = rev(c("c1", "c2", "c3", "c4"))))


p <- ggplot(plot_data) + 
				geom_tile(aes_string(x = "gene", y = "group", fill = "color"), color = "grey90", height = 0.9, width = 0.9, size = 0.1) + 
				theme_classic() +
				scale_fill_identity() + 
				scale_y_discrete(limits = c("c4", "c3", "c2", "c1")) +
				scale_x_discrete(limits = all_genes)+
				RotatedAxis() + 
	            theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 6), axis.line = element_blank(), legend.position = "bottom", axis.title = element_blank(), axis.ticks.x = element_blank())
pdf(paste0(outputdir, "SF3.conserved.markers.classification.anno.pdf"), width = 8, height = 2)
print(p)
dev.off()














