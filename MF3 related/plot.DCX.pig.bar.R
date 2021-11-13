## Plot the expression of DCX across all cluster in mouse, rhesus, Human DG and human EC cells.
source("../scripts/hip.fun.R")
library(ggpubr)



## Load the pig integration dataset
pig <- readRDS(file = paste0(inputdir, "Pig.017.DG.seu.rds")) 


## Get the data list for plot
cls_order <- c("RGL", "nIPC","NB","GC","MC","CA2-3","CA1 Sub","InN","OPC","Oligo","Astro","immune","Vas")
colors <- colorRampPalette(colors = c("grey95", "red3"))(100)[c(1, 5,35,100)]


plist <- lapply(c("DCX", "GAPDH"), function(x) {
	exp_data <- pig@meta.data %>%
				rownames_to_column("cell") %>%
				mutate(cluster = fig2cluster_new) %>%
				select(condition, cluster, samplename) %>%
				mutate(exp = pig$RNA@counts[x, ]) %>%
				mutate(cluster = factor(cluster, levels = cls_order)) %>%
				mutate(xaxis = pig$umap@cell.embeddings[, 1], 
						yaxis = pig$umap@cell.embeddings[, 2]) 

	plot_data <- exp_data %>%
					group_by(samplename, cluster, condition, .drop=FALSE) %>%
					summarize(ratio = mean(exp > 0) * 100) %>%
					mutate(condition = extract_field(samplename, 2, "_"))
	plot_data$ratio[is.na(plot_data$ratio)] <- 0

	plot_data <- plot_data %>%
					ungroup() %>%
					group_by(condition, cluster) %>%
					summarize(mratio = mean(ratio), rse = sd(ratio)/sqrt(n())) %>%
					ungroup() %>%
					mutate(condition = factor(as.character(condition), levels = c("h0", "h1", "h7")))

	p <- ggplot(plot_data, aes_string(x = "cluster", y = "mratio", fill = "condition")) +
              geom_bar(size = 0.2, stat="identity", color="black", position=position_dodge(0.45), width = 0.5) +
              geom_errorbar(aes(ymin=mratio-rse, ymax=mratio+rse), width=.1, position=position_dodge(0.45), size = 0.25) +
              theme_cowplot() +
			  RotatedAxis() +
              ##scale_fill_manual(values = gp_colors) +
              ##scale_x_discrete(limits = new_order, labels = extract_field(new_order, 2, ".")) +
              theme(legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
    p
	})

##ggpubr::ggarrange(plotlist = plist, nrow = 2, ncol = 1, align = "v", common.legend = TRUE, legend = "bottom") %>%
##		ggpubr::ggexport(filename = paste0(outputdir, "SF3.DCX.pig.bar.pdf"), width = 4, height = 3)




##  Bar plots showing the library-size normalized expression
pig@meta.data$newcls <- paste0(pig@meta.data$samplename, "|", pig@meta.data$fig2cluster_new)
Idents(pig) <- "newcls"
avgs <- log(AverageExpression(pig, features = c("DCX", "GAPDH"))$RNA + 1)

avg_data <- avgs %>%
			as.matrix() %>%
			reshape2::melt() %>%
			setNames(., c("gene", "cdt_cls", "exp")) %>%
			mutate(gene = as.character(gene)) %>%
			mutate(cluster = extract_field(as.character(cdt_cls), 2, "|")) %>%
			mutate(samplename = extract_field(as.character(cdt_cls), 1, "|")) %>%
			mutate(condition = extract_field(samplename, 2, "_")) %>%
			select(-c(cdt_cls, samplename)) %>%
			subset(gene == "DCX") %>%
			group_by(condition, cluster) %>%
			summarize(mexp = mean(exp), rse = sd(exp)/sqrt(n())) %>%
			ungroup() %>%
			mutate(condition = factor(as.character(condition), levels = c("h0", "h1", "h7"))) %>%
			mutate(cluster = factor(cluster, levels = cls_order))



p2 <- ggplot(avg_data, aes_string(x = "cluster", y = "mexp", fill = "condition")) +
              geom_bar(size = 0.2, stat="identity", color="black", position=position_dodge(0.45), width = 0.5) +
              geom_errorbar(aes(ymin=mexp-rse, ymax=mexp+rse), width=.1, position=position_dodge(0.45), size = 0.25) +
              theme_cowplot() +
			  RotatedAxis() +
              ##scale_fill_manual(values = gp_colors) +
              ##scale_x_discrete(limits = new_order, labels = extract_field(new_order, 2, ".")) +
              theme(legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1, size = 6))


p1 <- plist[[1]] +
		theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p_list2<- append(list(p1, p2), list(legend = patchwork::guide_area()), 2)


layout1<-"
A
B
C
"


pdf(paste0(outputdir, "SF3.DCX_exp.pig.bar.pdf"), width = 3, height = 3)
patchwork::wrap_plots(p_list2, guides = "collect", design = layout1)
dev.off() 








## Get the expression plots using schex
library(schex)
pig_list <- SplitObject(pig, split.by = "condition") %>%
			lapply(., function(x) make_hexbin(x, nbins = 100, dimension_reduction = "UMAP"))


gene <- "DCX"
p <- plot_hexbin_feature_byregion(sce = pig_list, assay = "RNA", slot = "data", features = gene, action = "mean", nrow = 1, ncol = 3, colors = c("#f0f0f0", "#d9d9d9", "#fdbb84", "#b30000","#7f0000"), file_name = paste0("SF3.Pig.PMI.", gene, ".hex"), output_dir = outputdir, pdf_size = c(12,4), region_order = c("h0", "h1", "h7"), legend = "bottom", return_rawp = TRUE, exp_ceiling = 2.5)[[gene]]
pdf(paste0(outputdir, "SF3.Pig.PMI.", gene, ".Exp.hex.pdf"), width = 15, height = 6)
plot_grid(plotlist = p, nrow = 1, ncol = 3) %>% print()
dev.off() 









