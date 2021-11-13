## Plot the expression of DCX across all cluster in mouse, rhesus, Human DG and human EC cells.
source("../scripts/hip.fun.R") 
library(ggpubr)


seulist <- readRDS(file = paste0(inputdir, "Subset.data.Cbn.all.rds"))
seulist$Human_EC@meta.data$fig2cluster <- gsub("InN CGE", "InN", seulist$Human_EC@meta.data$fig2cluster) %>%
								gsub("InN MGE", "InN", .)

order_list <- list(Mouse = c("RGL","nIPC","NB","GC","MC","CA2-3","InN","CR","OPC","Oligo","Astro","immune","Ependymal","Vas"), 
				Pig = c("RGL", "nIPC","NB","GC","MC","CA2-3","CA1 Sub","InN","OPC","Oligo","Astro","immune","Vas"),
			Rhesus = c("RGL", "nIPC","NB","GC","MC","CA2-3","CA1 Sub","InN","CR","OPC","Oligo","Astro","immune","Vas"), 
			Human = c("GC", "MC", "CA2-3", "InN", "CR", "OPC", "Oligo", "Astro", "immune", "Vas"), 
			Human_EC = c("EC L2-3", "EC L5-6", "EC L6", "InN", "CR", "OPC", "Oligo", "Astro", "immune", "Vas")) ## Merge two InN types in EC region
##Validate the vector length
for (ii in names(order_list)){
    print(length(order_list[[ii]]))
    print(sum(unique(seulist[[ii]]@meta.data$fig2cluster) %in% order_list[[ii]]))
}
colors <- colorRampPalette(colors = c("grey95", "red3"))(100)[c(1, 5,35,100)]


## Load the integration dataset for coordinates
hrpm <- readRDS(file = paste0(inputdir, "HRPM.all.seurat.expHVG.1500.slim.rds"))
ceb <- rbind(hrpm$umap@cell.embeddings, seulist$Human_EC$umap@cell.embeddings)


## Get the data list for plot
exp_data <- lapply(names(seulist), function(sp) {
	df <- data.frame(exp = seulist[[sp]]$RNA@counts["DCX", ], 
					species = sp,
					cluster = seulist[[sp]]@meta.data$fig2cluster,
					xaxis = ceb[colnames(seulist[[sp]]), 1], 
					yaxis = ceb[colnames(seulist[[sp]]), 2], 
					stringsAsFactors = FALSE)
	return(df)
	})  %>%
		do.call(rbind, .) %>%
		mutate(exp = MinMax(exp, min = 0, max = 3)) %>%
		mutate(exp = as.character(exp))


## Plot the expression [UMAP]
sub_exp <- subset(exp_data, species != "Human_EC")
bb <- 0
xrange <- c(min(sub_exp$xaxis)-bb, 
			max(sub_exp$xaxis) + bb)
yrange <- c(min(sub_exp$yaxis)- bb, 
			max(sub_exp$yaxis) + bb)


umap_list <- lapply(names(order_list), function(sp) {
			sub_data <- subset(exp_data, species == sp)
			sub_data <- sub_data[c(which(sub_data$exp == "0"), which(sub_data$exp != "0")), ]
	     	p <- ggplot(sub_data, aes(x = xaxis, y = yaxis, color = exp)) +
	      		geom_point(size = 0.2) +
	      		scale_color_manual(values = setNames(colors, c("0", "1", "2", "3"))) +
	        	theme_classic() + 
	        	theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")
	       	if (!sp %in% c("Human_EC")){
	       		p <- p + 
	       				xlim(limits = xrange) +
						ylim(limits = yrange)
	       	}
			return(p)
		})

umap_list <- umap_list %>%
          lapply(., function(x) x + theme(plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "inch")))

jpeg(paste0(outputdir, "MF2.DCX.UMAP.v2.jpeg"), width = 20, height = 4, units = "in", res = 300)
plot_grid(plotlist = umap_list, nrow = 1, ncol = 5) %>% print()
dev.off()


jpeg(paste0(outputdir, "MF2.DCX.UMAP.highreso.jpeg"), width = 5*9, height = 9, units = "in", res = 300)
plot_grid(plotlist = umap_list, nrow = 1, ncol = 5) %>% print()
dev.off()


## Plot the expression [Barplots]
maxy <- exp_data %>%
				group_by(species, cluster) %>%
				summarize(ratio = mean(exp!="0") * 100) %>%
				.$ratio %>% max()
bar_list <- lapply(names(order_list), function(sp) {
	data <- subset(exp_data, species == sp) %>%
				group_by(exp, cluster) %>%
				summarize(ncells = n()) %>%
				ungroup() %>%
				group_by(cluster) %>%
				mutate(ratio = ncells * 100/sum(ncells)) %>%
				subset(exp != "0")
	label_data <- data %>%
				summarize(sumcells = sum(ncells), sumratio = sum(ratio))
	new_range <- c(paste0("1"), paste0("2"), paste0(">=", ii))
	p <- ggplot(data) +
			geom_bar(aes(x = cluster, y = ratio, fill = exp), width = 1, stat = "identity") +
			geom_text(data = label_data, mapping = aes(x = cluster, y = sumratio, label = sumcells), size = 2, angle = 45, hjust = 0, vjust = 0) + 
			theme_cowplot() +
			RotatedAxis() +
			scale_fill_manual(values = setNames(colors[2:4], c("1", "2", "3")), labels = setNames(new_range, c("1", "2", "3"))) +
			scale_y_continuous(limits = c(0, maxy)) +
			scale_x_discrete(limits = order_list[[sp]]) +
			theme(legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
	if (sp != "Mouse"){
		p <- p + 
			theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
	}
	return(p)
	})

pdf(paste0(outputdir, "MF2.DCX.bar.pdf"), width = 8.5, height = 1.7)
patchwork::wrap_plots(bar_list, nrow = 1, ncol = 5) %>% print()
dev.off()



## Output a table summarize the number of cells expressing DCX across dataset 
sum_data <- exp_data %>%
				group_by(exp, cluster, species) %>%
				summarize(ncells = n()) %>%
				ungroup() %>%
				group_by(cluster, species) %>%
				mutate(ratio = ncells * 100/sum(ncells))
write.table(sum_data, file = paste0(inputdir, "DCX.exp.ratio.summarize.all.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


##--------------------------------------------------------------------
## Plot the library-size normalized expression [Barplots]
avgs <- lapply(names(seulist), function(sp) {
	xx <- seulist[[sp]]
	Idents(xx) <- "fig2cluster"
	avg <- log(AverageExpression(xx, features = c("DCX"))$RNA + 1)
	colnames(avg) <- paste0(sp, "|", colnames(avg))
	avg
	}) %>% 
		do.call(cbind, .)

avg_data <- avgs %>%
			as.matrix() %>%
			reshape2::melt() %>%
			setNames(., c("gene", "sp_cls", "exp")) %>%
			mutate(gene = as.character(gene)) %>%
			mutate(cluster = extract_field(as.character(sp_cls), 2, "|")) %>%
			mutate(species = extract_field(as.character(sp_cls), 1, "|")) %>%
			select(-sp_cls)


bar_list <- lapply(names(order_list), function(sp) {
	data <- subset(avg_data, species == sp)

	p <- ggplot(data) +
			geom_bar(aes(x = cluster, y = exp), fill = "white", color = "black", width = 1, size = 0.2, stat = "identity") +
			##geom_text(data = label_data, mapping = aes(x = cluster, y = sumratio, label = sumcells), size = 2, angle = 45, hjust = 0, vjust = 0) + 
			theme_cowplot() +
			RotatedAxis() +
			scale_y_continuous(limits = c(0, max(avg_data$exp))) +
			scale_x_discrete(limits = order_list[[sp]]) +
			theme(legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
	if (sp != "Mouse"){
		p <- p + 
			theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
	}
	return(p)
	})

pdf(paste0(outputdir, "MF2.DCX_normExp.bar.pdf"), width = 8.5, height = 1.7)
patchwork::wrap_plots(bar_list, nrow = 1, ncol = 5) %>% print()
dev.off()




