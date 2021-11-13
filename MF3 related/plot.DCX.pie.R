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


## Get the data list for plot
exp_data <- lapply(names(seulist), function(sp) {
	df <- data.frame(exp = seulist[[sp]]$RNA@counts["DCX", ], 
					species = sp,
					cluster = seulist[[sp]]@meta.data$fig2cluster,
					xaxis = seulist[[sp]]$umap@cell.embeddings[, 1], 
					yaxis = seulist[[sp]]$umap@cell.embeddings[, 2], 
					strinsAsFactors = FALSE)
	return(df)
	})  %>%
		do.call(rbind, .) ##%>%
		##mutate(exp = MinMax(exp, min = 0, max = 3)) %>%
		##mutate(exp = as.character(exp))


## Plot the expression [Barplots]
colors <- c("#e01134", "#d834b7", "#7dc6e8", "#08519c","#d8f796","#7bf28d","#92f7a1","#dbac55", "#d6821b", "#ef793e", "#b77305","#fcb5c5", "#f4f4f4", "#6b6b6b", "#3d3d3d", "#c1c1c1", "#f59847") %>% 
    setNames(., c("RGL","nIPC","NB", "GC","MC","CA2-3", "CA1 Sub", "InN","EC L2-3","EC L5-6", "EC L6","Astro", "OPC", "Oligo", "immune", "Vas", "CR"))
pplist <- lapply(c(1, 2), function(thre) {
	plot_data <- lapply(names(order_list)[1:4], function(sp){
		xx <- subset(exp_data, species == sp) %>%
				subset(exp >= thre) %>%
				group_by(cluster) %>%
				summarize(ncells = n())
		xx$cluster <- factor(xx$cluster, levels = names(colors)[names(colors) %in% xx$cluster])
		xx
		}) %>%
			setNames(., names(order_list)[1:4])


	pie_list <- lapply(names(plot_data), function(sp) {
		p <- ggplot(plot_data[[sp]], aes(x = "", y = ncells, fill = cluster)) +
				geom_bar(width = 1, stat = "identity") +
				coord_polar("y", start=0, direction = -1) +
				theme_classic() + 
				scale_fill_manual(values = colors[as.character(plot_data[[sp]]$cluster)]) +
				theme(legend.position = "right", axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), plot.margin = unit(rep(-0.025, 4), "inch"))
		return(p)
		})
	pie_list
	}) %>% do.call(c, .) %>%
		.[c(1,5,2,6,3,7,4,8)]




ggarrange(plotlist = pplist, nrow = 4, ncol = 2, common.legend = TRUE, legend = "right") %>%
		ggexport(filename = paste0(outputdir, "MF2.DCX.pie.pdf"), width = 5, height = 6)




## Output a table summarize the number of cells expressing DCX across dataset 
sum_data <- lapply(c(1, 2), function(thre) {
	plot_data <- lapply(names(order_list)[1:4], function(sp){
		xx <- subset(exp_data, species == sp) %>%
				subset(exp >= thre) %>%
				group_by(cluster) %>%
				summarize(ncells = n(), species = unique(species)) %>%
				mutate(threshold = thre, ratio = ncells * 100/sum(ncells))
		xx
		}) %>%
			do.call(rbind, .)
	plot_data
	}) %>%
		do.call(rbind, .)
write.table(sum_data, file = paste0(inputdir, "DCX.exp.pie.summarize.all.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)







