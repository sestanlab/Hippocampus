## Plot the integration UMAP figures
args <- commandArgs(trailingOnly = TRUE) 
source("../scripts/hip.fun.R")
library(ggpubr)


seu <- readRDS(paste0(inputdir, "Human.AdultFetal.inte.GC.seurat.2000.slim.rds"))
seu@meta.data$mres <- gsub("NB_GC$", "NB", seu@meta.data$mres) %>%
                        gsub("nIPC_cyc$", "nIPC", .) %>%
                        gsub("nIPC_GC$", "nIPC", .) %>%
                        gsub("RGC_cyc$", "RGL", .) %>%
                        gsub("RGC$", "RGL", .)
seulist <- SplitObject(seu, split.by = "species")


## Get the UMAP plot data
plot_data <- lapply(names(seulist), function(sp) {
	df <- data.frame(cluster = seulist[[sp]]@meta.data$mres,
					species = sp,
					xaxis = seulist[[sp]]$umap@cell.embeddings[, 1], 
					yaxis = seulist[[sp]]$umap@cell.embeddings[, 2], 
					stringsAsFactors = FALSE)
	return(df)
	})  %>%
		do.call(rbind, .)


##------------------------------------------------------------------------------------
## Plots for clusters
## Set the colors for clusters
colors <- c("#fcb5c5","#e01134", "#d834b7", "#7dc6e8", "#08519c") %>% 
              setNames(., c("Astro","RGL","nIPC","NB","GC"))

bb <- 0
xrange <- c(min(seu$umap@cell.embeddings[, 1])-bb, max(seu$umap@cell.embeddings[, 1]) +bb)
yrange <- c(min(seu$umap@cell.embeddings[, 2])-bb, max(seu$umap@cell.embeddings[, 2]) +bb)

sp_order <- c("HumanAdult", "HumanFetal")
annolist <- lapply(sp_order, function(sp) {
     	p <- ggplot(subset(plot_data, species == sp), aes(x = xaxis, y = yaxis, color = cluster)) +
      		geom_point(size = 0.01) +
      		scale_color_manual(values = colors[unique(plot_data$cluster[plot_data$species == sp])]) +
        	theme_classic() + 
        	labs(title = sum(plot_data$species == sp)) +
        	theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")  +
        	xlim(limits = xrange) +
			ylim(limits = yrange)
      return(p)
}) %>% 
    setNames(., sp_order)


##------------------------------------------------------------------------------------
## Plots for species
## Set the species colors
sp_cols <- setNames(paste0(c("#FF420E", "#8c510a"), ""), c("HumanAdult", "HumanFetal"))

p1 <- ggplot(plot_data, aes(x = xaxis, y = yaxis, color = species)) +
      	geom_point(size = 0.01) +
      	scale_color_manual(values = sp_cols) +
        theme_classic() + 
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")  +
        xlim(limits = xrange) +
		ylim(limits = yrange)

interval <- -0.025
newlist <- c(list(Species = p1), annolist) %>%
          lapply(., function(x) x + theme(plot.margin = unit(c(interval, interval, interval, interval), "inch")))
jpeg(paste0(outputdir, "SF2.Inte.HumanAdultFetal.GC.anno.jpeg"), width = 12, height = 4.5, units = "in", res = 300)
patchwork::wrap_plots(newlist, nrow = 1, ncol = 3) ##%>% print()
dev.off()


