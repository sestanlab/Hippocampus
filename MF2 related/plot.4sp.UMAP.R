## Plot the integration UMAP figures
args <- commandArgs(trailingOnly = TRUE) 
source("../scripts/hip.fun.R")
library(ggpubr)



## Plot all or only GC lineage
lev <- c("all", "GC")[as.numeric(args[1])]
hrpm <- switch(lev,
        all = readRDS(file = paste0(inputdir, "HRPM.all.seurat.expHVG.1500.slim.rds")),
        GC = readRDS(file = paste0(inputdir, "HRPM.GC.seurat.expHVG.1500.slim.rds")))


if (lev == "GC"){
	hrpm <- hrpm[, hrpm@meta.data$vis.batch != "Mouse young"]
}


## Get the UMAP plot data
plot_data <- data.frame(cluster = hrpm@meta.data$fig2cluster,
              species = hrpm@meta.data$species,
              xaxis = hrpm$umap@cell.embeddings[, 1], 
              yaxis = hrpm$umap@cell.embeddings[, 2],
              stringsAsFactors = FALSE)

## Mix the cells
set.seed(42)
plot_data <- plot_data[sample(1:nrow(plot_data), nrow(plot_data)), ]


##------------------------------------------------------------------------------------
## Plots for clusters
## Set the colors for clusters
colors <- switch(lev,
        all = c("#e01134", "#d834b7", "#7dc6e8", "#08519c","#d8f796","#7bf28d","#92f7a1","#dbac55", "#d6821b", "#ef793e", "#b77305","#fcb5c5", "#f4f4f4", "#6b6b6b", "#3d3d3d", "#bfa3a3", "#c1c1c1", "#f59847") %>% 
    setNames(., c("RGL","nIPC","NB", "GC","MC","CA2-3", "CA1 Sub", "InN","EC L2-3","EC L5-6", "EC L6","Astro", "OPC", "Oligo", "immune", "Ependymal", "Vas", "CR")),
        GC = c("#fcb5c5","#e01134", "#d834b7", "#7dc6e8", "#08519c") %>% 
              setNames(., c("Astro","RGL","nIPC","NB","GC")))

bb <- 0
xrange <- c(min(plot_data$xaxis)-bb, max(plot_data$xaxis) +bb)
yrange <- c(min(plot_data$yaxis)-bb, max(plot_data$yaxis) +bb)
print(xrange)
print(yrange)

sp_order <- c("Mouse", "Pig", "Rhesus", "Human")
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
sp_cols <- setNames(paste0(c("#ffa600", "#031a7f", "#89DA59", "#FF420E"), ""), c("Pig", "Mouse", "Rhesus", "Human"))

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
jpeg(paste0(outputdir, "MF2.Inte.", lev, ".anno.jpeg"), width = 20, height = 4.5, units = "in", res = 300)
patchwork::wrap_plots(newlist, nrow = 1, ncol = 5) ##%>% print()
dev.off()













