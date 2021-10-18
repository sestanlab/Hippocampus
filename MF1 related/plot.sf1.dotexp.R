## Plot the MainFigure and supplementary figure
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")

hipec <- readRDS(file = paste0(dataDir, "HIPEC.final.seu.01172021.rds"))
color_vec <- readRDS(file = paste0(dataDir, "cluster.color.rds"))


##------------------------------------------------------------------
## Plot the expression of selected genes
genes <- read.table(file = paste0(inputdir, "plot.genes.txt"), sep = "\t", header = TRUE, stringsAsFactor = FALSE) %>% .$gene %>% unique()


Idents(hipec) <- "fig1cluster"
p <- DotPlot(hipec, features = genes, cols = c("lightgrey", "red"), dot.min = 0.05, dot.scale = 2) +
			scale_x_discrete(limits = genes, position = "top") +
			scale_y_discrete(limits = names(color_vec) %>% rev()) +
			theme(legend.position = "bottom", axis.title = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text.x = element_text(size = 4.5, angle = 45,hjust = 0), axis.text.y = element_text(size = 4.5))

pdf(paste0(outputdir, "SF1.Exp.Dotplot.pdf"), width = 8, height = 6, useDingbats = FALSE)
print(p)
dev.off()



