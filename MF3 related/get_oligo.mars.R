source("../scripts/hip.fun.R")


## Get the Oligo markers prepared for calculating AUROC scores.

allgc <- readRDS(file = paste0(inputdir, "Konopka_HIP_seu.rds"))
Idents(allgc) <- "Cluster"
mm <- FindMarkers(allgc, ident.1 = c("Olig1", "Olig2", "Olig3", "Olig4", "Olig5"), max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.2)


mar <- mm %>%
		rownames_to_column("gene") %>%
		mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
		subset(pct.1 >= 0.4) %>%
		subset(ratio_fc >= 1.25) %>%
		mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
		subset(p_val_adj <= 0.01) %>%
		arrange(desc(avg_logFC)) %>%
		.$gene


## Dot Plots Visualizing the first 50 genes
all_cls <- c("Den.Gyr1", "Den.Gyr2", "Den.Gyr3", "Pyr1", "Pyr2", "In1", "In2", "In3", "Astro1", "Astro2", "Astro3","Endo", "Micro1", "Micro2", "Micro3", "OPC1", "OPC2", "OPC3", "OPC4", "Olig1", "Olig2", "Olig3", "Olig4", "Olig5")
allgc@meta.data$Cluster <- factor(as.character(allgc@meta.data$Cluster), levels = all_cls)
p <- DotPlot(object = allgc, features = mar[1:50], cols = c("lightgrey", "red"), dot.scale = 3, dot.min = 0.05) + 
        RotatedAxis() + 
        coord_flip() + 
        scale_y_discrete(limits = all_cls) +
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))


pdf(paste0(outputdir, "Konopka_Oligo_markers.pdf"), width = 6, height = 8)
print(p)
dev.off()


saveRDS(mar, file = paste0(inputdir, "Konopka_Oligo_mar.rds"))

