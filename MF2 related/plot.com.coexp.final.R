## Plot the expression of DCX & CALB2 across all cluster in mouse, rhesus, Human DG and human EC cells.
source("../scripts/hip.fun.R")
library(ggpubr)


seulist <- readRDS(file = paste0(inputdir, "HRPM.all.seurat.expHVG.1500.slim.rds")) %>%
                SplitObject(., split.by = "species")
seulist$Human_EC <- readRDS(paste0(inputdir, "Subset.data.Human_EC.rds"))
seulist$Human_EC@meta.data$fig2cluster <- gsub("InN CGE", "InN", seulist$Human_EC@meta.data$fig2cluster) %>%
                                gsub("InN MGE", "InN", .)

order_list <- list(Mouse = c("RGL","nIPC","NB","GC","CA2-3","InN","CR","OPC","Oligo","Astro","immune","Ependymal","Vas"), 
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
gene_list <- list(co1 = c("TOP2A", "CENPF", "SMC4", "MKI67", "EZH2", "NUSAP1", "CENPE"), 
                co2 = c("PROX1", "DCX", "CALB2", "DPYSL3"))
exp_list <- lapply(names(gene_list), function(tt) {
    exp_data <- lapply(names(seulist), function(sp) {
        thre <- switch(tt, co1 = 4, co2 = 4)
        df <- data.frame(coexp = as.numeric(colSums(seulist[[sp]]$RNA@counts[gene_list[[tt]], ]!=0) >=thre),
                        species = sp,
                        cluster = seulist[[sp]]@meta.data$fig2cluster,
                        xaxis = seulist[[sp]]$umap@cell.embeddings[, 1], 
                        yaxis = seulist[[sp]]$umap@cell.embeddings[, 2], 
                        strinsAsFactors = FALSE)
        if (tt == "co2"){
            df$coexp <- as.numeric(colSums(seulist[[sp]]$RNA@counts[gene_list[[tt]], ]!=0) >=thre)
        }
        return(df)
        })  %>%
            do.call(rbind, .)
    return(exp_data)
    }) %>% 
        setNames(., names(gene_list))


## Plot the expression [Barplots]
bar_list <- lapply(names(exp_list), function(tt) {
    exp_data <- exp_list[[tt]]
    maxy <- exp_data %>%
                group_by(species, cluster) %>%
                summarize(ratio = mean(coexp) * 100) %>%
                .$ratio %>% max()
    barp <- lapply(names(order_list), function(sp) {
        data <- subset(exp_data, species == sp) %>%
                    group_by(cluster) %>%
                    summarize(ncells = sum(coexp), ratio = mean(coexp) * 100) %>%
                    subset(ratio != 0)
        p2 <- ggplot(data) +
                geom_bar(aes(x = cluster, y = ratio), fill = "grey", alpha = 0.8, width = 1, stat = "identity") +
                geom_text(aes(x = cluster, y = ratio, label = ncells), size = 2, angle = 45, hjust = 0, vjust = 0) + 
                theme_cowplot() +
                RotatedAxis() +
                scale_x_discrete(limits = order_list[[sp]]) +
                scale_y_continuous(limits = c(0, maxy)) +
                theme(legend.position = "none", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
        if (sp != "Mouse"){
            p2 <- p2 + 
                theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
        return(p2)
    })
    return(barp)
    })



pdf(paste0(outputdir, "SF3.Conserved_coexp.bar.final.V2.pdf"), width = 8.5, height = 2)
pp <- bar_list[[1]] %>%
        lapply(., function(x) x + theme(axis.text.x = element_blank())) %>%
        patchwork::wrap_plots(., nrow = 1, ncol = 5)
qq <- patchwork::wrap_plots(bar_list[[2]], nrow = 1, ncol = 5)
patchwork::wrap_plots(pp, qq, nrow = 2, ncol = 1) %>% print()
dev.off()



































