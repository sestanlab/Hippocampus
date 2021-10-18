source("~/project/HIP/scripts/hip.fun.R")


hip <- readRDS("../data/HIPEC.final.seu.01172021.rds")


hip@meta.data$newregion <- gsub("eDG", "DG", hip@meta.data$region)
hip@meta.data$group[hip@meta.data$fig1cluster == "CR RELN NDNF"] <- "CR"

##------------------------------------------------------------------
## Plot the cell type proportions
gp_colors <- setNames(c("#C51B7D", "#4D9221", "#999999"), c("ExN", "InN", "NNC"))
size_data <- hip@meta.data %>%
              group_by(newregion) %>%
              mutate(newsp = as.character(as.numeric(as.factor(samplename)))) %>%
              ungroup() %>%
              group_by(newregion, group, newsp) %>%
              summarize(ncells = n()) %>%
              group_by(newregion) %>%
              mutate(ratio = ncells * 100/sum(ncells))

pp <- ggplot(size_data, aes_string(x = "newregion", y = "ratio", fill = "group")) +
              geom_bar(width = 0.8, size = 0.2, color = "black", stat = "identity", position="fill") +
              theme_bw() +
              RotatedAxis() +
              scale_fill_manual(values = gp_colors) +
              scale_x_discrete(limits = c("DG", "CA24", "CA1", "SUB", "EC")) +
              theme(legend.position = "right", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf(paste0(outputdir, "UMAP.HIP.proportion.pdf"), width = 5, height = 3)
print(pp)
dev.off()




ratio_data <- size_data %>%
        ungroup() %>%
        group_by(newregion, newsp) %>%
        mutate(ratio = ncells/sum(ncells)) %>%
        ungroup() %>%
        group_by(newregion, group) %>%
        summarize(mratio = mean(ratio), rse = sd(ratio)/sqrt(n())) %>%
        ungroup() %>%
        mutate(group = factor(as.character(group), levels = c("ExN", "InN", "NNC")))


qq <- ggplot(ratio_data, aes_string(x = "newregion", y = "mratio", fill = "group")) +
              geom_bar(size = 0.2, stat="identity", color="black", position=position_dodge()) +
              geom_errorbar(aes(ymin=mratio-rse, ymax=mratio+rse), width=.3, position=position_dodge(.9), size = 0.4) +
              theme_bw() +
              RotatedAxis() +
              scale_fill_manual(values = gp_colors) +
              scale_x_discrete(limits = c("DG", "CA24", "CA1", "SUB", "EC")) +
              theme(legend.position = "right", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.2))
pdf(paste0(outputdir, "UMAP.HIP.proportion.v2.pdf"), width = 5, height = 3)
print(qq) 
dev.off()

