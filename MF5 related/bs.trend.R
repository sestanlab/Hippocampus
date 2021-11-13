## For the HIP-specific genes, further filter ones that are elevated along the time-group
source("../scripts/hip.fun.R")
source("~/project/common_scripts/private_R_pack/plot_brainspan.R")
library(ggpubr)
library(lemon)
library(limma)


regmatch_list <- list(HIP=c("HIP"), AMY=c("AMY"), STR=c("STR", "VF", "MGE","CGE","LGE"), MD=c("DIE","DTH","MD"), CBC=c("URL","CBC"), NCX = c("FC", "OFC", "DFC", "VFC", "MFC", "M1C", "M1CS1C", "MSC", "PC", "S1C", "IPC", "TC", "A1C", "STC", "ITC", "OC", "V1C") %>% unique())
regmerge <- rep(names(regmatch_list), sapply(regmatch_list, length)) %>% 
				setNames(., unlist(regmatch_list, use.names = FALSE))


metadata$timegroup <- case_when(
			metadata$Period <= 7 ~ "prenatal", 
			metadata$Period > 7 & metadata$Period <= 12 ~ "postnatal",
			metadata$Period >= 13 ~ "adult")
metadata$region <- regmerge[metadata$Regioncode]



time <- factor(metadata$timegroup, levels = c("prenatal", "postnatal", "adult"))
region <- factor(metadata$region, levels = c("AMY", "CBC", "HIP", "MD", "NCX", "STR"))

design <- model.matrix(~ region + region:time)
fit <- lmFit(log2(bs_raw + 1), design)
fit <- eBayes(fit)
res <- topTable(fit, number = 1000000)


hip_genes <- readRDS(file = paste0(inputdir, "HIP-specific-genes.rds"))
match_genes <- map_gene(gene_names = rownames(bs_raw), input_genes = hip_genes)
gene_use <- match_genes

all_regions <- c("AMY", "CBC", "HIP", "MD", "NCX", "STR")

meanCOF <- lapply(all_regions, function(reg) 
	apply(res[gene_use, paste0(paste0("region", reg), c(".timepostnatal", ".timeadult"))], 1, mean)
	) %>% setNames(., all_regions) %>%
			as.data.frame()
meanCOF$diff <- meanCOF$HIP - apply(meanCOF[, c("AMY", "CBC", "MD", "NCX", "STR")], 1, max)
meanCOF <- meanCOF %>%
			rownames_to_column("gene") %>%
			arrange(desc(diff)) %>% 
			mutate(gene = extract_field(gene, 2, "|")) %>%
			column_to_rownames("gene") %>%
			as.matrix()



## Plot the bar plots
color_vec <- c("#2166ac", "#80cdc1", "#fdb863", "#4d4d4d", "#66bd63", "#d73027") %>%
				setNames(., c("NCX", "HIP", "AMY", "STR", "MD", "CBC"))
plot_data <- meanCOF[, 1:6] %>% 
			reshape2::melt() %>%
			setNames(., c("gene", "region", "COF")) %>%
			mutate(region = factor(as.character(region), levels = names(color_vec))) %>%
			mutate(gene = factor(as.character(gene), levels = rownames(meanCOF)))
p <- ggplot(plot_data) +
			geom_bar(aes_string(x = "gene", y = "COF", fill = "region"), color = NA, position = "dodge2", stat = "identity", width = 0.75, size = 1) +
			theme_cowplot() +
			coord_capped_cart(left='both') +
			RotatedAxis() +
			scale_fill_manual(values = color_vec) +
			labs(y = "time group-region coefficient") +
			theme(axis.text.y = element_text(size = 7), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 7)) +
			theme(axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25)) 

dif_data <- meanCOF[, "diff", drop = FALSE] %>%
			as.data.frame() %>%
			rownames_to_column("gene") %>%
			mutate(gene = factor(as.character(gene), levels = rownames(meanCOF)))
q <- ggplot(dif_data) + 
		geom_bar(aes_string(x = "gene", y = "diff"), fill = "lightgrey", color = "black", position = "dodge2", stat = "identity", width = 0.75, size = 0.15) + 
		coord_capped_cart(bottom='both', left='both') +
		theme_cowplot() +
		RotatedAxis() +
		labs(y = "coef diff [HIP - max(others)]") +
		theme(axis.text.y = element_text(size = 7), axis.title.x = element_blank(), axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5), axis.title.y = element_text(size = 7)) +
		theme(axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25)) 
par(mai = c(0,0,0,0), oma = c(0, 0,0,0))
ggarrange(p, q, nrow = 2, ncol = 1, heights = c(1, 1.2), common.legend = TRUE, legend = "right") %>%
		ggexport(., filename = paste0(outputdir, "MF4.HIP.ExN.genes.BrainSpan.coef.pdf"), width = 5.5, height = 2.5)







