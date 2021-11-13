## Organize the Developmental NHP bulktissue RNAseq dtat (Zhu et al, 2017, Science)
source("../scripts/hip.fun.R") 


rpkm <- read.table(file = paste0(inputdir, "nhp_development_RPKM_rmTechRep.txt"), sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE) %>%
			as.matrix()
meta <- read.table(file = paste0(inputdir, "evo_meta.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE) %>%
		column_to_rownames("Sample")
logrpkm <- log2(rpkm+1)


save(rpkm, logrpkm, meta, file = paste0(inputdir, "NHP.dev.Rdata"))


