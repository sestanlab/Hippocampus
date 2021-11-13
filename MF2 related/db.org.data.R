## Organize all the GC & the non-used cells 
source("../scripts/hip.fun.R")

fhdg <- readRDS(file = paste0("../MF2_neurogenesis/load_files/", "Human_DG_seu.rds"))


count1 <- readRDS(file = paste0(dataDir, "HIP.batch1st.count.rds"))
count1 <- count1[, grepl("_DG_", colnames(count1))]
count2 <- readRDS(file = paste0(dataDir, "HIP.batch2nd.human.count.rds"))

anno_new_cells <- lapply(c("HSB179", "HSB181", "HSB282"), function(sp) {
	cells <- readRDS(file = paste0("../DG/load_files/", "Integrate.hDG.", sp, ".slim.filtered.rds")) %>%
			subset(major == "new") %>%
			colnames()
	return(cells)
	})%>% 
		unlist()

all.raw <- cbind(count1, count2)

rmd_cells <- setdiff(colnames(all.raw), anno_new_cells) %>%
				setdiff(., colnames(fhdg))
gc_cells <- colnames(fhdg)[fhdg@meta.data$mres %in% c("Astro" , "GC")]

final.counts <- all.raw[, c(gc_cells, rmd_cells)]
dgdb <- seu_prepare(counts = final.counts, min.cells = 0, normalization.method = "LogNormalize", nfeatures = 2500, hvg.method = "vst", assay = "RNA")
dgdb@meta.data[gc_cells, "mres"] <- fhdg@meta.data[gc_cells, "mres"]
dgdb@meta.data$mres[is.na(dgdb@meta.data$mres)] <- "remove"
saveRDS(dgdb, file = paste0(inputdir, "hDG.GC.with.removed.cells.rds"))





