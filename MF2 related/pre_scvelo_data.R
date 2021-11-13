## conda activate R4
## R

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)


inputdir <- "./load_files/"
loomDir <- "../MF2_velocity/load_files/"
all_sps <- c("Human", "Rhesus", "Pig", "Mouse middle", "Mouse young")
hrpm_gc <- readRDS("/home/sm2726/project/HIP/MF2_final/load_files/HRPM.GC.seurat.expHVG.1500.slim.rds")
for (spe in all_sps[1:4]){
	subgc <- hrpm_gc %>%
				subset(vis.batch == spe)
	spe <- spe %>%
				gsub(" ", "-", .)
	cbn_ldat <- readRDS(file = paste0(loomDir, "Loom.rawmat.", spe, ".rds"))
	slim_ldat <- lapply(cbn_ldat, function(x) x[, colnames(subgc)])


	seu <- as.Seurat(x = slim_ldat)
	seu[["RNA"]] <- seu[["spliced"]]
	seu[["pca"]] <- CreateDimReducObject(embeddings = subgc$pca@cell.embeddings, key = "PC_", assay = "RNA")
	seu[["umap"]] <- CreateDimReducObject(embeddings = subgc$umap@cell.embeddings, key = "UMAP_", assay = "RNA")
	seu@meta.data$cluster <- subgc@meta.data$fig2cluster
	Idents(seu) <- "cluster"


	DefaultAssay(seu) <- "RNA"
	SaveH5Seurat(seu, filename = paste0(inputdir, "velo.", spe, ".h5Seurat"), overwrite = TRUE)
	Convert(paste0(inputdir, "velo.", spe, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)

}



















