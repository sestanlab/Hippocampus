## Integrate Human & Mouse dataset 
source("../scripts/hip.fun.R")
source("../MF2_final/inte.new.fun.R")
library(future)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 30*1000*1024^2)





## Integration
inte.hrm.file <- paste0(inputdir, "HM.Inte.GCDB.v2.slim.rds")
if (!file.exists(inte.hrm.file)){
	## Filtered GC lineage in each species
	mres_use <- c("Astro", "RGL", "nIPC", "NB", "GC")
	## Have consistent data version with MF2_final figures
	hdg_db <- readRDS(file = paste0(inputdir, "HDG.inte.GCdoublets.seurat.slim.F1.rds"))
	hdg_db <- hdg_db[, !hdg_db@meta.data$mres %in% c("GC", "Astro")]
	hdg_db@meta.data$fig2cluster <- hdg_db@meta.data$cluster
	hdg_gc <- readRDS(file = paste0("../MF2_final/load_files/", "Human_DG_seu.rds")) %>%
					subset(fig2cluster %in% c("GC", "Astro"))
	hdg <- merge(hdg_gc, hdg_db)
	hdg@meta.data$inte.batch <- extract_field(colnames(hdg), 1, "_")
	
	
	mdg_gc <- readRDS(file = paste0("../MF2_final/load_files/", "Mouse_DG_seu.rds")) %>%
					subset(mres %in% mres_use)


	## Load HVG
	load(file = paste0(inputdir, "HVG.HM.expHVG.1750.Rdata"))


	## Use the shared genes
	allgc <- merge(x = mdg_gc, y = hdg)
	allgc <- allgc[gene_use, ]
	gclist <- SplitObject(allgc, split.by = "inte.batch")
	rm(allgc)


	## Do the integration
	file_name <- "HM.Inte.GCDB.v2"
	seu <- Integratelist.seurat(obj.list = gclist, hvg = hm_hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:25, cluster.dims = 1:25, reference = which(names(gclist) %in% c("MMBold", "HSB237")))
	DimFig(seu, group.by = c("species", "vis.batch", "fig2cluster", "seurat_clusters"), file_name = file_name, plot.scale = 0.7)
	DimFig(seu, group.by = c("fig2cluster"), file_name = file_name, plot.scale = 0.7, split.by = "vis.batch", split.order = c("Mouse middle", "Human"))
}




