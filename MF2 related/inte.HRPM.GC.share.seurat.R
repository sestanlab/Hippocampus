## Integrate Rhesus & Mouse dataset
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")
source("./inte.new.fun.R")
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 20*1000*1024^2)


file_name <- paste0("HRPM.GC.seurat.expHVG.1500")

## Integration
inte.hrm.file <- paste0(inputdir, file_name, ".slim.rds")
if (!file.exists(inte.hrm.file)){
	## Prepare the dataset
	mres_use <- c("Astro", "RGL", "nIPC", "NB", "GC")
	allsps <- c("Mouse", "Pig", "Rhesus", "Human")
	data_list <- lapply(allsps[2:4], function(sp) {
		xx <- readRDS(file = paste0(inputdir, sp, "_DG_seu.rds")) %>%
					subset(fig2cluster %in% mres_use)
		xx
		}) %>%
			setNames(., allsps[2:4])
	data_list$Mouse <- readRDS(file = paste0(inputdir, "Mouse_DG_full_seu.rds")) %>%
					subset(fig2cluster %in% mres_use)


	## Load the HVG & gene use
	load(file = paste0(inputdir, "HVG.HRPM.expHVG.1500.Rdata"))


	## Use the shared genes
	allgc <- merge(x = data_list[[1]], y = data_list[2:4])
	allgc <- allgc[gene_use, ]
	gclist <- SplitObject(allgc, split.by = "inte.batch")
	rm(allgc)


	## Do the integration
	seu <- Integratelist.seurat(obj.list = gclist, hvg = hrpm_hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:22, cluster.dims = 1:22, reference = which(names(gclist) %in% c("MMByouth", "HSB231")))
}



