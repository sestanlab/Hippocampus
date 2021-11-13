## Integrate adult and fetal DG dataset
source("../scripts/hip.fun.R")
source("./inte.new.fun.R")
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 29*1000*1024^2)



## Load the dataset and find the Variable genes
mres_use <- c("Astro", "RGL", "nIPC", "NB", "GC")
hadult <- readRDS(file = paste0("../MF2_neurogenesis/load_files/", "Human_DG_seu.rds")) %>%
				subset(., mres %in% mres_use)
hadult@meta.data$species <- "HumanAdult"


hfetal <- readRDS(file = paste0(inputdir, "Human_fetal_DG.subset.rds"))
hfetal@meta.data$samplename <- "HSB1"
hfetal@meta.data$species <- "HumanFetal"
hfetal@meta.data$cluster <- hfetal@meta.data$label



## Get the highly variable genes
sharegenes <- intersect(rownames(hadult), rownames(hfetal))
adult_list <- SplitObject(hadult[sharegenes, ], split.by = "samplename")
file_name <- paste0("Human.AdultFetal.inte.GC.seurat.", "2000") 


hvg <- readRDS(file = paste0(inputdir, "Human.AF.hvg.", "2000", ".rds"))

seu <- Integratelist.seurat(obj.list = c(adult_list, list(HSB1= hfetal[sharegenes, ])), hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:30)
DimFig(seu, group.by = c("mres", "species", "seurat_clusters"), file_name = file_name, plot.scale = 0.7)
DimFig(seu, group.by = "mres", file_name = file_name, plot.scale = 0.7, split.by = "species", split.order = c("HumanAdult", "HumanFetal"))












