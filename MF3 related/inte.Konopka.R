## Integrate Rhesus & Mouse dataset
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")
source("../MF2_final/inte.new.fun.R")
library(future)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 80*1000*1024^2)


## Subset only the neurons & Astrocytes for integration
all_cls <- c("Astro1", "Astro2", "Astro3","Den.Gyr1", "Den.Gyr2", "Den.Gyr3", "Endo", "In1", "In2", "In3", "Micro1", "Micro2", "Micro3", "Olig1", "Olig2", "Olig3", "Olig4", "Olig5", "OPC1", "OPC2", "OPC3", "OPC4", "Pyr1", "Pyr2")
sel_cls <- c("Astro1", "Astro2", "Astro3","Den.Gyr1", "Den.Gyr2", "Den.Gyr3", "Olig1", "Olig2", "Olig3", "Olig4", "Olig5", "OPC1", "OPC2", "OPC3", "OPC4")


allgc <- readRDS(file = paste0(inputdir, "Konopka_HIP_seu.rds"))
allgc <- allgc[, !allgc@meta.data$Cluster %in% c("Micro1", "Micro2", "Micro3", "Endo")]


## Find HVG (slim the Oligo, otherwise there are too many Oligos)
set.seed(42)
sel_cells <- c(colnames(allgc)[grep("Olig", allgc@meta.data$Cluster, invert = TRUE)], 
				sample(colnames(allgc)[grep("Olig", allgc@meta.data$Cluster)], 5000))


hvg <- SplitObject(allgc[, sel_cells], split.by = "batch") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(object.list = ., nfeatures = 2500)
gclist <- SplitObject(allgc, split.by = "batch") ##%>%
			##lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500))


file_name <- paste0("Konopka.Inte.Seurat.2000.v2")
newseu <- Integratelist.seurat(obj.list = gclist, hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30)


## Plot some figures
DimFig(newseu, group.by = c("group", "seurat_clusters", "Cluster", "batch"), file_name = file_name, plot.scale = 0.8)

## Expression of selected genes
##newseu <- readRDS(file = paste0(inputdir, "Konopka.Inte.Seurat.2000.slim.rds"))
FeatureFig(newseu, features = c("LPAR1", "CALB2", "DCX", "PROX1", "SLC1A3", "AQP4", "SLC17A7"), file_name = paste0(file_name, "_sel.exp"), plot.scale = 0.6, ncol = 3)
FeatureFig(newseu, features = c("MOBP", "MBP", "PLP1", "SOX6", "PDGFRA", "OLIG2"), file_name = paste0(file_name, "_sel.exp.v2"), plot.scale = 0.6, ncol = 3)





