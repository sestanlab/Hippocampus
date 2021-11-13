## Integrate adult and fetal DG dataset
source("../scripts/hip.fun.R")
source("./inte.new.fun.R")



## Load the dataset and find the Variable genes
hadult <- readRDS(file = paste0("../MF2_neurogenesis/load_files/", "Human_DG_seu.rds"))
hadult@meta.data$mres[hadult@meta.data$mres %in% c("EC", "SUB")] <- "CA"
hadult@meta.data$species <- "HumanAdult"


hfetal <- readRDS(file = paste0(inputdir, "Human_fetal_DG.rds"))
hfetal@meta.data$samplename <- "HSB1"
hfetal@meta.data$species <- "HumanFetal"
hfetal@meta.data$cluster <- hfetal@meta.data$label



## Get the highly variable genes
nfeatures <- 2000
sharegenes <- intersect(rownames(hadult), rownames(hfetal))
adult_list <- SplitObject(hadult[sharegenes, ], split.by = "samplename") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 3000))
hvg_adult <- SelectIntegrationFeatures(object.list = adult_list[grep("HSB231|HSB237|HSB628", names(adult_list))], nfeatures = nfeatures)
hvg_fetal <- FindVariableFeatures(hfetal[sharegenes, ], nfeatures = nfeatures) %>% VariableFeatures() 
hvg <- union(hvg_adult, hvg_fetal)
saveRDS(hvg, file = paste0(inputdir, "Human.AF.hvg.", nfeatures, ".rds"))







