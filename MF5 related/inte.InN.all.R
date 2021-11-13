## Organize the dataset 
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")


##-----------------------------------------------------------------------
## Organize the dataset for integration
rawFile <- paste0(inputdir, "FourReg.InN.merged.raw.Rdata")
if (!file.exists(rawFile)){
	inn_4reg <- readRDS(file = paste0(dataDir, "FourReg.InN.final.01172021.rds"))
	inn_f4reg <- lapply(names(inn_4reg), function(reg) {
		yy <- inn_4reg[[reg]]
		yy@meta.data$merge_region <- reg
		yy@meta.data$inte_region <- ifelse(reg %in% c("HIP", "EC"), "HIPEC", reg)
		yy@meta.data$dataset <- paste0(yy@meta.data$samplename, "_", yy@meta.data$inte_region)
		return(yy)
		}) %>%
			setNames(., names(inn_4reg))
	rm(inn_4reg)


	## Subset to the shared genes
	share_genes <- Reduce("intersect", lapply(inn_f4reg, function(x) rownames(x)))
	seulist <- list(HIPEC = merge(x = inn_f4reg$HIP[share_genes, ], y = inn_f4reg$EC[share_genes, ])) %>%
					c(., list(MTG = inn_f4reg$MTG[share_genes, ], PFC = inn_f4reg$PFC[share_genes, ]))


	## Get the HVG
	hvglist <- lapply(names(seulist), function(reg) {
		if (reg == "MTG"){ ## Only one sample
			subhvg <- FindVariableFeatures(seulist[[reg]], nfeatures = 1500) %>%
						VariableFeatures() 
		} else if (reg == "HIPEC"){
			sublist <- SplitObject(seulist[[reg]], split.by = "dataset") %>%
								lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000))
			subhvg <- SelectIntegrationFeatures(object.list = sublist[grep("HSB231|HSB237|HSB628", names(sublist))], nfeatures = 1500)
		} else {
			sublist <- SplitObject(seulist[[reg]], split.by = "dataset") %>%
								lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000))
			subhvg <- SelectIntegrationFeatures(object.list = sublist, nfeatures = 1500)
		}
		
		message(paste0("Finish HVG selection for region: ", reg))
		return(subhvg)
		}) %>%
				setNames(., names(seulist))
	hvg <- unlist(hvglist) %>% unique()


	alld <- merge(x = seulist[[1]], y = seulist[2:length(seulist)])
	print(levels(as.factor(alld@meta.data$dataset)))
	save(alld, hvg, file = paste0(inputdir, "FourReg.InN.merged.raw.Rdata"))
}




library(harmony) 
Integrate.4reg.harmony <- function(object, split.by, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, theta = 2, lambda = 1, sigma = 0.1) {
	## Get the seurat object list and set the split.by="inteinfo"
	obj.list <- SplitObject(object, split.by = split.by)

    inte.slim.file <- paste0(input_dir, file_name, ".slim.rds")
    if (!file.exists(inte.slim.file)){
        ## Do the harmony integration
        object <- ScaleData(object, split.by = split.by, do.center = FALSE, features = hvg)%>%
                    RunPCA(., features = hvg, verbose = FALSE) %>%
                    RunHarmony(., group.by.vars = split.by, lambda = lambda, theta = theta, dims.use = inte.dims, sigma = sigma)
        object <- RunUMAP(object, dims = 1:35, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
        object <- FindNeighbors(object, dims =1:35, reduction = "harmony", k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        saveRDS(object, file = inte.slim.file)
    }else {
        object <- readRDS(file = inte.slim.file)
    }
    return(object)
}
        




load(file = paste0(inputdir, "FourReg.InN.merged.raw.Rdata"))

inte.method <- args[1]
file_name <- paste0("Inte.InN.4reg.", inte.method)
seu <- Integrate.4reg.harmony(object = alld, split.by = "dataset", hvg = hvg, file_name= file_name,  input_dir = inputdir, inte.dims = 1:30, theta = 2, lambda = 0.75, sigma = 0.1)


DimFig(seu, group.by = c("inte_region", "merge_region", "seurat_clusters"), file_name = file_name, plot.scale = 1)
DimFig(seu, group.by = "subtypes", file_name = file_name, split.by = "merge_region", split.order = c("HIP", "EC", "MTG", "PFC"), plot.scale = 0.7)
 



