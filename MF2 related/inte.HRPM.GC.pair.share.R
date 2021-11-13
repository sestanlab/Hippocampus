## Integrate Rhesus & Mouse dataset
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")
source("./inte.new.fun.R")


inte.method <- args[1] ##"seurat"
nfeatures <- 1500
pair <- c("Human.Mouse", "Human.Pig", "Human.Rhesus", "Rhesus.Pig", "Rhesus.Mouse", "Pig.Mouse")[as.numeric(args[2])]
file_name <- paste0("HRPM.pw.GC.", inte.method, ".expHVG.", nfeatures, ".", pair)


## Set parallel computint for seurat [have errors for harmony]
if (tolower(inte.method) == "seurat"){
	library(future)
	plan("multiprocess", workers = 4)
	options(future.globals.maxSize = 28*1000*1024^2)
}


## Integration
inte.hrm.file <- paste0(inputdir, file_name, ".slim.rds")
if (!file.exists(inte.hrm.file)){
	## Prepare the dataset
	mres_use <- c("Astro", "RGL", "nIPC", "NB", "GC")
	sp1 <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
	sp2 <- strsplit(pair, ".", fixed = TRUE)[[1]][2]


	## Get the data list
	data_list <- lapply(c(sp1, sp2), function(sp) {
		xx <- readRDS(file = paste0(inputdir, sp, "_DG_seu.rds"))
		if (sp == "Mouse"){
			xx <- FindVariableFeatures(xx, nfeatures = 2500)
			xx <- list(MMBold = xx)
		} else {
			xx <- SplitObject(xx, split.by = "inte.batch") %>%
					lapply(., function(yy) FindVariableFeatures(yy, nfeatures = 2500))
		}
		return(xx)
		}) 


	## Get the highly varaible genes
	hvg_list <- lapply(data_list, function(x) {
		if (length(x) > 1){
			sample_use <- grep("181|282|179", names(x), invert = TRUE, value = TRUE)
			print(sample_use)
			hh <- SelectIntegrationFeatures(object.list = x[sample_use], nfeatures = nfeatures)
		} else {
			hh <- FindVariableFeatures(x[[1]], nfeatures = nfeatures) %>%
					VariableFeatures()
		}
		return(hh)
		})
	hvg_use <- Reduce("union", hvg_list)


	## Get the genes to integration
	gc_list <- do.call(c, data_list)
	gene_use <- lapply(gc_list, rownames) %>%
					Reduce("intersect", .) %>%
					union(., hvg_use)


	## Use the shared genes 
	allgc <- merge(x = gc_list[[1]], y = gc_list[2:length(gc_list)])
	allgc <- subset(allgc, fig2cluster %in% mres_use)
	gc_list <- SplitObject(allgc[gene_use, ], split.by = "inte.batch")


	## Do the integration
	if (tolower(inte.method) == "seurat"){
		refs <- switch(pair, Human.Mouse = c("MMBold", "HSB237"), 
							Human.Pig = c("p41_h0", "HSB237"), 
							Human.Rhesus = c("RMB2_1", "HSB237"))
		seu <- Integratelist.seurat(obj.list = gc_list, hvg = hvg_use, file_name = file_name, input_dir = inputdir, inte.dims = 1:25, cluster.dims = 1:25, , reference = which(names(gc_list) %in% refs))
	} else if (tolower(inte.method) == "harmony"){
		seu <- Integratelist.harmony(obj.list = gc_list, split.by = "inte.batch", hvg = hvg_use, file_name = file_name, input_dir = inputdir, inte.dims = 1:25, theta = 2, lambda = 0.9, sigma = 0.1) 
	} else {
		stop("Unrecognize integration method")
	}
	


	## Plot some figures
	DimFig(seu, group.by = c("species", "fig2cluster", "seurat_clusters"), file_name = file_name, plot.scale = 0.7)
	DimFig(seu, group.by = c("fig2cluster"), file_name = file_name, plot.scale = 0.7, split.by = "species", split.order = c(sp1, sp2))
	FeatureFig(seu, features = c("MKI67", "TOP2A", "CENPF", "DCX", "PROX1", "SOX11", "CALB2", "DPYSL3"), file_name = paste0(file_name, ".exp.ANG"), plot.scale = 0.6, split.by = "species", split.order = c(sp1, sp2))
}



