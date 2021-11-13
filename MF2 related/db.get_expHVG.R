## Get the highly variable genes for all DG samples
source("../scripts/hip.fun.R")


## Human DG HVG
mdg <- readRDS(file = paste0("../MF2_final/load_files/", "Mouse_DG_full_seu.rds")) %>%
			SplitObject(., split.by = "inte.batch") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500))
hdg <- readRDS(file = paste0("../MF2_final/load_files/", "Human_DG_seu.rds")) %>%
			SplitObject(., split.by = "inte.batch") %>%
			.[c("HSB231", "HSB237", "HSB628")] %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500))


dglist <- list(Mouse = mdg, Human = hdg)
hm_share <- intersect(rownames(mdg[[1]]), rownames(hdg[[1]]))
print(length(hm_share))
print(names(dglist$Human))



## Union of H-R HVG genes
for (nfeatures in c(1200, 1500, 1750, 2000)[2:3]){
	hvglist <- lapply(dglist, function(x) {
		if (class(x) %in% "list"){
			hh <- SelectIntegrationFeatures(object.list = x, nfeatures = nfeatures)
		} else {
			hh <- FindVariableFeatures(x, nfeatures = nfeatures) %>%
					VariableFeatures()
		}
		return(hh)
		})
	hm_hvg <- Reduce('union', hvglist)
	gene_use <- union(hm_share, hm_hvg)
	save(gene_use, hm_hvg, file = paste0(inputdir, "HVG.HM.expHVG.", nfeatures, ".Rdata"))
}












