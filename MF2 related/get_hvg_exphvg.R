## Get the highly variable genes for all DG samples
source("../scripts/hip.fun.R")


## Human DG HVG
allsps <- c("Mouse", "Pig", "Rhesus", "Human")
dglist <- lapply(allsps[2:4], function(sp) {
	xx <- readRDS(file = paste0(inputdir, sp, "_DG_seu.rds")) %>%
			SplitObject(., split.by = "inte.batch") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500))
	xx
	}) %>% 
		setNames(., allsps[2:4])
dglist$Mouse <- readRDS(file = paste0(inputdir, "Mouse_DG_seu.rds"))


hrpm_share <- lapply(dglist[allsps[2:4]], function(x) rownames(x[[1]])) %>%
				Reduce("intersect", .) %>%
				intersect(., rownames(dglist$Mouse))
print(length(hrpm_share))

dglist$Human <- dglist$Human[grep("231|237|628", names(dglist$Human))]
print(names(dglist$Human))


## Union of H-R HVG genes
nfeatures <- 1500
hvglist <- lapply(dglist, function(x) {
		if (class(x) %in% "list"){
			hh <- SelectIntegrationFeatures(object.list = x, nfeatures = nfeatures)
		} else {
			hh <- FindVariableFeatures(x, nfeatures = nfeatures) %>%
					VariableFeatures()
		}
		return(hh)
		})
hrpm_hvg <- Reduce('union', hvglist)
gene_use <- union(hrpm_share, hrpm_hvg)
save(gene_use, hrpm_hvg, file = paste0(inputdir, "HVG.HRPM.expHVG.", nfeatures, ".Rdata"))













