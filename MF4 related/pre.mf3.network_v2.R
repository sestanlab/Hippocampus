## Get the marker similarity among the clusters across regions
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")
source("./network.fun.R")
library(igraph)


##------------------------------------------------------------------
## Markers organization
if (FALSE){
for (ctp in c("ExN", "InN", "NNC")){
	marDir <- paste0(inputdir, "rawMar.hres/")
	marFiles <- list.files(marDir) %>%
					grep("^Markers", ., value = TRUE) %>%
					grep(ctp, ., value = TRUE)

	all.mar <- lapply(marFiles, function(x) {
		reg <- strsplit(x, ".", fixed = TRUE)[[1]][3]
		mar <- readRDS(paste0(marDir, x)) %>%
					mutate(bg = "others") %>%
					mutate(region = reg, group = ctp) %>%
					mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01))
		message(paste0("Finish file: ", x))
		return(mar)
		}) %>% do.call(rbind, .)
	saveRDS(all.mar, file = paste0(inputdir, "Markers.4reg.", ctp, ".rds"))
}
}



##------------------------------------------------------------------
## Cluster similarity based on the markers

## Extract Markers
ctp <- args[1]

ctp <- c("ExN", "InN", "NNC")[3]
all.mar <- readRDS(file = paste0(inputdir, "Markers.4reg.", ctp, ".rds"))
all.reg <- c("HIP", "EC", "MTG", "PFC")


mars <- lapply(all.reg, function(reg) {
	sub.mar <- all.mar %>% 
				subset(region == reg)
	all.cls <- levels(as.factor(sub.mar$cluster))
	reg.mars <- lapply(all.cls, function(x) {
		xx <- sub.mar %>%
				subset(cluster == x) %>%
				subset(pct.1 >= 0.2 & avg_logFC >= log(2) & ratio_fc >= 1.25) %>%
				##group_by(cluster) %>%
				mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
				subset(p_val_adj <= 0.01)
		##yy <- setNames(xx$avg_logFC, xx$gene)
		yy <- setNames(xx$ratio_fc, xx$gene)
		return(yy)
		}) %>% setNames(., paste0(reg, "|", all.cls))
	return(reg.mars)
	}) %>% 
		do.call(c, .)


## No need to remove B cells cluster, as it will co-cluster with others
##if (ctp == "NNC"){
	##mars <- mars[setdiff(names(mars), "PFC|B EBF1 IGKC")]
##}
q_thre <- switch(ctp, ExN = 0.935, InN = 0.95, NNC = 0.935)
plot.simi.network(marker.list = mars, thre.quantile = q_thre, seed.use = 42, max.marker = 100, dot.size = 8, file_name = paste0("MF3.v2.marker.simi.network.", ctp))









##mars[["PFC|L2-3 CUX2 ARHGAP18"]] %>% length() 









