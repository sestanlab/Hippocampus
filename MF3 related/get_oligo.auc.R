source("../scripts/hip.fun.R")


allgc <- readRDS(file = paste0(inputdir, "Konopka_HIP_seu.rds"))
mar <- readRDS(file = paste0(inputdir, "Konopka_Oligo_mar.rds"))


### Get the AUROC scores of the Oligo markers
file_name <- "Konopka_Oligo"
auc_res <- GetModuleScore(assay.data = allgc$RNA@data, features = list(Oligo = mar[1:50]), seed = 42, method = "aucell", input_dir = inputdir, file_name = file_name, rethreshold_list = NULL, cellbin.size = 8000)$auc[colnames(allgc), ,drop = FALSE]

saveRDS(auc_res, file = paste0(inputdir, "Konopka_Oligo_auc_org.rds"))



##auc_res <- readRDS(file = paste0(inputdir, "Konopka_Oligo_auc_org.rds"))
##aa@meta.data$Oligo <- auc_res[colnames(aa), "Oligo"]

##FeatureFig(aa, features = "Oligo", file_name = paste0(file_name, "_Oligo_AUC"), plot.scale = 0.7)


##meta <- readRDS(file = paste0(inputdir, "Konopka_all_db_res.rds"))
##aa@meta.data <- cbind(aa@meta.data, meta[colnames(aa), ,drop = FALSE])
##FeatureFig(aa, features = c("doublet_scores", "doublet_assign"), file_name = paste0(file_name, "_doublets"), plot.scale = 0.7)










