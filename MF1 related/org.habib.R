## Read the mouse dataset 
source("../scripts/hip.fun.R")


## Prepare expression data
df <- read.table(file = paste0(inputdir, "GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

exp_data <- df[, 2:ncol(df)]
rownames(exp_data) <- df[, 1]
exp_data <- Matrix(as.matrix(exp_data), sparse = TRUE)



## Add cell annotation
cls_df <- read.table(file = paste0(inputdir, "GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.clusters.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
			column_to_rownames("V1")

seu <- seu_prepare(counts = exp_data, min.cells = 0, normalization.method = "LogNormalize", nfeatures = 1500, hvg.method = "vst", assay = "RNA")
seu@meta.data$seurat_clusters <- as.character(cls_df[colnames(seu), "V2"])
seu@meta.data$rawcluster <- get_ident(input_ident = seu, ident_col = "seurat_clusters", file_path = paste0(inputdir, "habib.anno.txt"), sample_name = "Habib", label_names = NULL)
seu@meta.data$cluster <- get_ident(input_ident = seu, ident_col = "seurat_clusters", file_path = paste0(inputdir, "habib.anno.txt"), sample_name = "Habib", label_names = "hres")



## Add tSNE
tsne_df <- read.table(file = paste0(inputdir, "GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.tsne.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
			column_to_rownames("X")
seu[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne_df[colnames(seu), ]), key = "tSNE_", assay = "RNA")



hip_sps <- c("hHP1","hHP2","hHP2a","hHP2b","hHP2c","hHP3b","HP2-A","HP2-B","HP3-A","HP3-B","PFC2-A1","PFC2-A2","PFC2-A3","PFC2-A5") ## The so called "PFC2-XX" were considered as HIP is that they have DG (cluster 7), CA (cluster 3, 4) and NSC (cluster 14) and have limited PFC EXN (cluster 1)
seu@meta.data$region <- ifelse(seu@meta.data$orig.ident %in% hip_sps, "HIP", "PFC")


saveRDS(seu, file = paste0(inputdir, "Habib_Human_DFC_HIP_seu.rds"))




