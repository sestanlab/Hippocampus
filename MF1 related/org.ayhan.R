source("~/project/HIP/scripts/hip.fun.R")


## Fast reading count matrix: assume the first column is the gene names
FastReadCounts <- function(file, showProgress = TRUE, header = TRUE, sep = "\t"){ 
    expmat <- data.table::fread(file, showProgress = showProgress, header = header, sep = sep)

    ## Save the gene names to another object and remove that column
    gene_names <- as.character(expmat[, 1][[1]])
    if (sum(duplicated(gene_names)) >= 1){
        stop("Gene names are duplicated, please check")
    }


    gcol <- colnames(expmat)[1]
    expmat[,(gcol):=NULL]


    ## Convert to Sparse marker
    to_sparse <- function(d_table){
      
      i_list <- lapply(d_table, function(x) which(x != 0))
      counts <- unlist(lapply(i_list, length), use.names = F)

      sparseMatrix(
        i = unlist(i_list, use.names = F),
        j = rep(1:ncol(d_table), counts),
        x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F),
        dims = dim(d_table),
        dimnames = list(NULL, names(d_table)))
    }

    X <- to_sparse(expmat)
    rownames(X) <- gene_names
    return(X)
}


## To seurat object
data <- FastReadCounts(file = "./load_files/GSE160189_Hippo_Counts.csv", showProgress = TRUE, header = TRUE, sep = ",")
seu <- seu_prepare(counts = data, min.cells = 0, normalization.method = "LogNormalize", nfeatures = 2500, hvg.method = NULL, assay = "RNA")


## Add meta data
meta <- read.table(file = "./load_files/meta.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
			column_to_rownames("Cell")


sh_cells <- intersect(rownames(meta), colnames(seu))
slim_seu <- seu[, sh_cells]
slim_seu@meta.data <- cbind(slim_seu@meta.data, meta[colnames(slim_seu), setdiff(colnames(meta), colnames(slim_seu@meta.data))])
slim_seu@meta.data$percent.mt <- meta[colnames(slim_seu), "percent.mt"]


umaps <- read.table(file = "./load_files/Seurat_umap.coords.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1) %>%
            setNames(., c("UMAP_1", "UMAP_2"))
slim_seu[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umaps)[colnames(slim_seu), ], key = "UMAP_", assay = "RNA")


saveRDS(slim_seu, file = paste0(inputdir, "Konopka_HIP_seu.rds"))










