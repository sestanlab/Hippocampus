Integratelist.seurat <- function(obj.list, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = TRUE) {
    if (length(obj.list) == 2){
        newseu <- merge(x = obj.list[[1]], y = obj.list[[2]])
    } else {
        newseu <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
    }


    inte.slim.file <- paste0(input_dir, file_name, ".slim.rds")
    if (!file.exists(inte.slim.file)){
        dg.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = inte.dims, assay = NULL, anchor.features = hvg, reference = reference)
        seuinte <- IntegrateData(anchorset = dg.anchors, dims = inte.dims)
        DefaultAssay(seuinte) <- "integrated"
        seuinte <- ScaleData(seuinte, verbose = FALSE) %>%
                                RunPCA(., npcs = 50, verbose = FALSE)
        
        newseu[["pca"]] <- CreateDimReducObject(embeddings = seuinte$pca@cell.embeddings[colnames(newseu), ], loadings = seuinte$pca@feature.loadings, stdev = seuinte$pca@stdev, key = "PC_", assay = "RNA")
        rm(seuinte)
        newseu <- RunUMAP(newseu, dims = cluster.dims, umap.method = "umap-learn", metric = "correlation")
        if (do.cluster){
            newseu <- FindNeighbors(newseu, dims = cluster.dims,  k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        } else {
            newseu@meta.data$seurat_clusters <- "empty"
        }
        saveRDS(newseu, file = inte.slim.file)
    }else {
        newseu <- readRDS(file = inte.slim.file)
    }
    return(newseu)
}




library(harmony) 
Integratelist.harmony <- function(obj.list, split.by, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, theta = 2, lambda = 1, sigma = 0.1) {
    if (length(obj.list) == 2){
        newseu <- merge(x = obj.list[[1]], y = obj.list[[2]])
    } else {
        newseu <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
    }
    rm(obj.list)

    inte.slim.file <- paste0(input_dir, file_name, ".slim.rds")
    if (!file.exists(inte.slim.file)){
        ## Do the LIGER integration
        ## Do the LIGER integration
        newseu <- ScaleData(newseu, split.by = split.by, do.center = FALSE, features = hvg)%>%
                    RunPCA(., features = hvg, verbose = FALSE) %>%
                    RunHarmony(., group.by.vars = split.by, lambda = lambda, theta = theta, dims.use = inte.dims, sigma = sigma)
        newseu <- RunUMAP(newseu, dims = 1:30, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
        newseu <- FindNeighbors(newseu, dims = 1:30, reduction = "harmony", k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        saveRDS(newseu, file = inte.slim.file)
    }else {
        newseu <- readRDS(file = inte.slim.file)
    }
    return(newseu)
}




plot.mvp <- function(object, file_name, output_dir = outputdir) {
    plot_data <- object$RNA@meta.features[, c("mvp.mean", "mvp.dispersion.scaled", "mvp.variable")] %>%
                    rownames_to_column("gene") %>%
                    mutate(color = ifelse(mvp.variable, "red", "grey"))

    top_var <- plot_data$gene[order(plot_data$mvp.dispersion.scaled, decreasing = TRUE)[1:50]]
    plot_data$label <- ifelse(plot_data$gene %in% top_var, plot_data$gene, "")
    p <- ggplot(plot_data, aes_string(x = "mvp.mean", y = "mvp.dispersion.scaled", color = "color")) +
            geom_point(size = 0.5) + 
            geom_text(aes_string(label = "label"), color = "black", size = 3)+
            scale_color_identity() + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 0.2), axis.ticks = element_line(size = 0.2), axis.line = element_line(size = 0.2))

    jpeg(paste0(output_dir, file_name, ".mvp.jpeg"), width = 8, height = 8, unit = "in", res = 300)
    print(p)
    dev.off()
}


FindVariableFeatures.mvp <- function(object, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf), file_name = "aa", output_dir = outputdir){
    object <- FindVariableFeatures(object, selection.method = "mean.var.plot", mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff)
    message(paste0("Number of HVG (via mvp method) was calculated: ", length(VariableFeatures(object))))
    plot.mvp(object = object, file_name = file_name, output_dir = output_dir)
    return(object)
}




## Do label transfer on the PCA dimensions
label_transfer <- function(object, reduction = "pca", dims.use = 40, transfer_cols = "label", k = 20){
    meta_use <- object@meta.data
    data.use <- object[[reduction]]@cell.embeddings[, 1:dims.use]

    
    ##A function to get the neighoring information
    get_most <- function(x){
        return(names(sort(table(x), decreasing = TRUE))[1])
    }


    ##Assume the NA values are the inquiry cell
    anno_list <- list()

    for (anno in transfer_cols){
        ref_cells <- rownames(meta_use)[!is.na(meta_use[, anno])]
        inquiry_cells <- setdiff(rownames(meta_use), ref_cells)


        ## Get the KNN cells for all the inquiry cells
        knn_cells <- FNN::get.knnx(data = data.use[ref_cells, ,drop = FALSE], query = data.use[inquiry_cells, ,drop = FALSE], k = k) %>%
                    .$nn.index


        ref_anno <- meta_use[ref_cells, anno] %>% setNames(., ref_cells)

        ## Do the annotation transfer
        if (class(ref_anno) %in% "character"){
            new_label <- apply(knn_cells, 1, function(x) get_most(ref_anno[x])) %>% 
                            setNames(., inquiry_cells) %>%
                            c(., ref_anno)
        } else if (class(ref_anno) %in% "numeric") {
            new_label <- apply(knn_cells, 1, function(x) median(ref_anno[x])) %>% 
                            setNames(., inquiry_cells) %>%
                            c(., ref_anno)
        } else {
            stop(paste0("column ", column, " has unsupported object type"))
        }

        anno_list[[anno]] <- new_label[rownames(meta_use)]
    }
    

    new_meta <- anno_list %>%  
                    as.data.frame(., stringsAsFactors = FALSE) %>% 
                    setNames(., paste0(names(anno_list), "_new"))
    object@meta.data <- cbind(object@meta.data[, setdiff(colnames(object@meta.data), colnames(new_meta))], new_meta)
    return(object)
}


