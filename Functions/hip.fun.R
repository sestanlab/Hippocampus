work_dir <- getwd()
inputdir <- paste0(work_dir, "/load_files/")
outputdir <- paste0(work_dir, "/report/")
dataDir <- "~/project/HIP/data/"


#Load all related packages 
library(ggplot2); library(ggrepel); library(cowplot); library(tibble);
library(dplyr); library(Matrix); 
library(parallel); library(viridis);library(batchelor)
#library(scran); 
if (sum(grepl("Seurat", (.packages()))) < 1){
    library(Seurat)#if source twice, the non-developmental version could be loaded and functions may change.
}


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


read_cellranger <- function(data.dir = NULL){ #The parameter gzip is removed in the new version

    ExtractField <- function(string, field = 1, delim = "_") {
        fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
        if (length(x = fields) == 1) {
            return(strsplit(x = string, split = delim)[[1]][field])
        }
        return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
    }

    full.data <- list()
    for (i in seq_along(data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(run)) {
            stop("Directory provided does not exist")
        }
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
            gzip <- ifelse(file.exists(paste0(run, "barcodes.tsv.gz")), TRUE, FALSE)
    
        if (gzip){
            gene.names <- readLines(gzfile(paste0(run, "features.tsv.gz")))
            cell.names <- readLines(gzfile(paste0(run, "barcodes.tsv.gz")))
            data <- readMM(gzfile(paste0(run, "matrix.mtx.gz")))
        } else {
            gene.names <- readLines(paste0(run, "genes.tsv"))
            cell.names <- readLines(paste0(run, "barcodes.tsv"))
            data <- readMM(paste0(run, "matrix.mtx"))
        }
        
        if (all(grepl(pattern = "\\-1$", x = cell.names))) {
            cell.names <- as.vector(x = as.character(x = sapply(X = cell.names,FUN = ExtractField, field = 1, delim = "-")))
        }
        rownames(x = data) <- sapply(gene.names, function(x) strsplit(x,split = "\t")[[1]][2]) %>%
                            setNames(., NULL)
       
        if (is.null(x = names(x = data.dir))) {
            stop("Supply names for the directory")
        } else {
            colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
        }
    
        full.data <- append(x = full.data, values = data)
    }
    full.data <- do.call(cbind, full.data)
    return(full.data)
}

#########################################################################################################

##Cluster the data 

#########################################################################################################
#This scirpt contains basic functions of single cell RNA-seq analysis (as well as many other technologies) 
map_gene <- function(gene_names, input_genes,ignore_case=TRUE){
    input_genes <- unique(input_genes)
  
    if (sum(grepl("\\|",gene_names))==length(gene_names)){
        if (sum(grepl("\\|",input_genes))==length(input_genes)){
              gene_id <- input_genes[input_genes %in% gene_names]
          }else{
            input_genes <- extract_field(input_genes=input_genes, 2, "|")
              gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
              gene_id <- gene_id[gene_id != "empty"]
          }
    } else if(sum(grepl("\\|",gene_names))==0){
        input_genes <- extract_field(input_genes=input_genes, 2, "|")
          gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
          gene_id <- gene_id[gene_id != "empty"]
    } else {
        stop("Inconsistent gene name format")
    }
    return(gene_id)
}

#Get the mito genes based on the inquiry genes
get_genes <- function(input_genes, gene_type = c("mito","ribo", "cc")[1], return_list = FALSE, revised = FALSE){
    gene_use <- list()
    if ("mito" %in% gene_type){
        mito.known <- map_gene(gene_names=input_genes, input_genes=c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")) #Refseq annotation 103 (Macaque)
        mito.combine <- grep(pattern = "\\|MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        mito.single <- grep(pattern = "^MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        gene_use[["mito"]] <- unique(c(mito.known, mito.combine, mito.single))
    }

    if ("ribo" %in% gene_type){
        ribo.combine <- c(grep("\\|RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        ribo.single <- c(grep("^RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        gene_use[["ribo"]] <- c(ribo.combine, ribo.single)
    }

    if ("cc" %in% gene_type){
        if (!revised){
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        } else {
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        }
        
        gene_use[["s"]] <- map_gene(gene_names=input_genes, input_genes=regev.s.genes,ignore_case=TRUE)
        gene_use[["g2m"]] <- map_gene(gene_names=input_genes, input_genes=regev.g2m.genes,ignore_case=TRUE)
    }

    if (return_list){
        return(gene_use)
    } else {
        all_genes <- setNames(unlist(gene_use), NULL)
        return(all_genes)
    }
}


seu_prepare <- function(counts, data = NULL, min.cells = 5, normalization.method = c("LogNormalize", "none","scran")[1], nfeatures = 2500, hvg.method = "vst", assay = "RNA") {
    #Check the normalization.method
    if (!is.null(data)){
        message("input-norm is provided and therefore normalization is not required")
        normalization.method <- "none"
    }

    if (normalization.method == "none"){
        message("normalization.method is none and therefore assume the counts has already been normlized")
    } else if (normalization.method == "LogNormalize"){
        message("normalization.method is LogNormalize and therefore will use seurat LogNormalize to transform the dataset")
    } else {
        stop("Unknown normalization.method")
    }

    inseu <- CreateSeuratObject(counts, meta.data = NULL, assay = "RNA", min.cells = min.cells, min.features = 0, names.field = 1, names.delim = "_"); #at least expressed in 5 cells #only RNA assay is supported

    #quality_gene list
    quality_genes <- get_genes(input_genes = rownames(inseu$RNA@data), gene_type = c("mito","ribo", "cc"), return_list = TRUE, revised = FALSE)
    inseu[["percent.mt"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["mito"]], col.name = NULL)
    inseu[["percent.ribo"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["ribo"]], col.name = NULL)


    #normalize the dataset if needed
    if (normalization.method == "none"){
        inseu[[assay]] <- CreateAssayObject(data = as(ifelse_check(is.null(data), inseu[["RNA"]]@data, data), "dgCMatrix"), min.cells = 0, min.features = 0)
        DefaultAssay(inseu) <- assay
    } else if (normalization.method == "LogNormalize"){
        inseu <- NormalizeData(inseu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    }
    
    if (!is.null(hvg.method)){
        inseu <- FindVariableFeatures(inseu, selection.method = hvg.method, nfeatures = nfeatures, verbose = FALSE)
    }
    return(inseu)
}

extract_field <- function(input_genes, field = 1, split_by = "_") {
    split_strings <- strsplit(input_genes, split_by, fixed=TRUE)
    if (is.numeric(field)){
        if (all(field > 0)){
            if (length(field) == 1){
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) split_strings[[x]][idx[x]])
            } else {
                idx <- sapply(split_strings, length) == 1
                output_strings <- sapply(1:length(split_strings), function(x) ifelse(idx[x], input_genes[x], paste(split_strings[[x]][field], collapse = split_by)))
            }
        } else {
            if (length(field) == 1) {
                field = as.integer(abs(field))
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) rev(split_strings[[x]])[idx[x]])
            } else {
                stop("currently doesnot support field with length >= 2")
            }
        }
    } else if (is.character(field)) {
        if (field == "rm_start"){
            idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
            output_strings <- sapply(1:length(split_strings), function(x) paste(split_strings[[x]][-1],collapse=split_by))
        } else if (field == "rm_end") {
            output_strings <- sapply(1:length(split_strings), function(x) paste(rev(rev(split_strings[[x]])[-1]),collapse=split_by))
        } else {
            stop("Currently only support rm_start, rm_end for the character field")
        }
    }
    output_strings <- output_strings %>% setNames(., input_genes)
    return(output_strings)
}


#Function in Seurat package, useful
MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}


remove_duplicates <- function(x, return_index = FALSE){
    dp_index <- !duplicated(x) & !rev(duplicated(rev(x)))
    if(return_index){
        return(dp_index)
    } else {
        return(x[dp_index])
    }
}


ifelse_check <- function(test, yes, no){if (test){return(yes)} else{ return(no) }}


find_hvg <- function(input_data = NULL, input_list = NULL, input_meta = NULL, group.by, nfeatures = 1500, ninte_features = 2000, min.cells = 5) {
    
    if (!is.null(input_list)){
        if (class(input_list[[1]]) == "Seurat"){
            seu_list <- lapply(input_list, function(seu) {
                seu <- seu  %>%
                    NormalizeData(., normalization.method = "LogNormalize", verbose = FALSE) %>%
                    FindVariableFeatures(., selection.method = "vst", nfeatures = nfeatures, verbose = FALSE) 
                seu
            })
            rm(input_list)
        } else if (grepl("matrix", class(input_list[[1]]), ignore.case = TRUE)){
            seu_list <- lapply(input_list, function(raw) {
                seu <- seu_prepare(input_raw = raw, input_norm = NULL, min.cells = min.cells, normalization.method = "LogNormalize", nfeatures = nfeatures, hvg_method = "vst", assay_name = "RNA")
                seu
                })
        }
    }


    if (class(input_data) %in% c("Seurat")){
        seu_use <- input_data; rm(input_data)
        all_samples <- levels(as.factor(seu_use@meta.data[, group.by]))
        seu_list <- lapply(all_samples, function(sample_name) {
            seu <- seu_use[, as.character(seu_use@meta.data[, group.by]) == sample_name]  %>%
                    NormalizeData(., normalization.method = "LogNormalize", verbose = FALSE) %>%
                    FindVariableFeatures(., selection.method = "vst", nfeatures = nfeatures, verbose = FALSE) 
            seu
            }) %>% setNames(., all_samples)
    } else if (grepl("matrix", class(input_data), ignore.case = TRUE)){
        if (is.null(input_meta)){
            stop("The input data is a matrix, please provide the input_meta")
        }

        all_samples <- levels(as.factor(input_meta[, group.by]))
        seu_list <- lapply(all_samples, function(sample_name){
            sample_cells <- rownames(input_meta)[as.character(input_meta[, group.by]) == sample_name]
            inseu <- seu_prepare(input_raw = input_data[, sample_cells], input_norm = NULL, min.cells = min.cells, normalization.method = "LogNormalize", nfeatures = nfeatures, hvg_method = "vst", assay_name = "RNA")
            return(inseu)
        }) %>% setNames(., all_samples)
    }


    ##Get the integration features
    feature_use <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = ninte_features)
    return(feature_use)
}       


mean.of.logs <- function (x, base = 2) {
    return(log(mean((base^x) - 1) + 1, base = base))
}



library(AUCell)
GetModuleScore <- function (assay.data, features, nbin = 24, ctrl = 100, k = FALSE, seed = 42, method = c("seurat","aucell")[2], input_dir = new_inputdir, file_name, output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000) {
     if (is.null(x = features)) {
        stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
        return(intersect(x = x, y = rownames(x = assay.data)))
    })
    cluster.length <- length(x = features) #number of feature list

    if (method == "seurat"){
        set.seed(seed = seed)
        pool <- rownames(x = assay.data)
        data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
        data.avg <- data.avg[order(data.avg)]
        data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
            n = nbin, labels = FALSE, right = FALSE)
        names(x = data.cut) <- names(x = data.avg)
        ctrl.use <- vector(mode = "list", length = cluster.length)
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            for (j in 1:length(x = features.use)) {
                ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                    data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
            }
        }
        ctrl.use <- lapply(X = ctrl.use, FUN = unique)
        ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
            ncol = ncol(x = assay.data))
        for (i in 1:length(ctrl.use)) {
            features.use <- ctrl.use[[i]]
            ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
               ])
        }
        features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
            ncol = ncol(x = assay.data))
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            data.use <- assay.data[features.use, , drop = FALSE]
            features.scores[i, ] <- Matrix::colMeans(x = data.use)
        }
        features.scores.use <- features.scores - ctrl.scores
        features.scores.use <- as.data.frame(x = t(x = features.scores.use))
        rownames(x = features.scores.use) <- colnames(x = assay.data)
        colnames(features.scores.use) <- names(features)

        return(features.scores.use)
    } else if (method == "aucell"){
        library(AUCell)
        if (!file.exists(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))){
            #Split the cells to bins[Sometimes, the matrix is too large for rankings]
            if (ncol(assay.data) < cellbin.size){
                cellbin.size <- ceiling(ncol(assay.data)/2)
            }
            bin.ind <- ceiling(c(1:ncol(assay.data))/cellbin.size)
            max.bin <- max(bin.ind)

            auc_list <- lapply(1:max.bin, function(tembin) {
                tem_matrix <- assay.data[, bin.ind == tembin]
                tem_rankings <- AUCell_buildRankings(tem_matrix, nCores=1, plotStats=FALSE) 
                tem_auc <- AUCell_calcAUC(features, tem_rankings)#, aucMaxRank = 500)
                tem_aucmatrix <- t(as.matrix(getAUC(tem_auc)))
                rm(tem_matrix, tem_rankings)
                return(tem_aucmatrix)
                })

            hauc_matrix <- do.call(rbind, auc_list)

            if (length(auc_list) == 1){ #When the input is not a named list, the colnames will be "Geneset"
                colnames(hauc_matrix) <- names(features)
            }
            saveRDS(hauc_matrix, file = paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        } else {
            hauc_matrix <- readRDS(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        }
        
        set.seed(seed)
        pdf(paste0(output_dir, file_name, "_modulescore_auc_auto_assignment.pdf"), paper="letter")
        par(mfrow=c(2,2)); cells_assignment <- AUCell_exploreThresholds(t(hauc_matrix), plotHist=TRUE, assign=TRUE) 
        dev.off()

        #Generate assignment matrix (rownames as cells)
        default_assign <- hauc_matrix * 0 #build an empty one
        for (gset in colnames(hauc_matrix)){
            default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
        }

        outdata <- list(auc = hauc_matrix, auto = default_assign, custom = default_assign)

        if (!is.null(rethreshold_list)){
            pdf(paste0(output_dir, file_name, "_modulescore_auc_custom_assignment.pdf"), paper="letter")
            for (j in names(rethreshold_list)){
                AUCell_plotHist(t(hauc_matrix)[j,,drop = FALSE], aucThr=rethreshold_list[[j]])
                abline(v=rethreshold_list[[j]])
                cells_assignment[[j]]$assignment <- rownames(hauc_matrix)[hauc_matrix[, j] >= rethreshold_list[[j]]]
            }
            dev.off()

            ##default_assign <- hauc_matrix * 0 #build an empty one
            for (gset in colnames(hauc_matrix)){
                default_assign[, gset] <- 0
                default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
            }
            outdata$custom <- default_assign
        }
        return(outdata)
    }
}



auto_levels <- function(object, split.by) { 
    all_levels <- levels(as.factor(object@meta.data[, split.by]))
    level_list <- list(c("FR", "MS", "TP", "OX"), 
                        c("Frontal", "MotorSensory", "Temporal", "Occipital"), 
                        c("Human", "Chimpanzee", "Macaque", "Marmoset"),
                        c("Human", "Chimpanzee", "Rhesus", "Marmoset"),
                        c("HSB", "PTB", "RMB", "MMB"),
                        c("E37", "E42-43", "E54", "E62-E64", "E77-78"))
    ## get the idx of reference names
    nshare <- sapply(level_list, function(x) sum(all_levels %in% x))
    idx <- which(idx == max(idx))[1]

    if (max(nshare) == 9){
        return(all_levels)
    } else {
        new_levels <- level_list[[idx]][level_list[[idx]] %in% all_levels]
        return(new_levels)
    }
}




FeatureFig <- function(object, features, cols = c("#f5f5dc", "#31a354","#253494"), split.by = NULL, plot.scale = 0.9, ncol = NULL, file_name, output_dir = outputdir, pt.size = "auto", split.order = "auto", output.ggplot = FALSE, ngradient = 10, ...) {
    
    extra_arguments <- list(...)


    if (pt.size == "auto"){
        point.size <- round(10/sqrt(ncol(object)), digits = 1) %>% MinMax(., min = 0.1, max = 1)
    } else {
        point.size <- pt.size
    }

    ##Set the levels of factors in the supervised mode
    if (!is.null(split.by)){
        ##Get the supervised order the split.by
        new_levels <- ifelse_check("auto" %in% split.order, 
                            auto_levels(object = object, split.by = split.by),
                            split.order)
        object@meta.data[, split.by] <- factor(as.character(object@meta.data[, split.by]), levels = new_levels) 
    }


    cols <- colorRampPalette(cols)(ngradient)

    ##Get the plot list
    if (length(features) == 1){
        plot_args <- c(list(object = object, cols = cols, features = features, pt.size = point.size, combine = FALSE, split.by = split.by), extra_arguments[names(extra_arguments) %in% c("order", "reduction", "min.cutoff", "max.cutoff", "shape.by", "slot", "blend","blend.threshold", "label", "label.size", "repel", "combine", "coord.fixed", "by.col", "sort.cell")])
        plist <- do.call(FeaturePlot, plot_args)
        if (!is.null(split.by)){
            names(plist) <- paste0(features, "-", rep(new_levels, length.out = length(plist)))
        } else {
            names(plist) <- features
        }
    } else {
        plot_args <- c(list(object = object, cols = cols, features = features, pt.size = point.size, combine = FALSE, split.by = split.by), extra_arguments[names(extra_arguments) %in% c("order", "reduction", "min.cutoff", "max.cutoff", "shape.by", "slot", "blend","blend.threshold", "label", "label.size", "repel", "combine", "coord.fixed", "by.col", "sort.cell")])
        plist <- do.call(FeaturePlot, plot_args)


        ##Set the names of the plot list 
        if (!is.null(split.by)){
            new_plist <- lapply(1:(length(plist)/length(new_levels)), function(x) {
                gidx <- seq(x, length(plist), length(plist)/length(new_levels))
                plist[gidx]
                }) %>% do.call(c, .)
            plist <- new_plist; rm(new_plist)

            gene_label <- sapply(plist, function(p) ggplot_build(p)$plot$plot_env$plot$layers[[1]]$mapping$colour[2] %>% as.character())
            gene_label <- paste0(gene_label, "-", rep(new_levels, length.out = length(plist)))
        } else {
            gene_label <- sapply(plist, function(p) ggplot_build(p)$plot$plot_env$plot$layers[[1]]$mapping$colour[2] %>% as.character())
        }

        names(plist) <- gene_label
    }


    ##Polish the figures by changing the theme
    plist <- lapply(names(plist), function(x) {
        p <- plist[[x]] + 
            coord_equal(ratio = 1) + 
                theme_classic() + 
                scale_color_gradientn(colors = cols, limits = c(1,ngradient)) + 
                theme(legend.position = "bottom",
                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))+ 
                labs(title =  text_wrapper(x, width=20))
        p
        })


    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist)
    } else {
        nlevel <- ifelse(is.null(split.by), 1, length(new_levels))
        ncol <- ifelse_check(is.null(ncol), nlevel, ncol)
        nrow <- ceiling(length(plist)/ncol)
        
        jpeg(paste0(output_dir, file_name, "_feature.jpeg"), width = 10 * plot.scale * ncol, height = 11 * plot.scale * nrow, units = "in", res = 300)
        plot_grid(plotlist = plist, nrow = nrow, ncol = ncol) %>% print() %>% print()
        dev.off()
    }
}



DimFig <- function(object, file_name, group.by = NULL, split.by = NULL, cols = NULL, plot.scale = 0.9, pt.size = "auto", label.size = 4, split.order = "auto", legend.position = "right", output_dir = outputdir, output.ggplot = FALSE, ...){
    p1 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = legend.position, label = TRUE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)
    p2 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = "none", label = TRUE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)
    p3 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = "none", label = FALSE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)


    plist <- lapply(1:length(group.by), function(x) list(p1[[x]], p2[[x]], p3[[x]]))

    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist %>% do.call(c, .))
    } else {

        nlevel <- ifelse_check(is.null(split.by), 1, ifelse_check("auto" %in% split.order, auto_levels(object = object, split.by = split.by),split.order) %>% length())
        nrow <- 3
        

        max_cls <- sapply(group.by, function(x) unique(object@meta.data[, x]) %>% length()) %>% max()
        pdf_width <- ifelse_check(legend.position == "none", 10 * plot.scale * nlevel, 10 * plot.scale * nlevel + ceiling(0.15 * length(max_cls)))

        for (xx in 1:length(group.by)){
            jpeg(paste0(output_dir, file_name, paste0("_", group.by[[xx]]), ifelse(is.null(split.by), "", paste0("_", split.by)), ".jpeg"), width = pdf_width, height = 10 * plot.scale * 3, units = "in", res = 300)
            plot_grid(plotlist = plist[[xx]], nrow = 3, ncol = 1) %>% print() %>% print()
            dev.off()
        }
    }
}



DimFig.default <- function(object, file_name, group.by = NULL, split.by = NULL, cols = NULL, plot.scale = 0.9, legend.position = "right", label = TRUE, pt.size = "auto", label.size = 4, split.order = "auto", output_dir = outputdir, output.ggplot = FALSE, ...){

    extra_arguments <- list(...)


    ##Set the point size
    if (pt.size == "auto"){
        point.size <- round(10/sqrt(ncol(object)), digits = 1) %>% MinMax(., min = 0.1, max = 1)
    } else {
        point.size <- pt.size
    }


    ##Set the levels of factors in the supervised mode
    if (!is.null(split.by)){
        ##Get the supervised order the split.by
        new_levels <- ifelse_check("auto" %in% split.order, 
                            auto_levels(object = object, split.by = split.by),
                            split.order)
        object@meta.data[, split.by] <- factor(as.character(object@meta.data[, split.by]), levels = new_levels) 
    }


    plot_args <- c(list(object = object, group.by = group.by, cols = cols, pt.size = point.size, combine = FALSE, split.by = split.by, label = label, label.size = label.size), extra_arguments[names(extra_arguments) %in% c("shape.by", "reduction", "cells", "repel", "cells.highlight" , "cols.highlight", "sizes.highlight", "na.value")])
    plist <- do.call(DimPlot, plot_args)


    ##Polish the figures by changing the theme
    plist <- lapply(plist, function(p) {
        p <- p + 
            #coord_equal(ratio = 1) + 
            #    theme_classic() + 
                theme(legend.position = legend.position)#,
            #        line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
            #        axis.text.x=element_blank(),axis.text.y=element_blank(), 
            #        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))+ 
            #    labs(title =  text_wrapper(x, width=20))
        p
        })


    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist)
    } else {
        nlevel <- ifelse(is.null(split.by), 1, length(new_levels))
        #ncol <- ifelse_check(is.null(ncol), nlevel, ncol)
        nrow <- length(plist)
        

        max_cls <- sapply(group.by, function(x) unique(object@meta.data[, x]) %>% length()) %>% max()
        pdf_width <- ifelse_check(legend.position == "none", 10 * plot.scale * nlevel, 10 * plot.scale * nlevel + ceiling(0.15 * length(max_cls)))

        jpeg(paste0(output_dir, file_name, ifelse(is.null(group.by), "", paste0("_", paste(group.by, collapse = "-"))), ifelse(is.null(split.by), "", paste0("_", split.by)), ifelse(label, "_labeled", ""), ".jpeg"), width = pdf_width, height = 10 * plot.scale * nrow, units = "in", res = 300)
        plot_grid(plotlist = plist, nrow = nrow, ncol = 1) %>% print() %>% print()
        dev.off()
    }
}


DotFig <- function(object, assay = "RNA", features, cols = c("white", "red"), dot.scale = 3, dot.min = 0, group.by = "hres", file_name, output_dir = outputdir) {
    p <- DotPlot(object = object, assay = assay, features = features, cols = cols, dot.scale = dot.scale,dot.min = dot.min, group.by = group.by) + 
        RotatedAxis() + 
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))


    ncls <- unique(object@meta.data[, group.by]) %>% length()
    plot_heights <- ceiling(0.2 * ncls + 2) %>% MinMax(., min = 5, max = 20)
    plot_width <- ceiling(0.6 * length(features) + 2) %>% MinMax(., min = 5, max = 20)
    jpeg(paste0(output_dir, file_name, "_markers.jpeg"), width = 4.3, height = 6, units = "in", res = 300)
    print(p)
    dev.off()
}



get_ident <- function(input_ident, ident_col = NULL, file_path, sample_name, label_names = NULL){
    ## Extract the raw identity
    if (grepl("Seurat", class(input_ident), ignore.case = TRUE)){
        tem_ident <- setNames(as.character(input_ident@meta.data[, ident_col]), colnames(input_ident))
    } else {
        tem_ident <- as.character(input_ident)
    }


    ## Set the ident names
    sub_file <- read.table(file_path, sep="\t", stringsAsFactors=FALSE, header=TRUE) %>%
                    mutate(cluster = as.character(cluster)) %>%
                    subset(sample == sample_name)
    ident_names <- setNames(sub_file[, ifelse(is.null(label_names), "label", label_names)], sub_file$cluster)


    ## check if all the identity has been annotated
    extras <- setdiff(unique(tem_ident), names(ident_names))
    if(length(extras) > 0){
        stop(paste0("The following idents are not annotated: ", paste(extras, collapse = ", ")))
    }


    output_ident <- setNames(ident_names[tem_ident], names(tem_ident))
    return(output_ident)
}


LinearRegress <- function (input_data, latent_data, vars_regress, genes_regress = NULL,display_progress = TRUE, num.cores = NULL){
    input_data <- as.matrix(input_data) #not implemented for the "Matrix" object, need to update
    genes_regress <- ifelse_check(is.null(genes_regress), rownames(input_data), genes_regress)
    genes_regress <- intersect(x = genes_regress, y = rownames(input_data))
    latent_data <- latent_data[colnames(input_data),vars_regress, drop=FALSE]
    bin.size <- 100
    bin.ind <- ceiling(c(1:length(x = genes_regress))/bin.size)
    max.bin <- max(bin.ind)
    if (display_progress) {
        message(paste("Regressing out:", paste(vars_regress,collapse = ", ")))
        pb <- txtProgressBar(min = 0, max = max.bin, style = 3,file = stderr())
    }
    data.resid <- c()
    
    print("prepare for parallel calculation")
    num.cores <- ifelse_check(is.null(num.cores), detectCores()/2, num.cores)
    cl <- parallel::makeCluster(num.cores)
    doSNOW::registerDoSNOW(cl)
    opts <- list()
    if (display_progress) {
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        time_elapsed <- Sys.time()
    }
    
    #print("Set regression variables")
    reg.mat.colnames <- c(colnames(x = latent_data), "GENE")
    fmla_str = paste0("GENE ", " ~ ", paste(vars_regress,collapse = "+"))  #linear regression formula
    regression.mat <- cbind(latent_data, input_data[1, ])
    colnames(regression.mat) <- reg.mat.colnames
    qr = lm(as.formula(fmla_str), data = regression.mat,qr = TRUE)$qr
    rm(regression.mat)

    print("Do parallel regression")
    data.resid <- foreach::foreach(i = 1:max.bin, .combine="c",.options.snow = opts) %dopar% {
            genes.bin.regress <- rownames(x = input_data)[bin.ind ==i]
            gene.expr <- as.matrix(x = input_data[genes.bin.regress, , drop = FALSE])
            empty_char <- "" # character(length = dim(gene.expr)[1])
            new.data <- sapply(X = genes.bin.regress, FUN = function(x) {
                resid <- qr.resid(qr, gene.expr[x, ])
                if (!is.list(resid)) {
                  resid <- list(resid = resid, mode = empty_char)
                }
                return(resid)
            })          
            #print(length(x = new.data))
            new.data.resid <- new.data[seq.int(from = 1, to = length(x = new.data),by = 2)]
            #print(head(new.data.resid, 5))
            new.data.resid <- matrix(unlist(new.data.resid), nrow = length(new.data.resid[[1]]))
            colnames(x = new.data.resid) <- genes.bin.regress
            new.data.mode <- unlist(x = new.data[seq.int(from = 2,to = length(x = new.data), by = 2)])
            names(x = new.data.mode) <- genes.bin.regress
            new.data <- list(resid = new.data.resid, mode = new.data.mode)
            return(new.data)
    }
    print(length(data.resid))
    if (display_progress) {
        time_elapsed <- Sys.time() - time_elapsed
        cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)))
        close(pb)
    }
    stopCluster(cl)
    modes <- unlist(x = data.resid[seq.int(from = 2, to = length(x = data.resid),by = 2)])
    modes <- modes[modes == "scale"]
    names(x = modes) <- gsub(pattern = "mode.", replacement = "", x = names(x = modes), fixed = TRUE)
    out_matrix <- data.resid[seq.int(from = 1, to = length(x = data.resid),by = 2)]
    out_matrix2 <- as.matrix(x = as.data.frame(x = out_matrix))
    out_matrix2 <- t(x = out_matrix2)
    print(dim(out_matrix2))
    if (length(x = modes)) {
        message("The following genes failed with glm.nb, and fell back to scale(log(y+1))\n\t",
            paste(names(x = modes), collapse = ", "))
    }
    rownames(out_matrix2) <- genes_regress
    colnames(out_matrix2) <- colnames(input_data)
    suppressWarnings(expr = gc(verbose = FALSE))
    return(out_matrix2)
}



#This scirpt wrap the long text to several lines, designed for the title of ggplot
text_wrapper <- function(x, width=15) {
  y <- sapply(x, function(x) paste(strsplit(x, "_", fixed=TRUE)[[1]], collapse=" "))
  return(paste(strwrap(y, width=width), collapse = "\n"))
}





#new_vln(object_use= hfint_snorm, input_genes=names(dex_max)[dex_max <= 0.04 & Matrix::rowSums(sel_matrix) < 5][1:20],input_identity= hfint_ident, cluster_color =HSB_colors, ident_level =HSB_order, logscale=FALSE, output_dir = outputdir, file_name="type2_dex")
new_vln <- function(object_use, input_genes,input_identity , cluster_color =NULL, ident_level, logscale=FALSE, output_dir = outputdir, file_name, pdf_width =8.5, pdf_height =11, nrow =1, ncol =2, plot_box = FALSE) {
    #Get the required data and Scale the data if necessary
    input_genes <- rev(input_genes)
    exp_data <- Extract_exp_data(input_object = object_use, vars_all = input_genes)
    colnames(exp_data) <- split_genename(input_genes = colnames(exp_data))
    exp_data <- exp_data[names(input_identity),]
    input_gene_number <- dim(exp_data)[2]
    genes <- colnames(exp_data)
    max_vals <- apply(exp_data, 2, max)
  
    for (i in 1:length(genes)) {
        gene <- genes[i]
        gene_max <- max_vals[i]
        if (logscale) {
            exp_data[gene] <- log10(exp_data[gene] + 1)/log10(gene_max + 1) * 0.99 + i
        } else {
            exp_data[gene] <- exp_data[gene]/gene_max * 0.99 + i
        }
    }
    exp_data <- as.data.frame(exp_data,check.names=TRUE)
    
    #Set the identity for the expression matrix, rows as cell, columns as genes, the last column is the identity
    if (!is.null(input_identity)){
        plot_identity <- ifelse_check(is.null(ident_level), as.factor(input_identity),factor(as.factor(input_identity), levels=ident_level))
    } else if(is.null(input_identity) & !is.null(object_use@ident)) {
        plot_identity <- ifelse_check(is.null(ident_level), object_use@ident, factor(object_use@ident, levels=ident_level))
    } else{
        cat("\nPlease provide the identity for the plot\n")
    }

    plot_data_infunction <- data.frame(exp_data, identity=plot_identity,check.names=TRUE, stringsAsFactors = FALSE)
    genes <- colnames(plot_data_infunction)[-dim(plot_data_infunction)[2]]
  
    inner_violn <- function(sel_genes = genes, x_order = ident_level, plot_data = plot_data_infunction, x_color = cluster_color, height_p = oo, height_q =pp){
        if (length(sel_genes) < 20){
            label_genes <- c(sel_genes, rep("", 20-length(sel_genes)))
            height_q <- 20
        } else {
            label_genes <- sel_genes
        }
        print(label_genes)

        p <- ggplot() + scale_y_continuous("", breaks = height_p:height_q + 0.45, labels = label_genes, expand = c(0,0))
        if(!is.null(x_order)){
            p <- p+ scale_x_discrete(x_order, expand = c(0, 0))
        }

        for (i in 1:length(sel_genes)) {
            if (!plot_box){
                p <- p + geom_violin(data = plot_data, aes_string(x = "identity", y = sel_genes[i], fill = "identity"), scale = "width", size = 0.1, adjust = 2,trim =TRUE) + stat_summary(data = plot_data, aes_string(x = "identity", y = sel_genes[i]), fun.y = "median", fun.ymin = "median", fun.ymax = "median", geom = "point", size = 0.1)
            } else {
                p <- p + geom_boxplot(data = plot_data, aes_string(x="identity", y=sel_genes[i], fill = "identity"), outlier.colour="white", outlier.shape=16, outlier.size=0.1)
            }
        }

        p <- p + scale_y_continuous("", breaks = height_p:height_q + 0.45, labels = label_genes, expand = c(0,0))
        if (!is.null(x_color)){
            cluster_levels <- ifelse_check(is.null(x_order), levels(as.factor(plot_data_infunction$identity)), x_order)
            plot_colors <- x_color[cluster_levels]   #extract the color based on the cluster level names
            p <- p + scale_fill_manual(values=plot_colors)
        }
        p <- p  +theme(legend.position = "none",axis.title.x = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),axis.title.y = element_text(angle = 0, size = 10,vjust=0.5),axis.text.x=element_text(angle = 270, size = 9, hjust=0, vjust = 0.5),axis.text.y=element_text(angle = 0, size = 9), axis.line = element_line(colour = "black"))
        return(p)
    }

    violin_list <- lapply(seq(1,length(seq(1,input_gene_number,20))), function(x) {
        oo <- seq(1,input_gene_number,20)[x]
        pp <- c(seq(min(20,input_gene_number),input_gene_number,20),input_gene_number)[x]
    
        tem_genes <- genes[oo:pp]
        tem_plot_data <- plot_data_infunction[, c(oo:pp, dim(plot_data_infunction)[2])]

        pp <- c(seq(min(20,input_gene_number),input_gene_number,20),40)[x]
        return(inner_violn(sel_genes = tem_genes, x_order = ident_level, plot_data = tem_plot_data, x_color = cluster_color, height_p = oo, height_q =pp))
    })
  
  #violin_list[23:40] <- lapply(1:18, function(x) NULL)
  file_type <- ifelse_check(plot_box, "_boxPlot.pdf", "_vlnPlot.pdf")
  pdf(paste0(output_dir,file_name, file_type),width = pdf_width, height = pdf_height)
  organize_figures(plot_list = violin_list, nrow = nrow, ncol = ncol, title_name = file_name)
  dev.off()
}


make_hexbin_byregion <- function(sce, nbins = 80, dimension_reduction = "UMAP", split.by = "species"){
    all_regions <- levels(as.factor(sce@meta.data[, split.by]))
    sce_list <- list()
    for (region in all_regions){
        sub_sce <- sce[, rownames(sce@meta.data)[sce@meta.data[, split.by] == region]]
        sce_list[[region]] <- make_hexbin(sub_sce, nbins = nbins, dimension_reduction = dimension_reduction)
    }
    return(sce_list)
} 



plot_hexbin_feature_byregion <- function (sce, assay, slot, features, action, nrow = 1, ncol = NULL, colors = viridis(3), file_name = NULL, output_dir, pdf_size = c(8,5), region_order = NULL, title_color = NULL, legend = "bottom", return_rawp = FALSE, exp_ceiling = NULL) {
    if (class(sce) != "list"){
        stop("The input sce should be a list output from make_hexbin_byregion")
    }
    if (!assay %in% names(sce[[1]])) {
        stop("Specify a valid assay.")
    }
    if (!slot %in% slotNames(GetAssay(sce[[1]], assay))) {
        stop("Specify a valid assay slot")
    }


    ##First prepare the data list
    all_regions <- ifelse_check(is.null(region_order), names(sce), region_order)
    features <- intersect(features, rownames(sce[[1]][[assay]]))
    if (length(features) == 0) {
            stop("No features can be matched to the dataset")
    }
    data_flist <- lapply(all_regions, function(region) GetAssayData(sce[[region]], assay = assay, slot = slot)[features, ,drop = FALSE]) %>% setNames(., all_regions)


    feature_plist <- list()
    allp_list <- list()
    for (feature in features){
        regp_list <- list()
        max_value <- 0
        exp_list <- lapply(data_flist, function(curdata) as.numeric(curdata[feature,])) %>% setNames(., names(data_flist))
        out_list <- list()
        for (region in all_regions){
            x <- exp_list[[region]]
            out <- sce[[region]]@misc$hexbin[[2]]
            cID <- sce[[region]]@misc$hexbin[[1]]
            hh <- schex:::.make_hexbin_function(x, action, cID)
            out <- as_tibble(out)
            if (grepl("^[[:digit:]]", feature)) { ##incase the feature start with a number
                feature <- paste0("F_", feature)
            }
            feature <- gsub("-", ".", feature) ##incase the "-" will be converted to "." in tibble
            col_hh <- paste0(feature, "_", action)
            eval(parse(text = paste0("out$", col_hh, " <- hh")))
            max_value <- max(max_value, max(hh))
            if (!is.null(exp_ceiling)){
                max_value <- exp_ceiling
                out[, col_hh] <- MinMax(out[, col_hh], min = 0, max = exp_ceiling)
            }
            out_list[[region]] <- out
        }
        max_value <- ifelse(max_value == 0, 1, max_value)
        

        regp_list <- list()
        for (region in all_regions){
            regp_list[[region]] <- schex:::.plot_hexbin(out_list[[region]], colour_by = paste0(feature, "_", action), title = paste0(feature, "  ",region), xlab = xlab, ylab = ylab) + 
                scale_fill_gradientn(limits = c(0,max_value),colours = colors) + 
                theme(legend.position = legend,line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
            if (!is.null(title_color)){
                regp_list[[region]] <- regp_list[[region]] + theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", color = title_color[region]))
            }
        }

        allp_list[[feature]] <- regp_list
        feature_plist[[feature]] <- plot_grid(plotlist = regp_list, nrow = nrow, ncol = ifelse(is.null(ncol), length(all_regions), ncol))
    }
    
    if (return_rawp){
        return(allp_list)
    }


    if (is.null(file_name)){
        return(feature_plist)
    } else {
        jpeg(paste0(output_dir, file_name, "_hexexpression.jpeg"),pdf_size[1], pdf_size[2], units = "in", res = 300)
        for (feature in names(feature_plist)){
            print(feature_plist[[feature]])
        }
        dev.off()
    }
}


plot_hex_cbn <- function(sce, feature, action = "mean", cols = viridis(3)){
    if (!feature %in% rownames(sce$RNA)) {
        stop("Specify a valid assay.")
    }

    ##First prepare the data list
    xx <- as.numeric(GetAssayData(sce, slot = "data")[feature, ])
    out <- sce@misc$hexbin[[2]]
    cID <- sce@misc$hexbin[[1]]
    plot.data <- data.frame(out, 
                exp = schex:::.make_hexbin_function(xx, action, cID), 
                stringsAsFactors = FALSE, check.names = FALSE)

    scale_exp <- function(x) {
        oo <- (x - mean(x))/sd(x)
        return(oo)
    }
    nbks <- 30
    plot.data <- plot.data %>%
                    mutate(scale.exp = scale_exp(exp)) %>%
                    mutate(scale.exp = MinMax(exp, -2.5, 2.5)) %>%
                    mutate(scale.exp = as.numeric(cut(scale.exp, nbks)))
    plot.data$scale.exp[is.na(plot.data$scale.exp)] <- 1


    p <- ggplot(data = plot.data) +
                    geom_hex(mapping = aes_string(x = "x", y = "y", fill = "scale.exp"), stat = "identity") + 
                    theme_classic() + 
                    labs(title = feature) + 
                    coord_fixed() + 
                    scale_fill_gradientn(colors = colorRampPalette(colors = cols)(nbks)) + 
                    theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    return(p)
}

plot_density.smooth <- function(cell_coordinates=NULL, dot_size=0.1, file_name="aa", output_dir=outputdir, output_ggplot = FALSE, pdf_size = c(9,9), region_infor = NULL, region_order = NULL,density_color = rev(rainbow(10, end = 4/6)), backgroud_color = "lightgrey", density_sf = 0.25, add_contour = FALSE, nbins = 100, legend = "bottom", add_point = TRUE, nrow = 1, ncol = NULL) {
    dr_name <- colnames(cell_coordinates) %>% gsub("_1", "", .)
    colnames(cell_coordinates) <- c("tSNE_1", "tSNE_2")
  
    region_infor <- region_infor[rownames(cell_coordinates)]
    plot_data_infunction <- data.frame(cell_coordinates, region = region_infor, check.names=FALSE, stringsAsFactors = FALSE)
    plot_data_infunction$dcolor <- densCols(plot_data_infunction[, 1], plot_data_infunction[, 2], colramp = colorRampPalette(density_color)) #rev(rainbow(10, end = 4/6))
    
    ##Plot all the cells as the backgroud cells
    p <- ggplot(data=plot_data_infunction, aes_string(x="tSNE_1", y="tSNE_2", color = "dcolor"))
   

    ##Plot the density across regions
    all_regions <- levels(as.factor(plot_data_infunction$region))
    region_order <- ifelse_check(is.null(region_order), all_regions, region_order)
    plist <- lapply(region_order, function(reg) {
        q <- p + geom_point(data = subset(plot_data_infunction, region != reg), aes_string(x="tSNE_1", y="tSNE_2"), color = backgroud_color, size = dot_size) +
            geom_point(data = subset(plot_data_infunction, region == reg), aes_string(x="tSNE_1", y="tSNE_2", color = "dcolor"), size = dot_size) + 
            scale_color_identity() + 
            theme_classic() +
            labs(title = reg) + 
            theme(legend.position = legend,line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
        return(q)
        })

    qlist <- lapply(region_order, function(reg) {
        p <- ggplot(data=subset(plot_data_infunction, region == reg), aes_string(x="tSNE_1", y="tSNE_2"))
        q <- p + stat_density2d(aes(fill = ..density..^density_sf), geom = "tile", contour = add_contour, n = nbins) +
            scale_fill_gradientn(colours = density_color) + 
            theme_classic() +
            labs(title = reg) + 
            theme(legend.position = legend,line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
        return(q)
        })


    if (output_ggplot){
        return(list(plist, qlist))
    } else {
        pdf(paste0(output_dir,file_name, "_", dr_name[1], "_density.smooth.pdf"),pdf_size[1], pdf_size[2])
        if (is.null(nrow)){
            ppp <- plot_grid(plotlist = plist, nrow = 1, ncol = length(region_order))
            qqq <- plot_grid(plotlist = qlist, nrow = 1, ncol = length(region_order))
        } else {
            ppp <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol)
            qqq <- plot_grid(plotlist = qlist, nrow = nrow, ncol = ncol)
        }
        
        print(ppp)
        print(qqq)
        dev.off()
        return()
    }
}




FastInte <- function(object, split.by, nfeatures = 2000, inte.dims = 20, cls.dims = 25) {
    ## Do the integration for the two species.
    seu_list <- SplitObject(object, split.by = split.by) %>%
                    lapply(., function(x) FindVariableFeatures(x, nfeatures = nfeatures))
    hvg_use <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = nfeatures)

    ## Do the CCA integration
    data.anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:inte.dims, assay = NULL, anchor.features = hvg_use)
    data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:inte.dims)
    DefaultAssay(data.integrated) <- "integrated"
    data.integrated <- ScaleData(data.integrated, verbose = FALSE) %>%
                        RunPCA(., npcs = 50, verbose = FALSE) %>%
                        RunUMAP(., dims = 1:cls.dims) %>%
                        FindNeighbors(., dims = 1:cls.dims) %>%
                        FindClusters(., resolution = 1.2)
    return(data.integrated)
}


###############################################################################################

## Dot Plot (modified version)

###############################################################################################
PointPlot <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, shape = 16, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (!is.null(x = group.by)) {
        Idents(object) <- group.by
    }

    data.features$id <- Idents(object = object)
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    ## Add a column for split
    data.plot$species <- extract_field(as.character(data.plot$id), -1, "_")
    data.plot$cluster <- extract_field(as.character(data.plot$id), "rm_end", "_")


    ## Define the order the species & cluster
    if (!is.null(species_order)){
        data.plot$species <- factor(data.plot$species, levels = species_order)
    } 
    if (!is.null(cluster_order)){
        data.plot$cluster <- factor(data.plot$cluster, levels = cluster_order)
    }
    

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
        y = "cluster")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by), shape = shape) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(PointPlot) <- environment(DotPlot)


## PointPlot2 is designed for stacked plots
##Change 1: data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features)), remove rev
##Change 2: ggplot, y axis is "species"
PointPlot2 <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, shape = 16, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (!is.null(x = group.by)) {
        Idents(object) <- group.by
    }

    data.features$id <- Idents(object = object)
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    ## Add a column for split
    data.plot$species <- extract_field(as.character(data.plot$id), -1, "_")
    data.plot$cluster <- extract_field(as.character(data.plot$id), "rm_end", "_")


    ## Define the order the species & cluster
    if (!is.null(species_order)){
        data.plot$species <- factor(data.plot$species, levels = species_order)
    } 
    if (!is.null(cluster_order)){
        data.plot$cluster <- factor(data.plot$cluster, levels = cluster_order)
    }
    

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "species",
        y = "cluster")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by), shape = shape) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(PointPlot2) <- environment(DotPlot)



SpeciesDot <- function(object, features, cols = color_list$sp, dot.scale = 5, dot.min = 0.05, group.by = "hres", split.by = "species", shape = 16, layout_type = c("v", "h", "hflip")[1], panel.space = 0, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), scale.by = "radius") {

    if (!is.null(cluster_order) && layout_type %in% c("v", "h")){
        cluster_order <- rev(cluster_order)
    }

    if (layout_type %in% c("h")){
        features <- rev(features)
    }

    p <- PointPlot(object, assay = "RNA", features = features, cols = cols, dot.scale = dot.scale, dot.min = dot.min, group.by = group.by, split.by = split.by, shape = shape, cluster_order = cluster_order, species_order = species_order, scale.by = scale.by)

    if (layout_type == "v"){
        p <- p + facet_wrap(vars(species), nrow = length(species_order), ncol = 1, scales = "free_y") + 
                    scale_x_discrete(position = "top") + 
                    theme(axis.text.x=element_text(size = 8, hjust = 0, angle = 45), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(panel.space, "in"))
    } else if (layout_type == "hflip") {
        p <- p + coord_flip() + 
                    RotatedAxis() + 
                    facet_wrap(vars(species), nrow = 1, ncol = length(species_order)) + 
                    theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(panel.space, "in"))
    } else if (layout_type == "h") {
        p <- p + RotatedAxis() + 
                    facet_wrap(vars(species), nrow = length(species_order), ncol = 1) + 
                    theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(panel.space, "in"))

    } else {
        stop("layout_type should be v, h or hflip")
    }
    return(p)
}
##, scales = "free_x"


SplitGenePointPlot <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, shape = 16, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (!is.null(x = group.by)) {
        Idents(object) <- group.by
    }

    data.features$id <- Idents(object = object)
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        if (is.null(names(cols))){
            names(x = cols) <- unique(x = splits)
        } 
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled

    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    ## Add a column for split
    data.plot$species <- extract_field(as.character(data.plot$id), -1, "_")
    data.plot$cluster <- extract_field(as.character(data.plot$id), "rm_end", "_")


    ## Set the new features
    all_sps <- levels(as.factor(object@meta.data[, split.by]))
    all_sps <- species_order[species_order %in% all_sps]
    new_features <- paste0(rep(all_sps, lengt.out = length(features) * length(all_sps)), "_", rep(features, each = length(all_sps)))
    print(new_features)
    print(levels(as.factor(data.plot$features.plot)))

    data.plot$features.plot <- factor(x = paste0(data.plot$species, "_", as.character(data.plot$features.plot)),
        levels = rev(x = new_features))

    ## Define the order the cluster
    if (!is.null(cluster_order)){
        data.plot$cluster <- factor(data.plot$cluster, levels = cluster_order)
    }
    

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
        y = "cluster")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by), shape = shape) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(SplitGenePointPlot) <- environment(DotPlot)

