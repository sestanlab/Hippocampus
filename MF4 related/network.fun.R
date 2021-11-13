
plot.simi.network <- function(marker.list, thre.quantile = 0.92, seed.use = 0, max.marker = 100, dot.size = 8, file_name = paste0("MF3.marker.simi.network.", ctp)){
	ji <- function(x, y) {
		gene1 <- names(x)[order(x, decreasing = TRUE)[1:min(length(x), max.marker)]]
		gene2 <- names(y)[order(y, decreasing = TRUE)[1:min(length(x), max.marker)]]
		return(length(intersect(gene1, gene2))/length(union(gene1, gene2)))
	}

	format_simi <- function(input_simi, nodes_meta = NULL, name_prefix = ""){
	    ## Organize the results in a table format
	    links <- as.matrix(input_simi) %>% 
	    			reshape2::melt() %>% 
	    			setNames(., c("from", "to", "weight")) %>%
	        		mutate(from = as.character(from), to = as.character(to)) %>%
	        		subset(!is.na(weight))


	    ## If nodes_meta is not provided, will generate one
	    new_cluster <- colnames(input_simi)[colnames(input_simi) %in% union(unique(links$from), unique(links$to))]
	    nodes_meta <- data.frame(cluster = paste0(name_prefix, new_cluster), stringsAsFactors = FALSE) %>% 
	                        mutate(cls_size = 5)
	    rownames(nodes_meta) <- nodes_meta$cluster

	    
	    simi_res <- list(link = links, node = nodes_meta, raw = as.matrix(input_simi))
	    return(simi_res)
	}


	## Get the jaccard index values
    simi_mat <- lapply(marker.list, function(m1)
        sapply(marker.list, function(m2) ji(m1, m2)) %>% setNames(., names(marker.list))) %>% as.data.frame(., check.names = FALSE) %>% as.matrix()
    diag(simi_mat) <- NA
    simi_thre <- quantile(simi_mat[grep("^PFC", colnames(simi_mat), value = TRUE), grep("^MTG", colnames(simi_mat), value = TRUE)], na.rm = TRUE, probs = thre.quantile)
    simi_mat[upper.tri(simi_mat, diag = TRUE)] <- NA


    ## Format the similarity matrix
    simi_res <- format_simi(input_simi = simi_mat, nodes_meta = NULL, name_prefix = "")


	reg_color <- c("#FF420E","#89DA59","#4CB5F5","#FFBB00") %>% 
					setNames(., c("HIP", "EC", "PFC", "MTG"))
	net <- graph_from_data_frame(d=simi_res$link, vertices=simi_res$node, directed=FALSE)

	## Vertices parameters
	V(net)$color <- extract_field(simi_res$node$cluster, 1, "|") %>% reg_color[.]
	V(net)$size <- dot.size

	## Edge parameters
	width_scale <- 8
	E(net)$width <- E(net)$weight * width_scale

	## Remove certain edges
	##simi_thre <- 0.1111111##0.11238001##0.09890110
	net <- delete.edges(net, which(E(net)$weight < simi_thre))
	set.seed(seed.use)
	lo <- do.call("layout_with_fr", list(net)) 


	pdf(paste0(outputdir, file_name, ".pdf"), width = 10, height = 10)
	plot(net, edge.curved=0, layout=lo, vertex.label=NA, vertex.label.cex= 0.5, add=FALSE, asp = 1, vertex.frame.color="#000000", edge.color = "#000000")
	legend(x = "bottomright", legend=as.character(c(0.1, 0.3, 0.5, 0.8)),lwd=c(0.1, 0.3, 0.5, 0.8) *width_scale , cex=1.25 * 1, col=rep(c("black"), 5), lty=rep(1, 5), title = "Similarity")
	plot(net, edge.curved=0, layout=lo, vertex.label=simi_res$node$cluster, vertex.label.cex= 0.7, add=FALSE, asp = 1, vertex.frame.color="#000000", edge.color = "#000000")
	dev.off()
} 


