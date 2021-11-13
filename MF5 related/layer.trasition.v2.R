## Generate the Plot data for the MF4   
source("../scripts/hip.fun.R")
library(igraph)


##------------------------------------------------------------------
## Markers organization
all.mar <- readRDS(file = paste0("../MF3_cellregion/load_files/", "Markers.4reg.ExN.rds"))
cls.order <- list(HIP = c("DG GC PROX1 PDLIM5","DG GC PROX1 SGCZ","DG MC ARHGAP24 DLC1","CA3 CFAP299 SYN3","CA2 CFAP299 HGF","CA1 dorsal GRIK1 GRM3","CA1 ventral ACVR1C SYT13","SUB proximal ROBO1 COL5A2","SUB proximal ROBO1 SEMA3E","SUB distal FN1 NTNG1"),
					EC = c("EC L2 CUX2 CALB1","EC L2 CUX2 IL1RAPL2","EC L2 CUX2 LAMA3","EC L2 CUX2 PDGFD","EC L2 RELN BCL11B","EC L2 RELN BMPR1B","EC L3 PCP4 ADARB2","EC L5 RORB TLL1","EC L5 RORB TPBG","EC L6 THEMIS CDH13","EC L6 THEMIS RGS12","EC L5 BCL11B ADRA1A","EC L56 TLE4 NXPH2","EC L6 TLE4 SULF1","EC L6b TLE4 CCN2"), 
					MTG = c("Exc L2 LAMP5 LTK", "Exc L2-3 LINC00507 FREM3", "Exc L2-4 LINC00507 GLP2R", "Exc L3-4 RORB CARM1P1", "Exc L3-5 RORB FILIP1L", "Exc L3-5 RORB TWIST2", "Exc L3-5 RORB COL22A1", "Exc L3-5 RORB ESR1", "Exc L4-6 RORB C1R", "Exc L5-6 RORB TTC12", "Exc L4-5 RORB FOLH1B", "Exc L4-5 RORB DAPK2", "Exc L4-6 RORB SEMA3E", "Exc L5-6 THEMIS C1QL3", "Exc L5-6 THEMIS CRABP1", "Exc L5-6 THEMIS DCSTAMP", "Exc L5-6 THEMIS FGF10", "Exc L4-5 FEZF2 SCN4B", "Exc L5-6 FEZF2 ABO", "Exc L4-6 FEZF2 IL26", "Exc L5-6 SLC17A7 IL15", "Exc L6 FEZF2 OR2T8", "Exc L6 FEZF2 SCUBE1", "Exc L5-6 FEZF2 EFTUD1P1"), 
					PFC = readRDS(file = "~/project/PFC/data/PFC.cluster.order.rds")[c(1:14, 16:22, 15, 23:40)])


##1. Get the markers for all clusters in each region
mars <- lapply(names(cls.order), function(reg) {
	reg.mars <- lapply(cls.order[[reg]], function(x) {
		xx <- all.mar %>%
				subset(region == reg & cluster == x) %>%
				subset(pct.1 >= 0.2 & avg_logFC >= log(2) & ratio_fc >= 1.1) %>%
				mutate(p_val_adj = p.adjust(p_val, method = "bonferroni")) %>%
				subset(p_val_adj <= 0.05)
		yy <- setNames(xx$avg_logFC, xx$gene)
		return(yy)
		}) %>% setNames(., cls.order[[reg]])
	return(reg.mars)
	}) %>% setNames(., names(cls.order))


ji <- function(x, y) {
	max.marker = 100
	gene1 <- names(x)[order(x, decreasing = TRUE)[1:min(length(x), max.marker)]]
	gene2 <- names(y)[order(y, decreasing = TRUE)[1:min(length(x), max.marker)]]
	return(length(intersect(gene1, gene2))/length(union(gene1, gene2)))
}




reg_pairs <- c("HIP|EC", "EC|MTG", "MTG|PFC")
links <- lapply(reg_pairs, function(x) {
	reg1 <- strsplit(x, "|", fixed = TRUE)[[1]][1]
	reg2 <- strsplit(x, "|", fixed = TRUE)[[1]][2]

	aa <- matrix(0, nrow = length(mars[[reg1]]), ncol = length(mars[[reg2]]), dimnames = list(names(mars[[reg1]]), names(mars[[reg2]])))
	for (i in rownames(aa)){
		for (j in colnames(aa)){
			aa[i, j] <- ji(mars[[reg1]][[i]], mars[[reg2]][[j]])
		}
	}
	##aa[upper.tri(aa, diag = TRUE)] <- NA
	df <- reshape2::melt(aa) %>%
				setNames(., c("from", "to", "weight")) %>%
				mutate(from = as.character(from), to = as.character(to)) %>%
				subset(!is.na(weight))
	##df$weight <- sapply(1:nrow(df), function(idx) ji(mars[[reg1]][[df$from[idx]]], mars[[reg2]][[df$to[idx]]]))
	df$pair <- x
	return(df)
	}) %>% 
			do.call(rbind, .)





## Plot the similarity network
width_scale <- 15

set.seed(42)

new_data <- lapply(c("HIP|EC", "EC|MTG", "MTG|PFC"), function(pp) {
	reg1 <- strsplit(pp, "|", fixed = TRUE)[[1]][1]
	reg2 <- strsplit(pp, "|", fixed = TRUE)[[1]][2]
	submeta <- data.frame(row.names = paste0(pp, ".", c(cls.order[[reg1]], cls.order[[reg2]])), 
					cluster = paste0(pp, ".", c(cls.order[[reg1]], cls.order[[reg2]])), 
					pair = pp, 
					region = c(rep(reg1, length(cls.order[[reg1]])), rep(reg2, length(cls.order[[reg2]]))), 
					stringsAsFactors = FALSE)
	sublink <- links %>%
				subset(pair == pp) %>%
				mutate(from = paste0(pp, ".", from), to = paste0(pp, ".", to))

	## Remove the 80th quantiles
	thre <- quantile(sublink$weight, 0.85)
	print(thre)
	sublink <- subset(sublink, weight >= thre)
	sublink$weight <- MinMax(sublink$weight, min = 0, max = 0.5)


	locv <- list(c(1, 1.7), c(4.3, 5), c(7.7, 8.4)) %>%
				setNames(., c("HIP|EC", "EC|MTG", "MTG|PFC"))  %>%
				.[[pp]]

	maxy <- ifelse(reg1 == "HIP", 26, 40)
	ll <- data.frame(row.names = submeta$cluster, 
					xloc = rep(locv, c(sum(submeta$region == reg1), sum(submeta$region == reg2))), 
					yloc = c(rev(seq(1, maxy, length.out = sum(submeta$region == reg1))), rev(seq(1, 40, length.out = sum(submeta$region == reg2)))),  
					stringsAsFactors = FALSE) %>%
			as.matrix()
	labels <- extract_field(submeta$cluster, 2, ".")
	colvec <- c("#FF420E","#89DA59","#4CB5F5","#FFBB00") %>% 
				setNames(., c("HIP", "EC", "PFC", "MTG")) %>% 
				.[submeta$region]
	sizevec <- ifelse(submeta$region == reg2, 0, 14)
	if (reg1 != "HIP"){
		labels[submeta$region == reg1] <- NA
	}
	return(list(node = submeta, link = sublink, loc = ll, labels = labels, color = colvec, size = sizevec))
	})

new_meta <- lapply(new_data, function(x) x$node) %>% do.call(rbind, .)
new_link <- lapply(new_data, function(x) x$link) %>% do.call(rbind, .)
new_loc <- lapply(new_data, function(x) x$loc) %>% do.call(rbind, .)
new_labels <- lapply(new_data, function(x) x$labels) %>% do.call(c, .)
new_color <- lapply(new_data, function(x) x$color) %>% do.call(c, .)
new_size <- lapply(new_data, function(x) x$size) %>% do.call(c, .)


gcolors <- colorRampPalette(colors = c("lightgrey", "black"))(100)[c(1:5, seq(60, 100, 1))]
edge_color <- gcolors[as.numeric(cut(new_link$weight, breaks = seq(0, max(new_link$weight) + 0.01, length.out = length(gcolors) + 2), include.lowest = TRUE))]
net <- graph_from_data_frame(d= new_link, vertices = new_meta, directed=FALSE)
E(net)$width <- E(net)$weight * width_scale

pdf(paste0(outputdir, "MF4.HIP.ExN.reg.similarity.v2.pdf"), width = 15, height = 6)
par(oma = c(0, 0, 0, 5)) 
plot(net, edge.arrow.size=0.5, edge.curved=0, layout=new_loc, main="MF4A Layer Transition", vertex.frame.color=NA,
			vertex.shape = "square", vertex.size = new_size, 
			vertex.color = new_color, 
			edge.color = edge_color,  
			vertex.label = NA, 
			vertex.label.dist= 0, vertex.label.degree=0, 
			vertex.label.cex = 0.8,vertex.label.family="Helvetica", rescale=FALSE, xlim = c(0, 12), ylim = c(0, 42), asp = 0.2, add = FALSE)
points(x = rep(11.5, 40), y = 1:40, pch = 22, cex = 1.2, col = NA, bg = "#4CB5F5")
text(x = new_loc[names(new_labels), "xloc"], y = new_loc[names(new_labels), "yloc"], labels = new_labels, cex = 0.7, pos = rep(c(2, 4), c(10, length(new_labels)-10)), adj = 0)
segments(x0 = c(1,4.3,7.7,11.5), y0 = 0, y1 = c(27, 41,41,41), lwd = 4, col = c("#FF420E","#89DA59","#FFBB00","#4CB5F5"))
lg_val <- c(0.05, seq(0.1, 0.5, 0.2))
lg_colors <- gcolors[as.numeric(cut(lg_val, breaks = seq(0, max(new_link$weight) + 0.01, length.out = length(gcolors) + 2), include.lowest = TRUE))]
lg_colors[is.na(lg_colors)] <- rev(gcolors)[1]
legend(-1.2, 41, legend=as.character(lg_val), lwd=lg_val *(width_scale) , cex=1.2, col= lg_colors, lty=1, title = "", box.lwd = 0,box.col = "white",bg = "white")

##Add titles
segments(x0 = c(-2, 8), y0 = c(43, 43), x1 = c(7.4, 11.3), y1 = c(43, 43), lwd = 1, col = "black")
text(x = c(2.8, 9.7), y = c(44, 44), labels = c("Temporal lobe: allo-, meso- and neocortex", "Prefrontal neocortex"), cex = 1)
text(x = c(-0.3, 3, 6.3, 9.8), y = rep(42, 4), labels = c("HIP", "EC", "MTG", "PFC"))

##Add dashed lines
segments(x0 = rep(1.85, 3), y0 = c(22, 19.2, 13.6), x1 = rep(3.9, 3), y1 = c(22, 19.2, 13.6), lwd = 2, col = "#89DA59", lty= 2)
segments(x0 = rep(5.2, 3), y0 = c(35.85, 22.4, 13.8), x1 = rep(7.3, 3), y1 = c(35.85, 22.4, 13.8), lwd = 2, col = "#FFBB00", lty= 2)
segments(x0 = rep(8.6, 3), y0 = c(32.6, 20.6, 14.6), x1 = rep(11.2, 3), y1 = c(32.6, 20.6, 14.6), lwd = 2, col = "#4CB5F5", lty= 2)
dev.off()




















