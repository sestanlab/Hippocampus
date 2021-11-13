## Plot DCX expression along the developmental periods
source("../scripts/hip.fun.R") 


load(file = paste0(inputdir, "NHP.dev.Rdata"))
#(rpkm, logrpkm, meta, )


gene <- "DCX" %>%
			paste0("\\|", ., "$") %>%
			grep(., rownames(logrpkm), value = TRUE)
plot_data <- data.frame(samplename = rownames(meta), 
						exp = logrpkm[gene, rownames(meta)],
						species = meta$Species, 
						age = meta$Predicted.age..PC.Days., 
						stringsAsFactors = FALSE) %>%
				mutate(age = log2(age))


tmp1=c(log2(c(69,111,132,167,447,1299,4648,7570)),max(plot_data$age));
tmp2=c(min(plot_data$age),log2(c(69,111,132,167,447,1299,4648,7570)));
text_data <- data.frame(age = (tmp1 + tmp2)/2, 
					label = 1:length(tmp1),
					exp = max(plot_data$exp),
					stringsAsFactors = FALSE)


p <- ggplot(data = plot_data, aes_string(x = "age", y = "exp")) + 
	  			##geom_point(size = 3, shape = 21, fill = "white", alpha = 0.75) +
	  			geom_smooth(aes_string(color = "species"), se = TRUE, method = "loess", span = 0.8, size = 1.2) + 
	  			geom_text(data = text_data, mapping = aes_string(x = "age", y = "exp", label = "label")) +
	  			scale_color_manual(name = "species", values=c(Human = "#FF420E", Macaque = "#89DA59")) + 
	  			theme_bw() + 
	  			labs(x = "Developmental Time", y = "log2(RPKM+ 1)") + 
	  			##RotatedAxis() + 
	  			##coord_capped_cart(bottom='both', left='both') +
	  			geom_vline(xintercept = log2(c(69.5,111.5,132.5,167.5,447.5,1299.5,4648.5,7570.5)), linetype = "dashed", size = 0.25) +
	  			theme(legend.position="bottom", panel.grid.major = element_blank(),
	                    panel.grid.minor = element_blank(), 
	                    axis.title = element_text(size = 10, hjust = 0.5), axis.text = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5),
	                    legend.title = element_blank(), legend.text = element_text(size = 8), axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25)) + 
	            guides(fill=guide_legend(ncol = 1))

pdf(paste0(outputdir, "SF2_DCX_developmental.HR.pdf"), width = 5, height = 5)
print(p)
dev.off()



