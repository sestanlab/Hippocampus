args <- commandArgs(trailingOnly = TRUE)
source("../scripts/hip.fun.R")



get_db <- function(object, file_name, input_dir = inputdir, output_dir = outputdir){
	##Execute the python. scrublet scripts to detect the doublets
	countm <- object$RNA@counts

	## Write the data for scrublet to read
	writeMM(countm, file = paste0(input_dir, file_name, "_matrix.mtx"))
	data.frame(ids = rownames(countm), name = rownames(countm), stringsAsFactors = FALSE) %>% 
		write.table(., file = paste0(input_dir, file_name, "_genes.tsv"), sep = "\t", col.names = FALSE, row.names = FALSE, quote= FALSE)

	## Run scrublet
	## scrublet is copied from: /home/sm2726/project/InN_evolution/scripts/scrublet.doublets.py
	paste0("python ./calc_doublet.py ", input_dir, file_name, "_matrix.mtx ", input_dir, file_name, "_genes.tsv ", output_dir, " ", file_name, " > ", input_dir, file_name, "_scrublet_out.txt") %>% 
	system(.)
	file.remove(paste0(input_dir, file_name, "_matrix.mtx"));
	file.remove(paste0(input_dir, file_name, "_genes.tsv"));


	##add the scrublet results to the seurat object
	object@meta.data[, c("doublet_scores", "doublet_assign")] <- read.table(paste0(output_dir, file_name, "_b_scrublet_out.tsv"), sep = "\t", header = FALSE) %>% t()

	return(object)
}


sp <- paste0("Batch", c(2,3,5,6))[as.numeric(args[1])]
allgc <- readRDS(file = paste0(inputdir, "Konopka_HIP_seu.rds"))
sp_seu <- subset(allgc, batch == sp) %>%
			get_db(object = ., file_name = paste0("Konopka_", sp), input_dir = inputdir)
meta_use <- sp_seu@meta.data[, c("doublet_scores", "doublet_assign")]
saveRDS(meta_use, file = paste0(inputdir, "Konopka_", sp, "_db_res.rds"))


## Summarize the db scores from all batches
## Only run once, manually 
if (FALSE){
	all_bch <- paste0("Batch", c(2,3,5,6))
	meta <- lapply(all_bch, function(sp) readRDS(file = paste0(inputdir, "Konopka_", sp, "_db_res.rds"))) %>%
					do.call(rbind,.)
	saveRDS(meta, file = paste0(inputdir, "Konopka_all_db_res.rds"))
}












