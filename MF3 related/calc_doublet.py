#matplotlib inline
#Format of the command line: 
#$python b.remove_doublets.py count_inputdir gene_inputdir outputdir file_name
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import gzip

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42


def load_gzip_genes(filename, delimiter='\t', column=0, skip_rows=0):
	gene_list = []
	gene_dict = {}

	with gzip.open(filename, "rt") as f:
		for iL in range(skip_rows):
			f.readline()
		for l in f:
			gene = l.strip('\n').split(delimiter)[column]
			if gene in gene_dict:
				gene_dict[gene] += 1
				gene_list.append(gene + '__' + str(gene_dict[gene]))
				if gene_dict[gene] == 2:
					i = gene_list.index(gene)
					gene_list[i] = gene + '__1'
			else: 
				gene_dict[gene] = 1
				gene_list.append(gene)
	return gene_list



#input_dir = './filtered_gene_bc_matrices/GRCh38/'
#counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
#genes = np.array(scr.load_genes(input_dir + 'genes.tsv', delimiter='\t', column=1))
if ".gz" in sys.argv[2]:
	genes = np.array(load_gzip_genes(sys.argv[2], delimiter='\t', column=1))
else:
	genes = np.array(scr.load_genes(sys.argv[2], delimiter='\t', column=1))
counts_matrix = scipy.io.mmread(sys.argv[1]).T.tocsc()







print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06, sim_doublet_ratio = 2) #sim_doublet_ratio, n_neighbors = 0.5*sqrt(ncells) are default values

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

scrub.plot_histogram();
plt.savefig(sys.argv[3] + sys.argv[4] + "_b_scrublet_plot.pdf")

#np.savetxt("doublet_score.tsv", doublet_scores, delimiter="\t")
#np.savetxt("doublet_assignment.tsv", predicted_doublets, delimiter="\t")

#Only two rows with first row as the score and the second row as assignment of doublets. 
new_scores = np.concatenate((np.array([doublet_scores]), np.array([predicted_doublets])), axis = 0)
np.savetxt(sys.argv[3] + sys.argv[4] + "_b_scrublet_out.tsv", new_scores, delimiter="\t")


