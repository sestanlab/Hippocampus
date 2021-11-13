'''import required modules'''
import scvelo as scv
import pandas as pd
out_h5ad = "./load_files/velo_Human_ana.h5ad"
adata = scv.read(out_h5ad)


xran = [-5.910893, 17.496210]
yran = [-4.438446, 10.984851]
##scv.pl.velocity_embedding_stream(adata, basis="umap", save = "scVelo_final_Human_stream.png", dpi=300, color = "cluster", xlim = xran, ylim = yran, size = 2.5, fontsize = 0, title = "", alpha = 1, legend_loc='none', linewidth = 0.5, density = 0.8, figsize = (7,7))
##scv.pl.velocity_embedding(adata, basis='umap', save = "scVelo_final_Human_arrow.png", dpi=300, arrow_length=2, arrow_size=1.5, color = "cluster", xlim = xran, ylim = yran, fontsize = 0, title = "", alpha = 1, legend_loc='none', size = 2.5, figsize = (7,7))
scv.pl.velocity_embedding_grid(adata, basis='umap',save = "scVelo_final_Human_grid.png", dpi=300,  arrow_length = 1.5, arrow_size = 1.5, color = "cluster", xlim = xran, ylim = yran, fontsize = 0, title = "", alpha = 1, legend_loc='none', size = 5, density = 0.8, figsize = (7,7), min_mass = 5, linewidth = 0.35)

