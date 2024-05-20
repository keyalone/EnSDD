# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:12:14 2024

@author: lihs
"""

import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


def visualization_clustering_shiny(counts_path, meta_path, img_path, cluster_num = 10, 
                                   spot_diameter_fullres = 150, spot_size = 1.5, 
                                   visual_base = True):
    adata = ad.read_text(counts_path).T
    meta_data = pd.read_csv(meta_path, sep="\t")
    
    for col in meta_data.columns.values:
        meta_data[col] = meta_data[col].astype("category")
    
    adata.obs = meta_data
    
    adata.obsm["spatial"] = meta_data.loc[:,["pixel_y", "pixel_x"]].values
    
    
    spatial_key = "spatial"
    library_id = 'HBC'
    adata.uns.setdefault(spatial_key, {}).setdefault(library_id, {})
    adata.uns[spatial_key][library_id]["images"] = {"hires": None}
    adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 1.,
                                                          "spot_diameter_fullres": spot_diameter_fullres}
    
    
    
    #### visualization of ground-truth
    plt.rcParams["figure.figsize"] = (8, 4)
    sc.pl.spatial(adata, img_key = "hires", color = 'EnSDD', size = spot_size, show=False, 
                  save = False, legend_loc='on data', legend_fontsize = 8, return_fig = True,
                  legend_fontweight = 'bold')
    plt.savefig("EnSDD.png", bbox_inches='tight', dpi=300)
    plt.close()

    

    
    #plt.show()
    #### Sometimes need close more than one times.
    #plt.close()
    
    if visual_base:
        methods_name_inte = meta_data.columns.values
        all_methods_name = np.array(['BayesSpace', 'DR.SC', 'SpaGCN', 'stLearn', 'GraphST', 'STAGATE'])
        methods_20 = np.intersect1d(methods_name_inte, all_methods_name)
        fig, axs = plt.subplots(1, len(methods_20), figsize = (10,3))
        for col_idx, mtd in enumerate(methods_20):
            ax = axs[col_idx]
            sc.pl.spatial(adata, img_key = "hires", color = mtd, size = spot_size, 
                          title = methods_20[col_idx], legend_loc='on data', legend_fontsize = 5,
                          legend_fontweight = 'bold', ax = ax, save=False, show = False)
            ax.set(xlabel=None)
            ax.set(ylabel = None)
            ax.set_title(methods_20[col_idx], size = 5)
        plt.subplots_adjust(wspace=0.05)
        plt.tight_layout()
        fig.savefig("base_methods.png", dpi = 300)  
        plt.close()
        
        
            
            
        
        
    
    
