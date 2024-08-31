# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:12:14 2024

@author: lihs
"""

import anndata as ad
import pandas as pd
import scanpy as sc
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt


def visualization_clustering(counts_path, meta_path, img_path, cluster_num = 10, 
                           spot_diameter_fullres = 150, spot_size = 1.5):
    adata = ad.read_text(counts_path).T
    meta_data = pd.read_csv(meta_path, sep="\t")
    
    for col in meta_data.columns.values:
        meta_data[col] = meta_data[col].astype("category")
    meta_data['ground_truth'] = meta_data['ground_truth'].astype("category")
    
    adata.obs = meta_data
    
    adata.obsm["spatial"] = meta_data.loc[:,["pixel_y", "pixel_x"]].values
    
    
    # Image.MAX_IMAGE_PIXELS=None
    # img = Image.open(img_path)
    # img = np.array(img)
    
    spatial_key = "spatial"
    library_id = 'HBC'
    adata.uns.setdefault(spatial_key, {}).setdefault(library_id, {})
    adata.uns[spatial_key][library_id]["images"] = {"hires": None}
    adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 1.,
                                                          "spot_diameter_fullres": spot_diameter_fullres}
    
    #### no pic plot
    ##### no pic plot
    
    
    
    #### visualization of ground-truth
    plt.rcParams["figure.figsize"] = (8, 4)
    # file_name_gt = os.path.join(path_no_pic, f'{col_gt}_no_pic.pdf')
    sc.pl.spatial(adata, img_key = "hires", color = 'EnSDD', size = spot_size, show=False, 
                  return_fig = True, legend_loc='on data', legend_fontsize = 8,
                  legend_fontweight = 'bold')
    plt.show()
    #### Sometimes need close more than one times.
    plt.close()
    
    
        
    # GT_file = os.path.join(folder_path, 'GT.pdf')
    
    # methods_name_inte = ['BayesSpace', 'DR.SC', 'GraphST', 'SpaGCN', 'stLearn', 'STAGATE', 'SiGra', 'spaVAE', 'EnSDD']
    # 
    # if cluster_num == 20:
    #     methods_20 = methods_name_inte[0:8]
    #     fig, axs = plt.subplots(1, 8, figsize = (10,3))
    #     for col_idx, mtd in enumerate(methods_20):
    #         ax = axs[col_idx]
    #         sc.pl.spatial(adata, img_key = "hires", color = mtd, size = spot_size, 
    #                       title = methods_20[col_idx], legend_loc= None,ax = ax, save=False, show = False)
    #         ax.set(xlabel=None)
    #         ax.set(ylabel = None)
    #         ax.set_title(methods_name_inte[col_idx], size = 5)
    #     plt.subplots_adjust(wspace=0.1)
    #     plt.tight_layout()
    #     plt.show()
    #     
    #                 
    # else:
    #     fig, axs = plt.subplots(1, 9, figsize = (9.5,3))
    #     for col_idx, mtd in enumerate(methods_name_inte):
    #         ax = axs[col_idx]
    #         if col_idx == 6:
    #             # ax = fig.add_subplot(fig_left[0,0])
    #             sc.pl.spatial(adata, img_key = "hires", color = mtd, size = spot_size, 
    #                           ax=ax, save=False)
    #             ax.set(xlabel=None)
    #             ax.set(ylabel = None)
    #             ax.set_title(methods_name_inte[col_idx], size = 5)
    #             ax.legend(fontsize = 4.5, markerscale = 0.5, ncol = 2,loc='center left', bbox_to_anchor=(1, 0.5))
    #             # ax.legend(fontsize = 4.5, markerscale = 0.5)
    #         else:
    #             # ax = fig.add_subplot(fig_right[0, 0])
    #             sc.pl.spatial(adata, img_key = "hires", color = mtd, size = spot_size, 
    #                           title = methods_name_inte[col_idx], legend_loc= None,ax = ax, save=False)
    #             ax.set(xlabel=None)
    #             ax.set(ylabel = None)
    #             ax.set_title(methods_name_inte[col_idx], size = 5)
    #     plt.subplots_adjust(wspace=0.05)
    #     plt.tight_layout()
