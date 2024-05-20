# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 15:34:16 2023

@author: lihs
"""

import pandas as pd
import stlearn as st
import scanpy as sc
import os
import STAGATE
from sklearn.decomposition import PCA
import warnings
import tensorflow as tf
warnings.filterwarnings("ignore")



#lamda1 : float, optional Weight factor to control the influence of reconstruction loss in mapping matrix learning. 
# The default is 10.
# lambda2 : float, optional Weight factor to control the influence of contrastive loss in mapping matrix learning. 
# The default is 1.
# tools: 'leiden', 'louvain', 'mclust'
def run_STAGATE (counts_path, meta_path, img_path, alpha = 0):
    
    
    adata = sc.read_text(counts_path).T
    meta_data = pd.read_csv(meta_path, sep="\t")
    adata.var_names_make_unique()
    adata.obs = meta_data
    
    
    imgcol = meta_data.loc[:,"pixel_x"]
    imgrow = meta_data.loc[:,"pixel_y"]
    adata.obs["imagecol"] = imgrow
    adata.obs["imagerow"] = imgcol
    adata.obs["array_row"] = meta_data.loc[:,"x"]
    adata.obs["array_col"] = meta_data.loc[:,"y"]
    st.add.image(adata,library_id="151673",quality="fulres", 
                 imgpath=img_path,scale=1,spot_diameter_fullres=150)
    
    #Normalization
    st.pp.filter_genes(adata,min_cells=3)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    
    STAGATE.Cal_Spatial_Net(adata, rad_cutoff = 150)
    #STAGATE.Stats_Spatial_Net(adata)
     
    tf.compat.v1.disable_eager_execution()
    adata = STAGATE.train_STAGATE(adata, alpha = alpha)
    
    sc.pp.neighbors(adata, use_rep='STAGATE')
    # sc.tl.umap(adata) 
    temp_np = adata.obsm['STAGATE']
    var_names = ["STAGATE_" + str(i) for i in range(1, temp_np.shape[1] + 1,1)]
    df = pd.DataFrame(temp_np, columns = var_names, index = adata.obs_names.values)
    LABEL_PATH = os.path.dirname(img_path) + "/STAGATE_res"
    os.makedirs(LABEL_PATH, exist_ok=True)
    LABEL_FILE = LABEL_PATH + "/STAGATE_mclust.txt"
    df.to_csv(LABEL_FILE, sep = "\t")
    
    
    
        
    
    