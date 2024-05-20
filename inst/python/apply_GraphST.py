# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 15:34:16 2023

@author: lihs
"""

import pandas as pd
import ot
import stlearn as st
import scanpy as sc
import os
from GraphST import GraphST 
import torch
from GraphST.utils import clustering
from sklearn.decomposition import PCA

#lamda1 : float, optional Weight factor to control the influence of reconstruction loss in mapping matrix learning. 
# The default is 10.
# lambda2 : float, optional Weight factor to control the influence of contrastive loss in mapping matrix learning. 
# The default is 1.
# tools: 'leiden', 'louvain', 'mclust'
def run_GraphST(counts_path, meta_path, img_path, n_setting = 7,lambda1 = 10, lambda2 = 1, 
                tool = "louvain", radius = 50, n_PCs = 20, refinement = False):
    # get the device id use the torch.cuda.device_count()
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
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
    
    st.pp.filter_genes(adata,min_cells=3)
    GraphST.preprocess(adata)
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()
    
    if tool in ['louvain', 'leiden']:
        clustering(adata, 7, radius=radius, method=tool, refinement=refinement)
        temp = adata.obs["domain"].cat.codes
        LABEL_PATH = os.path.dirname(img_path) + "/GraphST_res"
        os.makedirs(LABEL_PATH, exist_ok=True)
        LABEL_FILE = LABEL_PATH + "/GraphST_louvain.txt"
        temp.to_csv(LABEL_FILE, header = None, sep = "\t")
    elif tool == 'mclust':
        pca = PCA(n_components=n_PCs, random_state=42) 
        embedding = pca.fit_transform(adata.obsm['emb'].copy())
        adata.obsm['emb_pca'] = embedding
        var_names = ["PC_" + str(i) for i in range(1,n_PCs+1,1)]
        df = pd.DataFrame(embedding, columns = var_names, index=adata.obs_names.values)
        LABEL_PATH = os.path.dirname(img_path) + "/GraphST_res"
        os.makedirs(LABEL_PATH, exist_ok=True)
        LABEL_FILE = LABEL_PATH + "/GraphST_mclust.txt"
        df.to_csv(LABEL_FILE, sep = "\t")
    
    
        
    
    
