# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:33:23 2023

@author: lihs
"""


import pandas as pd
import scanpy as sc
import SpaGCN as spg
import random, torch
import numpy as np
import cv2
import os

## b: represents the area of each spot when extracting color intensity.
## a: determines the weight given to histology when calculating Euclidean distance between every two spots.
## ‘a = 1’ means that the histology pixel intensity value has the same scale variance as the (x,y) coordinates,
## whereas higher value of ‘a’ indicates higher scale variance, hence,higher weight to histology, when 
## calculating the Euclidean distance.
## p: represents the percentage of total expression contributed by neighborhoods.
def run_spaGCN(counts_path, meta_path, img_path, platform = "Visium", res_setting_by_user = None, 
                b = 49, a = 1, p = 0.5, n_clusters = 7, lr = 0.05, max_epochs = 200):
    adata = sc.read_text(counts_path).T
    meta_data = pd.read_csv(meta_path, sep="\t")
    adata.var_names_make_unique()
    adata.obs = meta_data
    
    adata.var_names_make_unique()
    spg.prefilter_genes(adata, min_cells = 3)
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    
    adata.obs["x_array"] = adata.obs["x"]
    adata.obs["y_array"] = adata.obs["y"]
    adata.obs["x_pixel"] = adata.obs["pixel_x"]
    adata.obs["y_pixel"] = adata.obs["pixel_y"]
    
    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()
    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()


    img = cv2.imread(img_path)
    
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel,
                                 image=img, beta=b, alpha=a, histology=True)
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    
    r_seed=t_seed=n_seed=100
    if res_setting_by_user is None:
        print("The user need input the n_clusters for the determination of resolution if the res_setting_by_user is None...")
        res=spg.search_res(adata, adj, l, n_clusters, start = 0.7, step=0.05, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    else:
        print("The user don't need set the n_clusters if set the res_setting_by_user...")
        res = res_setting_by_user
    
    
    ### 4.3 Run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=lr, max_epochs=max_epochs)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    
    if platform == "Visium":
        shape = "hexagon"
        
    if platform == "ST":
        shape = "square"
    
    #Do cluster refinement(optional)
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape=shape)
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    
    temp = adata.obs["refined_pred"].cat.codes
    LABEL_PATH = os.path.dirname(img_path) + "/spaGCN_res"
    os.makedirs(LABEL_PATH, exist_ok=True)
    LABEL_FILE = LABEL_PATH + "/spaGCN_label.txt"
    temp.to_csv(LABEL_FILE, header = None, sep = "\t")
    
    
    
    
    
    
    
    
    
    