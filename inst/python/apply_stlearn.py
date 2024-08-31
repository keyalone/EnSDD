# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 22:31:05 2023

@author: lihs
"""
import pandas as pd
import stlearn as st
import scanpy as sc
import os
import shutil

# normalize_type
#'weights_matrix_all', 'weights_matrix_pd_gd', 
#'weights_matrix_pd_md', 'weights_matrix_gd_md', 
#'gene_expression_correlation', 'physical_distance', 
#'morphological_distance'

## Weighting matrix for imputation.
# if `weights_matrix_all`, matrix combined all information from spatial location (S),
#tissue morphological feature (M) and gene expression (E)
#if `weights_matrix_pd_md`, matrix combined information from spatial location (S),
#tissue morphological feature (M)

def run_stlearn(counts_path, meta_path, img_path, n_comps = 30, n_setting = 7, 
                normalize_type = "physical_distance"):
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
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)
    TILL_PATH = os.path.dirname(img_path) + "/till"
    os.makedirs(TILL_PATH, exist_ok=True)
    st.pp.tiling(adata, TILL_PATH)
    st.pp.extract_feature(adata)
    shutil.rmtree(TILL_PATH) 
    
    st.em.run_pca(adata,n_comps=n_comps)
    
    st.spatial.SME.SME_normalize(adata, use_data="raw", weights = normalize_type)
    adata_SME = adata.copy()
    # apply stSME to normalise log transformed data
    adata_SME.X = adata_SME.obsm['raw_SME_normalized']
    st.pp.scale(adata_SME)
    st.em.run_pca(adata_SME,n_comps=n_comps)
    st.tl.clustering.kmeans(adata_SME,n_clusters=n_setting, use_data="X_pca", key_added="X_pca_kmeans", algorithm = 'lloyd')
    temp = adata_SME.obs["X_pca_kmeans"].cat.codes
    LABEL_PATH = os.path.dirname(img_path) + "/stLearn_res"
    os.makedirs(LABEL_PATH, exist_ok=True)
    LABEL_FILE = LABEL_PATH + "/stLearn_kmeans.txt"
    temp.to_csv(LABEL_FILE, header = None, sep = "\t")
    

    
    
