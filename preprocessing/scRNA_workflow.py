#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################################
# 
# Copyright (c) 2023 biomap.com, Inc. All Rights Reserved
# 
########################################################################

import scanpy as sc
import time
from scipy.sparse import csc_matrix
import h5py
import pandas as pd
import numpy as np
from tqdm import tqdm
import math
from scipy import sparse

sc.settings.set_figure_params(dpi=80, facecolor='white')

def BasicFilter(adata,qc_min_genes=200,qc_min_cells=0,plot_show=False):
    """
    do basic filtering
    :param adata: adata object
    :param qc_min_genes: filter cell, which expresses genes number less than this paramters
    :param qc_min_cells: filter gene, which expressed in cells number less than this paramters
    :param plot_show: show the plot not not
    :return: adata object, and plot stored in figures sub-folder
    """

    sc.pl.highest_expr_genes(adata, n_top=20, save='_20TopGene.png' ,show=plot_show)

    print('Before filter, %d Cells, %d Genes' % (adata.shape))
    sc.pp.filter_cells(adata, min_genes=qc_min_genes)
    sc.pp.filter_genes(adata, min_cells=qc_min_cells)
    print('After filter, %d Cells, %d Genes' % (adata.shape))

    return adata

def QC_Metrics_info(adata,doublet_removal=False,plot_show=False):
    """
    display quality control plot, for QC parameters selection
    :param adata: adata object
    :param plot_show: show the plot not not
    :return: adata object, and plot stored in figures sub-folder
    """
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True,save='_QC_guide.png',show=plot_show)

    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',save='_pct_counts_mt.png',show=plot_show)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',save='_n_genes_by_counts.png',show=plot_show)
    return adata

def save_adata_h5ad(adata,out_path='None',shrink_to_sparse=False):
    """
    save adata with h5ad format
    :param adata: adata object to save
    :param out_path: path of saved adata
    :return: None
    """
    if shrink_to_sparse:

        if adata.raw is None:
            # convert adata X only
            adata.X = sparse.csr_matrix(adata.X)  
        else:
            # convert adata X
            adata.X = sparse.csr_matrix(adata.X)
            # convert adata.raw X
            adata_raw_tmp=adata.raw.to_adata()
            adata_raw_tmp.X=sparse.csr_matrix(adata_raw_tmp.X)
            adata.raw=adata_raw_tmp

        adata.write_h5ad(out_path)
    else:
        adata.write_h5ad(out_path)
    print('Current data saved')

def read_adata_h5ad(h5ad_path='None', expand_from_sparse=False, N=None):
    """
    Legge un file h5ad e opzionalmente seleziona solo le prime N righe.
    
    :param h5ad_path: Percorso del file h5ad.
    :param expand_from_sparse: Se True, converte da sparse a array denso.
    :param N: Numero di righe da selezionare (None per caricare tutto).
    :return: Oggetto adata ridotto con le prime N righe (se specificato).
    """
    
    adata = sc.read_h5ad(h5ad_path)

    if N is not None:
        adata = adata[:N].copy()  # Prendi solo le prime N righe

    if expand_from_sparse:
        if adata.raw is None:
            adata.X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        else:
            adata.X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
            adata_raw_tmp = adata.raw.to_adata()
            adata_raw_tmp.X = adata_raw_tmp.X.toarray() if sparse.issparse(adata_raw_tmp.X) else adata_raw_tmp.X
            adata.raw = adata_raw_tmp

    return adata

def main_gene_selection(X_df, gene_list):
    """
    Describe:
        rebuild the input adata to select target genes encode protein 
    Parameters:
        adata->`~anndata.AnnData` object: adata with var index_name by gene symbol
        gene_list->list: wanted target gene 
    Returns:
        adata_new->`~anndata.AnnData` object
        to_fill_columns->list: zero padding gene
    """
    to_fill_columns = list(set(gene_list) - set(X_df.columns))
    padding_df = pd.DataFrame(np.zeros((X_df.shape[0], len(to_fill_columns))), 
                              columns=to_fill_columns, 
                              index=X_df.index)
    X_df = pd.DataFrame(np.concatenate([df.values for df in [X_df, padding_df]], axis=1), 
                        index=X_df.index, 
                        columns=list(X_df.columns) + list(padding_df.columns))
    X_df = X_df[gene_list]
    
    var = pd.DataFrame(index=X_df.columns)
    var['mask'] = [1 if i in to_fill_columns else 0 for i in list(var.index)]
    return X_df, to_fill_columns,var

def main():
    data_path = r"D:\Repositories\scFoundationModule-main\data\raw\data_yuto.h5ad"
    N=100000
    save_path = r'D:\Repositories\scFoundationModule-main\data\interim\data_yuto.h5ad'

    adata = read_adata_h5ad(data_path, N=N)

    gene_list_df = pd.read_csv(r'D:\Repositories\scFoundationModule-main\preprocessing\OS_scRNA_gene_index.19264.tsv', header=0, delimiter='\t')
    gene_list = list(gene_list_df['gene_name'])

    X_df= pd.DataFrame(sparse.csr_matrix.toarray(adata.X),index=adata.obs.index.tolist(),columns=adata.var.index.tolist())
    X_df, to_fill_columns, var = main_gene_selection(X_df, gene_list)
    adata_uni = sc.AnnData(X_df)  
    adata_uni = BasicFilter(adata_uni,qc_min_genes=200,qc_min_cells=0) # filter cell and gene by lower limit
    
    save_adata_h5ad(adata_uni,save_path)

if __name__=='__main__':
    main()