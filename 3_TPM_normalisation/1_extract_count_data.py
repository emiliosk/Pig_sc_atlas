import scanpy as sc

import numpy as np
import pandas as pd
import scipy.sparse as sp


def load_file_mgi(acc, directory = './', use_filtered = True, min_cells = 0, min_genes = 0):
    if use_filtered == True:
        location= 'Solo.out/Gene/filtered/'
    elif use_filtered == False:
        location = 'Solo.out/Gene/raw/'
    
    filename = directory + acc + '/' + location
    print("loading .... " + filename)
    adata = sc.read_10x_mtx(filename, var_names='gene_ids', cache=False)
    print(adata.shape)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
  #  adata_temp.obs['type'] = tissue
    adata.obs['sample'] = acc
    if min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        print(adata.shape)
    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(adata.shape)
    return(adata)


def integrate_files_mgi(tissue, ls_acc, min_cells = 0, min_genes = 0, results_path = './', directory = './', plot = True, use_filtered = True):
    adatas = []
    
    for acc in ls_acc:
        adata_temp = load_file_mgi(acc = acc, directory = directory, use_filtered = use_filtered, min_cells = min_cells, min_genes = min_genes)
        #adata_temp = individual_qc(adata_temp,acc = acc, results_path = results_path,  plot = plot)
        adatas.append(adata_temp.copy())
    
    adata = adatas[0].concatenate(adatas[1:])
    adata.obs["ori_barcode"]=[i.split(sep="-")[0] for i in adata.obs_names]
    adata.obs['type'] = tissue
    
    del(adatas, adata_temp)
    print("***" + tissue + "***")
    print("Final shape: ")
    print(adata.shape)
    print("value counts per sample")
    print(adata.obs['sample'].value_counts())
    print("*********************************")
    return adata
def load_file_10x(acc, directory = './', use_filtered = True, use_soupx = False,  min_cells = 0, min_genes = 0):
    if use_filtered == True:
        location= 'outs/filtered_feature_bc_matrix'
    elif use_filtered == False:
        location = 'outs/raw_feature_bc_matrix'
    if use_soupx == True:
        location= 'outs/soupx/'

    filename = directory + acc + '/' + location
    print("loading .... " + filename)
    adata = sc.read_10x_mtx(filename, var_names='gene_ids', cache=False)
    print(adata.shape)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
  #  adata_temp.obs['type'] = tissue
    adata.obs['sample'] = acc

    if min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        print(adata.shape)
    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(adata.shape)
    adata.obs['type'] = acc
    adata.obs["ori_barcode"]=[i.split(sep="-")[0] for i in adata.obs_names]
    return(adata)




single_nuclei = {
    #MGI samples
   'Area Postrema': ["Area_Postrema_1","Area_Postrema_2","Area_Postrema_3"],
                'Cerebellum': ["Cerebellum_1","Cerebellum_2","Cerebellum_3","Cerebellum_4","Cerebellum_5","Cerebellum_6","Cerebellum_7","Cerebellum_8"],
                'Heart': ["Heart_1","Heart_2","Heart_3","Heart_4","Heart_5"] ,
                'Kidney': ["Kidney_1","Kidney_2","Kidney_3","Kidney_4"],
                 'Liver': ["Liver_1","Liver_2","Liver_3","Liver_4"],
                 'OVoLT': ["OVoLT_1","OVoLT_2","OVoLT_3","OVoLT_4"],
                 'Retina': ["Retina_1","Retina_2"],
                 'Spleen': ["Spleen_1","Spleen_2","Spleen_3"],
                 'Subfornical organ': ["SubfomicalOrgan_1","SubfomicalOrgan_2","SubfomicalOrgan_3"],
                 'Lung':['Lung_1', 'Lung_2', 'Lung_3', 'Lung_4', 'Lung_5', 'Lung_6', 'Lung_7', 'Lung_8'],
    #10x vs samples
   'Frontal lobe':['frontal'],
    'Hypothalamus':['hypothalamus1', 'hypothalamus2'],
    'Occipital lobe':['occipital'],
    'Parietal lobe':['parietal'],
    'Temporal lobe':['temporal'],
    
                }


method = 'snRNA-seq'
directory_mgi = '../../single_nuclei/'
directory_10x = '../../single_nuclei_10x_star/'



for tissue, ls_acc in single_nuclei.items():
    try: 
        adata = integrate_files_mgi(tissue = tissue, 
                                    ls_acc = ls_acc, 
                                    min_cells = 0, min_genes = 0, 
                                    results_path = '.', 
                                    directory =  directory_mgi, 
                                    plot = False, use_filtered = True)
    except:
        adata = integrate_files_mgi(tissue = tissue, 
                                    ls_acc = ls_acc, 
                                    min_cells = 0, min_genes = 0, 
                                    results_path = '.', 
                                    directory =  directory_10x, 
                                    plot = False, use_filtered = True)
    #barcodes = adata.obs_names.tolist()
    obs_df = pd.read_csv('_'.join([tissue, method, 'cluster_labels.csv']), index_col = 0)

    index_to_keep = list(obs_df.index)

    adata_filtered = adata[adata.obs.index.isin(index_to_keep),:]

    counts=adata_filtered.to_df()
    counts=counts.T
    if len(index_to_keep) ==counts.shape[1]:
            counts.astype(int).to_csv('_'.join([tissue, method, 'counts.tsv']),sep="\t")
            print ('done!')
    else:
            print ('ERROR ' + tissue + ' is missing genes!' )
            
            
tissues = ["brain", "Adipose-S", "Adipose-V", "Intestine", "liver", "Lung", "PBMC", "retina", "spleen"]

method = 'scRNA-seq'

for tissue in tissues:
    adata = load_file_10x(acc = tissue, 
                          directory = "/proj/rnaatlas/nobackup/private/EmilioTemp/pig_sc_109/single_cell/", 
                          use_filtered = True,  min_cells = 0, min_genes = 0)
    #barcodes = adata.obs_names.tolist()
    obs_df = pd.read_csv('_'.join([tissue, method, 'cluster_labels.csv']), index_col = 0)

    index_to_keep = list(obs_df.index)

    adata_filtered = adata[adata.obs.index.isin(index_to_keep),:]

    counts=adata_filtered.to_df()
    counts=counts.T
    if len(index_to_keep) ==counts.shape[1]:
            counts.astype(int).to_csv('_'.join([tissue, method, 'counts.tsv']),sep="\t")
            print ('done!')
    else:
            print ('ERROR ' + tissue + ' is missing genes!' )