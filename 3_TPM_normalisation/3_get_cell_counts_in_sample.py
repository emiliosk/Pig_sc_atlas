


import numpy as np
import pandas as pd


cluster_labels = [
   "Adipose-S_scRNA-seq_cluster_labels.csv",
   "Adipose-V_scRNA-seq_cluster_labels.csv",
   "Area Postrema_snRNA-seq_cluster_labels.csv",
   "brain_scRNA-seq_cluster_labels.csv",
   "Cerebellum_snRNA-seq_cluster_labels.csv",
   "Frontal lobe_snRNA-seq_cluster_labels.csv",
   "Heart_snRNA-seq_cluster_labels.csv",
   "Hypothalamus_snRNA-seq_cluster_labels.csv",
   "Intestine_scRNA-seq_cluster_labels.csv",
   "Kidney_snRNA-seq_cluster_labels.csv",
   "liver_scRNA-seq_cluster_labels.csv",
   "Liver_snRNA-seq_cluster_labels.csv",
   "Lung_scRNA-seq_cluster_labels.csv",
"Lung_snRNA-seq_cluster_labels.csv",#error
   "Occipital lobe_snRNA-seq_cluster_labels.csv",
   "OVoLT_snRNA-seq_cluster_labels.csv",
   "Parietal lobe_snRNA-seq_cluster_labels.csv",
   "PBMC_scRNA-seq_cluster_labels.csv", #['mixture', 'confidence'] not in index"
    "retina_scRNA-seq_cluster_labels.csv",
    "Retina_snRNA-seq_cluster_labels.csv",
    "spleen_scRNA-seq_cluster_labels.csv",
    "Spleen_snRNA-seq_cluster_labels.csv",
    "Subfornical organ_snRNA-seq_cluster_labels.csv",
    "Temporal lobe_snRNA-seq_cluster_labels.csv",
]
drop_mix = False
exclude_low_confidence = True
exclusivly_hq = False

#input_dir = 'nTPM_output_unweighted/'
output_dir = 'metadata/wMix/'

for metadata_name in cluster_labels:
    print(metadata_name + ' reading ....')
  #  pTPM = pd.read_csv(metadata_name, index_col = 0)
    metadata = pd.read_csv(metadata_name, index_col = 0)
    #nTPM = pd.read_csv(input_dir +'_'.join(metadata_name.split('_')[:2])+'_nTPM.tsv', index_col = 0, sep = '\t')
    
    metadata['celltype_cluster'] = metadata['manual_celltype_annotation'] + '_' + metadata['leiden'].astype(str)
    
    metadata = metadata.dropna()
    
    if drop_mix == True:
        metadata = metadata[metadata.mixture != True]

    if exclude_low_confidence == True:
        metadata = metadata[metadata.confidence != 'low']

    if exclusivly_hq == True:
        metadata = metadata[metadata.confidence == 'high']

    #map weights (counts) to celltype_cluster
    cell_counts = metadata['celltype_cluster'].value_counts()
    metadata["cell_count"] = metadata.celltype_cluster.map(cell_counts)
    total_counts = cell_counts.sum()
    metadata["total_counts"] = total_counts
    
    #nTPM_temp1 = nTPM.multiply(cell_counts.to_dict)
    #nTPM_temp1 = nTPM.multiply(cell_counts.to_dict(), axis = 1)
    metadata.to_csv(output_dir +'_'.join(metadata_name.split('_')[:2])+'.csv')
    

    
    