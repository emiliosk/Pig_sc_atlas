import scanpy as sc

import numpy as np
import pandas as pd
import scipy.sparse as sp


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
    "Lung_snRNA-seq_cluster_labels.csv",
   "Occipital lobe_snRNA-seq_cluster_labels.csv",
   "OVoLT_snRNA-seq_cluster_labels.csv",
   "Parietal lobe_snRNA-seq_cluster_labels.csv",
   "PBMC_scRNA-seq_cluster_labels.csv", 
    "retina_scRNA-seq_cluster_labels.csv",
    "Retina_snRNA-seq_cluster_labels.csv",
    "spleen_scRNA-seq_cluster_labels.csv",
    "Spleen_snRNA-seq_cluster_labels.csv",
    "Subfornical organ_snRNA-seq_cluster_labels.csv",
    "Temporal lobe_snRNA-seq_cluster_labels.csv",
]




###### based on protein coding genes of ensemble plus 23 genes defined by the HPA as protein coding genes:

gtf = pd.read_csv('/proj/rnaatlas/nobackup/private/EmilioTemp/ref/Sus_scrofa.Sscrofa11.1.109.gtf', sep='\t', comment='#', header=None, usecols=[2, 8], names=['type', 'attributes'])

gtf['gene_id'] = gtf['attributes'].str.extract('gene_id "([^"]+)"', expand=False)
gtf['gene_type'] = gtf['attributes'].str.extract('gene_biotype "([^"]+)"', expand=False)
gtf['gene_name'] = gtf['attributes'].str.extract('gene_name "([^"]+)"', expand=False)

#rename with gene_id if NaN
#gtf.loc[gtf['gene_name'].isna(), 'gene_name'] = gtf['gene_id']
                #check for correct change 
                #gtf[gtf.gene_name.isna()]
                #gtf[gtf.gene_name.str.startswith('ENSS')]

#load transcriptlength info
from pybiomart import Server
server = Server(host='http://www.ensembl.org')
dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                     .datasets['sscrofa_gene_ensembl'])

gene_length = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','transcript_length'])
gene_length.index = gene_length['Gene stable ID']


#add transcript length info
gtf['transcript_length'] = gtf.gene_id.map(
    gene_length['Transcript length (including UTRs and CDS)'].to_dict()
)
#nicer version of the gtf

gene_meta = gtf[gtf['gene_type'].isin([ 
    'protein_coding', 'IG_V_gene', 'TR_V_gene', 'TR_J_gene', 'IG_C_gene'
    ])][['gene_id', 'gene_type', 'gene_name', 'transcript_length']]
gene_meta = gene_meta.drop_duplicates()

gene_meta.set_index('gene_id', inplace=True)


#list of protein coding genes

protein_coding_genes = gtf.loc[gtf['gene_type'].isin([ 'protein_coding', 'IG_V_gene', 'TR_V_gene', 'TR_J_gene', 'IG_C_gene']), 'gene_id']
protein_coding_genes =  list(protein_coding_genes.unique())




drop_mix = False
exclude_low_confidence = True
exclusivly_hq = False


for metadata_name in cluster_labels:
    print(metadata_name + ' reading ....')
    metadata = pd.read_csv(metadata_name, index_col = 0)
    raw_counts = pd.read_csv('_'.join(metadata_name.split('_')[:2])+'_counts.tsv', index_col = 0, sep = '\t')

    #generate celltype cluster laballing
    metadata['celltype_cluster'] = metadata['manual_celltype_annotation'] + '_' + metadata['leiden'].astype(str)

    #map weights (counts) to celltype_cluster
    cell_counts = metadata['celltype_cluster'].value_counts().to_dict()
    metadata["cell_count"] = metadata.celltype_cluster.map(cell_counts)
    metadata["sample_cell_count"] = metadata['celltype_cluster'].value_counts().sum()

    #filter MIxtures, NAs and low quality labellings
    pseudobulk = pd.merge(raw_counts.T, metadata[['celltype_cluster','mixture','confidence' ]], left_index=True, right_index=True)

    pseudobulk = pseudobulk.dropna()

    if drop_mix == True:
        pseudobulk = pseudobulk[pseudobulk.mixture != True]

    if exclude_low_confidence == True:
        pseudobulk = pseudobulk[pseudobulk.confidence != 'low']

    if exclusivly_hq == True:
        pseudobulk = pseudobulk[pseudobulk.confidence == 'high']

    #generate pseudulk through sum
    pseudobulk = pseudobulk.groupby('celltype_cluster').sum()


    #generate protein coding filtered pseudobulk
    pseudobulk_filtered = pseudobulk[pseudobulk.columns[pseudobulk.columns.isin(protein_coding_genes)]]
    pseudobulk_filtered = pseudobulk_filtered.T

    #calculate RPK: divide gene count for each transcript by it's length
    RPK = pseudobulk_filtered.divide(gene_meta['transcript_length'].to_dict(), axis=0)

    #calculate TPM scaling factor, sum of RPK for every celltype/million
    TPM_scaling = RPK.sum(axis=0)/1000000

    #get TPM by divideing through the saling factor. pTPM; since only including proteins
    pTPM = RPK.divide(TPM_scaling.to_dict(), axis = 1)
    
    
    filename = '_'.join(metadata_name.split('_')[:2])+'_pTPM.tsv'
    
    if drop_mix == False:
        filename = '_'.join(metadata_name.split('_')[:2])+'_wMix_pTPM.tsv'
    
    pTPM.to_csv(filename, sep="\t")
    
    
    
    
    

  