import scanpy as sc
import scrublet as scr
from pybiomart import Server
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns
import pandas as pd
import os
import pickle
import warnings
from matplotlib.pyplot import rc_context

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

sc.set_figure_params(scanpy=True, dpi=150)


def load_file_mgi(acc, directory = './', use_filtered = True, min_cells = 0, min_genes = 0):
    if use_filtered == True:
        location= 'Solo.out/Gene/filtered/'
    elif use_filtered == False:
        location = 'Solo.out/Gene/raw/'

    filename = directory + acc + '/' + location
    print("loading .... " + filename)
    adata = sc.read_10x_mtx(filename, var_names='gene_symbols', cache=False)
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

def individual_qc(adata, acc, results_path = './',  plot = True):
    server = Server(host='http://www.ensembl.org')
    dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                     .datasets['sscrofa_gene_ensembl'])
    mito_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                  filters={'chromosome_name': ['MT']})

    ## ribosomal GO:0022626
    ribo_genes = dataset.query(attributes=['ensembl_gene_id','external_gene_name'],
                 filters={'link_go_closure': ['GO:0022626']})
    ## Extracellular Matrix: GO:0031012
    exc_genes = dataset.query(attributes=['ensembl_gene_id','external_gene_name'],
                 filters={'link_go_closure': ['GO:0031012']})
    ## Cytosol: GO:0005829
    cyt_genes = dataset.query(attributes=['ensembl_gene_id','external_gene_name'],
                 filters={'link_go_closure': ['GO:0005829']})
    ## Membrane: GO:0016020
    mem_genes = dataset.query(attributes=['ensembl_gene_id','external_gene_name'],
                 filters={'link_go_closure': ['GO:0016020']})
    ## PcG: GO:0031519
    pcg_genes = dataset.query(attributes=['ensembl_gene_id','external_gene_name'],
                 filters={'link_go_closure': ['GO:0031519']})
    ## Transcription factors: GO:0006355
            ##actually: regulation of DNA-templated transcription
    tf_genes = dataset.query(attributes=['ensembl_gene_id','external_gene_name'],
                 filters={'link_go_closure': ['GO:0006355']})
    ## hemoglobin complex
    hb_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                 filters={'link_go_closure': ['GO:0005833']})
    ## platelet alpha granule
    pl_genes = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'],
                 filters={'link_go_closure': ['GO:0031091']})

    adata.var['mt'] = adata.var['gene_ids'].isin(mito_genes['Gene stable ID'])
    adata.var['ribo'] = adata.var['gene_ids'].isin(ribo_genes['Gene stable ID'])
    adata.var['exc'] = adata.var['gene_ids'].isin(exc_genes['Gene stable ID'])
    adata.var['cyt'] = adata.var['gene_ids'].isin(cyt_genes['Gene stable ID'])
    adata.var['mem'] = adata.var['gene_ids'].isin(mem_genes['Gene stable ID'])
    adata.var['pcg'] = adata.var['gene_ids'].isin(pcg_genes['Gene stable ID'])
    adata.var['tf'] = adata.var['gene_ids'].isin(tf_genes['Gene stable ID'])
    adata.var['hb'] = adata.var['gene_ids'].isin(hb_genes['Gene stable ID'])
    adata.var['pl'] = adata.var['gene_ids'].isin(pl_genes['Gene stable ID'])

    ## Cell Cycle genes as defined by cc.genes in Seurat_4.3.0
    hs_s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
        "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
        "RAD51AP1","GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
        "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2","USP1", "CLSPN",
        "POLA1", "CHAF1B", "BRIP1", "E2F8"]
    hs_g2m_genes = [ "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
        "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2",
        "AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20",
        "TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR",
        "AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA" ]

    hs_homologs = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name',
        'hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name',
        'hsapiens_homolog_orthology_type'])

    s_genes = list(hs_homologs[hs_homologs['Human gene name'].isin(hs_s_genes)]['Gene stable ID'])
    g2m_genes = list(hs_homologs[hs_homologs['Human gene name'].isin(hs_g2m_genes)]['Gene stable ID'])

    adata.var['s_genes'] = adata.var_names.isin(adata.var[adata.var['gene_ids'].isin(s_genes)].index)
    adata.var['g2m_genes'] = adata.var_names.isin(adata.var[adata.var['gene_ids'].isin(g2m_genes)].index)

    #QC
    if plot == True:
        pdf=PdfPages(results_path+ "/" + acc + "_figures_qc.pdf")

    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=['mt', 'ribo', 'exc', 'cyt', 'mem', 'pcg', 'tf', 'hb', 'pl'],
                               percent_top=None, log1p=False, inplace=True)
    if plot == True:
        sns.displot(adata.obs["total_counts"], bins=100, kde=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.violin(adata,
                    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt',
                        'pct_counts_ribo', 'pct_counts_hb'],
                    jitter=0.4,
                    groupby = 'sample',
                    rotation= 45,
                    show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.violin(adata,
                    ['pct_counts_pl',
                        'pct_counts_exc', 'pct_counts_mem', 'pct_counts_pcg',
                        'pct_counts_tf'],
                    jitter=0.4,
                    groupby = 'sample',
                    rotation= 45,
                    show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.highest_expr_genes(adata,n_top=20, show=False)
        pdf.savefig(bbox_inches='tight')

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05, random_state=4)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
    if plot == True:
        scrub.plot_histogram()
        pdf.savefig(bbox_inches='tight')

    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

    if plot == True:
        sc.pl.violin(adata, 'n_genes_by_counts',jitter=0.4, groupby = 'doublet_info', rotation=45, show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="predicted_doublets", show=False)
        pdf.savefig(bbox_inches='tight')


    print("adata shape: ")
    print(adata.shape)
    print("value counts per sample")
    print(adata.obs['sample'].value_counts())

    print("after doublet cutoff")
    adata = adata[adata.obs['doublet_info'] == 'False',:]
    print("adata shape: ")
    print(adata.shape)
    print("value counts per sample")
    print(adata.obs['sample'].value_counts())

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
    pdf.savefig(bbox_inches='tight')


    pdf.close()
    return adata


def integrate_files_mgi(tissue, ls_acc, min_cells = 0, min_genes = 0, results_path = './', directory = './', plot = True, use_filtered = True):
    adatas = []

    for acc in ls_acc:
        adata_temp = load_file_mgi(acc = acc, directory = directory, use_filtered = use_filtered, min_cells = min_cells, min_genes = min_genes)
        adata_temp = individual_qc(adata_temp,acc = acc, results_path = results_path,  plot = plot)
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



def integrated_qc(adata, tissue, results_path, plot = True):
    ### Plots initial QC plots on adata
    # Load ensembl data on genes

    ##################changed log1p to true!!!! + èrcent_tp`
    sc.pp.calculate_qc_metrics(adata,
                               percent_top=[20], log1p=True, inplace=True)



    if plot == True:

        pdf=PdfPages(results_path+"/" + tissue + "_figures_qc.pdf")

        sns.displot(adata.obs["total_counts"], bins=100, kde=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.violin(adata,
                    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt',
                        'pct_counts_ribo', 'pct_counts_hb'],
                    jitter=0.4,
                    groupby = 'sample',
                    rotation= 45,
                    show=False)
        pdf.savefig(bbox_inches='tight')
        sc.pl.violin(adata,
                    ['pct_counts_pl',
                        'pct_counts_exc', 'pct_counts_mem', 'pct_counts_pcg',
                        'pct_counts_tf'],
                    jitter=0.4,
                    groupby = 'sample',
                    rotation= 45,
                    show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "pct_counts_in_top_20_genes", "log1p_n_genes_by_counts",color="pct_counts_mt",  show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "pct_counts_in_top_20_genes", "log1p_total_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')


        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="doublet_info", show=False)
        pdf.savefig(bbox_inches='tight')


        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')

        pdf.close()


    print("***" + tissue + "***")
    print("Final shape, after mito and doublet filter: ")
    print(adata.shape)
    print("value counts per sample")
    print(adata.obs['sample'].value_counts())
    print("*********************************")

    return adata


def is_outlier(adata, metric: str, nmads: int):
    ### taken from: https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * M.mad()) | (
        np.median(M) + nmads * M.mad() < M
    )
    return outlier



def preprocessing(adata, tissue,  res_path, filter_genes = True, filter_cells = True, plot = True):
    print("before preprocessing:")
    print("adata shape:")
    print(adata.shape)

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    print("outlier value counts:")
    print(adata.obs.outlier.value_counts())

    if plot == True:
        pdf=PdfPages(res_path+ "/"  + tissue + "_figures_pre_process.pdf")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')

    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3)
    adata.obs.mt_outlier.value_counts()

    adata = adata[adata.obs['mt_outlier'] == False,:]
    if plot == True:
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')



    if plot == True:


        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="outlier", show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "log1p_total_counts", "log1p_n_genes_by_counts", color="outlier", show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "pct_counts_in_top_20_genes", "log1p_n_genes_by_counts",color="outlier",  show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "pct_counts_in_top_20_genes", "log1p_total_counts", color="outlier", show=False)
        pdf.savefig(bbox_inches='tight')

    adata = adata[adata.obs['outlier'] == False,:]

    if filter_genes == True:
       # min_cells = min(adata.n_obs*0.1, 1000)
        min_cells = 20
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print("after gene cutuff:")
        print("adata shape:")
        print(adata.shape)
    if filter_cells == True:
        min_genes = 200
        #max_genes = 5000
        sc.pp.filter_cells(adata, min_genes=min_genes)
       # sc.pp.filter_cells(adata, max_genes=max_genes)
        print("after cell cutoff:")
        print("adata shape:")
        print(adata.shape)

    if plot == True:
        sc.pl.violin(adata,
                    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt',
                        'pct_counts_ribo', 'pct_counts_hb'],
                    jitter=0.4,
                    groupby = 'sample',
                    rotation= 45,
                    show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.violin(adata,
                    ['pct_counts_pl',
                        'pct_counts_exc', 'pct_counts_mem', 'pct_counts_pcg',
                        'pct_counts_tf'],
                    jitter=0.4,
                    groupby = 'sample',
                    rotation= 45,
                    show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "pct_counts_in_top_20_genes", "log1p_n_genes_by_counts",color="pct_counts_ribo",  show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.scatter(adata, "pct_counts_in_top_20_genes", "log1p_total_counts", color="pct_counts_ribo", show=False)
        pdf.savefig(bbox_inches='tight')

    sc.pp.normalize_total(adata, target_sum=1e4)

    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata

    if plot == True:
        sc.pl.highly_variable_genes(adata, show=False)
        pdf.savefig(bbox_inches='tight')
        sc.pl.highest_expr_genes(adata,n_top=20, show=False)
        pdf.savefig(bbox_inches='tight')


   # sc.pp.scale(adata, max_value=10)
    sc.pp.regress_out(adata, [ 'pct_counts_mt'])

    sc.tl.pca(adata, svd_solver='arpack')

    if plot == True:
        sc.pl.pca(adata, color='total_counts', title=tissue, show=False)
        pdf.savefig(bbox_inches='tight')
        sc.pl.pca_variance_ratio(adata, log=True, show=False)
        pdf.savefig(bbox_inches='tight')


    sc.pp.neighbors(adata)
    sc.tl.umap(adata, random_state=4)
    sc.tl.louvain(adata)
    sc.tl.leiden(adata)

    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')


    results = adata.uns['rank_genes_groups']
    out = np.array([[0,0,0,0,0]])
    for group in results['names'].dtype.names:
        out = np.vstack((out, np.vstack((results['names'][group],
                                         results['scores'][group],
                                         results['pvals_adj'][group],
                                         results['logfoldchanges'][group],
                                         np.array([group] * len(results['names'][group])).astype('object'))).T))
    markers = pd.DataFrame(out[1:], columns = ['Gene', 'scores', 'pval_adj', 'lfc', 'cluster'])
    #adata.uns['markers'] = markers.to_dict() #save marker df to uns
    markers[(markers.pval_adj < 0.05) & (markers.lfc > 1.5)].to_csv(results_path+ "/"  + tissue + "_diffex_genes_per_cluster.csv")

    if plot == True:
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.rank_genes_groups(adata, file="top20 Diff_Genes",n_genes=20, show=False)
        pdf.savefig(bbox_inches='tight')

        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="leiden", show=False)
        pdf.savefig(bbox_inches='tight')


    sc.tl.dendrogram(adata, groupby = "leiden")
    if plot == True:
        sc.pl.dendrogram(adata, groupby = "leiden", show=False)
        pdf.savefig(bbox_inches='tight')

        adata_obs_names=['louvain', 'leiden', "sample",'n_genes_by_counts','total_counts', 'pct_counts_mt',
                         'pct_counts_ribo', 'pct_counts_hb', 'pct_counts_pl',
                         'pct_counts_exc', 'pct_counts_mem',
                         'pct_counts_pcg', 'pct_counts_tf', "S_score", "G2M_score", "phase"]
        for i in adata_obs_names:
            if i in adata.obs.columns:
                x=sc.pl.umap(adata,legend_fontsize=9, title=tissue + ' UMAP based on '+i,legend_fontweight='bold',color=[i],show=False)
                pdf.savefig(bbox_inches='tight')


        pdf.close()

    return adata



def plot_markers(adata, tissue, prep_method, results_dir):
    pdf=PdfPages(results_path+ "/"  + tissue + "_markers.pdf")
    #open
    with open('cell_markers_dict_dk.pkl', 'rb') as f:
        cell_markers_dict_dk = pickle.load(f)
    sample = tissue.title() + '_' + prep_method
    cell_markers_dict_dk = cell_markers_dict_dk[sample]
    cell_markers_dict_dk = get_markers_in_adata(adata, cell_markers_dict_dk)

    sc.pl.dotplot(
        adata,
        title = 'Cell markers on tissue used by Wang et al 2022 for this tissue',
        groupby="leiden",
        var_names=cell_markers_dict_dk,
        standard_scale="var", dendrogram = True, show=False # standamarker_genes_in_datard scale: normalize each gene to range from 0 to 1
    )
    pdf.savefig(bbox_inches='tight')

    dp = sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=cell_markers_dict_dk,
    standard_scale="var", dendrogram = True , return_fig=True# standamarker_genes_in_datard scale: normalize each gene to range from 0 to 1
    ).add_totals().make_figure()
    pdf.savefig(bbox_inches='tight')

    sc.pl.matrixplot(adata, cell_markers_dict_dk, groupby='leiden',standard_scale = 'var',
                     #cmap='Blues',
                     cmap=sns.color_palette("blend:#F1ED6F,#9C2964,#000000", as_cmap=True),
                     colorbar_title = 'Mean expression\nin group\nscaled to max var expression',
                     dendrogram=True, return_fig=True
    ).add_totals().make_figure()#show=False
    pdf.savefig(bbox_inches='tight')


    sc.pl.correlation_matrix(adata, groupby = 'leiden', figsize=(8,6.5),
                         cmap = 'magma_r', show = False)
    pdf.savefig(bbox_inches='tight')
    sc.pl.correlation_matrix(adata, groupby = 'leiden',vmin = 0, figsize=(8,6.5),
                         cmap = 'magma_r',show = False)
    pdf.savefig(bbox_inches='tight')

    sc.tl.dendrogram(adata, groupby='leiden')
    sc.pl.dendrogram(adata, 'leiden', show = False)
    pdf.savefig(bbox_inches='tight')

    all_markers_df = pd.read_csv('celltype_markers_pooled.csv', index_col = 0)
    markers_dk = all_markers_df[all_markers_df['author'] == 'dk'][['celltype', "marker"]].groupby('celltype')['marker'].apply(list).to_dict()
    markers_dk = get_markers_in_adata(adata, markers_dk)
    markers_hpa = all_markers_df[all_markers_df['author'] == 'hpa'][['celltype', "marker"]].groupby('celltype')['marker'].apply(list).to_dict()
    markers_hpa = get_markers_in_adata(adata, markers_hpa)

    sc.pl.dotplot(
        adata,
        groupby="leiden",
        title = 'All Wang et al cell markers',
        var_names=markers_dk,
        standard_scale="var",
        dendrogram = True,
        return_fig=True
    ).add_totals().make_figure()#.show()
    pdf.savefig(bbox_inches='tight')

    sc.pl.matrixplot(
        adata,
        groupby="leiden",
        standard_scale = 'var',
        var_names=markers_dk, #standamarker_genes_in_datard scale: normalize each gene to range from 0 to 1
        dendrogram = True,
        cmap=sns.color_palette("blend:#F1ED6F,#9C2964,#000000", as_cmap=True),
        return_fig=True
    ).add_totals().make_figure()#.show()
    pdf.savefig(bbox_inches='tight')



    dp = sc.pl.dotplot(
        adata,
        groupby="leiden",
        title = 'All HPA cell markers',
        var_names=markers_hpa,
        standard_scale="var", dendrogram = True , return_fig=True, show=False# standamarker_genes_in_datard scale: normalize each gene to range from 0 to 1
    ).add_totals().make_figure()#.show()
    pdf.savefig(bbox_inches='tight')

    sc.pl.matrixplot(
        adata,
        groupby="leiden",
        standard_scale = 'var',
        var_names=markers_hpa, #standamarker_genes_in_datard scale: normalize each gene to range from 0 to 1
        dendrogram = True,
        cmap=sns.color_palette("blend:#F1ED6F,#9C2964,#000000", as_cmap=True),
        return_fig=True
    ).add_totals().make_figure()#.show()
    pdf.savefig(bbox_inches='tight')

    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
    pdf.savefig(bbox_inches='tight')

    sc.pl.rank_genes_groups(adata, file="top20 Diff_Genes",n_genes=20, show=False)
    pdf.savefig(bbox_inches='tight')

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="leiden", show=False)
    pdf.savefig(bbox_inches='tight')

    sc.pl.rank_genes_groups_matrixplot(adata,
                                       n_genes=5,
                                       groupby="leiden",
                                       standard_scale = 'var',
                                       cmap=sns.color_palette("blend:#F1ED6F,#9C2964,#000000", as_cmap=True),
                                       return_fig=True ).add_totals().make_figure()#.show()
    pdf.savefig(bbox_inches='tight')

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, values_to_plot='logfoldchanges', min_logfoldchange=3,vmax=7, vmin=-7, cmap='bwr', show=False)
    pdf.savefig(bbox_inches='tight')

    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, values_to_plot='logfoldchanges', min_logfoldchange=3,vmax=7, vmin=-7, cmap='bwr', show=False)
    pdf.savefig(bbox_inches='tight')



    pdf.close()
    return adata


def get_markers_in_adata(adata, marker_genes):
    marker_genes_in_data = dict()
    for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found
    marker_genes_in_data = {k: v for k, v in marker_genes_in_data.items() if v}
    return marker_genes_in_data



single_nuclei = {'Area Postrema': ["Area_Postrema_1","Area_Postrema_2","Area_Postrema_3"],
                 'Cerebellum': ["Cerebellum_1","Cerebellum_2","Cerebellum_3","Cerebellum_4","Cerebellum_5","Cerebellum_6","Cerebellum_7","Cerebellum_8"],
                 'Heart': ["Heart_1","Heart_2","Heart_3","Heart_4","Heart_5"] ,
                 'Kidney': ["Kidney_1","Kidney_2","Kidney_3","Kidney_4"],
                 'Liver': ["Liver_1","Liver_2","Liver_3","Liver_4"],
                 'OVoLT': ["OVoLT_1","OVoLT_2","OVoLT_3","OVoLT_4"],
                 'Retina': ["Retina_1","Retina_2"],
                 'Spleen': ["Spleen_1","Spleen_2","Spleen_3"],
                 'Subfornical Organ': ["SubfomicalOrgan_1","SubfomicalOrgan_2","SubfomicalOrgan_3"]
                }

for tissue, ls_acc in single_nuclei.items():
    try:
        os.mkdir('../results_qc/single_nuclei_MGI/'+ tissue)
        results_path = '../results_qc/single_nuclei_MGI/'+ tissue
        adata = integrate_files_mgi(tissue = tissue,
                            ls_acc = ls_acc,
                            min_cells = 0, min_genes = 0,
                            results_path = results_path,
                            directory =  "/proj/rnaatlas/nobackup/private/EmilioTemp/pig_sc_109/single_nuclei/", 
                            plot = True, use_filtered = True)
        adata =  integrated_qc(adata, tissue = tissue, results_path = results_path)
        adata = preprocessing(adata, tissue = tissue, res_path = results_path, filter_genes = True, filter_cells = True, plot = True)
        adata.write_h5ad(results_path  + '/'+ tissue  + '_preprocessed.h5ad')
        adata = plot_markers(adata, tissue = tissue, prep_method = "snRNA-seq", results_dir = results_path )
    except:
        pass
