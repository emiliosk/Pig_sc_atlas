
library(tidyverse)
library(NOISeq)
library(patchwork)
library(ggplot2)




calc_tmm_normfactors <- 
  function (object, method = c("TMM", "quantile"), refColumn = NULL, 
            logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
            Acutoff = -1e+10, quantile = 0.75) {
    method <- match.arg(method)
    if (is.matrix(object)) {
      if (is.null(refColumn)) 
        refColumn <- 1
      data <- object
      libsize <- colSums(data)
    } else {
      stop("calcNormFactors() only operates on 'matrix' objects")
    }
    
    if(refColumn == "median") {
      ref <- 
        apply(data, MARGIN = 1, median)
    } else {
      ref <- data[, refColumn]
    }
    
    f <- switch(method, TMM = apply(data, 2, NOISeq:::.calcFactorWeighted, 
                                    ref = ref, logratioTrim = logratioTrim, 
                                    sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff), 
                quantile = NOISeq:::.calcFactorQuantile(data, libsize, q = quantile))
    f <- f/exp(mean(log(f)))
    return(f)
  }

input_dir  <- 'pTPM_output/wMix'
output_dir <- 'nTPM_output_unweighted'
metadata <- read_csv('pooled_data/celltype_name_table.csv')

pTPM_files = list.files(input_dir, pattern = 'pTPM.tsv')





####TMM Normalize the pTPM to nTPM

pal_tbl <- metadata %>% select(celltype_cluster, color) %>% distinct()
pal <- pal_tbl %>% pull(color)
pal <- setNames(pal, pal_tbl %>% pull(celltype_cluster))


all_tmm_factors <- tibble()

for (pTPM_filename in pTPM_files){
  pTPM <- read_tsv(paste0(input_dir, '/',pTPM_filename)) %>% rename('genes' = '...1')
  
  tmm_factors <- pTPM %>% column_to_rownames("genes") %>% as.matrix() %>% 
    calc_tmm_normfactors(method = "TMM", 
                         
                         # Use median column as reference distribution:
                         refColumn = "median",
                         
                         # Trim parameters:
                         logratioTrim = 0.3, 
                         sumTrim = 0.3, 
                         
                         # Weighting should be done only if count data:
                         doWeighting = F) %>% 
    enframe("celltype_cluster", "tmm_factor")
  
  TMM <- t(t(column_to_rownames(pTPM, "genes")) / tmm_factors$tmm_factor) 
  names_sample <- colnames(TMM)
  names_genes <- rownames(TMM)
  
  TMM <- TMM %>%
    as_tibble(rownames = "genes")
  
  
  name <- list(c(strsplit(pTPM_filename, '_')[[1]][1:2], 'nTMM'))
  tissue <- strsplit(pTPM_filename, '_')[[1]][1]
  method <-  strsplit(pTPM_filename, '_')[[1]][2]
  filename = paste(tissue, method,'nTPM.tsv', sep = '_')
  
  write_delim(TMM, paste0(output_dir, '/', filename), delim = '\t')
  
  plot_data <-pTPM %>% gather(cluster,pTPM, -1) %>% left_join(TMM %>% gather(cluster,nTPM, -1), by = c('genes', 'cluster')) %>% 
    gather(expression_type, expression, -1, -2) %>% mutate_if(is.numeric, function(x) {
      log10(x + 1)
    }) %>% mutate(expression_type = factor(expression_type, c('pTPM', 'nTPM')))
  
  ggplot(data = plot_data, aes(x = expression, y =cluster, fill = cluster)) +
    geom_boxplot(draw_quantiles = 0.5, outlier.size = 0.5, outlier.alpha = 0.3)+
    facet_wrap(~expression_type)+
    scale_fill_manual(values = pal) +
    theme(legend.position = "none", 
          axis.title = element_blank())
  ggsave(paste0('cluster_plots/nTPM_pTPM_expression/', tissue, '_', method, '.pdf'))
  
  tmm_factors['sample']  = paste(tissue, method, sep = '_')
  
  all_tmm_factors <- bind_rows(all_tmm_factors, tmm_factors)
  
}


all_tmm_factors %>% 
  ggplot(aes(1/tmm_factor)) +
  geom_histogram()
ggsave('cluster_plots/histogram_tmm_factors.pdf')

all_tmm_factors %>% arrange(tmm_factor)




##Pool unweighted nTPMs together
nTPM_input_dir  <- 'nTPM_output_unweighted/'
output_dir <- 'pooled_data/'


nTPM_files = list.files(nTPM_input_dir, pattern = 'nTPM.tsv')

nTPM_unweighted <- tibble()

for (nTPM_filename in nTPM_files){
  
  tissue = strsplit(nTPM_filename, split = '_')[[1]][1]
  method = strsplit(nTPM_filename, split = '_')[[1]][2]
  sample = paste(tissue, method, sep = '_')
  nTPM <- read_delim(paste0(nTPM_input_dir,nTPM_filename), delim = '\t') 
  nTPM <- nTPM %>% gather(celltype_cluster, nTPM, -1)
  nTPM$sample <- sample
  nTPM <-nTPM %>%  select(genes, celltype_cluster, sample, nTPM)
  nTPM_unweighted <- bind_rows(nTPM_unweighted, nTPM)
  
}

#nTPM_unweighted <- nTPM_unweighted %>% unite(cell_id, celltype_cluster, sample, sep = '_', remove = FALSE)

nTPM_unweighted %>% write_csv(paste0(output_dir,'pooled_nTPM_unweighted.csv'))




######## Generating weighted nTPM - with mixed and low q clusters
###get weighted nTPM
nTPM_unweighted_directory <- 'pooled_data/'
metadata_input_dir <- 'pooled_data/'
weighted_nTPM_output <- 'pooled_data/'

nTPM_unweighted <-  read_csv(paste0(nTPM_unweighted_directory,'pooled_nTPM_unweighted.csv'))

naming_scheme <- read_csv(paste0(metadata_input_dir, 'celltype_name_table.csv')) %>% 
  select(celltype_cluster, manual_celltype_annotation, sample_celltype_consensus, cell_type_name, group_name, color, mixture, pooled_sample)

cluster_count <- read_csv(paste0(metadata_input_dir,'wMix/metadata_nTPM_pooling_wMix.csv' )) %>% select(celltype_cluster, cell_count, pooled_sample) %>% distinct()

metadata <- naming_scheme %>% left_join(cluster_count, by = c('celltype_cluster', 'pooled_sample'))

nTPM_unweighted <- nTPM_unweighted %>% left_join(metadata, by = join_by(celltype_cluster == celltype_cluster, sample == pooled_sample)) 

nTPM_unweighted %>% select(celltype_cluster,sample_celltype_consensus, cell_count, sample) %>% distinct()

#if want to keep it more specfic: change cell_type_name for sample_celltype_consensus
#nTPM_unweighted <- nTPM_unweighted %>% left_join(naming_scheme %>% select(manual_celltype_annotation,cell_type_name ) %>% unique(), by = #join_by(manual_celltype_annotation))


weighted_nTPM_per_sample = tibble()
#sample_name = "Frontal lobe_snRNA-seq"

for (sample_name in nTPM_unweighted$sample %>% unique() ){
  #celltype_cell_count <- nTPM_unweighted %>%  filter(sample == sample_name) %>% select(cell_type_name, cell_count) %>% unique() %>% group_by(cell_type_name) %>% summarise(total_cells_in_celltype = sum(cell_count))
  
  celltype_cell_count <- nTPM_unweighted %>%  filter(sample == sample_name) %>%
    select(-genes, -nTPM) %>% distinct() %>%  
    select(sample_celltype_consensus, cell_count) %>% group_by(sample_celltype_consensus) %>% summarise(total_cells_in_celltype = sum(cell_count))
  
  weighted_nTPM <- nTPM_unweighted %>%  filter(sample == sample_name) %>% 
    left_join(celltype_cell_count, by = join_by(sample_celltype_consensus == sample_celltype_consensus)) %>% 
    mutate(weighted_ntpm = nTPM * cell_count / total_cells_in_celltype) %>% 
    group_by(sample_celltype_consensus, genes) %>% 
    summarise(nTPM = sum(weighted_ntpm))
  weighted_nTPM['sample_name'] = sample_name
  
  weighted_nTPM_per_sample <- rbind(weighted_nTPM_per_sample, weighted_nTPM)
  
}

weighted_nTPM_per_sample %>% group_by(sample_celltype_consensus, sample_name) %>% summarise(sum(nTPM))
weighted_nTPM_per_sample %>% group_by(sample_celltype_consensus, sample_name) %>% count()


weighted_nTPM_per_sample<- weighted_nTPM_per_sample %>% unite(cell_id, sample_celltype_consensus, sample_name, sep = '_' , remove = FALSE)


weighted_nTPM_per_sample %>% select(genes, sample_celltype_consensus, sample_name, nTPM, cell_id ) %>% write_csv(paste0(weighted_nTPM_output, 'weighted_nTPM_per_sample_celltype_consensus.csv'))




######## Generating weighted nTPM - WITHOUT mixed and low q clusters
##get weighted
nTPM_unweighted_directory <- 'pooled_data/'
metadata_input_dir <- 'pooled_data/'
weighted_nTPM_output <- 'pooled_data/'

nTPM_unweighted <-  read_csv(paste0(nTPM_unweighted_directory,'pooled_nTPM_unweighted.csv'))

naming_scheme <- read_csv(paste0(metadata_input_dir, 'celltype_name_table.csv')) %>% 
  select(celltype_cluster, manual_celltype_annotation, sample_celltype_consensus, cell_type_name, group_name, color, mixture, pooled_sample)

cluster_count <- read_csv(paste0(metadata_input_dir,'wMix/metadata_nTPM_pooling_wMix.csv' )) %>% select(celltype_cluster, cell_count, pooled_sample) %>% distinct()

metadata <- naming_scheme %>% left_join(cluster_count, by = c('celltype_cluster', 'pooled_sample'))

nTPM_unweighted <- nTPM_unweighted %>% left_join(metadata, by = join_by(celltype_cluster == celltype_cluster, sample == pooled_sample)) 


##FILTER FOR WIHTOUT MIXES
nTPM_unweighted <- nTPM_unweighted %>% filter(mixture == FALSE)

nTPM_unweighted %>% select(celltype_cluster,cell_type_name, cell_count, sample) %>% distinct()

#if want to keep it more specfic: change cell_type_name for sample_celltype_consensus
#nTPM_unweighted <- nTPM_unweighted %>% left_join(naming_scheme %>% select(manual_celltype_annotation,cell_type_name ) %>% unique(), by = #join_by(manual_celltype_annotation))


weighted_nTPM_per_sample = tibble()
#sample_name = "Frontal lobe_snRNA-seq"

for (sample_name in nTPM_unweighted$sample %>% unique() ){
  #celltype_cell_count <- nTPM_unweighted %>%  filter(sample == sample_name) %>% select(cell_type_name, cell_count) %>% unique() %>% group_by(cell_type_name) %>% summarise(total_cells_in_celltype = sum(cell_count))
  
  celltype_cell_count <- nTPM_unweighted %>%  filter(sample == sample_name) %>%
    select(-genes, -nTPM) %>% distinct() %>%  
    select(cell_type_name, cell_count) %>% group_by(cell_type_name) %>% summarise(total_cells_in_celltype = sum(cell_count))
  
  weighted_nTPM <- nTPM_unweighted %>%  filter(sample == sample_name) %>% 
    left_join(celltype_cell_count, by = join_by(cell_type_name == cell_type_name)) %>% 
    mutate(weighted_ntpm = nTPM * cell_count / total_cells_in_celltype) %>% 
    group_by(cell_type_name, genes) %>% 
    summarise(nTPM = sum(weighted_ntpm))
  weighted_nTPM['sample_name'] = sample_name
  
  weighted_nTPM_per_sample <- rbind(weighted_nTPM_per_sample, weighted_nTPM)
  
}
#weighted_nTPM%>% group_by(cell_type_name, sample_name)%>% count()# %>% summarise(sum(nTPM))
#weighted_nTPM_per_sample %>% group_by(cell_type_name, sample_name) %>% summarise(sum(nTPM))
#weighted_nTPM_per_sample %>% group_by(cell_type_name, sample_name) %>% count()


weighted_nTPM_per_sample<- weighted_nTPM_per_sample %>% unite(cell_id, cell_type_name, sample_name, sep = '_' , remove = FALSE)


weighted_nTPM_per_sample %>% select(genes, cell_type_name, sample_name, nTPM, cell_id ) %>% write_csv(paste0(weighted_nTPM_output, 'weighted_nTPM_per_sample.csv'))







#####Get global weighted


weighted_nTPM_input <- 'pooled_data/'
weighted_nTPM_per_sample <- read_csv(paste0(weighted_nTPM_input, 'weighted_nTPM_per_sample.csv'))
weighted_nTPM_output <- 'pooled_data/'

#weighted_nTPM_per_sample %>% group_by(sample_name, cell_type_name) %>% count()
weighted_nTPM_global <- weighted_nTPM_per_sample %>% group_by(cell_type_name, genes) %>% 
  summarise(nTPM = mean(nTPM))


weighted_nTPM_global %>% write_csv(paste0(weighted_nTPM_output, 'weighted_nTPM_global.csv'))



