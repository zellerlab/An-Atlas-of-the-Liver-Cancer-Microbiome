#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Association between HCC tumor tissue microbial compositons (from bulk RNA-Seq cohorts) with host-gene expression profiles.
# by Fabian Springer
########

library(tidyverse)
library(here)
library(parallel)
library(pbapply)
library(msigdbr)
library(decoupleR)


source(here("src","analysis","functions_analysis.R"))
#source(here("src","plotting","functions_plotting.R"))

#* Import corrected rel. abundance matrix ----
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))

#* Import host gene expression data ----
host_gene_expression_mat <- readRDS(here("data","raw","hostGeneExp_allTumors.rds"))
# The host gene expression matrix contains Trimmed-means of M-values (TMM)-normalized gene counts per million (CPM)
# For host gene expression analysis only tumor sampples were considered with at least 1e6 total gene counts of protein-coding genes.

# 19 Tumor samples of the DKFZ cohort did not meet this threshold and will not be considered for downstream analyses of microbiome vs host gene expression
meta_combined_df %>%
  filter(Tissue_type == "Tumor") %>%
  filter(!str_detect(Dataset, "16S")) %>% 
  filter(!(Sample_ID %in% colnames(host_gene_expression_mat))) %>% 
  dplyr::select(Sample_ID,Dataset)

#* Add bacterial counts per million sequencing reads to the metadata  ----
meta_test_df <- meta_combined_df %>%
  filter(Sample_ID %in% colnames(host_gene_expression_mat)) %>%
  dplyr::select(Sample_ID, Dataset, Batch, FastQ_libsize, Tissue_type,Inflammation_status) %>% 
  mutate(Dataset = ifelse(str_detect(Dataset,"DKFZ"),"DKFZ",Dataset)) # remove _RNA-Seq suffix from DKFZ since 5R16S not analysed here

# Compute bacterial counts per million sequencing reads and merge with metadata (from raw bacterial counts)
genus_counts_raw_mat <- readRDS(here("data","raw","rawScores_bulkRNAseq_genus.rds"))
total_bact_counts_df <- genus_counts_raw_mat[,meta_test_df$Sample_ID] %>% colSums() %>% enframe(name = "Sample_ID",value = "Bacteria_total")
meta_test_df <- meta_test_df %>%
  left_join(., total_bact_counts_df) %>%
  mutate(Bacteria_CPM = Bacteria_total / FastQ_libsize * 1e6) #fix later: for one Barcelona sample the FastQ_libsize is missing
rm(genus_counts_raw_mat)

#* Prepare the testing: ----
samples_for_testing <- meta_test_df$Sample_ID
# 1. Gene expression vs tumor inflammation status (inflamed vs non-inflamed)
# 2. Gene expression vs total bacterial load (bacterial CPM)
# 3. Gene expression vs relative abundance of individual genera

# The functions used for the testing are identical to the ones in the HCC-meta analyses:
# Linear (mixed) models for every gene vs every comparison (e.g. inflammation status, bacterial load) with dataset/batch as random effect (e.g. Gene1~TotalBactCPM + (1|Batch))
# For the continuous comparisons (bacteria CPM and genus rel. abudannces) bacterial features will be used as predictor.
# All tests will be performed with log10-transformed features.

# The linear models are implemented in a testing function that expects two matrices as input and builds models
# for every row in matrix 1 vs every row in matrix 2.

# Matrix with total bacterial load (bact. CPM and Inflammation status)
meta_test_df$Dataset %>% unique()
test_mat1 <- meta_test_df %>%  
  dplyr::select(Sample_ID, Dataset, Inflammation_status, Bacteria_CPM) %>%
  mutate(Bacteria_CPM = log10(Bacteria_CPM + 1e-5)) %>% # log-transform bacterial counts per million sequencing reads

  mutate(
    # Inflammation status with all datasets combined and each dataset individually
    test_Inflammed_vs_nonInflamed_All = Inflammation_status,
    test_Inflammed_vs_nonInflamed_DKFZ = ifelse(Dataset == "DKFZ", Inflammation_status, NA),
    test_Inflammed_vs_nonInflamed_IDIBAPS = ifelse(Dataset == "ISMMS-IDIBAPS", Inflammation_status, NA),
    test_Inflammed_vs_nonInflamed_TCGA = ifelse(Dataset == "TCGA", Inflammation_status, NA),

    # Bacterial load with all datasets combined and each dataset individually
    test_total_Bacteria_All = Bacteria_CPM,
    test_total_Bacteria_DKFZ = ifelse(Dataset == "DKFZ", Bacteria_CPM, NA),
    test_total_Bacteria_IDIBAPS = ifelse(Dataset == "ISMMS-IDIBAPS", Bacteria_CPM, NA),
    test_total_Bacteria_TCGA = ifelse(Dataset == "TCGA", Bacteria_CPM, NA),
    test_total_Bacteria_INSERM = ifelse(Dataset == "INSERM", Bacteria_CPM, NA)
  ) %>%
  dplyr::select(-Inflammation_status, -Bacteria_CPM, -Dataset) %>%   
  gather(-Sample_ID, key = "feat", value = "val") %>%
  pivot_wider(names_from = Sample_ID, values_from = val) %>%
  column_to_rownames("feat") %>%
  as.matrix()
rowSums(!is.na(test_mat1)) #For the INSERM cohort the inflammation status is not available
# Reorder columns (Sample_IDs)
test_mat1 <- test_mat1[,samples_for_testing]

# Matrix with bacterial relative abundances on genus level
test_mat_genus <- genus_RelAbundance_corrected_mat[,samples_for_testing]
# log-transform
test_mat_genus_l10 <- log10(test_mat_genus+1e-5)

# generate dataframe with random-effects
meta_randomEff_df <-
  meta_test_df %>%
  dplyr::select(Sample_ID, Dataset, Batch) %>%
  mutate(rowN = Sample_ID) %>%
  column_to_rownames("rowN") %>%
  as.data.frame()
meta_randomEff_df <- meta_randomEff_df[samples_for_testing,]


# log-transform and re-order gene expression matrix
host_gene_expression_mat_l10 <- log10(host_gene_expression_mat[,samples_for_testing]+1e-2)
#hist(host_gene_expression_mat_l10,breaks = 100)
dim(host_gene_expression_mat_l10)

stopifnot(colnames(host_gene_expression_mat_l10) == colnames(test_mat1))
stopifnot(colnames(host_gene_expression_mat_l10) == colnames(test_mat_genus_l10))
stopifnot(rownames(meta_randomEff_df) == colnames(test_mat_genus_l10))

#* Run the testing for total bacteria and inflammation status ----
# If possible, run on a compute cluster since this consumes a lot of resources: 1 linear (mixed) model per gene and comparison
nrow(test_mat1) * nrow(host_gene_expression_mat_l10)
n_cores = 16

# Generate vector that indicates which rows in test_mat1 are continuous and which are categorical
cont_or_cat_vec <- case_when(
  str_detect(rownames(test_mat1), "Inflammed") ~ "categorical",
  str_detect(rownames(test_mat1), "total") ~ "continuous"
)

test_res_df <-
  f_run_linear_models_parallel(
    dset_name = "HostGenes_vs_Bacteria",
    mat1 = test_mat1,
    mat2 = host_gene_expression_mat_l10,
    meta = meta_randomEff_df,
    random_effect_variable = "Batch", 
    threshold_for_prev = -2, # Set to pseudocount, prevalence will always be one (prevalence will be counted for genes which does not make sense anyways)
    n_cores_max = n_cores,
    compute_CI = FALSE,
    cont_or_cat_vec = cont_or_cat_vec
  ) %>%
  group_by(feat1) %>% # Group by comparison (e.g. Inflammation_Status_All) before p-value correction)
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

datasets <- meta_test_df$Dataset %>% unique()

# Clean up the results
test_res_sampleMetrics_clean_df <-
  test_res_df %>%
  glimpse() %>%
  mutate(feat1 = str_replace(feat1, "IDIBAPS", "ISMMS-IDIBAPS")) %>%
  mutate(
    Dataset = case_when(
      str_detect(feat1, datasets[1]) ~ datasets[1],
      str_detect(feat1, datasets[2]) ~ datasets[2],
      str_detect(feat1, datasets[3]) ~ datasets[3],
      str_detect(feat1, datasets[4]) ~ datasets[4],
      str_detect(feat1, datasets[5]) ~ datasets[5],
      TRUE ~ "All"
    )
  ) %>%
  mutate(
    testing = str_remove(feat1, paste0("_", Dataset)),
    testing = str_remove(testing, "test_"),
    testing = str_remove(testing, paste0("_", Group2))
  ) %>%
  transmute(
    comparison = testing,
    gene = feat2, Group1, Group2,
    effect.size = effect_size, t_value, p.val = p_value, p.val_adj = p_adj,
    N_samples = N_Samples,N_Group1, N_Group2, 
    Dataset
  ) %>%
  mutate(comparison_Dataset = paste(Dataset, comparison, sep = "_")) %>%
  glimpse()
# save as RDS
saveRDS(test_res_sampleMetrics_clean_df, here("data","results","genes_vs_sampleMetrics_df.rds"))

#*  Run the testing for every combination of genus abundance and gene expression ----

nrow(test_mat_genus_l10) * nrow(host_gene_expression_mat_l10)
n_cores = 16

# all genus rel. abundances are continuous
cont_or_cat_vec_genus <- rep("continuous", nrow(test_mat_genus_l10))

test_res_genus_df <-
  f_run_linear_models_parallel(
    dset_name = "HostGenes_vs_Genera",
    mat1 = test_mat_genus_l10[,,drop=F],
    mat2 = host_gene_expression_mat_l10,
    meta = meta_randomEff_df,
    random_effect_variable = "Batch", 
    threshold_for_prev = -2, # Set to pseudocount, prevalence will always be one (prevalence will be counted for genes which does not make sense anyways)
    n_cores_max = n_cores,
    compute_CI = FALSE,
    cont_or_cat_vec = cont_or_cat_vec_genus
  ) %>%  
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Clean up the results
test_res_genus_clean_df <-
  test_res_genus_df %>%
  glimpse() %>%
  mutate(feat1 = str_replace(feat1, "IDIBAPS", "ISMMS-IDIBAPS")) %>%
  mutate(
    Dataset = case_when(
      str_detect(feat1, datasets[1]) ~ datasets[1],
      str_detect(feat1, datasets[2]) ~ datasets[2],
      str_detect(feat1, datasets[3]) ~ datasets[3],
      str_detect(feat1, datasets[4]) ~ datasets[4],
      str_detect(feat1, datasets[5]) ~ datasets[5],
      TRUE ~ "All"
    )
  ) %>%
  mutate(
    testing = str_remove(feat1, paste0("_", Dataset)),
    testing = str_remove(testing, "test_")   
  ) %>%
  transmute(
    genus = testing,
    gene = feat2,
    effect.size = effect_size, t_value, p.val = p_value, p.val_adj = p_adj,
    N_samples = N_Samples,
    Dataset
  ) %>%
  mutate(genus_Dataset = paste(Dataset, genus, sep = "_")) %>%
  glimpse()

# save as RDS
saveRDS(test_res_genus_clean_df, here("data","results","genes_vs_bacteria_df.rds"))


#*#############################
#* Part II: Infer pathway activities based on microbiome / gene associations ----
#*#############################
# Pathway activities are inferred based on the association of gene expression with genus abundances, total bacterial load or inflammation status.
# We use the decoupleR tool that computes linaer models to predict the pathway activities of every hallmark pathway based on the provided gene 'scores'.
# We provide the t-values of the linear models from the testings performed above as scores since they contain information about the directionality as well the significance
# of each association between gene expression and genus abundances / total bacterial load / inflammation status.
# For more detailed explanation see: https://doi.org/10.1093/bioadv/vbac016 or https://saezlab.github.io/decoupleR/articles/pw_bk.html

# Import results of the individual gene testings
test_res_sampleMetrics_clean_df <- readRDS(here("data","results","genes_vs_sampleMetrics_df.rds"))
test_res_genus_clean_df <- readRDS(here("data","results","genes_vs_bacteria_df.rds"))

# import msigDB Hallmarks gene sets
hallmarks_msigdb <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
hm_network <- hallmarks_msigdb %>%
    mutate(mor = 1) %>%
    dplyr::select(gs_name, gene_symbol,mor) %>%
    distinct()

f_run_association <- function(df, meta_var_to_split) {
  datasets <- df$Dataset %>% unique()
  # compute the pathway activities based on every comparison individually (e.g. Gene vs Total bacteria in all datasets, in DKFZ only,...)
  all_features <- sort(unique(df[[meta_var_to_split]]))
  res_df <- tibble()
  feat <- all_features[1]
  # loop over the input dataset, generate a matrix as input for the decoupleR tool.
  for (feat in all_features) {
    c_df <- df %>%
      filter(!!as.symbol(meta_var_to_split) == feat) %>%
      filter(!is.na(p.val)) # make sure to remove incomplete cases
    if (nrow(c_df) < 1) {
      message(feat, " is empty - skipping")
      next
    }
    message("-----", feat, "----")

    # generate matrix with 1 column (every gene is a row, the current feature (e.g. Inflamed vs non-inflamed) is the column, values represent t-values from the linear models)
    tmp_mat <-
      c_df %>%
      dplyr::select(!!as.symbol(meta_var_to_split), gene, t_value) %>%
      pivot_wider(names_from = !!as.symbol(meta_var_to_split), values_from = t_value, values_fill = 0) %>%
      column_to_rownames("gene") %>%
      as.matrix()

    # run decoupleR against the hallmarks gene set
    gsea_hallmarks_res <-
      decoupleR::run_ulm(
        mat = tmp_mat,
        network = hm_network,
        .source = gs_name,
        .target = gene_symbol,
        .mor = mor
      ) %>%
      mutate(gene_set_name = "MSigDB_Hallmark")

    # combine with other DFs
    res_df <- bind_rows(res_df, gsea_hallmarks_res)
  }

  res_df <- res_df %>%
    mutate(
      dataset = case_when(
        str_detect(condition, datasets[1]) ~ datasets[1],
        str_detect(condition, datasets[2]) ~ datasets[2],
        str_detect(condition, datasets[3]) ~ datasets[3],
        str_detect(condition, datasets[4]) ~ datasets[4],
        str_detect(condition, datasets[5]) ~ datasets[5],
        TRUE ~ "All"
      )
    ) %>%
    mutate(feature = str_remove(condition, paste0(dataset, "_"))) %>%
    dplyr::relocate(source, feature, dataset, score, p_value) %>%
    dplyr::rename(p.val = p_value)

  return(res_df)
}

# Compute the pathway activities based on the gene vs sample metrics associations: split by comparison and Dataset
pw_activities_sampleMetrics_df <- f_run_association(test_res_sampleMetrics_clean_df, "comparison_Dataset") %>%
  group_by(feature) %>% # correct for multiple testing when comparing the same feature across datasets (more conservative than grouping by feature and dataset)
  mutate(p.val_adj = p.adjust(p.val, method = "fdr")) %>% 
  ungroup() %>% 
  dplyr::relocate(source,feature,dataset,score,p.val,p.val_adj)

pw_activities_genus_df <- f_run_association(test_res_genus_clean_df, "genus_Dataset") %>%
  mutate(p.val_adj = p.adjust(p.val, method = "fdr")) %>% 
  dplyr::relocate(source,feature,dataset,score,p.val,p.val_adj)

saveRDS(pw_activities_sampleMetrics_df, here("data","results","pathwayActivities_bySampleMetrics_df.rds"))
saveRDS(pw_activities_genus_df, here("data","results","pathwayActivities_byGenus_df.rds"))

