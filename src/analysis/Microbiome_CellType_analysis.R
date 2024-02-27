#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Association between HCC tumor tissue microbial compositons (from bulk RNA-Seq cohorts) with host-cell type abundances.
# by Fabian Springer
########

library(tidyverse)
library(here)
library(parallel)
library(pbapply)

source(here("src","analysis","functions_analysis.R"))

#* Import corrected rel. abundance matrix ----
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))

#* Import host gene expression data ----
host_celltype_proportion_mat <- readRDS(here("data","raw","celltypeProportions_allTumors.rds"))

#* Add bacterial counts per million sequencing reads, shannon diversity and genus richness to the metadata  ----
meta_test_df <- meta_combined_df %>%
  filter(Sample_ID %in% colnames(host_celltype_proportion_mat)) %>%
  dplyr::select(Sample_ID, Dataset, Batch, FastQ_libsize) %>% 
  mutate(Dataset = ifelse(str_detect(Dataset,"DKFZ"),"DKFZ",Dataset)) # remove _RNA-Seq suffix from DKFZ since 5R16S not analysed here

# Import raw PathSeq scores for the RNA-Seq samples
genus_counts_raw_mat <- readRDS(here("data","raw","rawScores_bulkRNAseq_genus.rds"))
# Get total bacterial scores and compute bacterial counts per million sequencing reads
total_bact_counts_df <- genus_counts_raw_mat[, meta_test_df$Sample_ID] %>%
  colSums() %>%
  enframe(name = "Sample_ID", value = "Bacteria_total") %>%
  left_join(., meta_test_df %>% dplyr::select(Sample_ID, FastQ_libsize)) %>%
  mutate(`Total bacteria (CPM)` = Bacteria_total / FastQ_libsize * 1e6)

set.seed(1)
# Perform rarefaction and calculate diversity and richness
count_rar <- as.matrix(t(vegan::rrarefy(t(round(genus_counts_raw_mat[,meta_test_df$Sample_ID],0)),1000)))
# There is 1 sample with a rrarefied total number of bacteria smaller than 1000 (979).
# This is because the sum of the PathSeq scores in this sample becomes less than 1000 when rounded to integers.
count_rar_rel <- prop.table(count_rar,2)
div <- vegan::diversity(t(count_rar_rel), index = "shannon") %>% enframe(name = "Sample_ID",value = "Shannon diversity")
rich <- colSums(count_rar_rel > 0) %>% enframe(name = "Sample_ID",value = "Genus richness")
diversity_richness_df <- full_join(div,rich,by = "Sample_ID")

additional_metrics_df <- left_join(diversity_richness_df,total_bact_counts_df %>% dplyr::select(Sample_ID,Bacteria_CPM))

#join with the bacterial relative abundance matrix
bac_feature_mat <- bind_rows(
  genus_RelAbundance_corrected_mat[, samples_for_testing] %>%
    as_tibble(rownames = "feature") %>%
    gather(-feature, key = "Sample_ID", value = "value"),
  additional_metrics_df %>%
    gather(-Sample_ID, key = "feature", value = "value")
) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("feature") %>% 
  as.matrix()
    
# log-transform genus relative abundances (and shannon/richness/bactCPM metrics)
l10_bac_feature_mat <- log10(bac_feature_mat + 1e-5)

# log transform celltype proportions
l10_celltype_proportion_mat <- log10(host_celltype_proportion_mat + 1e-5)

#reorder matrices and generate dataframe with random effects
samples_for_testing <- meta_test_df$Sample_ID
l10_bac_feature_mat <- l10_bac_feature_mat[,samples_for_testing]
l10_celltype_proportion_mat <- l10_celltype_proportion_mat[,samples_for_testing]

meta_randomEff_df <-
  meta_test_df %>%
  dplyr::select(Sample_ID, Dataset, Batch) %>%
  mutate(rowN = Sample_ID) %>%
  column_to_rownames("rowN") %>%
  as.data.frame()
meta_randomEff_df <- meta_randomEff_df[samples_for_testing,]

stopifnot(colnames(l10_bac_feature_mat) == colnames(l10_celltype_proportion_mat))
stopifnot(colnames(l10_bac_feature_mat) == rownames(meta_randomEff_df))

nrow(l10_bac_feature_mat) * nrow(l10_celltype_proportion_mat)
n_cores = 6

# Generate vector that indicates which rows in test_mat1 are continuous and which are categorical
cont_or_cat_vec <- rep("continuous", nrow(l10_bac_feature_mat))
# define whether a feature is continuous or categorical
# Run linear models to predict celltype abundance based on the microbial features
test_res_df <-
  f_run_linear_models_parallel(
    dset_name = "CelltypeProportions_vs_Bacteria",
    mat1 = l10_bac_feature_mat,
    mat2 = l10_celltype_proportion_mat,
    meta = meta_randomEff_df,
    random_effect_variable = "Batch", 
    threshold_for_prev = -5, # Set to pseudocount, prevalence will always be one (prevalence will be counted for genes which does not make sense anyways)
    n_cores_max = n_cores,
    compute_CI = FALSE,
    cont_or_cat_vec = cont_or_cat_vec
  ) %>%  
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
  transmute(
    genus = feat1,
    celltype = feat2,
    effect.size = effect_size, t_value, p.val = p_value, p.val_adj = p_adj,
    N_samples = N_Samples,
    Dataset
  ) %>%
  glimpse()
# save as RDS
saveRDS(test_res_sampleMetrics_clean_df, here("data","results","celltypeProportion_vs_bacteria.rds"))









