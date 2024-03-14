#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Analysis of differential abundances on species level and associations with metadata
# by Fabian Springer
########

library(tidyverse)
library(here)
library(parallel)
library(pbapply)
library(progress)

source(here("src","analysis","functions_analysis.R"))
source(here("src","plotting","functions_plotting.R"))


#* Import 16S species level data and metadata ----
meta_all_df <- read_tsv(here("data","metadata","meta_5R16S_allSamples.tsv"))
species_relAbundance_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_allSamples.rds"))
log10_species_relAbundance_mat <- log10(species_relAbundance_mat + 1e-5)
log10_species_relAbundance_mat[1:5,1:5]

#* Also import the immunotheryapy cohort
meta_immunotherapy_df <- read_tsv(here("data","metadata","meta_5R16S_ImmunotherapyCohort.tsv"))
species_relAbundance_immunotherapy_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_ImmunotherapyCohort.rds"))


#* Do the following things:
# 1) Test differential species abundance across HCC etiologies and tumor vs adj. non-tumor tissues
# 2) Test associations of species with clinical parameters in HCC, CRLM and iCCA
# 3) Test associations of species with response to immunotherapy
# 4) Compute alpha and beta-diversities across etiologies and within immunotherapy cohort


table(meta_all_df$Etiology)
HCC_types <- c("HBV_HCC","HCV_HCC","ALD/ASH_HCC","MAFLD/MASH_HCC","other_HCC")
summary(HCC_types %in% meta_all_df$Etiology)

tmp_df <- meta_all_df %>%
  dplyr::select(Sample_ID, Dataset, Etiology) %>%
  mutate(
    # Cancer comparisons
    test_HCC_vs_iCCA = case_when(Etiology %in% HCC_types ~ "HCC", Etiology == "iCCA" ~ "iCCA", TRUE ~ NA_character_),
    test_HCC_vs_CRLM = case_when(Etiology %in% HCC_types ~ "HCC", Etiology == "CRLM" ~ "CRLM", TRUE ~ NA_character_),
    test_iCCA_vs_CRLM = case_when(Etiology == "iCCA" ~ "iCCA", Etiology == "CRLM" ~ "CRLM", TRUE ~ NA_character_),
    test_HCC_iCCA_vs_CRLM = case_when(Etiology %in% c(HCC_types, "iCCA") ~ "HCC_iCCA", Etiology == "CRLM" ~ "CRLM", TRUE ~ NA_character_),

    # HCC subtype comparisons
    test_Viral_vs_nonViral_HCC = case_when(Etiology %in% c("HBV_HCC", "HCV_HCC") ~ "viral_HCC", Etiology %in% c("ALD/ASH_HCC", "MAFLD/MASH_HCC") ~ "non-viral_HCC", TRUE ~ NA_character_),
    test_ALD_vs_MAFLD = case_when(Etiology %in% c("ALD/ASH_HCC", "MAFLD/MASH_HCC") ~ Etiology, TRUE ~ NA_character_),
    test_HBV_vs_HCV = case_when(Etiology %in% c("HBV_HCC", "HCV_HCC") ~ Etiology, TRUE ~ NA_character_),
    test_Tumor_vs_Adj_HCC =
      case_when(Etiology %in% HCC_types ~ "HCC", Etiology == "Adj. non-tumor_HCC" ~ "Adj. non-tumor_HCC", TRUE ~ NA_character_),

    # Fibrosis vs HCC
    test_Fibrosis_early_vs_HCC = case_when(Etiology == "Liver fibrosis early stage" ~ "EarlyFib", Etiology %in% HCC_types ~ "HCC", TRUE ~ NA_character_),
    test_Fibrosis_late_vs_HCC = case_when(Etiology == "Liver fibrosis late stage" ~ "LateFib", Etiology %in% HCC_types ~ "HCC", TRUE ~ NA_character_),
    test_Fibrosis_all_vs_HCC = case_when(str_detect(Etiology, "fibrosis") ~ "Fibrosis", Etiology %in% HCC_types ~ "HCC", TRUE ~ NA_character_),
    test_Fibrosis_late_vs_early = case_when(Etiology == "Liver fibrosis late stage" ~ "LateFib", Etiology == "Liver fibrosis early stage" ~ "EarlyFib", TRUE ~ NA_character_),

    # BTC subtypes
    test_iCCA_vs_phCCA_dCCA = case_when(Etiology == "iCCA" ~ "iCCA", Etiology %in% c("phCCA", "dCCA") ~ "phCCA_dCCA", TRUE ~ NA_character_),
    test_iCCA_vs_GBC = case_when(Etiology == "iCCA" ~ "iCCA", Etiology == "GBC" ~ "GBC", TRUE ~ NA_character_),
    test_GBC_vs_phCCA_dCCA = case_when(Etiology == "GBC" ~ "GBC", Etiology %in% c("phCCA", "dCCA") ~ "phCCA_dCCA", TRUE ~ NA_character_),

    # Tumor vs normal of adj
    test_Tumor_vs_Adj_iCCA = case_when(Etiology == "iCCA" ~ "iCCA", Etiology == "Adj. non-tumor_CCC" ~ "Adj. non-tumor_CCC", TRUE ~ NA_character_),
    test_Tumor_vs_Adj_CRLM = case_when(Etiology == "CRLM" ~ "CRLM", Etiology == "Adj. non-tumor_CRLM" ~ "Adj. non-tumor_CRLM", TRUE ~ NA_character_)
  ) %>% 
  dplyr::select(-Etiology, -Dataset)

# Check the levels in each column
lapply(tmp_df[,-1], function(x) table(x, useNA = "always"))

  # Looks good, prepare matrix for testing
meta_test_mat <-
  tmp_df %>%
  gather(-Sample_ID, key = "test", value = "group") %>%
  pivot_wider(names_from = Sample_ID, values_from = group, values_fill = NA) %>%
  column_to_rownames("test") %>%
  as.matrix()



###
# Now start the testing for the different entities/etiologies
###

# Make sure that orders of samples in metadata, meta_test_mat and relAB matrix are identical
Sample_IDs <- colnames(species_relAbundance_mat)
meta_test_mat <- meta_test_mat[, Sample_IDs]
meta <- meta_all_df %>% select(Sample_ID) %>% filter(Sample_ID %in% Sample_IDs) %>% mutate(rowN = Sample_ID) %>% column_to_rownames("rowN") %>% mutate(dummy = "dummy")
meta <- meta[Sample_IDs,]
log10_relAB_mat <- log10(species_relAbundance_mat[, Sample_IDs] + 1e-5)

# Define number of CPUs to be used
num_cores <- parallel::detectCores()
n_cores <- 4 #Define cores to use
stopifnot(n_cores < num_cores)

# Run the testing
# Run linear models
test_res_lm_df <-
  f_run_linear_models_parallel(
    dset_name = "5R16S",
    mat1 = meta_test_mat,
    mat2 = log10_relAB_mat,
    meta = meta,
    random_effect_variable = "dummy", #select this to run without random effect
    threshold_for_prev = -5, # Relative abundance threshold (log10-transformed) to consider a genus "prevalent"
    prevalence_threshold = 0.05, #threshold for prevalence applied to the genera before testing
    n_cores_max = n_cores,
    compute_CI = FALSE) %>% # Only use true if you really need it - it slows down the computation massively.
  group_by(feat1) %>% # group by comparison (e.g. HCC vs non-tumor HCC) prior to p-value correction
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Fisher's test for comparison of prevalences
test_res_fisher_df <-
  f_run_fisher_test_parallel(
    dset_name = "5R16S",
    mat1 = meta_test_mat,
    mat2 = log10_relAB_mat,    
    threshold_for_prev = -5, # Relative abundance threshold (log10-transformed) to consider a genus "prevalent"
    prevalence_threshold = 0.05, #threshold for prevalence applied to the genera before testing
    n_cores_max = n_cores) %>% # Only use true if you really need it - it slows down the computation massively.
  group_by(feat1) %>% # group by comparison (e.g. HCC vs non-tumor HCC) prior to p-value correction
  mutate(p_adj = p.adjust(p.val_fisher, method = "fdr")) %>%
  ungroup()

# Combine Fisher's and linear model results
Etiology_combined_test_res_df <-
  full_join(
    test_res_lm_df %>%
      glimpse() %>%
      transmute(feat1, feat2, Group1, Group2,
        p.val_lm = p_value, p.val_adj_lm = p_adj, effect.size = effect_size,
        N_Group1, N_Group2, Prev_Group1, Prev_Group2
      ),
    test_res_fisher_df %>%
      glimpse() %>%
      transmute(feat1, feat2, Group1, Group2,
        p.val_fisher,
        p.val_adj_fisher = p_adj, odds.ratio = odds_ratio
      )
  ) %>%
  mutate(comparison = str_remove(feat1, "test_"),tax = feat2) %>% 
  mutate(comparison = str_remove(comparison,paste0("_",Group2))) %>% 
  dplyr::select(-feat1, -feat2) %>% 
  dplyr::relocate(-N_Group1, -N_Group2, -Prev_Group1, -Prev_Group2) %>% 
  dplyr::relocate(comparison,tax)


#* Run the same testing for the immunotherapy cohort ----
meta_test_immuno_mat <- meta_immunotherapy_df %>%
  dplyr::select(Sample_ID, Response_to_immunotherapy) %>%
  pivot_wider(names_from = Sample_ID, values_from = Response_to_immunotherapy) %>%
  mutate(feat = "Response_to_immunotherapy") %>%
  column_to_rownames("feat") %>%
  as.matrix(drop = F)
log10_relAB_ICI_mat <- log10(species_relAbundance_immunotherapy_mat[,colnames(meta_test_immuno_mat)] + 1e-5)
meta_ICI <- meta_immunotherapy_df %>%
  mutate(Sample_ID = factor(Sample_ID, levels = colnames(meta_test_immuno_mat))) %>%
  arrange(Sample_ID) %>%
  dplyr::select(Sample_ID) %>%
  mutate(rowN = Sample_ID, dummy = "dummy") %>% 
  column_to_rownames("rowN")
test_res_ICI_lm_df <-
  f_run_linear_models_parallel(
    dset_name = "5R16S_ICI",
    mat1 = meta_test_immuno_mat,
    mat2 = log10_relAB_ICI_mat,
    meta = meta_ICI,
    random_effect_variable = "dummy", 
    threshold_for_prev = -5, 
    prevalence_threshold = 0.05,
    n_cores_max = n_cores,
    compute_CI = TRUE) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# Fisher's test for comparison of prevalences
test_res_ICI_fisher_df <-
  f_run_fisher_test_parallel(
    dset_name = "5R16S_ICI",
    mat1 = meta_test_immuno_mat,
    mat2 = log10_relAB_ICI_mat,
    threshold_for_prev = -5, 
    prevalence_threshold = 0.05, 
    n_cores_max = n_cores) %>%   
  mutate(p_adj = p.adjust(p.val_fisher, method = "fdr")) %>%
  ungroup()
ICI_combined_test_res_df <-
  full_join(
    test_res_ICI_lm_df %>%
      glimpse() %>%
      transmute(feat1, feat2, Group1, Group2,
        p.val_lm = p_value, p.val_adj_lm = p_adj, effect.size = effect_size,
        N_Group1, N_Group2, Prev_Group1, Prev_Group2
      ),
    test_res_ICI_fisher_df %>%
      glimpse() %>%
      transmute(feat1, feat2, Group1, Group2,
        p.val_fisher,
        p.val_adj_fisher = p_adj, odds.ratio = odds_ratio
      )
  ) %>%
  mutate(comparison = str_remove(feat1, "test_"),tax = feat2) %>% 
  mutate(comparison = str_remove(comparison,paste0("_",Group2))) %>% 
  dplyr::select(-feat1, -feat2) %>% 
  dplyr::relocate(-N_Group1, -N_Group2, -Prev_Group1, -Prev_Group2) %>% 
  dplyr::relocate(comparison,tax)

# save testing results for etiological comparisons and immunotherapy response
saveRDS(Etiology_combined_test_res_df,file = here("data","results","5R16S_EtiologyTesting_result_df.rds"))
saveRDS(ICI_combined_test_res_df,file = here("data","results","5R16S_ImmunotherapyCohortTesting_result_df.rds"))


#*#########
#* Run the clinical testing ----
#*#########

# Consider only tumor samples of the HCC, CRLM and iCCA cohorts
clinical_all_cancers_df <- meta_all_df %>%
  mutate(cancer = case_when(Etiology %in% HCC_types ~ "HCC", Etiology == "iCCA" ~ "iCCA", Etiology == "CRLM" ~ "CRLM", TRUE ~ NA_character_)) %>%
  filter(!is.na(cancer)) %>%
  dplyr::select(-Dataset, -Patient_ID, -Batch, -Medical_center, -Etiology) %>%
  dplyr::relocate(cancer)

cancer_type  <- "HCC"
f_run_clinical_testing <- function(dset_name = "all", cancer_type, clinical_all_cancers_df, log10_relAB_mat) {
  # Run clinical association testing with selected cancer type
  n_cores <- 4
  # Subset dataframe
  clinical_test_mat <- clinical_all_cancers_df %>%
    filter(cancer == cancer_type) %>%
    dplyr::select(-cancer) %>%
    gather(-Sample_ID, key = "feat", value = "val") %>%
    pivot_wider(names_from = Sample_ID, values_from = val) %>%
    column_to_rownames("feat") %>%
    as.matrix()

  # Split continuous and categorical features; 
  # continuous features will be tested using Spearman correlations, categorical using linear-models and Fisher's exact test
  numeric_features <- names(clinical_all_cancers_df)[sapply(clinical_all_cancers_df, is.numeric)]
  categorical_features <- colnames(clinical_all_cancers_df)[!colnames(clinical_all_cancers_df) %in% c(numeric_features,"Sample_ID","cancer")]

  # Split continuous and categorical features
  cont_feat_mat <- clinical_test_mat[numeric_features, ]
  cat_feat_mat <- clinical_test_mat[categorical_features, ]

  # Subset rel. Abundance matrix
  Sample_IDs <- colnames(cont_feat_mat)
  tmp_log10_relAB_mat <- log10_species_relAbundance_mat[, Sample_IDs]
  meta <- clinical_all_cancers_df %>% select(Sample_ID) %>% filter(Sample_ID %in% Sample_IDs) %>% mutate(rowN = Sample_ID) %>% column_to_rownames("rowN") %>% mutate(dummy = "dummy")

  stopifnot(colnames(cont_feat_mat)==colnames(cat_feat_mat))

  # run spearman correlations for continuous features
  message("Running Spearman correlations for continuous features")
  spearman_res_df <-
    f_run_spearman(
      mat1 = cont_feat_mat[,,drop=F],
      mat2 = tmp_log10_relAB_mat[,,drop=F],
      prevalence_threshold = 0.05, # threshold for prevalence applied to the genera before testing
      threshold_for_prev = -5 # Relative abundance threshold (log10-transformed) to consider a genus "prevalent
    )  %>% group_by(feat1) %>% #* Important: Define how pvalues are controlled!
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    ungroup()


  # run linear-models for categorical features
  message("Running linear models for categorical features")
  test_res_lm_df <-
    f_run_linear_models_parallel(
      dset_name = "Meta-analysis",
      mat1 = cat_feat_mat,
      mat2 = tmp_log10_relAB_mat,
      meta = meta,
      random_effect_variable = "dummy", # run without a random effect (basic linear models)
      threshold_for_prev = -5,
      prevalence_threshold = 0.05,
      n_cores_max = n_cores,
      compute_CI = FALSE) %>%  # Only use true if you really need it - it slows down the computation massively.     
    group_by(feat1_group) %>% # control for every comparison group to correctly account for tests with multiple levels (e.g. Child Pugh C vs all, Child Pugh B vs all, etc.)
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
  
  # Run Fisher's exact test for categorical features as well
  message("Running Fisher's exact test for categorical features")
  test_res_fisher_df <-
    f_run_fisher_test_parallel(
      mat1 = cat_feat_mat,
      mat2 = tmp_log10_relAB_mat,
      threshold_for_prev = -5,
      prevalence_threshold = 0.05,
      n_cores_max = n_cores) %>% # Only use true if you really need it - it slows down the computation massively.
    group_by(feat1_group) %>% # control for every comparison group to correctly account for tests with multiple levels (e.g. Child Pugh C vs all, Child Pugh B vs all, etc.)
    mutate(p_adj = p.adjust(p.val_fisher, method = "fdr")) %>%
    ungroup()
  
  # Combine results
  combined_test_res_df <-
    full_join(
      test_res_lm_df %>%
        glimpse() %>%
        transmute(feat1_group, feat1, feat2, Group1, Group2,
          p.val_lm = p_value, p.val_adj_lm = p_adj, effect.size = effect_size,
          N_Group1, N_Group2, Prev_Group1, Prev_Group2
        ),
      test_res_fisher_df %>%
        glimpse() %>%
        transmute(feat1_group, feat1, feat2, Group1, Group2,
          p.val_fisher,
          p.val_adj_fisher = p_adj, odds.ratio = odds_ratio
        )
    ) %>%
    mutate(feature_type = "categorical", N_Samples = N_Group1 + N_Group2) %>%
    full_join(., spearman_res_df %>%
      transmute(feat1, feat2,
        p.val_spearman = p_value, p.val_adj_spearman = p_adj, Spearman.R = effect_size,
        N_Samples = as.numeric(N_samples), Prevalence = as.numeric(Prevalence)
      ) %>%
      mutate(feature_type = "continuous")) %>%
    dplyr::rename(comparison = feat1, tax = feat2,comparison_group = feat1_group) %>%
    dplyr::relocate(-N_Group1, -N_Group2, -Prevalence, -Prev_Group1, -Prev_Group2, -feature_type) %>%
    dplyr::relocate(comparison, tax, Group1, Group2, N_Samples) %>% 
    dplyr::relocate(-feature_type,-comparison_group)
  
  return(combined_test_res_df)
}

clin_features_HCC_df <- f_run_clinical_testing(cancer_type = "HCC", clinical_all_cancers_df = clinical_all_cancers_df, log10_relAB_mat = log10_relAB_mat)
clin_features_iCCA_df <- f_run_clinical_testing(cancer_type = "iCCA", clinical_all_cancers_df = clinical_all_cancers_df, log10_relAB_mat = log10_relAB_mat)
clin_features_CRLM_df <- f_run_clinical_testing(cancer_type = "CRLM", clinical_all_cancers_df = clinical_all_cancers_df, log10_relAB_mat = log10_relAB_mat)

saveRDS(clin_features_HCC_df,file = here("data","results","5R16S_ClinicalFeaturesAssociation_HCC_df.rds"))
saveRDS(clin_features_iCCA_df,file = here("data","results","5R16S_ClinicalFeaturesAssociation_iCCA_df.rds"))
saveRDS(clin_features_CRLM_df,file = here("data","results","5R16S_ClinicalFeaturesAssociation_CRLM_df.rds"))


#* Shannon, Richness and Beta-Diversity testing ----
tmp_df <- meta_all_df %>%
  dplyr::select(Sample_ID, Dataset, Etiology) %>%
  mutate(
    # Cancer comparisons
    test_HCC_vs_iCCA_vs_CRLM = case_when(Etiology %in% HCC_types ~ "HCC", Etiology == "iCCA" ~ "iCCA",Etiology == "CRLM" ~ "CRLM", TRUE ~ NA_character_),    

    # ALD vs MAFLD vs HBV vs HCV    
    test_ALD_vs_MAFLD_vs_HBV_vs_HCV = case_when(Etiology %in% c("ALD/ASH_HCC", "MAFLD/MASH_HCC","HBV_HCC", "HCV_HCC") ~ Etiology, TRUE ~ NA_character_),

    # iCCA vs GBC vs phCCA_dCC
    test_iCCA_GBC_vs_phCCA_dCCA = case_when(Etiology %in% c("GBC","iCCA") ~ Etiology, Etiology %in% c("phCCA", "dCCA") ~ "phCCA_dCCA", TRUE ~ NA_character_),

    # iCCA and CRLM vs Adj
    test_Tumor_vs_Adj_iCCA = case_when(Etiology == "iCCA" ~ "iCCA", Etiology == "Adj. non-tumor_CCC" ~ "Adj. non-tumor_CCC", TRUE ~ NA_character_),
    test_Tumor_vs_Adj_CRLM = case_when(Etiology == "CRLM" ~ "CRLM", Etiology == "Adj. non-tumor_CRLM" ~ "Adj. non-tumor_CRLM", TRUE ~ NA_character_),
    
    # Fibrosis --> HCC
    test_Fibrosis_to_HCC = case_when(str_detect(Etiology, "fibrosis|non-tumor_HCC") ~ Etiology, Etiology %in% c(HCC_types) ~ "HCC", TRUE ~ NA_character_)
    ) %>% 
  dplyr::select(-Etiology, -Dataset)

lapply(tmp_df[,-1], function(x) table(x, useNA = "always"))

# Looks good, prepare matrix for testing
meta_test_shannon_richness_mat <-
  tmp_df %>%
  gather(-Sample_ID, key = "test", value = "group") %>%
  pivot_wider(names_from = Sample_ID, values_from = group, values_fill = NA) %>%
  column_to_rownames("test") %>%
  as.matrix()

# Compute Shannon diversity
div <- vegan::diversity(t(species_relAbundance_mat[,meta_all_df$Sample_ID]),index = "shannon") %>% enframe(name = "Sample_ID",value = "ShannonDiv")

# Compute Richness
threshold_for_prev <- 0
rich <- colSums(species_relAbundance_mat>threshold_for_prev) %>% enframe(name = "Sample_ID",value = "SpeciesRichness")

# Combine diversity metrics
diversity_richness_df <- full_join(div,rich)


alpha_beta_diversity_df <- tibble()
i <- 2
for(i in seq(1,nrow(meta_test_shannon_richness_mat))){
  # Select sample IDs and their respective groups for the comparison, merge with diversity metrics
  message(comparison <- rownames(meta_test_shannon_richness_mat)[i])

  df <- meta_test_shannon_richness_mat[i, ] %>%
    enframe(name = "Sample_ID", value = "group") %>%
    filter(!is.na(group)) %>%
    left_join(., diversity_richness_df) %>% 
    dplyr::select(-Sample_ID) %>% 
    gather(-group,key = "metric",value = "value")
  
  # Compute wilcoxon test
  res_df <- f_pairwise_wilcoxon_tests(df) %>% mutate(comparison = comparison)
  res_df <- res_df %>%
    pivot_wider(names_from = metric, values_from = p.val_adj_wilcox) %>%
    dplyr::rename("p.val_adj_wilcox_richness_species" = "SpeciesRichness", "p.val_adj_wilcox_shannon_diversity" = "ShannonDiv")

  # Compute median and sample number
  median_df <- df %>%
    group_by(metric, group) %>%
    summarise(median = median(value), N_Samples = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = metric, values_from = c(median)) %>%
    dplyr::rename(median_richness_species = SpeciesRichness, median_shannon_diversity = ShannonDiv) %>%
    mutate(
      Group1 = group, Group2 = group,
      N_Group1 = N_Samples, N_Group2 = N_Samples,
      median_shannon_diversity_Group1 = median_shannon_diversity,
      median_shannon_diversity_Group2 = median_shannon_diversity,
      median_species_richness_Group1 = median_richness_species,
      median_species_richness_Group2 = median_richness_species
    )
  
  res_df <- res_df %>% 
    left_join(., median_df %>% select(c(contains("Group1")))) %>% 
    left_join(., median_df %>% select(c(contains("Group2"))))

  #* Compute Beta diversity (on Bray-Curtis and Euclidean distances)
  message("Computing beta-diversities")
  df <- meta_test_shannon_richness_mat[i, ] %>%
    enframe(name = "Sample_ID", value = "group") %>%
    filter(!is.na(group))

  permanova_res_df <- f_compute_distance_metrics(df,species_relAbundance_mat) %>% transmute(Group1,Group2,p.val_adj_permanova = p.adj,metric = paste(distance_metric,"distance",sep = " "))
  
  permanova_res_df <- permanova_res_df %>%
    pivot_wider(names_from = metric, values_from = p.val_adj_permanova) %>%
    dplyr::rename(p.val_adj_permanova_Bray = `Bray-Curtis distance`, p.val_adj_permanova_Euc = `Euclidean distance`)
  # Duplicate permanova result with opposite group names to match all cases with shannon and richness df
  permanova_res_df_full <-
    bind_rows(
      permanova_res_df,
      permanova_res_df %>%
        mutate(Group1c = Group2, Group2c = Group1) %>%
        transmute(Group2 = Group2c, Group1 = Group1c, p.val_adj_permanova_Bray, p.val_adj_permanova_Euc)
    )
  

  res_combined_df <- left_join(res_df,permanova_res_df_full)
  alpha_beta_diversity_df <- bind_rows(alpha_beta_diversity_df,res_combined_df)

}
alpha_beta_diversity_clean_df <- alpha_beta_diversity_df %>% 
  glimpse() %>% 
  mutate(comparison = str_remove(comparison, "test_")) %>%
  dplyr::relocate(
    comparison,Group1,Group2,N_Group1,N_Group2,
    median_species_richness_Group1,median_species_richness_Group2,p.val_adj_wilcox_richness_species,
    median_shannon_diversity_Group1,median_shannon_diversity_Group2,p.val_adj_wilcox_shannon_diversity,
    p.val_adj_permanova_Bray,p.val_adj_permanova_Euc
  )

#* Immunotherapy Cohort: Shannon, Richness and Beta-Diversity testing ----
meta_immunotherapy_df <- read_tsv(here("data","metadata","meta_5R16S_ImmunotherapyCohort.tsv"))
species_relAbundance_immunotherapy_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_ImmunotherapyCohort.rds"))

# Compute diversity
div <- vegan::diversity(t(species_relAbundance_immunotherapy_mat[,meta_immunotherapy_df$Sample_ID]),index = "shannon") %>% enframe(name = "Sample_ID",value = "ShannonDiv")
# Compute Richness
threshold_for_prev <- 0
rich <- colSums(species_relAbundance_immunotherapy_mat>threshold_for_prev) %>% enframe(name = "Sample_ID",value = "SpeciesRichness")
diversity_richness_df <- full_join(div,rich)

immuno_shannon_richness_df <- diversity_richness_df %>%
  gather(-Sample_ID, key = "metric", value = "value") %>%
  left_join(., meta_immunotherapy_df %>% dplyr::select(Sample_ID, Response_to_immunotherapy))
    
immuno_median_df <-
  immuno_shannon_richness_df %>%
  group_by(metric, Response_to_immunotherapy) %>%
  summarise(median = median(value), N_Samples = n()) %>%
  pivot_wider(names_from = metric, values_from = median)

df <- immuno_shannon_richness_df %>% dplyr::rename(group = Response_to_immunotherapy)
res_df <- f_pairwise_wilcoxon_tests(df) %>%
  pivot_wider(names_from = metric, values_from = p.val_adj_wilcox) %>%
  dplyr::rename("p.val_adj_wilcox_richness_species" = "SpeciesRichness", "p.val_adj_wilcox_shannon_diversity" = "ShannonDiv")
median_df <- df %>%
  group_by(metric, group) %>%
  summarise(median = median(value), N_Samples = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = metric, values_from = c(median)) %>%
  dplyr::rename(median_richness_species = SpeciesRichness, median_shannon_diversity = ShannonDiv) %>%
  mutate(
    Group1 = group, Group2 = group,
    N_Group1 = N_Samples, N_Group2 = N_Samples,
    median_shannon_diversity_Group1 = median_shannon_diversity,
    median_shannon_diversity_Group2 = median_shannon_diversity,
    median_species_richness_Group1 = median_richness_species,
    median_species_richness_Group2 = median_richness_species
  )
res_df <- res_df %>%
  left_join(., median_df %>% select(c(contains("Group1")))) %>%
  left_join(., median_df %>% select(c(contains("Group2"))))

#* Compute Bray-Curtis and Euclidean distances
df <- meta_immunotherapy_df %>%
  dplyr::select(Sample_ID, Response_to_immunotherapy) %>%
  dplyr::rename(group = Response_to_immunotherapy)
permanova_res_df <- f_compute_distance_metrics(
  df, relAB_mat = species_relAbundance_immunotherapy_mat) %>% transmute(Group1, Group2, p.val_adj_permanova = p.adj, metric = paste(distance_metric, "distance", sep = " "),
  threshold_for_prevalence = 0, prevalence_threshold = 0.05
)
permanova_res_df %>% glimpse()
permanova_res_df <- permanova_res_df %>%
  pivot_wider(names_from = metric, values_from = p.val_adj_permanova) %>%
  dplyr::rename(p.val_adj_permanova_Bray = `Bray-Curtis distance`, p.val_adj_permanova_Euc = `Euclidean distance`)
# Duplicate permanova result with opposite group names to match all cases with shannon and richness df
permanova_res_df_full <-
  bind_rows(
    permanova_res_df,
    permanova_res_df %>%
      mutate(Group1c = Group2, Group2c = Group1) %>%
      transmute(Group2 = Group2c, Group1 = Group1c, p.val_adj_permanova_Bray, p.val_adj_permanova_Euc)
  )
res_combined_df <- left_join(res_df, permanova_res_df_full)

alpha_beta_diversity_immuno_clean_df  <- 
  res_combined_df %>% 
  mutate(comparison = "Responder_vs_NonResponder") %>%
  glimpse() %>%   
  dplyr::relocate(
    comparison,Group1,Group2,N_Group1,N_Group2,
    median_species_richness_Group1,median_species_richness_Group2,p.val_adj_wilcox_richness_species,
    median_shannon_diversity_Group1,median_shannon_diversity_Group2,p.val_adj_wilcox_shannon_diversity,
    p.val_adj_permanova_Bray,p.val_adj_permanova_Euc
  )

# Save difersity comparisons
write_tsv(alpha_beta_diversity_clean_df,file = here("data","results","5R16S_diversity_comparisons_df.tsv"))
write_tsv(alpha_beta_diversity_immuno_clean_df,file = here("data","results","5R16S_ImmunotherapyCohort_diversity_comparisons_df.tsv"))