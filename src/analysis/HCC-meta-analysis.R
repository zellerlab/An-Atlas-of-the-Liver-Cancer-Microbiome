#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Perform meta-analysis of tumor tissue microbiomes across bulk-RNAseq and 5R16S cohorts
# by Fabian Springer
########

library(tidyverse)
library(here)
library(parallel)
library(pbapply)

source(here("src","analysis","functions_analysis.R"))
source(here("src","plotting","functions_plotting.R"))


#* Import corrected rel. abundance matrix ----
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))

table(meta_combined_df$Dataset)
table(meta_combined_df$Tissue_type)

# We do the following things: 
#1) Compute mean prevalence and mean abundance for every genus in every dataset
#2) Perform testing of differential abudnance in tumor vs. adj. non-tumor and in different etiologies
# Testing is performed using linear (mixed effects) models with the Dataset/Batch variable as random effect. 
colnames(meta_combined_df)
table(meta_combined_df$Dataset,meta_combined_df$Tissue_type)

#* Compute mean prevalence and mean abundance ----
threshold_for_prev_RNAseq <- 1e-3 #relative abundance in the RNA-Seq cohort for which a genus is considered prevalent
threshold_for_prev_5R16S <- 0 #same in the 5R16S cohort - In the preprocessing small values (<1e-4) were floored, so the threshold can be 0

tumor_abundance_prevalence_df <-
  genus_RelAbundance_corrected_mat %>%
  as_tibble(rownames = "genus") %>%
  gather(-genus, key = "Sample_ID", value = "rel") %>%
  left_join(., meta_combined_df %>% dplyr::select(Sample_ID, Tissue_type, Dataset)) %>%
  filter(Tissue_type == "Tumor") %>%
  mutate(thresh_for_prev = ifelse(str_detect(Dataset, "16S"), threshold_for_prev_5R16S, threshold_for_prev_RNAseq)) %>%
  group_by(Dataset) %>%
  mutate(N_samples = length(unique(Sample_ID))) %>%
  group_by(Dataset, genus) %>%
  mutate(
    prevalence = (sum(rel > thresh_for_prev)) / N_samples,
    mean_abundance = mean(rel)
  ) %>%
  ungroup() %>%
  dplyr::select(genus, Dataset, N_samples, prevalence, mean_abundance) %>%
  distinct() %>%
  arrange(genus)
 
tumor_abundance_prevalence_df %>%
  dplyr::select(-N_samples, -mean_abundance) %>% 
  pivot_wider(names_from = Dataset,values_from = prevalence) %>% 
  arrange(-TCGA)

tumor_abundance_prevalence_df %>%
  dplyr::select(-N_samples, -prevalence) %>% 
  pivot_wider(names_from = Dataset,values_from = mean_abundance) %>% 
  arrange(-TCGA)


#* Run the meta-analysis ----
# Define a dataframe with one row per comparison (e.g. Tumor vs Adj. Normal or Viral vs non-Viral HCC) and the Sample_IDs in columns.
# Take every row of the genus relative abundance matrix and build a linear (mixed) model for every row in the comparison matrix. 
colnames(meta_combined_df)

# Following tests will be performed: 
# Tumor vs Adj. non-tumor
# Inflamed vs non-inflamed 
# Viral vs non-viral
# MAFLD vs ALD
# HBV vs HCV
# This will be performed once in the meta-analysis setting and once within every Dataset individually
meta_combined_df$Etiology %>% table()
meta_combined_df$Dataset %>% table()
meta_combined_df %>% filter(Dataset == "ISMMS-IDIBAPS") %>% glimpse()
tmp_df <-
  meta_combined_df %>%
  dplyr::select(Sample_ID, Dataset, Tissue_type, Etiology, Inflammation_status) %>%
  mutate(
    test_Tumor_vs_Adj = Tissue_type,
    test_Inflamed_vs_nonInflamed_HCC = case_when(Tissue_type == "Tumor" ~ Inflammation_status, TRUE ~ NA_character_),
    test_Viral_vs_nonViral_HCC = case_when(
      Tissue_type == "Tumor" & Etiology %in% c("HBV_HCC", "HCV_HCC") ~ "viral_HCC",
      Tissue_type == "Tumor" & Etiology %in% c("ALD/ASH_HCC", "MAFLD/MASH_HCC") ~ "non-viral_HCC", TRUE ~ NA_character_
    ),
    test_ALD_vs_MAFLD = case_when(Tissue_type == "Tumor" & Etiology %in% c("ALD/ASH_HCC", "MAFLD/MASH_HCC") ~ Etiology, TRUE ~ NA_character_),
    test_HBV_vs_HCV = case_when(Tissue_type == "Tumor" & Etiology %in% c("HBV_HCC", "HCV_HCC") ~ Etiology, TRUE ~ NA_character_)
  ) %>%
  # Do the same but for the individual cohorts
  mutate(
    # DKFZ RNAseq
    test_Tumor_vs_Adj_DKFZ_RNAseq = case_when(Dataset == "DKFZ_RNAseq" ~ test_Tumor_vs_Adj,TRUE ~ NA_character_),
    test_Inflamed_vs_nonInflamed_HCC_DKFZ_RNAseq = case_when(Dataset == "DKFZ_RNAseq" ~ test_Inflamed_vs_nonInflamed_HCC,TRUE ~ NA_character_),
    test_Viral_vs_nonViral_HCC_DKFZ_RNAseq = case_when(Dataset == "DKFZ_RNAseq" ~ test_Viral_vs_nonViral_HCC,TRUE ~ NA_character_),
    test_ALD_vs_MAFLD_DKFZ_RNAseq = case_when(Dataset == "DKFZ_RNAseq" ~ test_ALD_vs_MAFLD,TRUE ~ NA_character_), 
    test_HBV_vs_HCV_DKFZ_RNAseq = case_when(Dataset == "DKFZ_RNAseq" ~ test_HBV_vs_HCV,TRUE ~ NA_character_),
    
    # DKFZ 5R16S
    test_Tumor_vs_Adj_DKFZ_5R16S = case_when(Dataset == "DKFZ_5R16S" ~ test_Tumor_vs_Adj,TRUE ~ NA_character_),
    test_Inflamed_vs_nonInflamed_HCC_DKFZ_5R16S = case_when(Dataset == "DKFZ_5R16S" ~ test_Inflamed_vs_nonInflamed_HCC,TRUE ~ NA_character_), 
    test_Viral_vs_nonViral_HCC_DKFZ_5R16S = case_when(Dataset == "DKFZ_5R16S" ~ test_Viral_vs_nonViral_HCC,TRUE ~ NA_character_),
    test_ALD_vs_MAFLD_DKFZ_5R16S = case_when(Dataset == "DKFZ_5R16S" ~ test_ALD_vs_MAFLD,TRUE ~ NA_character_), 
    test_HBV_vs_HCV_DKFZ_5R16S = case_when(Dataset == "DKFZ_5R16S" ~ test_HBV_vs_HCV,TRUE ~ NA_character_),

    # TCGA
    test_Tumor_vs_Adj_TCGA = case_when(Dataset == "TCGA" ~ test_Tumor_vs_Adj,TRUE ~ NA_character_), 
    test_Inflamed_vs_nonInflamed_HCC_TCGA = case_when(Dataset == "TCGA" ~ test_Inflamed_vs_nonInflamed_HCC,TRUE ~ NA_character_),
    test_Viral_vs_nonViral_HCC_TCGA = case_when(Dataset == "TCGA" ~ test_Viral_vs_nonViral_HCC,TRUE ~ NA_character_),
    test_ALD_vs_MAFLD_TCGA = case_when(Dataset == "TCGA" ~ test_ALD_vs_MAFLD,TRUE ~ NA_character_),
    test_HBV_vs_HCV_TCGA = case_when(Dataset == "TCGA" ~ test_HBV_vs_HCV,TRUE ~ NA_character_),

    # IDIBAPS
    test_Tumor_vs_Adj_IDIBAPS = case_when(Dataset == "ISMMS-IDIBAPS" ~ test_Tumor_vs_Adj,TRUE ~ NA_character_), 
    test_Inflamed_vs_nonInflamed_HCC_IDIBAPS = case_when(Dataset == "ISMMS-IDIBAPS" ~ test_Inflamed_vs_nonInflamed_HCC,TRUE ~ NA_character_),
    test_Viral_vs_nonInflassed_HCC_IDIBAPS = case_when(Dataset == "ISMMS-IDIBAPS" ~ test_Viral_vs_nonViral_HCC,TRUE ~ NA_character_),
    test_ALD_vs_NAFLD_IDIBAPS = case_when(Dataset == "ISMMS-IDIBAPS" ~ test_ALD_vs_MAFLD,TRUE ~ NA_character_),
    test_HBV_vs_HCV_IDIBAPS = case_when(Dataset == "ISMMS-IDIBAPS" ~ test_HBV_vs_HCV,TRUE ~ NA_character_),

    # INSERM
    test_Tumor_vs_Adj_INSERM = case_when(Dataset == "INSERM" ~ test_Tumor_vs_Adj,TRUE ~ NA_character_),
    test_Inflamed_vs_nonInflamed_HCC_INSERM = case_when(Dataset == "INSERM" ~ test_Inflamed_vs_nonInflamed_HCC,TRUE ~ NA_character_),
    test_Viral_vs_nonViral_HCC_INSERM = case_when(Dataset == "INSERM" ~ test_Viral_vs_nonViral_HCC,TRUE ~ NA_character_),
    test_ALD_vs_MAFLD_INSERM = case_when(Dataset == "INSERM" ~ test_ALD_vs_MAFLD,TRUE ~ NA_character_),
    test_HBV_vs_HCV_INSERM = case_when(Dataset == "INSERM" ~ test_HBV_vs_HCV,TRUE ~ NA_character_)
  ) %>%   
  dplyr::select(-Tissue_type, -Etiology, -Inflammation_status, -Dataset)
# Inspect whether every row contains the desired levels
#lapply(tmp_df[,-1], function(x) table(x, useNA = "always"))

meta_test_mat  <- 
  tmp_df %>% 
  gather(-Sample_ID, key = "test", value = "group") %>%
  pivot_wider(names_from = Sample_ID, values_from = group, values_fill = NA) %>%
  column_to_rownames("test") %>%
  as.matrix()
# Define a meta_df with the random effects that should be tested
meta_randomEff_df <-
  meta_combined_df %>%
  dplyr::select(Sample_ID, Dataset, Batch) %>%
  mutate(rowN = Sample_ID) %>%
  column_to_rownames("rowN") %>%
  as.data.frame()

# Make sure that orders of samples in metadata, meta_test_mat and relAB matrix are identical
Sample_IDs <- colnames(genus_RelAbundance_corrected_mat)
meta_test_mat <- meta_test_mat[, Sample_IDs]
meta_randomEff_df <- meta_randomEff_df[Sample_IDs, ]
log10_relAb_corrected_mat <- log10(genus_RelAbundance_corrected_mat[, Sample_IDs] + 1e-5)

# Define number of CPUs to be used
num_cores <- parallel::detectCores()
n_cores <- 4 #Define cores to use
stopifnot(n_cores < num_cores)

# This function does the following: It takes every row of mat1 (the meta-data matrix) and builds a vector of the respective non-NA samples with their labels
# then it loops over every row in mat2 (log10-transformed genus-level relative abudnance matrix) and performes a linear-mixed model of the form:
# log10(bac-rel) ~ group + (1|random_effect_variable) where random_effect_variable represents the Dataset/Batch variable

# Run the testing
test_res_df <-
  f_run_linear_models_parallel(
    dset_name = "Meta-analysis",
    mat1 = meta_test_mat,    
    mat2 = log10_relAb_corrected_mat,
    meta = meta_randomEff_df,
    random_effect_variable = "Batch", 
    threshold_for_prev = -3, # Relative abundance threshold (log10-transformed) to consider a genus "prevalent"
    n_cores_max = n_cores,
    compute_CI = TRUE # Only use true if you really need it - it slows down the computation massively.
  ) %>%
  group_by(feat1) %>% # Group by meta-variable (e.g. Viral- vs non-viral or Tumor vs Adj.non-tumor before p-value correction)
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

datasets <- meta_combined_df$Dataset %>% unique()
test_res_df$feat1 %>% unique()

# Restructure the meta-analysis output dataframe
# Add a column to test_res_df that searches for any dataset name in feat1
tmp_df <- test_res_df %>%
  mutate(feat1 = str_replace(feat1,"IDIBAPS","ISMMS-IDIBAPS")) %>% 
  mutate(
    Dataset = case_when(
      str_detect(feat1, datasets[1]) ~ datasets[1],
      str_detect(feat1, datasets[2]) ~ datasets[2],
      str_detect(feat1, datasets[3]) ~ datasets[3],
      str_detect(feat1, datasets[4]) ~ datasets[4],
      str_detect(feat1, datasets[5]) ~ datasets[5],
      TRUE ~ "Meta-analysis"
    )
  ) %>%
  mutate(
    testing = str_remove(feat1, paste0("_", Dataset)),
    testing = str_remove(testing, "test_"),
    testing = str_remove(testing, paste0("_", Group2))
  ) %>%
  transmute(
    comparison = testing,
    bacteria = feat2, Group1, Group2,
    effect.size = effect_size,lower95CI,upper95CI, t_value, p.val = p_value, p.val_adj = p_adj,
    N_samples = N_Group1+N_Group2,
    N_Group1, N_Group2,Prev_Group1, Prev_Group2,
    Dataset
  ) %>% 
  mutate(comparison_Dataset = paste(Dataset,comparison,sep = "_"))

#* Store the results of the HCC meta-analysis ----
tmp_df$comparison_Dataset %>% unique() %>% sort()
saveRDS(tmp_df,file = "./data/results/HCC-meta-analysis_result_df.rds")

#* Compare shannon / richness / CPM /  across bulk RNA-Seq cohorts - from RAW PathSeq scores without batch-effect correction ----
# meta_combined_df <- read_tsv(here("data","metadata","meta_combined_bulkRNAseq_5R16S.tsv"))
# genus_RelAbundance_corrected_mat <- readRDS(here("data/processed/genus_relAbundance_BatchEffectCorrected.rds"))
genus_counts_raw_mat <- readRDS(here("data","raw","rawScores_bulkRNAseq_genus.rds"))
dim(genus_counts_raw_mat)
meta_combined_df$Sample_ID %>% unique() %>% length()
dim(meta_combined_df)

# Select all RNA-Seq samples for the meta-analysis and calculate bacterial counts and bacterial CPM
RNAseq_IDs <- meta_combined_df %>%
  dplyr::select(Sample_ID, Dataset, Tissue_type, FastQ_libsize) %>%
  filter(!str_detect(Dataset, "16S")) %>% 
  pull(Sample_ID)

total_bact_counts_df <- genus_counts_raw_mat[,RNAseq_IDs] %>% colSums() %>% enframe(name = "Sample_ID",value = "Bacteria_total")
dim(total_bact_counts_df)

# Compute shannon and richness based on rrarefied samples
set.seed(1)
count_rar <- as.matrix(t(vegan::rrarefy(t(round(genus_counts_raw_mat[,RNAseq_IDs],0)),1000)))
# There is 1 sample with a rrarefied total number of bacteria smaller than 1000.
# This is because the sum of the PathSeq scores in this sample becomes less than 1000 when rounded to integers. 

count_rar_rel <- prop.table(count_rar,2)
div <- vegan::diversity(t(count_rar_rel), index = "shannon") %>% enframe(name = "Sample_ID",value = "Shannon")
rich <- colSums(count_rar_rel > 0) %>% enframe(name = "Sample_ID",value = "Richness")
diversity_richness_df <- full_join(div,rich,by = "Sample_ID")


# Combine with metadata 
diversity_richness_Dataset_df <- meta_combined_df %>%
  dplyr::select(Sample_ID, Dataset, FastQ_libsize) %>% 
  filter(Sample_ID %in% RNAseq_IDs) %>% 
  full_join(.,total_bact_counts_df) %>%
  full_join(.,diversity_richness_df) %>% 
  mutate(Bacteria_CPM = Bacteria_total/FastQ_libsize*1e6) %>% 
  mutate(Dataset = ifelse(str_detect(Dataset,"DKFZ"), "DKFZ RNA-Seq", Dataset))

# Do the same per "Batch"
meta_combined_df$Batch %>% table()
diversity_richness_Batch_df <- meta_combined_df %>%
  mutate(Dataset = Batch) %>%
  dplyr::select(Sample_ID, Dataset, FastQ_libsize) %>% 
  filter(Sample_ID %in% RNAseq_IDs) %>% 
  full_join(.,total_bact_counts_df) %>%
  full_join(.,diversity_richness_df) %>% 
  mutate(Bacteria_CPM = Bacteria_total/FastQ_libsize*1e6)
table(diversity_richness_Batch_df$Dataset)

# One barcelona sample without libsize #ignore for now - only for supplementary figure relevant

#* Compute alpha and beta diversity ----
f_compute_diversity <- function(diversity_richness_df, count_rar_rel) {
  # Compute wilcoxon test
  df <- diversity_richness_df %>%
    pivot_longer(cols = c(FastQ_libsize, Shannon, Richness, Bacteria_CPM, Bacteria_total), names_to = "metric", values_to = "value") %>%
    dplyr::rename(group = Dataset)
  res_df <- f_pairwise_wilcoxon_tests(df)

  median_df <- df %>%
    group_by(group, metric) %>%
    summarise(N = n(), median = median(value, na.rm = T)) %>%
    ungroup() %>%
    mutate(Group1 = group, Group2 = group, Median_Group1 = median, Median_Group2 = median, N_Group1 = N, N_Group2 = N)

  res_df <- res_df %>%
    left_join(., median_df %>% select(c(metric, contains("Group1")))) %>%
    left_join(., median_df %>% select(c(metric, contains("Group2"))))

  #* Compute Beta diversity (on Bray-Curtis and Euclidean distances)
  message("Computing beta-diversities")
  df <- diversity_richness_df %>% transmute(Sample_ID, group = Dataset)

  permanova_res_df <- f_compute_distance_metrics(df,count_rar_rel) %>% transmute(Group1, Group2, p.val_adj_permanova = p.adj, metric = paste(distance_metric, "distance", sep = " "))

  res_combined_df <-
    res_df %>%
    bind_rows(., permanova_res_df) %>%
    mutate(Metric = case_when(
      str_detect(metric, "Shannon") ~ "Shannon diversity",
      str_detect(metric, "Richness") ~ "Genus richness",
      str_detect(metric, "FastQ") ~ "Number of raw reads",
      str_detect(metric, "CPM") ~ "Bacterial CPM",
      str_detect(metric, "total") ~ "Bacterial counts",
      TRUE ~ metric
    ))
  # * Batch
  alpha_beta_diversity_clean_df <-
    res_combined_df %>%


    mutate(comparison = "02_Bulk RNA-Seq comparison by Batch") %>%
    glimpse() %>%
    transmute(
      comparison, Metric, Group1, Group2, N_Group1, N_Group2, Median_Group1, Median_Group2, p.val_adj_wilcox, p.val_adj_permanova
    )
}
dataset_diversity_res_df <- f_compute_diversity(diversity_richness_df = diversity_richness_Dataset_df,count_rar_rel = count_rar_rel)
batch_diversity_res_df <- f_compute_diversity(diversity_richness_df = diversity_richness_Batch_df,count_rar_rel = count_rar_rel)

table(dataset_diversity_res_df$Group1)
table(batch_diversity_res_df$Group2)

#save output of diversity analysis as tsv
write_tsv(dataset_diversity_res_df,here("data","results","HCC-meta-analysis_diversity_comparisons_byDataset.tsv"))  
write_tsv(batch_diversity_res_df,here("data","results","HCC-meta-analysis_diversity_comparisons_byBatch.tsv"))  

