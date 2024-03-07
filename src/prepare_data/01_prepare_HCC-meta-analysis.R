#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Script 01: Import microbial profiles from liver tissue bulk RNA-Seq datasets and combine with
#   genus level profiles from the 5R-16S rDNA-Seq cohort to perform a meta-analysis of the liver cancer tissue microbiome.
# by Fabian Springer
########

library(here)
library(tidyverse)

#* Import microbial counts and metadata  ----
meta_bulkRNAseq_df <- read_tsv(here("data","metadata","meta_bulkRNAseq.tsv"))
meta_5R16S_df <- read_tsv(here("data","metadata","meta_5R16S_HCC-meta-analysis.tsv"))

genus_counts_bulkRNAseq_mat <- readRDS(here("data","raw","rawScores_bulkRNAseq_genus.rds"))
genus_relAbundance_5R16S_mat <- readRDS(here("data","raw","relAbundance_5R16S_genus_HCC-meta-analysis.rds"))

genus_counts_bulkRNAseq_mat[1:5,1:5]

#* Prepare profiles from bulk RNA-Seq data ----
# Filter samples with <1000 total bacterial counts

bacterial_counts_bulkRNAseq_df <- colSums(genus_counts_bulkRNAseq_mat) %>% enframe(name = "Sample_ID",value = "Bacteria_total")
# Total bacterial counts are not integers since PathSeq defines scores based also on multi-mappers and therefore distributes counts across several taxa 
# (which is more accurate than only assigning integer counts to unique mappers)
Sample_IDs_passed_bacterial_threshold <- bacterial_counts_bulkRNAseq_df %>% filter(Bacteria_total >= 1000) %>% pull(Sample_ID)
#bacterial_counts_bulkRNAseq_df %>% filter(Sample_ID %in% Sample_IDs_passed_bacterial_threshold) %>% pull(Bacteria_total) %>% range()

# Compute relative abundances
genus_relAbundance_bulkRNAseq_mat <- prop.table(genus_counts_bulkRNAseq_mat[,Sample_IDs_passed_bacterial_threshold],2) # Divide by colSums

# Filter genera based on prevalence criteria in individual cohorts
threshold_for_prevalence <- 1e-3 #a genus is considered "prevalent" in a RNA-seq sample if its relative abundance is > 1e-3
prevalence_cutoff <- 0.1 #we consider only genera with at least 10% prevalence per dataset (ignoring the "Tumor" or "Adj. non-tumor label")

genus_prevalence_df <-
  genus_relAbundance_bulkRNAseq_mat %>%
  as_tibble(rownames = "genus") %>%
  gather(-genus, key = "Sample_ID", value = "rel") %>%
  left_join(., meta_bulkRNAseq_df %>% dplyr::select(Sample_ID, Dataset)) %>%
  group_by(Dataset) %>%
  mutate(N_samples = length(unique(Sample_ID))) %>%
  ungroup() %>%
  group_by(Dataset, genus) %>%
  mutate(prevalence = sum(rel > threshold_for_prevalence) / N_samples) %>%
  ungroup() %>% 
  dplyr::select(genus,Dataset,prevalence) %>% 
  distinct()

# Next: Filter out genera with prevalence < cutoff and keep only genera that fulfull the prevalence cutoff in all four cohorts
genera_passed_threshold_bulkRNAseq <- genus_prevalence_df %>% 
  filter(prevalence > prevalence_cutoff) %>%   
  group_by(genus) %>% 
  summarise(Passed_cutoff_in_N_cohorts = n()) %>% 
  ungroup() %>% 
  filter(Passed_cutoff_in_N_cohorts == length(unique(genus_prevalence_df$Dataset))) %>% 
  pull(genus)

#* Filter genera based on prevalence in the 5R-16S Dataset  ----

# The 5R 16S data was profiled, decontaminated and converted to relative abundances as described in:
# Nejman, D. et al. The human tumor microbiome is composed of tumor type–specific intracellular bacteria. Science 368, 973–980 (2020).
prevalence_cutoff_16S <- 0.01 #The prevalence cutoff is only 1% since contamination removal (on species level) and flooring of rel. abundances is performed during preprocessing (see Nejman et al.)

genus_prev_5R16S <- rowSums(genus_relAbundance_5R16S_mat>0) / ncol(genus_relAbundance_5R16S_mat)
genera_passed_threshold_5R16S <- names(genus_prev_5R16S[genus_prev_5R16S > prevalence_cutoff_16S])

#* Define final testset and combine bulkRNA-seq profiles with 16S profiles ----
# 63 genera pass these thresholds
genera_kept_for_analysis <- intersect(genera_passed_threshold_bulkRNAseq,genera_passed_threshold_5R16S)
length(genera_kept_for_analysis)
# Generate counts matrix combining the relative abundances of the profiles from bulkRNA-seq and 5R16S of the selected genera
genus_RelAbundance_mat <- cbind(
  genus_relAbundance_bulkRNAseq_mat[genera_kept_for_analysis, ],
  genus_relAbundance_5R16S_mat[genera_kept_for_analysis, ]
) %>% 
as.matrix()

# Combine metadata of bulk-RNAseq samples and 5R16S samples
meta_combined_df <- bind_rows(
  meta_bulkRNAseq_df %>%
    filter(Sample_ID %in% Sample_IDs_passed_bacterial_threshold),
  meta_5R16S_df
) %>%
  mutate(Etiology = case_when(
    str_detect(Etiology, "Adj. non-tumor") ~ "Adj. non-tumor",    
    TRUE ~ Etiology
  ))

#* Next step: Take filtered relative abundance matrix and perform batchwise linear transformation ----
### Install the BEATLE package from this GitLab Repository: https://git.embl.de/grp-zeller/batch_effect_correction
library(BEATLE)

# We use an R package

# Define TCGA as reference dataset
reference_data_mat <- genus_RelAbundance_mat[,meta_combined_df %>% filter(Dataset == "TCGA") %>% pull(Sample_ID)]


# Function to apply bacth effect correction using the BEATLE package
f_correct_batch <- function(target_data_mat,reference_data_mat){
  diagnostic_plots_dir <- file.path("./data/processed/batch_effect_correction_diagnostic_plots/")
  if(!dir.exists(diagnostic_plots_dir)){dir.create(diagnostic_plots_dir,recursive = T)}
  
  res <- correct.batch(target.data = target_data_mat,  # Target dataset to be transformed
                    reference.data = reference_data_mat,        # Reference dataset
                    log.n0=1e-08,          # Pseudo-count to be added
                    ref.ctr.idx=NULL,      # Indices of control instances in the reference dataset
                    trg.ctr.idx=NULL,      # Indices of control instances in the target dataset
                    diagnostic.plot=file.path(diagnostic_plots_dir,paste0(dset,"_vs_TCGA.pdf")),  # Filename of diagnostic plot
                    verbose=1)
  
  return(res)
}

# Loop over all Datasets that should be corrected (all except TCGA) and apply the batchwise linear transformation
datasets_to_transform <- meta_combined_df %>% filter(!str_detect(Dataset,"TCGA")) %>% pull(Dataset) %>% unique()
dset <- datasets_to_transform[1]
genus_RelAbundance_corrected_mat <- matrix(nrow = nrow(genus_RelAbundance_mat),ncol = 0)
rownames(genus_RelAbundance_corrected_mat) <- rownames(genus_RelAbundance_mat)
for(dset in datasets_to_transform){
  res <- f_correct_batch(
    target_data_mat = genus_RelAbundance_mat[, meta_combined_df %>%
                                                  filter(Dataset == dset) %>%
                                                  pull(Sample_ID)],
    reference_data_mat = reference_data_mat
  )
  #Join corrected matrix of current batch to the main matrix
  genus_RelAbundance_corrected_mat <- cbind(genus_RelAbundance_corrected_mat,res$transformed.data)
}
# Also add the TCGA samples
genus_RelAbundance_corrected_mat <- cbind(genus_RelAbundance_corrected_mat,res$reference.data)
# Revert the log transformation performed during batch effect correction
genus_RelAbundance_corrected_mat <- (10^genus_RelAbundance_corrected_mat)-1e-8

# For consistence: Order columns in transformed and untransformed matrices equally
genus_RelAbundance_corrected_mat <- genus_RelAbundance_corrected_mat[,colnames(genus_RelAbundance_mat)]

# Export transformed and untransformed relative abundance matrix4
write_tsv(meta_combined_df,file = here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
saveRDS(genus_RelAbundance_mat,file = here("data/processed/relAbundance_combined_genus_HCC-meta-analysis_raw.rds"))
saveRDS(genus_RelAbundance_corrected_mat,file = here("data/processed/relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))
