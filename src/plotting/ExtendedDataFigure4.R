#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 4
# by Fabian Springer
########

library(tidyverse)
library(here)
require(yaml)
library(GGally)
library(patchwork)

source(here("src","plotting","functions_plotting.R"))

# Define folder to store plots
save_fig_folder <- here("figures","ExtendedDataFigures","ExtendedDataFigure4")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

#* Import data ----
# Import 5R16S profiles
meta_16S_all_samples_df <- read_tsv(here("data","metadata","meta_5R16S_allSamples.tsv"))
relAB_phylum_5R16S_mat <- readRDS(here("data","raw","relAbundance_5R16S_phylum_allSamples.rds"))
# renormalize phylum relative abundances in the 5R 16S cohort for better visualization
relAB_phylum_5R16S_mat <- prop.table(relAB_phylum_5R16S_mat,2)

# Import bulk RNAseq profiles
meta_bulkRNAseq_all_samples_df <- read_tsv(here("data","metadata","meta_bulkRNAseq.tsv"))  
relAB_phylum_bulkRNAseq_mat <- readRDS(here("data","raw","relAbundance_DKFZ_bulkRNAseq_phylum.rds"))

# generate meta_df with only shared samples
meta_shared_df <- inner_join(
  meta_bulkRNAseq_all_samples_df %>%
  filter(Tissue_type == "Tumor") %>%
    transmute(Sample_ID_RNAseq = Sample_ID, Patient_ID),
  meta_16S_all_samples_df %>%
    transmute(Sample_ID_5R16S = Sample_ID, Patient_ID)
)

#* Define phylum colors ----
phylum_colors_vec <- c("#009F4D", "#E40046","#307FE2", "#FFA300", "#8246AF", "#FFCD00", "#A8C700", "#BE5400", "#A8A99E")
names <- c(
  "Proteobacteria", "Actinobacteria", "Firmicutes", "Bacteroidetes",
  "Cyanobacteria", "Fusobacteria", "Acidobacteria", "Verrucomicrobia", "other"
)
names(phylum_colors_vec) <- names

#* Panel A: Sample by sample comparison of bulk RNA-seq derived microbial profiles with 5R16S-Seq in the DKFZ cohort ----
# take all phyla from 16S and RNA-Seq to build consensus on "what is matchable"
taxa_16S <- str_remove_all(rownames(relAB_phylum_5R16S_mat),".*;")
taxa_16S <- str_remove_all(taxa_16S,"[a-z]__")
rownames(relAB_phylum_5R16S_mat) <- taxa_16S
shared_taxa <- intersect(taxa_16S,rownames(relAB_phylum_bulkRNAseq_mat))
message("Shared phyla: ",length(shared_taxa))

# Since not all phyla are shared, check the fraction of relative abundance explained by shared phyla
relAB_5R16S_mat <- relAB_phylum_5R16S_mat
relAB_RNAseq_mat <- relAB_phylum_bulkRNAseq_mat

fract_relAB_byShared_taxa_df <- left_join(
  # Compute fraction of rel. Abundance described by shared Taxa in 16S 
  left_join(
    colSums(relAB_5R16S_mat[, meta_shared_df$Sample_ID_5R16S]) %>%
      enframe(name = "Sample_ID", value = "sum_relTaxa"),
    colSums(relAB_5R16S_mat[shared_taxa, meta_shared_df$Sample_ID_5R16S]) %>%
      enframe(name = "Sample_ID", value = "sum_relTaxa_shared")
  ) %>%
    mutate(fract_abundace_shared = sum_relTaxa_shared / sum_relTaxa) %>%
    left_join(., meta_shared_df %>% transmute(Sample_ID = Sample_ID_5R16S, Patient_ID)) %>%
    transmute(Patient_ID, fract_relAB_5R16S = fract_abundace_shared),
  
  # Compute fraction of rel. Abundance explained by shared Taxa in RNA-Seq
  left_join(
    colSums(relAB_RNAseq_mat[, meta_shared_df$Sample_ID_RNAseq]) %>%
      enframe(name = "Sample_ID", value = "sum_relTaxa"),
    colSums(relAB_RNAseq_mat[shared_taxa, meta_shared_df$Sample_ID_RNAseq]) %>%
      enframe(name = "Sample_ID", value = "sum_relTaxa_shared")
  ) %>%
    mutate(fract_abundace_shared = sum_relTaxa_shared / sum_relTaxa) %>%
    left_join(., meta_shared_df %>% transmute(Sample_ID = Sample_ID_RNAseq, Patient_ID)) %>%
    transmute(Patient_ID, fract_relAB_RNAseq = fract_abundace_shared)
)
# At least 78% of rel abudnance are epxlained by shared phyla
range(fract_relAB_byShared_taxa_df$fract_relAB_5R16S)
range(fract_relAB_byShared_taxa_df$fract_relAB_RNAseq)

# Its valid to consider only shared taxa in order to compare homogeneity.
# Combine relative abundance matrices

relAb_shared_df <-
  bind_rows(
    relAB_5R16S_mat[shared_taxa, ] %>%
      as_tibble(rownames = "tax") %>%
      gather(-tax, key = "Sample_ID_5R16S", value = "relAB") %>%
      inner_join(., meta_shared_df) %>%
      dplyr::select(-Sample_ID_5R16S, -Sample_ID_RNAseq) %>%
      mutate(Assay = "5R16S"),
    relAB_RNAseq_mat[shared_taxa, ] %>%
      as_tibble(rownames = "tax") %>%
      gather(-tax, key = "Sample_ID_RNAseq", value = "relAB") %>%
      inner_join(., meta_shared_df) %>%
      dplyr::select(-Sample_ID_5R16S, -Sample_ID_RNAseq) %>%
      mutate(Assay = "RNAseq")
  ) %>%  
  pivot_wider(names_from = Assay, values_from = relAB) %>% 
  group_by(tax) %>% 
  # Compute prevalence
  mutate(prev_5R16S = sum(`5R16S` > 0)/n(),
         prev_RNAseq = sum(RNAseq > 0)/n()) %>% 
  ungroup() %>% 
  mutate(
    l10_5R16S = log10(`5R16S` + 1e-4),
    l10_RNAseq = log10(RNAseq + 1e-4)
  )

# Filter only phyla with a prevalence > 5% in both cohorts
prev_threshold <- 0.05
relAb_shared_filtered_df <- relAb_shared_df %>%
  filter(prev_5R16S > prev_threshold, prev_RNAseq > prev_threshold) %>% 
  mutate(phylum = tax)

# Order phyla based on mean relative abundance
phylum_levels <- relAb_shared_filtered_df %>%
  group_by(tax) %>% 
  summarise(mean_rel = mean(`5R16S`)) %>% 
  dplyr::select(tax, mean_rel) %>%
  distinct() %>%
  arrange(-mean_rel) %>%
  pull(tax)

#* Generate the barplot plot. Basic idea: Generate two barplots (for 5R16S and RNAseq) separately and arrange them in a grid next to each other ----
xBreaks <- seq(0, 1, 0.25)
xLabels <- 100 * xBreaks

# Arrange the Patients based on "Proteobacteria" rel abundance in both cohorts
Patient_ID_levels <- relAb_shared_df %>%
  filter(tax == "Proteobacteria") %>%
  transmute(Patient_ID,rel16S = `5R16S`,relRNAseq = RNAseq) %>% 
  mutate(
    rank_16S = rank(-rel16S),
    rank_RNAseq = rank(-relRNAseq)) %>% 
  mutate(mean_rank = (rank_16S+rank_RNAseq) / 2) %>% 
  arrange(-mean_rank) %>% 
  pull(Patient_ID)

### left plot: 5R16S
bar_16S_df <-
  relAb_shared_df %>%
  transmute(tax, Patient_ID, rel = `5R16S`) %>%
  mutate(phylum = ifelse(tax %in% names(phylum_colors_vec), tax, "other")) %>%
  group_by(phylum, Patient_ID) %>%
  summarise(rel = sum(rel)) %>%
  ungroup()

# Since we consider only shared genera, the sum of relative abudnance explained is not 1. Add "other" up to one to account for this
bar_16S_df <- bind_rows(
  bar_16S_df,
  relAb_shared_df %>%
    transmute(tax, Patient_ID, rel = `5R16S`) %>%
    group_by(Patient_ID) %>%
    mutate(sum_rel_explained = sum(rel)) %>%
    dplyr::select(Patient_ID, sum_rel_explained) %>%
    distinct() %>%
    ungroup() %>% 
    transmute(rel = 1 - sum_rel_explained, phylum = "other",Patient_ID)
) %>% 
mutate(phylum = factor(phylum, levels = rev(c(phylum_levels, "other"))))

pt_16S <- bar_16S_df %>% 
  mutate(Patient_ID = factor(Patient_ID, levels = Patient_ID_levels)) %>%
  ggplot(aes(x = rel, y = Patient_ID, fill = phylum)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  scale_fill_manual(values = phylum_colors_vec) +
  # invert x axis
  scale_x_reverse(expand = c(0, 0), breaks = xBreaks, labels = xLabels,name = "Abundance [%]")+
  theme_paper+
  ggtitle("DKFZ 5R 16S (N=76)")+
  scale_y_discrete(position = "right")+
  guides(fill = guide_legend(reverse = TRUE))

### right plot: Bulk RNAseq
bar_RNAseq_df <-
  relAb_shared_df %>%
  transmute(tax, Patient_ID, rel = RNAseq) %>%
  mutate(phylum = ifelse(tax %in% names(phylum_colors_vec), tax, "other")) %>%
  group_by(phylum, Patient_ID) %>%
  summarise(rel = sum(rel)) %>%
  ungroup()

# Since we consider only shared genera, the sum of relative abudnance explained is not 1. Add "other" up to one to account for this
bar_RNAseq_df <- bind_rows(
  bar_RNAseq_df,
  relAb_shared_df %>%
    transmute(tax, Patient_ID, rel = RNAseq) %>%
    group_by(Patient_ID) %>%
    mutate(sum_rel_explained = sum(rel)) %>%
    dplyr::select(Patient_ID, sum_rel_explained) %>%
    distinct() %>%
    ungroup() %>% 
    transmute(rel = 1 - sum_rel_explained, phylum = "other",Patient_ID)
) %>% 
mutate(phylum = factor(phylum, levels = rev(c(phylum_levels, "other"))))

pt_RNAseq <- bar_RNAseq_df %>%
  mutate(Patient_ID = factor(Patient_ID, levels = Patient_ID_levels)) %>%
  ggplot(aes(x = rel, y = Patient_ID, fill = phylum)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  scale_fill_manual(values = phylum_colors_vec) +
  # invert x axis
  scale_x_continuous(expand = c(0, 0), breaks = xBreaks, labels = xLabels, name = "Abundance [%]") +
  theme_paper+
  ggtitle("DKFZ RNA-Seq (N=76)")+
  guides(fill = guide_legend(reverse = TRUE))


### Combine the two barplots, remove y-labels and y-ticks
pt_bar_combined <- pt_16S + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  pt_RNAseq + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(pt_bar_combined, filename = file.path(save_fig_folder,"A_Barplot_5R16SvsBulkRNAseq_DKFZ.pdf"), width = 7, height = 10)

#* Panel B: Scatterplot of relative abundances of matched samples at phylum level ----
# Compute pearson and spearman correlation
pearson_R <- cor(relAb_shared_filtered_df$l10_5R16S,relAb_shared_filtered_df$l10_RNAseq,method = "pearson")
spearman_R <- cor(relAb_shared_filtered_df$l10_5R16S,relAb_shared_filtered_df$l10_RNAseq,method = "spearman")

pt_scatter <- relAb_shared_filtered_df %>% 
  mutate(phylum = factor(phylum,levels = phylum_levels)) %>% 
  ggplot(aes(y = l10_5R16S, x = l10_RNAseq)) +
  geom_point(pch = 21,size = 2,alpha = 0.75,aes(fill = phylum)) +
  geom_abline(slope = 1, intercept = 0) +
  tune::coord_obs_pred()+
  theme_paper+
  scale_fill_manual(values = phylum_colors_vec)+
  annotate("text", x = -4, y = 0, hjust = 0, vjust = 0.75,size = 3,  
           label = paste("Pearson =", signif(pearson_R, 2), "\nSpearman =", signif(spearman_R, 2)))+
  ylab("Rel. abundance 5R16S-Seq [log10]")+
  xlab("Rel. abundance RNA-Seq [log10]")+
  labs(fill = "Phylum")+
  ggtitle("Phylum level relative abundances\nin matched DKFZ HCC samples (N=76)")

ggsave(pt_scatter ,filename = file.path(save_fig_folder,"B_Scatter_5R16SvsBulkRNAseq_DKFZ.pdf"), width = 5, height = 5)

#* Panel C: Boxplots of matched samples grouped by phylum ----
yBreaks <- seq(-4, 0, 1)
# Generate boxplot for 5R16S
box_5R16S_df <-
  relAb_shared_filtered_df %>%
  mutate(
    tax = factor(tax, levels = phylum_levels), Group = tax, y = l10_5R16S
  )
pt_box_16S <-
  f_boxplot_by_group(box_5R16S_df, xlab = "", ylab = "Rel. abundance [log10]") +
  scale_fill_manual(values = phylum_colors_vec) +
  scale_color_manual(values = phylum_colors_vec) +
  theme(legend.position = "none", axis.title.x = element_blank())+
  ggtitle("DKFZ 16S-Seq (N=76)")+
  scale_y_continuous(breaks = yBreaks, labels = yBreaks,limits = c(-4,0))

# Generate boxplot for RNA-Seq
box_RNAseq_df <-
  relAb_shared_filtered_df %>%
  mutate(
    tax = factor(tax, levels = phylum_levels), Group = tax, y = l10_RNAseq
  )
pt_box_RNAseq <-
  f_boxplot_by_group(box_RNAseq_df, xlab = "", ylab = "Rel. abundance [log10]") +
  scale_fill_manual(values = phylum_colors_vec) +
  scale_color_manual(values = phylum_colors_vec) +
  theme(legend.position = "none", axis.title.x = element_blank())+
  ggtitle("DKFZ RNA-Seq (N=76)")+
  scale_y_continuous(breaks = yBreaks, labels = yBreaks,limits = c(-4,0))

# Combine the two boxplots
pt_box_combined <- (pt_box_16S+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())) / pt_box_RNAseq
ggsave(pt_box_combined ,filename = file.path(save_fig_folder,"C_Boxplot_5R16SvsBulkRNAseq_DKFZ.pdf"), width = 12, height = 8)


#* Panel D: Pairsplot of genus level relative abundances for all datasets against each other ----
# Import the batch-effect corrected genus level profiles and metadata of all samples used in the HCC meta-analysis 
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))

# Deinfe the threshold for a genus to be considered prevalent
threshold_for_prev_RNAseq <- 1e-3 #relative abundance in the RNA-Seq cohort for which a genus is considered prevalent
threshold_for_prev_5R16S <- 0 #same in the 5R16S cohort - In the preprocessing small values (<1e-4) were floored, so the threshold can be 0

# Get prevalence and abundance by Dataset (tumor samples only)
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
  arrange(genus) %>% 
  mutate(mean_abundance = log10(mean_abundance+1e-5)) # log transform mean abundances

# Pairsplot of Prevalences
pairs_prev_df <-
  tumor_abundance_prevalence_df %>%
  mutate(Dataset_lab = paste0(Dataset," (N=",N_samples,")")) %>% 
  dplyr::select(genus, Dataset_lab, prevalence) %>%
  pivot_wider(names_from = Dataset_lab, values_from = prevalence)

pt_pair_prev <- ggpairs(pairs_prev_df, columns = 2:ncol(pairs_prev_df)) +
  theme_paper +
  xlab("Prevalence") +
  ylab("Prevalence") +
  ggtitle("Comparison of genus prevalences between datasets")+  
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

# Pairsplot of mean abundances
pairs_abundance_df <-
  tumor_abundance_prevalence_df %>%
  mutate(Dataset_lab = paste0(Dataset," (N=",N_samples,")")) %>% 
  dplyr::select(genus, Dataset_lab, mean_abundance) %>%
  pivot_wider(names_from = Dataset_lab, values_from = mean_abundance)

pt_pair_mean <- ggpairs(pairs_abundance_df, columns = 2:ncol(pairs_prev_df)) +
  theme_paper +
  xlab("Mean abundance [log10]") +
  ylab("Mean abundance [log10]") +
  ggtitle("Comparison of genus abundances between datasets")+
  #scale_y_log10() +
  #scale_x_log10()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

# save plots
ggsave(pt_pair_mean,filename = file.path(save_fig_folder,"D_Pairsplot_MeanGenusAbundances.pdf"),width = 10,height = 10)
ggsave(pt_pair_prev,filename = file.path(save_fig_folder,"E_Pairsplot_GenusPrevalences.pdf"),width = 10,height = 10)



