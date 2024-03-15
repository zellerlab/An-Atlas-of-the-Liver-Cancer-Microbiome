#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 6
# by Fabian Springer
########

library(tidyverse)
library(here)
require(yaml)
library(patchwork)

source(here("src","plotting","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","ExtendedDataFigures","ExtendedDataFigure6")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import data
HCC_types <- c("HBV_HCC","HCV_HCC","ALD/ASH_HCC","MAFLD/MASH_HCC","other_HCC")

# Import the metadata and define the Groups for the individual comparisons in the figure:
# 1) Fibrosis -> HCC
# 2) HCC subtypes
meta_all_df <- read_tsv(here("data", "metadata", "meta_5R16S_allSamples.tsv")) %>%
  mutate(
    Group_FibrosisComparison = case_when( # Summarise phCCA and dCCA as one group
      str_detect(Etiology,"fibrosis") |  Etiology == "Adj. non-tumor_HCC" ~ Etiology,
      Etiology %in% HCC_types ~ "HCC",
      TRUE ~ NA_character_      
    ),
    Group_HCCsubtypes = case_when(
      str_detect(Etiology,"HBV|HCV|ALD|MAFLD") ~ Etiology,
      TRUE ~ NA_character_      
    )
  )
table(meta_all_df$Etiology,meta_all_df$Group_FibrosisComparison)
table(meta_all_df$Etiology,meta_all_df$Group_HCCsubtypes)

# Import the species profiles
species_relAbundance_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_allSamples.rds"))


#* Shannon diversities, richness and PCoAs ----
# Compute diversity
div <- vegan::diversity(t(species_relAbundance_mat),index = "shannon") %>% enframe(name = "Sample_ID",value = "ShannonDiv")
# Compute Richness
threshold_for_prev <- 0
rich <- colSums(species_relAbundance_mat>threshold_for_prev) %>% enframe(name = "Sample_ID",value = "SpeciesRichness")
diversity_richness_df <- full_join(div,rich)
sample_diversity_df <- diversity_richness_df %>%
  gather(-Sample_ID, key = "metric", value = "value") %>%
  left_join(., meta_all_df %>% dplyr::select(Sample_ID, Group_FibrosisComparison,Group_HCCsubtypes))

# IMport test result of alpha end beta diversity
alpha_beta_diversity_clean_df <-
  read_tsv(here("data", "results", "5R16S_diversity_comparisons_df.tsv"))

# Since this is now always the same, generate a helper function
table(meta_all_df$Etiology)
comparison_list <- list(
  c("Liver fibrosis early stage", "Liver fibrosis late stage", "Adj. non-tumor_HCC", "HCC"),
  c("ALD/ASH_HCC", "MAFLD/MASH_HCC", "HBV_HCC", "HCV_HCC")
)

# define how the comparisons are called in the dataframe with the alpha- and beta-diversity results
alpha_beta_diversity_clean_df$comparison %>% unique()
comparison_name <- c("Fibrosis_to_HCC","ALD_vs_MAFLD_vs_HBV_vs_HCV")

# define corresponding column in the metadata DF
meta_group_columns <- c("Group_FibrosisComparison","Group_HCCsubtypes")

f_helper_plot_diversity <- function(c_comparison, c_name, c_meta_group_column) {
  # takes a comparison vector (e.g. c("early Fibrosis","late fibrosis,..")) and a comparison name (e.g. "Fibrosis_to_HCC"), 
  # generates a dataframe with richness and shannon diversities and computes PCoA and pulls the p-values from the alpha_beta_diversity_clean_df
  # returns plots of shannon diversity, richness and PCoA

  message("Plotting ", c_name, " with comparison ", paste(c_comparison,sep = "  :"))

  tmp_df <- sample_diversity_df %>%
    dplyr::rename(Group = !!c_meta_group_column) %>% 
    filter(Group %in% c_comparison) %>%
    filter(metric == "SpeciesRichness") %>%
    mutate(
      Group = factor(Group, levels = c_comparison),
      y = value
    )

  # Get the p-values from the Wilcoxon testing of the species richness and shannon diversity
  # Get list of comparisons and the p-values:
  df <- alpha_beta_diversity_clean_df %>% filter(comparison == c_name)
  c_comparison_list <- lapply(1:nrow(df), function(i) c(df$Group1[i], df$Group2[i]))
  p_adj_vec_richness <- df$p.val_adj_wilcox_richness_species
  
  pt_richness <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Species richness", corral.width = 0.49) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors)

  # For visualization: Ceil the y-axis at 100
  if (max(pt_richness$data$y) > 100) {
    # For better visualization: Manually fix y-values > 100 to defined values in order to decrease the plotting window
    larger_y <- pt_richness$data %>% filter(y > 101) %>% pull(y) %>% sort()
    pt_richness$data <- pt_richness$data %>%
      mutate(y = case_when(y > 101 & y < 300 ~ 110, y > 300 ~ 120, TRUE ~ y))
    pt_richness <- pt_richness +
      scale_y_continuous(
        breaks = c(seq(0, 100, 25), 110,120), labels = as.character(c(seq(0, 100, 25), larger_y)),
        limits = c(0, 121)
      )
  }
  
  ### Add matrix indicating significance of pairwise comparison
  # extract pvalues and the respective comparison
  p_adj_richness_vec <- alpha_beta_diversity_clean_df %>%
    filter(comparison == c_name) %>%
    transmute(pair_comp = paste0(Group1,"_vs_",Group2), p.val_adj_wilcox_richness_species) %>% 
    deframe()
  
  # generate matrix with p-values
  tri_mat <- f_create_upper_triangle_matrix(p_adj_richness_vec, condition_levels = c_comparison)
  
  # For visualization: Update names in tri_mat wit line breaks
  rownames(tri_mat) <- str_replace_all(rownames(tri_mat), "fibrosis ","fibrosis\n")
  rownames(tri_mat) <- str_replace_all(rownames(tri_mat), "non-tumor","non-\ntumor")
  colnames(tri_mat) <- str_replace_all(colnames(tri_mat), "fibrosis ","fibrosis\n")
  colnames(tri_mat) <- str_replace_all(colnames(tri_mat), "non-tumor","non-\ntumor")

  # plot significance matrix
  pt_sig_mat <- f_plot_signif_matrix(tri_mat)
  
  # merge with richness
  pt_richness_merged <- (pt_richness + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())) /
    plot_spacer() / pt_sig_mat + plot_layout(heights = c(2, -0.25, 1), widths = c(1, 1, 1))


  #* Shannon diversity ----
  tmp_df <- sample_diversity_df %>%
    dplyr::rename(Group = !!c_meta_group_column) %>% 
    filter(Group %in% c_comparison) %>%
    filter(metric == "ShannonDiv") %>%
    mutate(
      Group = factor(Group, levels = c_comparison),
      y = value
    )  
  
  pt_shannon <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Shannon diversity", corral.width = 0.49) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors)

  # Add matrix indicating significance of pairwise comparison
  # extract pvalues and the respective comparison
  p_adj_shannon_vec <- alpha_beta_diversity_clean_df %>%
    filter(comparison == c_name) %>%
    transmute(pair_comp = paste0(Group1,"_vs_",Group2), p.val_adj_wilcox_shannon_diversity) %>% 
    deframe()
  
  # generate matrix with p-values
  tri_mat <- f_create_upper_triangle_matrix(p_adj_shannon_vec, condition_levels = c_comparison)
  
  # For visualization: Update names in tri_mat wit line breaks
  rownames(tri_mat) <- str_replace_all(rownames(tri_mat), "fibrosis ","fibrosis\n")
  rownames(tri_mat) <- str_replace_all(rownames(tri_mat), "non-tumor","non-\ntumor")
  colnames(tri_mat) <- str_replace_all(colnames(tri_mat), "fibrosis ","fibrosis\n")
  colnames(tri_mat) <- str_replace_all(colnames(tri_mat), "non-tumor","non-\ntumor")

  # plot significance matrix
  pt_sig_mat <- f_plot_signif_matrix(tri_mat)
  
  # merge with richness
  pt_shannon_merged <- (pt_shannon + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())) /
    plot_spacer() / pt_sig_mat + plot_layout(heights = c(2, -0.25, 1), widths = c(1, 1, 1))


  #* PCoA ----
  meta <- tmp_df %>% transmute(Sample_ID,Group)
  
  pt_pcoa <- f_plot_PCoA(
    meta_df = meta,
    mat = species_relAbundance_mat[,meta$Sample_ID], method = "bray"
  ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    theme(legend.position = "none")+
    theme(aspect.ratio = 3/4) # fix aspect ratio before merging with significance matrix 
  
  # Also add the triangular matrix with significance levels
  p_adj_bray_vec <- alpha_beta_diversity_clean_df %>%
    filter(comparison == c_name) %>%
    transmute(pair_comp = paste0(Group1,"_vs_",Group2), p.val_adj_permanova_Bray) %>% 
    deframe()
  
  # generate matrix with p-values
  tri_mat <- f_create_upper_triangle_matrix(p_adj_bray_vec, condition_levels = c_comparison)
  
  # For visualization: Update names in tri_mat wit line breaks
  rownames(tri_mat) <- str_replace_all(rownames(tri_mat), "fibrosis ","fibrosis\n")
  rownames(tri_mat) <- str_replace_all(rownames(tri_mat), "non-tumor","non-\ntumor")
  colnames(tri_mat) <- str_replace_all(colnames(tri_mat), "fibrosis ","fibrosis\n")
  colnames(tri_mat) <- str_replace_all(colnames(tri_mat), "non-tumor","non-\ntumor")

  # plot significance matrix
  pt_sig_mat <- f_plot_signif_matrix(tri_mat)

  pt_pcoa_merged <- pt_pcoa / plot_spacer() / pt_sig_mat +
    plot_layout(heights = c(2, -0.4, 1), widths = c(1, 1, 1))


  return(list(richness = pt_richness_merged, shannon = pt_shannon_merged, pcoa = pt_pcoa_merged))
}
comparison_list
figure_indices <- c("A","B","C","G","H","I")
meta_group_columns
i <- 1
w <- 4.5
h <- 4.5
id <- 1
for(i in seq(1,length(comparison_list))){
  # loop through the comparisons and generate the plots
  res <- f_helper_plot_diversity(
    c_comparison = comparison_list[[i]],
    c_name = comparison_name[i],
    c_meta_group_column = meta_group_columns[i]
  )  
  # save the plots
  ggsave(paste0(save_fig_folder,"/",figure_indices[id],"_SpeciesRichness.pdf"), res$richness, width = w, height = h)
  id <- id + 1
  ggsave(paste0(save_fig_folder,"/",figure_indices[id],"_ShannonDiv.pdf"), res$shannon, width = w, height = h)
  id <- id + 1
  ggsave(paste0(save_fig_folder,"/",figure_indices[id],"_PCoA.pdf"), res$pcoa, width = w, height = h)
  id <- id + 1
}

#* Volcano plots of pairwise comparisons ----
# load test result
all_test_results_df <- readRDS(here("data","results","5R16S_EtiologyTesting_result_df.rds"))
all_test_results_df$comparison %>% table()
size_definition <- scale_size_continuous(
  name = "Prevalence",
  range = c(1, 4),
  limits = c(0.01, 0.5),
  breaks = c(0.05, 0.1, 0.2, 0.4)  
)

### Starti with Fibrosis -> HCC comparison
man_y_breaks <- -log10(c(0.05,0.01,0.001,0.0001,1e-5,1e-6))

all_test_results_df$comparison %>% unique()

# Early fibrosis vs HCC
plot_df <- all_test_results_df %>% filter(comparison == "Fibrosis_early_vs_HCC")
plot_df$Group1 <- str_replace(plot_df$Group1,"EarlyFib","Fibrosis early stage")
plot_df$Group2 <- str_replace(plot_df$Group2,"EarlyFib","Fibrosis early stage")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.4,0.3,0.1),1)
xLims <- range(xBreaks)
pt_D <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors,drop = F)+
  size_definition

# Late fibrosis vs HCC
plot_df <- all_test_results_df %>% filter(comparison == "Fibrosis_late_vs_HCC")
plot_df$Group1 <- str_replace(plot_df$Group1,"LateFib","Fibrosis late stage")
plot_df$Group2 <- str_replace(plot_df$Group2,"LateFib","Fibrosis late stage")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.8,0.4,0.2),1)
xLims <- range(xBreaks)
pt_E <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors,drop = F)+
  size_definition

# Adjacent non-tumor vs HCC
plot_df <- all_test_results_df %>% filter(comparison == "Tumor_vs_Adj_HCC")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.3,0.3,0.1),1)
xLims <- range(-0.33,0.33)
pt_F <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors,drop = F)+
  size_definition


### Compare HCC etiologies
man_y_breaks <- -log10(c(0.05,0.01,0.001,1e-4))
all_test_results_df$comparison %>% table()

# Viral vs Non-Viral
plot_df <- all_test_results_df %>% filter(comparison == "Viral_vs_nonViral_HCC")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.4,0.2,0.2),1)
xLims <- range(-0.42,0.23)
pt_J <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors,drop = F)+
  size_definition

# ALD vs MAFLD #! Direction swapped compared to manuscript
plot_df <- all_test_results_df %>% filter(comparison == "ALD_vs_MAFLD")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.4,0.4,0.2),1)
xLims <- range(-0.44,0.44)
pt_K <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors,drop = F)+
  size_definition

# HBV vs HCV #! Direction swapped compared to manuscript
plot_df <- all_test_results_df %>% filter(comparison == "HBV_vs_HCV")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.3,0.3,0.15),2)
xLims <- range(-0.34,0.34)
pt_L <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors,drop = F)+
  size_definition

# Save volcano plots
width_volcano <- 4.5
height_volcano <- 5
ggsave(paste0(save_fig_folder,"/D_Volcano_earlyFib_vs_HCC.pdf"), pt_D, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/E_Volcano_lateFib_vs_HCC.pdf"), pt_E, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/F_Volcano_HCC_vs_AdjNonTumor.pdf"), pt_F, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/J_Volcano_Viral_vs_NonViral.pdf"), pt_J, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/K_Volcano_ALD_vs_MAFLD.pdf"), pt_K, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/L_Volcano_HBV_vs_HCV.pdf"), pt_L, width = width_volcano, height = height_volcano)

