#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 2
# by Fabian Springer
########

library(tidyverse)
library(here)
require(yaml)

source(here("src","plotting","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","ExtendedDataFigures","ExtendedDataFigure2")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import data
meta_all_df <- read_tsv(here("data", "metadata", "meta_5R16S_allSamples.tsv")) %>%
  mutate(
    Group = case_when( # Summarise phCCA and dCCA as one group
      Etiology == "phCCA" ~ "phCCA/dCCA",
      Etiology == "dCCA" ~ "phCCA/dCCA",
      Etiology == "Adj. non-tumor_CCC" ~ "Adj. non-tumor_iCCA",
      TRUE ~ Etiology
    )
  )
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
  left_join(., meta_all_df %>% dplyr::select(Sample_ID, Group))

alpha_beta_diversity_clean_df <-
  read_tsv(here("data", "results", "5R16S_diversity_comparisons_df.tsv")) %>%
  mutate(    
    Group1 = case_when(Group1 == "phCCA_dCCA" ~ "phCCA/dCCA", Group1 == "Adj. non-tumor_CCC" ~ "Adj. non-tumor_iCCA", TRUE ~ Group1),
    Group2 = case_when(Group2 == "phCCA_dCCA" ~ "phCCA/dCCA", Group2 == "Adj. non-tumor_CCC" ~ "Adj. non-tumor_iCCA", TRUE ~ Group2)
  )

# Since this is now always the same, generate a helper function
comparison_list <- list(
  c("Adj. non-tumor_iCCA","iCCA"),
  c("Adj. non-tumor_CRLM","CRLM"),
  c("iCCA","GBC",'phCCA/dCCA')
  )

# define how the comparisons are called in the dataframe with the testing result
comparison_name <- c("Tumor_vs_Adj_iCCA","Tumor_vs_Adj_CRLM","iCCA_GBC_vs_phCCA_dCCA")


i <- 3
c_comparison <- comparison_list[[i]]
c_name <- comparison_name[i]

f_helper_plot_diversity <- function(c_comparison, c_name) {
  # takes a comparison vector (e.g. c(Adj. non-tumor iCCA, iCCA)) and a comparison name (e.g. "Tumor_vs_Adj_iCCA"), 
  # generates a dataframe with richness and shannon diversities and computes PCoA and pulls the p-values from the alpha_beta_diversity_clean_df
  # returns plots of shannon diversity, richness and PCoA

  message("Plotting ", c_name, " with comparison ", paste(c_comparison,sep = "  :"))

  tmp_df <- sample_diversity_df %>%
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
    scale_color_manual(values = group_colors) +
    # Just use the ggsignif function for plotting the p-value symbols
    ggsignif::geom_signif(
      comparisons = c_comparison_list,
      annotations = p_to_symbol(p_adj_vec_richness),
      vjust = -0, step_increase = 0.075, textsize = 3
    )

  # Shannon diversity
  tmp_df <- sample_diversity_df %>%
    filter(Group %in% c_comparison) %>%
    filter(metric == "ShannonDiv") %>%
    mutate(
      Group = factor(Group, levels = c_comparison),
      y = value
    )  
  p_adj_vec_shannon <- df$p.val_adj_wilcox_shannon_diversity

  pt_shannon <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Shannon diversity", corral.width = 0.49) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    # Just use the ggsignif function for plotting the p-value symbols
    ggsignif::geom_signif(
      comparisons = c_comparison_list,
      annotations = p_to_symbol(p_adj_vec_shannon),
      vjust = -0, step_increase = 0.075, textsize = 3
    )

  # PCoA  
  meta <- tmp_df %>% transmute(Sample_ID,Group)
  # get the corrected p-value for all comparisons from the PERMANOVA testing of bray-curtis dissimilarities
  p_text <- paste(df %>% mutate(p_text =  paste0(Group1," vs ",Group2,": ",signif(p.val_adj_permanova_Bray,1))) %>% pull(p_text),collapse = "\n")  
  
  pt_pcoa <- f_plot_PCoA(
    meta_df = meta,
    mat = species_relAbundance_mat[,meta$Sample_ID], method = "bray"
  ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    theme(legend.position = "none")+
    # add pvalue to bottom left
    annotate("text", x = -Inf, y = -Inf, label = p_text, size = 3, color = "black", hjust = -0.05, vjust = -0.5)
  

  return(list(richness = pt_richness, shannon = pt_shannon, pcoa = pt_pcoa))
}

i <- 1
comparison_list
figure_indices <- c("A","B","C","E","F","G","I","J","K")
id <- 1
w <- 3.5
h <- 4.5
for(i in seq(1,length(comparison_list))){
  # loop through the comparisons and generate the plots
  res <- f_helper_plot_diversity(c_comparison = comparison_list[[i]], c_name = comparison_name[i])

  # save the plots
  ggsave(paste0(save_fig_folder,"/",figure_indices[id],"_SpeciesRichness.pdf"), res$richness, width = w, height = h)
  id <- id + 1
  ggsave(paste0(save_fig_folder,"/",figure_indices[id],"_ShannonDiv.pdf"), res$shannon, width = w, height = h)
  id <- id + 1
  ggsave(paste0(save_fig_folder,"/",figure_indices[id],"_PCoA.pdf"), res$pcoa, width = 6, height = 5)
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

### iCCA and CRLM vs adjacent
man_y_breaks <- -log10(c(0.05,0.01,0.001))

# iCCA vs non-tumor iCCA
plot_df <- all_test_results_df %>% filter(comparison == "Tumor_vs_Adj_iCCA") %>% mutate(Group1 = str_replace(Group1,"CCC","iCCA"))
range(plot_df$effect.size)
xBreaks <- round(seq(-0.4,0.4,0.2),1)
xLims <- range(xBreaks)
pt_D <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors, drop = F)+
  size_definition

# CRLM vs non-TUmor CRLM
plot_df <- all_test_results_df %>% filter(comparison == "Tumor_vs_Adj_CRLM")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.5,0.5,0.25),1)
xLims <- c(-0.5,0.52)
pt_H <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors, drop = F)+
  size_definition


# iCCA vs GBC and phCCA/dCCA
man_y_breaks <- -log10(c(0.05,0.01,1e-4,1e-6,1e-8,1e-10))

# iCCA vs phCCA/dCCA #! Note: Labels are swapped compard to the manuscript figure (iCCA is on the left in the manuscript). Values are all the same with swapped sign.
all_test_results_df$comparison %>% table()
plot_df <- all_test_results_df %>% filter(comparison == "iCCA_vs_phCCA_dCCA") %>% mutate(Group1 = str_replace(Group1,"_","/"))
range(plot_df$effect.size)
xBreaks <- round(seq(-0.5,0.25,0.25),1)
xLims <- c(-0.5,0.26)
pt_L <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors, drop = F)+
  size_definition

# iCCA vs GBC #! Note: Labels are swapped compard to the manuscript figure (iCCA is on the left in the manuscript). Values are all the same with swapped sign.
plot_df <- all_test_results_df %>% filter(comparison == "iCCA_vs_GBC")
range(plot_df$effect.size)
xBreaks <- round(seq(-1.5,0.5,0.5),1)
xLims <- c(-1.5,0.5)
pt_M <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors, drop = F)+
  size_definition

# GBC vs phCCA/dCCA #! Note: Labels are swapped compard to the manuscript figure (iCCA is on the left in the manuscript). Values are all the same with swapped sign.
plot_df <- all_test_results_df %>% filter(comparison == "GBC_vs_phCCA_dCCA") %>% mutate(Group1 = str_replace(Group1,"_","/")) %>% filter(!is.na(p.val_lm))
range(plot_df$effect.size)
xBreaks <- round(seq(-0.5,1,0.25),1)
xLims <- c(-0.51,1.1)
pt_N <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks) +
  scale_fill_manual(values = group_colors, drop = F)+
  size_definition


# Save volcanos
width_volcano <- 4.5
height_volcano <- 5
ggsave(paste0(save_fig_folder,"/D_Volcano_iCCA_vs_Adj_iCCA.pdf"), pt_D, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/H_Volcano_CRLM_vs_Adj_CRLM.pdf"), pt_H, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/L_Volcano_iCCA_vs_phCCA_dCCA.pdf"), pt_L, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/M_Volcano_iCCA_vs_GBC.pdf"), pt_M, width = width_volcano, height = height_volcano)
ggsave(paste0(save_fig_folder,"/N_Volcano_GBC_vs_phCCA_dCCA.pdf"), pt_N, width = width_volcano, height = height_volcano)