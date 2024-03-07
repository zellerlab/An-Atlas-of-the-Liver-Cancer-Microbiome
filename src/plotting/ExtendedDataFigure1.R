#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 1
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
save_fig_folder <- here("figures","ExtendedDataFigures")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import data
meta_all_df <- read_tsv(here("data","metadata","meta_5R16S_allSamples.tsv"))
species_relAbundance_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_allSamples.rds"))
HCC_types <- c("HBV_HCC","HCV_HCC","ALD/ASH_HCC","MAFLD/MASH_HCC","other_HCC")

# Subset metadata to consider only cancer samples
table(meta_all_df$Etiology)
meta_cancer_df <- meta_all_df %>%
  transmute(Sample_ID, Etiology) %>%
  mutate(cancer = case_when(
    Etiology %in% HCC_types ~ "HCC",
    Etiology %in% c("CRLM", "iCCA") ~ Etiology,
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(cancer)) %>% 
  dplyr::select(-Etiology)
table(meta_cancer_df$cancer)

# Factor levels for cancer types
cancer_levels <- c("HCC", "iCCA", "CRLM")


#* H: Volcano plots of rel. abundance testings between HCC/iCCA vs CRLM ----
# Load the testing result
test_res_cancer_df <- readRDS(here("data", "results", "5R16S_EtiologyTesting_result_df.rds")) %>%
  filter(comparison %in% c("HCC_iCCA_vs_CRLM")) %>% 
  mutate(Group1 = ifelse(Group1 == "HCC_iCCA", "HCC/iCCA", Group1))
table(test_res_cancer_df$comparison)
test_res_cancer_df$Group1
# Define size and ranges of the points in the volcano plots
size_definition <- scale_size_continuous(
  name = "Prevalence",
  range = c(1, 4),
  limits = c(0, 0.3),
  breaks = c(0.05, 0.10, 0.15, 0.2),
  guide = guide_legend(reverse = TRUE)
)
# Define y-range for all volcano plots
man_y_breaks <- c(-log10(c(0.05,0.01,0.001,1e-4,1e-5)))

plot_df <- test_res_cancer_df
range(plot_df$effect.size)
xBreaks <- c(seq(-0.3,0,0.1),seq(0.1,0.4,0.1))
xLims <- c(-0.31,0.4)
range(plot_df$p.val_lm)
pt_H <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors) +
  size_definition

w <- 4.5
h <- 5 
ggsave(pt_H, filename = file.path(save_fig_folder,"ExtDataFig1H_Volcano_HCCiCCAvsCRLM.pdf"), width = w,height = h)

#* G: PCoA of different cancer samples ----
alpha_beta_diversity_clean_df <- read_tsv(here("data","results","5R16S_diversity_comparisons_df.tsv")) %>% filter(comparison == "HCC_vs_iCCA_vs_CRLM")
alpha_beta_diversity_clean_df %>% glimpse()

df <- df %>%
  mutate(p_text = paste0(Group1, " vs ", Group2, ": ", signif(p.val_adj_permanova_Bray,1))) %>% 
  dplyr::relocate(p_text)
p_text <- paste(df$p_text, collapse="\n")

pt_G <- f_plot_PCoA(
  meta_df = meta_cancer_df %>% mutate(Group = cancer),
  mat = species_relAbundance_mat, method = "bray"
)+
scale_fill_manual(values = group_colors)+
scale_color_manual(values = group_colors)+
theme(legend.position = "none")+
# add pvalue to bottom left
annotate("text", x = -Inf, y = -Inf, label = p_text, size = 3, color = "black",hjust = -0.1, vjust = -0.5)
rowSums(species_relAbundance_mat > 0) %>% summary()

ggsave(pt_H,filename = file.path(save_fig_folder,"ExtDataFig1G_PCoA_HCCvsiCCAvsCRLM.pdf"), width = 6, height = 5)




