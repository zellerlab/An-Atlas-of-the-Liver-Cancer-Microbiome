#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Figure 1
# by Fabian Springer
########

library(tidyverse)
library(here)
require(yaml)
require(vegan)
require(labdsv)
require(patchwork)

source(here("src","plotting","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","Figure1")
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


#* Panel i and j: Barplots of Shannon Diversity and genus richness in cancer samples ----
# Load data and metadata
# Load the testing results of alpha- and beta-diversity
alpha_beta_diversity_clean_df <- read_tsv(here("data","results","5R16S_diversity_comparisons_df.tsv")) %>% filter(comparison == "HCC_vs_iCCA_vs_CRLM")

# Compute diversity
div <- vegan::diversity(t(species_relAbundance_mat[,meta_cancer_df$Sample_ID]),index = "shannon") %>% enframe(name = "Sample_ID",value = "ShannonDiv")
# Compute Richness
threshold_for_prev <- 0
rich <- colSums(species_relAbundance_mat>threshold_for_prev) %>% enframe(name = "Sample_ID",value = "SpeciesRichness")
diversity_richness_df <- full_join(div,rich)

cancer_richness_df <- diversity_richness_df %>%
  gather(-Sample_ID, key = "metric", value = "value") %>%
  inner_join(., meta_cancer_df)


# K: Species richness by cancer type
tmp_df <- cancer_richness_df %>%
  filter(metric == "SpeciesRichness") %>%
  mutate(
    Group = factor(cancer, levels = cancer_levels),
    y = value
  )
# Get the p-values from the Wilcoxon testing of the species richness and shannon diversity
# Get list of comparisons and the p-values:
df <- alpha_beta_diversity_clean_df
comparisons_list <- lapply(1:nrow(df), function(i) c(df$Group1[i], df$Group2[i])) %>% rev() # reverse order for better visualization
p_adj_vec_richness <- df$p.val_adj_wilcox_richness_species %>% rev()

pt_tmp <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Species richness", corral.width = 0.49) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # Just use the ggsignif function for plotting the p-value symbols
  ggsignif::geom_signif(
    comparisons = comparisons_list,
    annotations = p_to_symbol(p_adj_vec_richness),
    vjust =-0,step_increase = 0.075,textsize = 5
  )
# For better visualization: Manually set y-values > 100 to 110 in order to decrease the plotting window
orig_max <- max(pt_tmp$data$y)
pt_tmp$data <- pt_tmp$data %>%
    mutate(y = ifelse(y > 100, 110, y))
pt_I <- pt_tmp+
  scale_y_continuous(
    breaks = c(seq(0,100,25),110),labels = as.character(c(seq(0,100,25),orig_max)),
    limits = c(0,135))

# L: Shannon diversity of responders vs non-responders
tmp_df <- cancer_richness_df %>%
  filter(metric == "ShannonDiv") %>%
  mutate(
    Group = factor(cancer, levels = cancer_levels),
    y = value
  )
# Get the p-value from Wilcoxon-Tests of Shannon diversity
alpha_beta_diversity_clean_df %>% colnames()
p_adj_vec_shannon <- alpha_beta_diversity_clean_df$p.val_adj_wilcox_shannon_diversity %>% rev()

pt_J <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Shannon diversity", corral.width = 0.49) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # Just use the ggsignif function for plotting the p-value symbols
  ggsignif::geom_signif(
    comparisons = comparisons_list,
    annotations = p_to_symbol(p_adj_vec_shannon),
    vjust =-0,step_increase = 0.075,textsize = 5
  )+
  ylim(c(0,5))

# Save shannon and Richness plot
w <- 3.5
h <- 4.4
ggsave(pt_I,filename = file.path(save_fig_folder,"I_SpeciesRichness.pdf"),width = w,height = h)
ggsave(pt_J,filename = file.path(save_fig_folder,"J_ShannonDiversity.pdf"),width = w,height = h)


#* Panel l: Prevalence of 10 highest prevalent species in tumor samples ----

# get highest prevalence species across cancer types
species_rank_df <- species_relAbundance_mat %>%
  as_tibble(rownames = "tax") %>%
  gather(Sample_ID, relAbundance, -tax) %>%
  inner_join(., meta_cancer_df) %>% 
  group_by(tax, cancer) %>%
  summarise(prev = sum(relAbundance > 0) / n()) %>% 
  ungroup() %>% 
  group_by(cancer) %>% 
  mutate(prev_rank = rank(-prev)) %>% 
  ungroup()

top10_species_vec <-   
  species_rank_df %>% 
  filter(prev_rank <= 10) %>% 
  arrange(prev_rank) %>% 
  pull(tax) %>% unique()

# Get factor levels: based on prevalence in HCC samples
species_levels <- species_rank_df %>%
  filter(cancer == "HCC",tax %in% top10_species_vec) %>% 
  arrange(prev_rank) %>% 
  pull(tax)

# Perform Fisher's exact tests
# generate df with presence/absence of the selected species in the cancer samples
presence_df <- species_relAbundance_mat[species_levels, meta_cancer_df$Sample_ID] %>%
  as_tibble(rownames = "tax") %>%
  gather(Sample_ID, relAbundance, -tax) %>% 
  mutate(presence = relAbundance > 0) %>% 
  inner_join(., meta_cancer_df) %>% 
  dplyr::select(tax,Sample_ID,presence,cancer)

# Generate contingency tables for each tax
tables <- with(presence_df, table(tax, cancer, presence))
fisher_res <-
  apply(tables, 1, function(t) {
    tbl <- matrix(t, ncol = 3, byrow = TRUE)
    return(fisher.test(tbl)$p.value)
  })
fisher_res_df <- fisher_res %>%
  enframe() %>%
  dplyr::rename(tax = name, p.val_fisher = value) %>%
  mutate(p.val_adj_fisher = p.adjust(p.val_fisher, method = "fdr")) %>% # correct pvalues
  mutate(fdr_symbol = sapply(p.val_adj_fisher, p_to_symbol)) %>%  # convert pvalues to symbols
  mutate(fdr_symbol = ifelse(fdr_symbol == "n.s.", "", fdr_symbol))

# Generate df for plotting with prevalences, pvalues and clean names.
plot_df <- species_rank_df %>%
  filter(tax %in% top10_species_vec) %>%
  left_join(., fisher_res_df) %>% 
  dplyr::relocate(cancer,prev,p.val_adj_fisher,fdr_symbol) %>% 
  mutate(
    tax = factor(tax, levels = species_levels),
    tax_clean = f_create_label(tax)
  ) %>% 
  arrange(tax) %>% 
  mutate(tax_clean = as_factor(tax_clean)) %>% 
  mutate(cancer = factor(cancer,levels = c("HCC","iCCA","CRLM")))

# Generate the barplot of prevalences
p1 <- plot_df %>%
  ggplot(aes(y = fct_rev(tax_clean), x = prev * 100, fill = fct_rev(cancer))) +
  theme_paper +
  theme(
    plot.margin = margin(1, 0, 1, 1, "cm"),
    legend.position = "left"
  ) +
  annotate(
    geom = "segment",
    y = Inf,
    yend = Inf,
    x = -Inf,
    xend = Inf, size = 1.5
  ) +
  theme(panel.border = element_blank(), axis.line = element_line()) + # some black magic to have the right side open
  labs(
    title = NULL,
    y = "",
    x = "Prevalence [%]", fill = "Cancer Type"
  ) +
  ggforestplot::geom_stripes(odd = "lightgrey", even = "white") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.y = element_text(face = "bold.italic"))

# Add a right plot with the symbols for pvalues
fdr_df <- plot_df %>% group_by(tax_clean,fdr_symbol,p.val_adj_fisher) %>% summarise(prev = max(prev))
p2 <- fdr_df %>% 
  ggplot(aes(y = fct_rev(tax_clean), x = 1)) +
   theme_void()+
    theme(
    plot.margin = margin(1,1,1,0,"cm"))+
    #ggforestplot::geom_stripes(odd = "lightgrey",even = "white")+
    geom_text(aes(label = fdr_symbol,hjust = 0.125,vjust = 0.8),size=8)    # adjust left and right margin as needed    
  #xlab("FDR <0.1")
pt_l <- (p1+theme(legend.position="none"))+p2+plot_layout(ncol=2,widths = c(4, 0.2))
pt_l

ggsave(pt_l,filename = file.path(save_fig_folder,"L_Barplot_MostPrevalentSpecies.pdf"),width = 10,height = 8)



#* m-o: Volcano plots of rel. abundance testings between cancer samples ----
# Load the testing result
test_res_cancer_df <- readRDS(here("data", "results", "5R16S_EtiologyTesting_result_df.rds")) %>%
  filter(comparison %in% c("HCC_vs_CRLM", "HCC_vs_iCCA", "iCCA_vs_CRLM"))
table(test_res_cancer_df$comparison)

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

# m: HCC vs iCCA
#! Important note: In the manuscript the order of the levels in the volcano plot are swapped:
#! HCC-enriched species are on the left (Enrichment effect size <1), iCCA enriched species are on the right (Enrichment effect size >1).
#! Absolute values and p-values are identical so this is just a matter of visualization.

# Check effect sizes and p-value ranges; define x- and y-limits for plotting
plot_df <- test_res_cancer_df %>% filter(comparison == "HCC_vs_iCCA")
range(plot_df$effect.size)
xBreaks <- c(-0.2,-0.1,0,0,0.1,0.2)
xLims <- c(-0.25,0.25)
range(plot_df$p.val_lm)
pt_M <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors) +
  size_definition

# HCC vs CRLM
plot_df <- test_res_cancer_df %>% filter(comparison == "HCC_vs_CRLM")
range(plot_df$effect.size)
xBreaks <- c(seq(-0.3,0,0.1),seq(0.1,0.4,0.1))
xLims <- range(xBreaks)
pt_N <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors) +
  size_definition

# iCCA vs CRLM
plot_df <- test_res_cancer_df %>% filter(comparison == "iCCA_vs_CRLM")
range(plot_df$effect.size)
xBreaks <- c(seq(-0.4,0,0.2),seq(0.2,0.6,0.2))
xLims <- range(xBreaks)
pt_O <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors) +
  size_definition

# save volcano plots
w <- 4.5
h <- 5 
ggsave(pt_M, filename = file.path(save_fig_folder,"M_Volcano_HCCvsiCCA.pdf"), width = w,height = h)
ggsave(pt_N, filename = file.path(save_fig_folder,"N_Volcano_HCCvsCRLM.pdf"), width = w,height = h)
ggsave(pt_O, filename = file.path(save_fig_folder,"O_Volcano_iCCAvsCRLM.pdf"), width = w,height = h)

#* qPCR Plotting ----
# Import the qPCR result table
qPCR_res_df <- read_tsv(here("data","results","qPCR_results_df.tsv"))
qPCR_concentrations_df <- read_tsv(here("data","raw","qPCR_concentrations_df.tsv"))


yBreaks <- c(0.1, 1, 10, 100)
yLabs  <-  as.character(yBreaks)

plot_df <- qPCR_concentrations_df %>%
  left_join( # add sample numbers
    .,
    qPCR_concentrations_df %>%
      group_by(Group) %>%
      summarise(N = n())
  ) %>%
  mutate(Group_N = paste0(Group, "\n(N=", N, ")")) %>%
  mutate(Group = factor(Group, levels = c("CTRL", "HCC", "BTC", "CRLM"))) %>% # define factor levels
  arrange(Group) %>%
  mutate(Group_N = as_factor(Group_N))

xLabs <- plot_df %>% select(Group,Group_N) %>% distinct() %>% deframe()


# Get p-values from the pairwise testing
df <- qPCR_res_df
comparisons_list <- lapply(1:nrow(df), function(i) c(df$Group1[i], df$Group2[i]))
p_adj_vec_qPCR <- df$p.val_adj_wildox

pt_K <- f_boxplot_by_group(plot_df %>% mutate(y = bacterial_cells_per_40ngDNA + 0.1), xlab = "", ylab = "#Bacterial cells/40ng DNA", corral.width = 0.49) +
  scale_y_log10(breaks = yBreaks,labels = yLabs)+
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(breaks = names(xLabs),labels = as.character(xLabs))+
  theme(axis.title.x = element_blank(),legend.position = 'none')+
  # Just use the ggsignif function for plotting the p-value symbols
  ggsignif::geom_signif(
    comparisons = wilcox_tests_list,
    annotations = p_to_symbol(p_adj_vec_qPCR),
    vjust =-0,step_increase = 0.075,textsize = 5
  )
ggsave(pt_K,filename = file.path(save_fig_folder,"K_qPCR_BacterialLoad.pdf"), width = 6,height = 5)










