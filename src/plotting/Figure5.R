#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Figure 5
# by Fabian Springer
########

library(tidyverse)
library(here)
require(yaml)
require(vegan)
require(labdsv)

source(here("src","figures","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","Figure5")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

#* Panel A-D: Selected Kaplan-Meier Plots ----
# Load the survival results
survRes_HCC_list <- readRDS(here("data","results","survival_resObjects_HCC.rds"))

# Define the species for which Kaplan-Meier plots should be generated
to_plot <- c("Roseomonas mucosa","Methylobacterium oryzae","Acinetobacter schindleri","Leuconostoc gelidum")
panel_number <- c("A","B","C","D")

# select survival result object of the given taxon and plot the Kaplan-Meier curve
p <- 1
for(species in to_plot){
  message(species)

  # get index of species in list
  idx <- which(str_detect(names(survRes_HCC_list), species))
  stopifnot(length(idx) == 1)
  pt <- f_kaplan_meier(survRes_HCC_list[[idx]]$survFit,species)
  
  # save the plot
  out_name <- file.path(save_fig_folder,paste0(panel_number[p],"_KaplanMeier_",species,".pdf"))
  pdf(out_name, width = 4, height = 5)
  print(pt, newpage = FALSE)
  dev.off()
  p <- p+1
}

#* Panel E-J: Selected clinical comparisons ----
# Load data and metadata
meta_all_df <- read_tsv(here("data","metadata","meta_5R16S_allSamples.tsv"))
species_relAbundance_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_allSamples.rds"))
HCC_types <- c("HBV_HCC","HCV_HCC","ALD/ASH_HCC","MAFLD/MASH_HCC","other_HCC")

# Subset metadata to consider only HCC samples
meta_HCC_df <- meta_all_df %>% filter(Etiology %in% HCC_types)

# Define the species for which comparisons should be plotted
to_plot <- c(
  "s__Lactobacillus casei",
  "s__Sphingomonas yanoikuyae",
  "s__Pseudomonas oleovorans",
  "g__Acinetobacter;s__Unknown species264",
  "s__Massilia timonae"
)
# generate clean species labels
to_plot_clean <- f_create_label(to_plot)

# define the clinical features to be plotted (order same as for species)
to_plot_clinical <- c(
  "AFP_ug_ul",
  "pre_op_alt_u_l",
  "Age_years",
  "multifocal_disease_status",
  "multifocal_disease_status"
)

meta_df <- meta_HCC_df
relAB_mat <- species_relAbundance_mat
tax <- to_plot[4]
clin_feature <- to_plot_clinical[4]

f_prepare_boxplot_plot_df <- function(tax,clin_feature,relAB_mat,meta_df){        
    # Takes metadata and rel. abundance matrix. 
    # Generates a dataframe for a given taxon and a given clinical continuous parameter that can be used for plotting. 
    idx <- which(str_detect(rownames(relAB_mat),tax))
    stopifnot(length(idx) == 1)
    plot_df <-
        relAB_mat[idx, , drop = F] %>%
        as_tibble(rownames = "bac") %>%
        gather(-bac, key = "Sample_ID", value = "rel") %>%
        left_join(., meta_df %>% transmute(Sample_ID, clin_feat := !!as.symbol(clin_feature))) %>%
        filter(!is.na(clin_feat)) %>% 
        mutate(isPrev = ifelse(rel > 0,"Present","Absent")) %>% 
        mutate(Group = as_factor(isPrev))
    
    # add sample numbers
    plot_df <- 
      left_join(plot_df, 
      plot_df %>% dplyr::select(isPrev, Sample_ID) %>%
        group_by(isPrev) %>%
        summarise(N = n())) %>%
      mutate(bac_lab = paste0(isPrev, "\n(N=", N, ")"))
    
    return(plot_df)
}
f_prepare_barplot_plot_df <- function(tax, clin_feature, relAB_mat, meta_df) {
  # Takes metadata and rel. abundance matrix.
  # Generates a dataframe for a given taxon and a given clinical categorical parameter that can be used for plotting.
  # Computes Fisher's exact test for the given taxon prevalence and clinical feature.

  idx <- which(str_detect(rownames(relAB_mat), tax))
  stopifnot(length(idx) == 1)
  plot_df <-
    relAB_mat[idx, , drop = F] %>%
    as_tibble(rownames = "bac") %>%
    gather(-bac, key = "Sample_ID", value = "rel") %>%
    left_join(., meta_df %>% transmute(Sample_ID, clin_feat := !!as.symbol(clin_feature))) %>%
    filter(!is.na(clin_feat)) %>%
    mutate(isPrev = ifelse(rel > 0, 1, 0))

  # Add sample numbers
  plot_df <-
    plot_df %>%
    left_join(plot_df %>%
      dplyr::select(clin_feat, Sample_ID) %>%
      group_by(clin_feat) %>%
      summarise(N = n())) %>%
    mutate(clin_feat_lab = paste0(clin_feat, "\n(N=", N, ")"))

  # Compute Fisher's exact test
  tables <- with(plot_df, table(bac, clin_feat_lab, isPrev))
  fisher_res <- apply(tables, 1, function(t) {
    tbl <- matrix(t, ncol = 2, byrow = TRUE)
    return(fisher.test(tbl)$p.value)
  })
  
  # Compute prevalence of the taxon per clinical category
  bar_df <-
    plot_df %>%
    group_by(clin_feat, isPrev, N, clin_feat_lab) %>%
    summarise(n_category = n()) %>%
    mutate(n_category_perc = n_category / N) %>%
    filter(isPrev == 1) %>%
    mutate(n_category_perc = n_category_perc * 100) %>%
    mutate(p.val_fisher = as.numeric(fisher_res))

  # If there are samples with one clinical condition, manually add the second one and set it to 0
  if (nrow(bar_df) == 1) {
    missing_feat <- unique(plot_df$clin_feat)[!(unique(plot_df$clin_feat) %in% bar_df$clin_feat)]
    bar_df <- bind_rows(
      bar_df,
      tibble(
        clin_feat = missing_feat,
        isPrev = 1,
        N = plot_df %>% filter(clin_feat == missing_feat) %>% nrow(),
        n_category = 0,
        mutate(p.val_fisher = as.numeric(fisher_res)),
        n_category_perc = 0
      ) %>%
        mutate(clin_feat_lab = paste0(clin_feat, "\n(N=", N, ")"))
    )
  }
  
  return(bar_df %>% mutate(clin_feat_lab = as_factor(clin_feat_lab)))
}




col_vec <- c("Absent" = "#D95F02","Present" = "#1B9E77")

# E: Lactobacillus casei vs AFP levels
tax_idx <- 1
plot_df <- f_prepare_boxplot_plot_df(
  tax = to_plot[tax_idx], clin_feature = to_plot_clinical[tax_idx],
  relAB_mat = species_relAbundance_mat,
  meta_df = meta_HCC_df
)
pt_E <- f_boxplot_by_group(plot_df %>% mutate(y = clin_feat + 0.1), xlab = "", ylab = "AFP [ug/ul]", corral.width = 0.49) +
  scale_y_log10(breaks = c(0.1, 10, 1000, 100000), limits = c(0.08, 250000), labels = c("0.1", "10", "1000", "100,000")) +
  # Use ggsignif to compute a wilcoxon test
  ggsignif::geom_signif(test = "wilcox.test", comparisons = list(c(unique(plot_df$isPrev)[1], unique(plot_df$isPrev)[2])), tip_length = 0) +
  theme_paper +
  scale_x_discrete(labels = c(unique(plot_df$bac_lab)[1], unique(plot_df$bac_lab)[2])) +
  scale_fill_manual(values = col_vec, labels = plot_df %>% dplyr::select(isPrev, bac_lab) %>% deframe() %>% unique()) +
  scale_color_manual(values = col_vec, labels = plot_df %>% dplyr::select(isPrev, bac_lab) %>% deframe() %>% unique()) +
  ggtitle(to_plot_clean[tax_idx]) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold.italic")
  )

# F: Sphingomonas yanoikuyae vs pre-op ALT levels
tax_idx <- 2
plot_df <- f_prepare_boxplot_plot_df(
  tax = to_plot[tax_idx], clin_feature = to_plot_clinical[tax_idx],
  relAB_mat = species_relAbundance_mat,
  meta_df = meta_HCC_df
)
pt_F <- f_boxplot_by_group(plot_df %>% mutate(y = clin_feat + 0.1), xlab = "", ylab = "ALT [U/L]", corral.width = 0.49) +
  scale_y_log10(limits = c(8,3000))+
  # Use ggsignif to compute a wilcoxon test
  ggsignif::geom_signif(test = "wilcox.test", comparisons = list(c(unique(plot_df$isPrev)[1], unique(plot_df$isPrev)[2])), tip_length = 0) +
  theme_paper +
  scale_x_discrete(labels = c(unique(plot_df$bac_lab)[1], unique(plot_df$bac_lab)[2])) +
  scale_fill_manual(values = col_vec, labels = plot_df %>% dplyr::select(isPrev, bac_lab) %>% deframe() %>% unique()) +
  scale_color_manual(values = col_vec, labels = plot_df %>% dplyr::select(isPrev, bac_lab) %>% deframe() %>% unique()) +
  ggtitle(to_plot_clean[tax_idx]) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold.italic")
  )

# G: Pseudomonas oleovorans vs Age
tax_idx <- 3
plot_df <- f_prepare_boxplot_plot_df(
  tax = to_plot[tax_idx], clin_feature = to_plot_clinical[tax_idx],
  relAB_mat = species_relAbundance_mat,
  meta_df = meta_HCC_df
)
pt_G <- f_boxplot_by_group(plot_df %>% mutate(y = clin_feat), xlab = "", ylab = "Age [years]", corral.width = 0.49) +  
  # Use ggsignif to compute a wilcoxon test
  ggsignif::geom_signif(test = "wilcox.test", comparisons = list(c(unique(plot_df$isPrev)[1], unique(plot_df$isPrev)[2])), tip_length = 0) +
  theme_paper +
  scale_x_discrete(labels = c(unique(plot_df$bac_lab)[1], unique(plot_df$bac_lab)[2])) +
  scale_fill_manual(values = col_vec, labels = plot_df %>% dplyr::select(isPrev, bac_lab) %>% deframe() %>% unique()) +
  scale_color_manual(values = col_vec, labels = plot_df %>% dplyr::select(isPrev, bac_lab) %>% deframe() %>% unique()) +
  ggtitle(to_plot_clean[tax_idx]) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold.italic")
  )+
  ylim(c(45,91))

# H: Acinetobacter US264 vs multifocal disease status
tax_idx <- 4
plot_df <- f_prepare_barplot_plot_df(
  tax = to_plot[tax_idx], clin_feature = to_plot_clinical[tax_idx],
  relAB_mat = species_relAbundance_mat,
  meta_df = meta_HCC_df
)
plot_df$clin_feat_lab %>% levels()
pt_H <- plot_df %>%
  ggplot(aes(x = fct_rev(clin_feat_lab), y = n_category_perc, fill = fct_rev(clin_feat_lab))) +
  theme_paper +
  geom_bar(stat = "identity", color = "black") +
  # geom_text(aes(label = paste0("N=", n_category), y = n_category_perc), vjust = -0.3) +
  ylab("Prevalence [%]") +
  # Use the signif package just for annotation of the Fisher test p-value
  ggsignif::geom_signif(
    comparisons = list(c(unique(as.character(plot_df$clin_feat_lab)[1]), unique(as.character(plot_df$clin_feat_lab)[2]))),
    annotations = signif(unique(plot_df$p.val_fisher), 1), tip_length = 0, y_position = 20
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 23)) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4")) +
  ggtitle(to_plot_clean[tax_idx]) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold.italic")
  )

# I: Masilla timonae vs multifocal disease status
tax_idx <- 5
plot_df <- f_prepare_barplot_plot_df(
  tax = to_plot[tax_idx], clin_feature = to_plot_clinical[tax_idx],
  relAB_mat = species_relAbundance_mat,
  meta_df = meta_HCC_df
)
plot_df$clin_feat_lab %>% levels()
pt_I <- plot_df %>%
  ggplot(aes(x = fct_rev(clin_feat_lab), y = n_category_perc, fill = fct_rev(clin_feat_lab))) +
  theme_paper +
  geom_bar(stat = "identity", color = "black") +
  # geom_text(aes(label = paste0("N=", n_category), y = n_category_perc), vjust = -0.3) +
  ylab("Prevalence [%]") +
  # Use the signif package just for annotation of the Fisher test p-value
  ggsignif::geom_signif(
    comparisons = list(c(unique(as.character(plot_df$clin_feat_lab)[1]), unique(as.character(plot_df$clin_feat_lab)[2]))),
    annotations = signif(unique(plot_df$p.val_fisher), 1), tip_length = 0, y_position = 18
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 21)) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4")) +
  ggtitle(to_plot_clean[tax_idx]) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold.italic")
  )

# Save boxplots and barplots of selected clinical comparisons: 
w <- 3.5
h <- 4
ggsave(file.path(save_fig_folder,"E_ClinicalComparison_LactobacillusCasei.pdf"), pt_E, width = w, height = h)
ggsave(file.path(save_fig_folder,"F_ClinicalComparison_SphingomonasYanoikuyae.pdf"), pt_F, width = w, height = h)
ggsave(file.path(save_fig_folder,"G_ClinicalComparison_PseudomonasOleovorans.pdf"), pt_G, width = w, height = h)
ggsave(file.path(save_fig_folder,"H_ClinicalComparison_AcinetobacterUS264.pdf"), pt_H, width = w, height = h)
ggsave(file.path(save_fig_folder,"I_ClinicalComparison_MassiliaTimonae.pdf"), pt_I, width = w, height = h)


#* Compute Shannon diversity and species richness of immunotherapy responders vs non-responders ----
# Load data and metadata
meta_immunotherapy_df <- read_tsv(here("data","metadata","meta_5R16S_ImmunotherapyCohort.tsv"))
species_relAbundance_immunotherapy_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_ImmunotherapyCohort.rds"))
# Load the testing results of alpha- and beta-diversity
alpha_beta_diversity_immuno_clean_df <- read_tsv(here("data","results","5R16S_ImmunotherapyCohort_diversity_comparisons_df.tsv"))
alpha_beta_diversity_immuno_clean_df %>% 

# Compute diversity
div <- vegan::diversity(t(species_relAbundance_immunotherapy_mat[,meta_immunotherapy_df$Sample_ID]),index = "shannon") %>% enframe(name = "Sample_ID",value = "ShannonDiv")
# Compute Richness
threshold_for_prev <- 0
rich <- colSums(species_relAbundance_immunotherapy_mat>threshold_for_prev) %>% enframe(name = "Sample_ID",value = "SpeciesRichness")
diversity_richness_df <- full_join(div,rich)

immuno_shannon_richness_df <- diversity_richness_df %>%
  gather(-Sample_ID, key = "metric", value = "value") %>%
  left_join(., meta_immunotherapy_df %>% dplyr::select(Sample_ID, Response_to_immunotherapy))


# K: Species richness of responders vs non-responders
tmp_df <- immuno_shannon_richness_df %>%
  filter(metric == "SpeciesRichness") %>%
  mutate(
    Group = factor(Response_to_immunotherapy, levels = c("non-responder", "responder")),
    y = value
  )
# Get the p-value for the Wilcoxon test
p_val <- alpha_beta_diversity_immuno_clean_df %>% pull(p.val_adj_wilcox_richness_species) # It pulls the adjusted p-value but since there is only one comparison, it is the same as the unadjusted p-value
pt_K <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Species richness", corral.width = 0.49) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # Just use the ggsignif function for plotting the p-value symbols
  ggsignif::geom_signif(
    comparisons = list(c("non-responder", "responder")),
    annotations = p_to_symbol(p_val)
  )

# L: Shannon diversity of responders vs non-responders
tmp_df <- immuno_shannon_richness_df %>%
  filter(metric == "ShannonDiv") %>%
  mutate(
    Group = factor(Response_to_immunotherapy, levels = c("non-responder", "responder")),
    y = value
  )
# Get the p-value for the Wilcoxon test
p_val <- alpha_beta_diversity_immuno_clean_df %>% pull(p.val_adj_wilcox_shannon_diversity) # It pulls the adjusted p-value but since there is only one comparison, it is the same as the unadjusted p-value
pt_L <- f_boxplot_by_group(tmp_df, xlab = "", ylab = "Shannon diversity", corral.width = 0.49) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # Just use the ggsignif function for plotting the p-value symbols
  ggsignif::geom_signif(
    comparisons = list(c("non-responder", "responder")),
    annotations = p_to_symbol(p_val)
  )

# M: PCoA of immunotherapy responders vs non-responders
# Compute PCoA of Bray-Curtis distances at species level
p_val <- alpha_beta_diversity_immuno_clean_df$p.val_adj_permanova_Bray
# p_val_text
p_text <- paste("NR vs. R:", signif(p_val,1))
pt_M <- f_plot_PCoA(
  meta_df = meta_immunotherapy_df %>% mutate(Group = Response_to_immunotherapy),
  mat = species_relAbundance_immunotherapy_mat, method = "bray"
)+
scale_fill_manual(values = group_colors)+
scale_color_manual(values = group_colors)+
theme(legend.position = "none")+
# add pvalue to bottom left
annotate("text", x = -Inf, y = -Inf, label = p_text, size = 3, color = "black",hjust = -0.1, vjust = -0.5)

w <- 3.5
h <- 4.5
ggsave(pt_K,filename = file.path(save_fig_folder,"K_SpeciesRichness_ImmunotherapyCohort.pdf"), width = w, height = h)
ggsave(pt_L,filename = file.path(save_fig_folder,"L_ShannonDiv_ImmunotherapyCohort.pdf"), width = w, height = h)
ggsave(pt_M,filename = file.path(save_fig_folder,"M_PCoA_ImmunotherapyCohort.pdf"), width = 6, height = 5)

#* Volcano plot of responder vs non-responder ----
# Import testing results
test_res_immonotherapy_df <- readRDS(here("data","results","5R16S_ImmunotherapyCohortTesting_result_df.rds"))

# Check effect sizes and p-value ranges; define x- and y-limits for plotting
test_res_immonotherapy_df %>%
  glimpse()
plot_df <- test_res_immonotherapy_df
range(plot_df$effect.size)
xBreaks <- c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
xLims <- range(xBreaks)
range(plot_df$p.val_lm)
man_y_breaks <- c(-log10(c(0.05,0.01,0.001)))

# plot the volcano of immunotherapy responders vs non-responders
pt_N <- f_plot_volcano(
  plot_df = test_res_immonotherapy_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks
) +
  scale_fill_manual(values = group_colors) +
  scale_size_continuous(
    name = "Prevalence",
    range = c(1, 6),
    limits = c(0, 1),
    breaks = c(0.01, 0.2, 0.4, 0.6)
  )
width_volcano <- 4.5
height_volcano <- 5 
ggsave(pt_N, filename = file.path(save_fig_folder,"N_Volcano_ResponderVsNonResponder.pdf"), width = width_volcano,height = height_volcano)


