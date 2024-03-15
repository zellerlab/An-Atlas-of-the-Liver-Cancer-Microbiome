#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 10
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
save_fig_folder <- here("figures","ExtendedDataFigures","ExtendedDataFigure10")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import raw data of the immunotherapy cohort
meta_immunotherapy_df <- read_tsv(here("data","metadata","meta_5R16S_ImmunotherapyCohort.tsv"))
species_relAbundance_immunotherapy_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_ImmunotherapyCohort.rds"))

# Import testing results
test_res_immonotherapy_df <- readRDS(here("data","results","5R16S_ImmunotherapyCohortTesting_result_df.rds"))

# select all species with p-value < 0.05 in the linear models (the same shown as colored dots in the Volcano plot in Figure 5)
to_plot <- test_res_immonotherapy_df %>%
  mutate(lab_clean = f_create_label(tax)) %>%
  dplyr::relocate(lab_clean) %>%
  filter(p.val_lm < 0.05) %>%
  arrange(-Prev_Group2,-Prev_Group1) %>% # arrange by prevalence in Group2  
  pull(tax)

# generate clean species labels
to_plot_clean <- f_create_label(to_plot)
# replace the last space with a linebreak for better visualization
to_plot_clean <- sub(" (?!.* )", "\n", to_plot_clean, perl = TRUE)

f_generate_barplot <- function(tax,tax_clean) {
  # Just a wrapper to generate the barplot for a given taxon
  bar_df <- suppressMessages(f_prepare_barplot_plot_df(
    tax = tax,
    clin_feature = "Response_to_immunotherapy",
    relAB_mat = species_relAbundance_immunotherapy_mat,
    meta_df = meta_immunotherapy_df
  )) %>%
    mutate(
      n_category_perc_orig = n_category_perc, # save original values
      n_category_perc = ifelse(n_category_perc == 0, 0.1, n_category_perc) # for visualization: add small prevalence when real prevalence is zero
    )

  pt_bar <- bar_df %>%
    ggplot(aes(x = clin_feat_lab, y = n_category_perc, fill = clin_feat)) +
    theme_paper +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = paste0(round(n_category_perc_orig, 0), "%"), y = n_category_perc + 0.1), vjust = -0.3,size = 3) +
    ylab("Prevalence [%]") +
    # Use the signif package just for annotation of the Fisher test p-value
    ggsignif::geom_signif(
      comparisons = list(c(unique(as.character(bar_df$clin_feat_lab)[1]), unique(as.character(bar_df$clin_feat_lab)[2]))),
      annotations = paste0("P=", signif(unique(bar_df$p.val_fisher), 2)),
      tip_length = 0, y_position = 25.5,margin_top = 0,textsize = 3
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 27), breaks = seq(0,25,5),labels = seq(0,25,5)) +
    scale_fill_manual(values = group_colors) +
    ggtitle(tax_clean) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5, face = "bold.italic")
    )
  return(pt_bar)
}
i <- 1
pt_list <- list()
for(i in seq(1,length(to_plot))){  
  tax <- to_plot[i]
  tax_clean <- to_plot_clean[i]
  message(tax_clean)
  p <- f_generate_barplot(tax = tax,tax_clean = tax_clean)+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())
  if(i > 1){
    p <- p + theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())
  }
  pt_list[[i]] <- p
}
pt_A <- wrap_plots(pt_list,nrow = 1)
ggsave(pt_A,filename = file.path(save_fig_folder,"A_Barplot_PrevalenceRespVsNonResp.pdf"),width = 14,height = 6)


#* B: Forest plot of Hazard ratios for the influence of species abundances on overall survival ----
# Import survival analysis result
survival_result_immunotherapy_df <- read_tsv(here("data","results","survival_results_ImmunotherapyCohort_df.tsv"))
p_threshold <- 0.05
fdr_threshold <- 0.2

plot_df <-
  survival_result_immunotherapy_df %>%
  # Filter to show only species with an nominally significant p-value
  filter(p.val_Cox < p_threshold) %>%
  dplyr::select(tax, HR, HR_lower_95_CI, HR_upper_95_CI, p.val_Cox, p.val_adj_Cox) %>%
  # assign significance level
  mutate(sig_type = case_when(
    p.val_adj_Cox < fdr_threshold & HR > 1 ~ "sigPos",
    p.val_adj_Cox < fdr_threshold & HR < 1 ~ "sigNeg",
    TRUE ~ "n.s."
  )) %>%
  arrange(HR) %>%
  mutate(
    tax_label = as_factor(f_create_label(tax)) # assign clean species names and convert to factor
  ) %>% 
  relocate(tax_label) %>%
  mutate(HR_upper_95_CI = case_when( # if upper CI of HR is infinite, replace with HR
    is.infinite(HR_upper_95_CI) ~ HR,
    TRUE ~ HR_upper_95_CI
  )) %>%

  # floor and ceil all values between a HR of 0.1 and 10
  mutate(
    HR_mod = ifelse(HR < 0.1,0.1,HR),
    HR_upper = ifelse(HR_upper_95_CI > 10, 10, HR_upper_95_CI),
    HR_lower = ifelse(HR_lower_95_CI < 0.1, 0.1, HR_lower_95_CI)
  ) 

# Generate the plot
xBreaks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
xLabels <- c("\u2264 0.1", 0.2, 0.5, 1, 2, 5, "\u2265 10")
sig_color_vec <- c("sigPos" = "#E40046", "sigNeg" = "#307FE2")
sig_color_label <- c("sigPos" = "HR > 1", "sigNeg" = "HR < 1")

pt_B <- plot_df %>%
  ggplot(aes(y = fct_rev(tax_label))) +
  theme_paper +
  geom_vline(xintercept = 1, color = "black", lwd = 0.5) +
  geom_point(aes(x = HR_mod, fill = sig_type, color = sig_type), size = 1.5, show.legend = T) +
  geom_linerange(aes(
    xmin = HR_lower,
    xmax = HR_upper, color = sig_type,
  )) +
  scale_x_log10(
    breaks = xBreaks,
    labels = xLabels,
    limits = c(0.09, 11)
  ) +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_color_manual(values = sig_color_vec, labels = sig_color_label) +
  scale_fill_manual(values = sig_color_vec, labels = sig_color_label) +
  labs("color" = "FDR<0.2") +
  guides(fill = FALSE) +
  ylab("") +
  xlab("Hazard ratio for Death") +
  theme(
    legend.position = c(0, 0), legend.justification = c(-0.1, -0.1),
    legend.background = element_rect(color = "black")
  )

# Note: No species has an FDR-significant HR in the immunotherapy cohort
ggsave(pt_B,filename = file.path(save_fig_folder,"B_ForestPlot_HazardRatioImmunotherapyCohort.pdf"),width = 6,height = 10,device=cairo_pdf)
