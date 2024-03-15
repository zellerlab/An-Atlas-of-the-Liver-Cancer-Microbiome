#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 8 and 9
# by Fabian Springer
########

library(tidyverse)
library(here)
require(yaml)
library(ComplexHeatmap)
library(circlize)

source(here("src","plotting","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","ExtendedDataFigures","ExtendedDataFigure8and9")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

#* Generate heatmaps of associations between species abundances and clinical characteristics in HCC, iCCA and CRLM patients ----
# Import test resluts of the 5R16S dataset
test_results_clinicalAssociations_HCC_df <- readRDS(here("data","results","5R16S_ClinicalFeaturesAssociation_HCC_df.rds"))
test_results_clinicalAssociations_CRLM_df <- readRDS(here("data","results","5R16S_ClinicalFeaturesAssociation_CRLM_df.rds"))
test_results_clinicalAssociations_iCCA_df <- readRDS(here("data","results","5R16S_ClinicalFeaturesAssociation_iCCA_df.rds"))

# Select features representative of cancer stage, liver health and demography
cancer_stage_features <- c(
    "AFP_ug_ul", "tumor_size_max_cm", "macrovasvular_invasion_yes",
    "microvascular_invasion_yes", "multifocal_disease_status_multinodular",
    "T_T1", "T_T2", "T_T3", "T_T4", "N_N1", "M_M1", "surgical_margin_R0", "surgical_margin_R1",
    "surgical_margin_R2", "BCLC_stage_0", "BCLC_stage_A", "BCLC_stage_B", "BCLC_stage_C"
)
liver_health_features <- c(
    "Child_Pugh_score_A", "Child_Pugh_score_B", "Child_Pugh_score_C", "pre_op_bilirubin_mg_dl",
    "pre_op_blood_urea_nitrogen", "pre_op_gfr_ml_min", "pre_op_albumin_g_l", "pre_op_alp_u_l",
    "pre_op_ast_u_l", "pre_op_alt_u_l", "pre_op_erythrocytes_x_10e9_l", "pre_op_ggt_u_l",
    "pre_op_hemoglobin_10_9_l", "pre_op_inr", "pre_op_leukocytes_x_10e9_l", "pre_op_platelet_count_10_9_l",
    "pre_op_prothrombin_time_quick_percent","diabetes_yes"
)
demographic_features <- c("Age_years", "Gender_male", "BMI")

# Generate tibble
selected_features_df <- tibble(original_feature = c(cancer_stage_features, liver_health_features, demographic_features))

# Generate clean names of the clinical features selected
selected_features_df$feature_clean <- c(
    "AFP", "Tumor size", "Macrovasc. invasion", "Microvasc. invasion", "Multifocal disease",
    "T1", "T2", "T3", "T4", "N1", "M1", "R0 margin", "R1 margin", "R2 margin",
    "BCLC 0", "BCLC A", "BCLC B", "BCLC C",
    "Child-Pugh A", "Child-Pugh B", "Child-Pugh C", "Bilirubin", "Urea",
    "GFR", "Albumin", "ALP", "AST", "ALT", "Erythrocytes", "GGT",
    "Hemoglobin", "INR", "Leukocytes", "Platelets", "Prothrombin","Diabetes",
    
    "Age", "Male sex", "BMI"
    #demographic_features
)

# Assign feature categories for heatmap spit
selected_features_df <- selected_features_df %>%    
    mutate(original_feature = factor(original_feature, levels = c(cancer_stage_features, liver_health_features, demographic_features))) %>%
    arrange(original_feature) %>%
    mutate(feature_clean = as_factor(feature_clean)) %>%
    mutate(feat_group = case_when(
        original_feature %in% cancer_stage_features ~ "Cancer stage",
        original_feature %in% liver_health_features ~ "Liver health",
        original_feature %in% demographic_features ~ "Demography"
    ))


# Function to generate heatmap 
f_clin_assoc_heatmap <- function(test_res_df,clin_features_clean_df = selected_features_df) {
  # Takes a dataframe of testing results between species abundances and clinical features (categorical and continuous) and plots a heatmap
  # Privided clinical features are further filtered: At least 20 samples for continuous features and 5 samples in each group for categorical features
  # Taxa are filtered: At least 1 FDR-significant association. If less than 5 taxa are FDR-significant, show all taxa with p < 0.05

  # Generate plot_df by subsetting the provided testing results to keep only clinical features of interest.
  # Also filter the clinical features to keep only those with at least 20 samples for continuous features and 5 samples in each group for categorical features
  # Select p-value and effect size measure: Spearman correlations for continuous features and linear models for categorical features
  plot_df <- test_res_df %>%
    inner_join(., clin_features_clean_df %>% dplyr::rename(comparison = original_feature)) %>%
    dplyr::relocate(feature_type) %>%
    # Select pvalues and correlation coefficients from spearman correlations for continuous features
    # and from linear models for categorical features
    mutate(
      p_val = ifelse(feature_type == "continuous", p.val_spearman, p.val_lm),
      p_val_adj = ifelse(feature_type == "continuous", p.val_adj_spearman, p.val_adj_lm),
      eff_size = ifelse(feature_type == "continuous", Spearman.R, effect.size)
    ) %>%
    # keep only clinical features with at elast 20 samples for continuous features and 5 samples in each group for categorical features
    filter((feature_type == "continuous" & N_Samples > 20) | (feature_type == "categorical" & N_Group1 > 5 & N_Group2 > 5)) %>%
    mutate(feature_clean_label = paste0(feature_clean, " (", N_Samples, ")")) %>%
    dplyr::select(comparison, feature_clean, feature_clean_label, tax, feat_group, feature_type, p_val, p_val_adj, eff_size)

  plot_df$comparison %>% unique()

  # Define thresholds for significance
  p_threshold <- 0.05
  fdr_threshold <- 0.2

  # Now generate the heatmap
  # Generate matrices of effect sizes, pvalues and p-adj values
  mat_effect_size <- plot_df %>%
    dplyr::select(tax, feature_clean_label, eff_size) %>%
    pivot_wider(names_from = feature_clean_label, values_from = eff_size, values_fill = 0) %>%
    column_to_rownames("tax") %>%
    as.matrix()

  mat_p <- plot_df %>%
    dplyr::select(tax, feature_clean_label, p_val) %>%
    pivot_wider(names_from = feature_clean_label, values_from = p_val, values_fill = 1) %>%
    column_to_rownames("tax") %>%
    as.matrix()

  mat_p_adj <- plot_df %>%
    dplyr::select(tax, feature_clean_label, p_val_adj) %>%
    pivot_wider(names_from = feature_clean_label, values_from = p_val_adj, values_fill = 1) %>%
    column_to_rownames("tax") %>%
    as.matrix()

  ### Additional filtering of taxa and clinical features:
  # 1) Filter for taxa:
  # At least 1 FDR-significant association.
  # If less than 5 taxa are FDR-significant, show all taxa with p < 0.05
  row_levels <- rownames(mat_p_adj)[rowSums(mat_p_adj < fdr_threshold) > 0] %>% sort()
  if (length(row_levels) < 5) {
    # when less than 5 taxa are FDR-significant, show all taxa with p < 0.05
    row_levels <- sort(rownames(mat_p)[rowSums(mat_p < p_threshold) > 2])
  }

  # 2) Filter for clinical features:
  # At least one nominally significant association (p < 0.05)
  column_levels <- colnames(mat_p)[colSums(mat_p[row_levels, ] < p_threshold) > 0] %>% sort()


  ### Generate a split vector based on the categories defined above
  split_vector <- plot_df %>%
    dplyr::select(feature_clean_label, feat_group) %>%
    distinct() %>%
    deframe()
  split_vector <- split_vector[column_levels]


  ### Define clinical features to be plotted next to each other
  # The code is a bit messy but basically it checks which clinincal features are present that should be plotted next to each other
  # and re-orders the column_levels accordingly.
  column_clean_labels <- plot_df %>%
    dplyr::select(feature_clean_label, feature_clean) %>%
    mutate(feature_clean = as.character(feature_clean)) %>%
    deframe()
  unique(column_clean_labels)
  ref_names_ordered <- column_clean_labels[column_levels]
  columns_to_group <- c("T1", "T2", "T3", "T4", "M1", "N1", "AFP", "Tumor size", "R0 margin", "Albumin", "Bilirubin", "ALT", "AST", "GGT", "ALP", "Urea", "GFR")
  which(ref_names_ordered %in% columns_to_group)

  present_columns <- columns_to_group[columns_to_group %in% ref_names_ordered]
  new_col_labels <- names(ref_names_ordered)[match(present_columns, ref_names_ordered)]
  new_col_labels <- c(new_col_labels, names(ref_names_ordered)[!(ref_names_ordered %in% columns_to_group)])

  stopifnot(all(column_levels %in% new_col_labels))

  # Subset matrices with selected row_levels and column_levels
  mat_effect_size <- mat_effect_size[row_levels, new_col_labels]
  mat_p <- mat_p[row_levels, new_col_labels]
  mat_p_adj <- mat_p_adj[row_levels, new_col_labels]
  # reorder the split vector
  split_vector <- split_vector[new_col_labels]

  colnames(mat_effect_size)

  stopifnot(colnames(mat_effect_size) == colnames(mat_p) & colnames(mat_effect_size) == colnames(mat_p_adj))

  # Normalize the score matrix for better visualization
  mat_effect_size_rescaled <- f_rescale_effect_sizes(mat_effect_size)

  # Define color mapping for the heatmap
  range(mat_effect_size_rescaled)
  max_val <- signif(max(abs(range(na.omit(mat_effect_size_rescaled)))), 1)
  col_range <- c(-max_val, 0, max_val)
  color_mapping <- circlize::colorRamp2(col_range, c("dodgerblue", "white", "tomato"))

  # Define legend for the heatmap
  leg_title <- "Effect size\n(scaled)"
  lgd_main <- Legend(
    col_fun = color_mapping,
    title = leg_title,
    legend_height = unit(3, "cm"),
    border = "black"
  )
  l2 <- Legend(labels = c(paste0("p<", p_threshold)), title = "", type = "points", pch = 20, legend_gp = gpar(col = 1), size = unit(1, "mm"), background = NA)
  l3 <- Legend(labels = c(paste0("q<", fdr_threshold)), title = "", type = "points", pch = 8, legend_gp = gpar(col = 1), size = unit(1, "mm"), background = NA)
  leg_combined <- packLegend(lgd_main, l2, l3, direction = "vertical")

  # Modify appearanece of row and column labels
  colLabs <- as.expression(lapply(colnames(mat_effect_size_rescaled), function(a) bquote(bold(.(a)))))
  rowLabs <- as.expression(lapply(f_create_label(rownames(mat_effect_size_rescaled)), function(a) bquote(bolditalic(.(a)))))

  hm <- ComplexHeatmap::Heatmap(mat_effect_size_rescaled,
    cluster_rows = T,
    cluster_columns = F,

    # Column style
    column_names_rot = 45, column_names_side = "bottom", column_names_centered = F,
    column_names_gp = gpar(fontsize = 8),
    column_labels = colLabs,

    # Add the column split by category
    column_split = factor(split_vector, levels = c("Cancer stage", "Liver health", "Demography")),
    column_title_gp = gpar(fill = c("#E58F9E", "#BE5400", "#CBA3D8"), fontsize = 8),

    # Row style
    row_dend_side = "right",
    row_names_gp = gpar(fontsize = 8),
    row_names_side = c("left"),
    row_gap = unit(c(1.5), "mm"),
    row_labels = rowLabs,

    # general style
    border = T,
    rect_gp = gpar(col = "black", lwd = 1),
    col = color_mapping,
    gap = c(3, 3),
    show_heatmap_legend = FALSE,


    # Indicate p-values < 0.05 and p.val_adj < 0.001
    layer_fun = function(j, i, x, y, width, height, fill) {
      # Show border around significant tiles
      sig_lvls_fdr <- pindex(mat_p_adj, i, j)
      sig_lvls_p <- pindex(mat_p, i, j)
      l <- sig_lvls_fdr < fdr_threshold
      p <- (sig_lvls_p < p_threshold & sig_lvls_fdr >= fdr_threshold)
      if (any(l)) {
        grid.points(x[l], y[l], pch = 8, size = unit(1.5, "mm"))
      }
      if (any(p)) {
        grid.points(x[p], y[p], pch = 16, size = unit(1, "mm"))
      }
    }
  )
  return(draw(hm, annotation_legend_list = leg_combined))
}

# Generate Heatmap for HCC
hmap_HCC <- f_clin_assoc_heatmap(
  test_res_df = test_results_clinicalAssociations_HCC_df %>% filter(!str_detect(comparison, "T_T1|T_T2|Child_Pugh_score_A|stage_0")), # Remove comparisons of categorical variables in order to keep only the "most extreme" per group
  clin_features_clean_df = selected_features_df)

# Generate Heatmap for iCCA
hmap_iCCA <- f_clin_assoc_heatmap(
  test_res_df = test_results_clinicalAssociations_iCCA_df %>% filter(!str_detect(comparison,"T_T1|T_T2|stage_0")),
  clin_features_clean_df = selected_features_df)

# Generate Heatmap for CRLM
hmap_CRLM <- f_clin_assoc_heatmap(
  test_res_df = test_results_clinicalAssociations_CRLM_df,
  clin_features_clean_df = selected_features_df)

# Save the heatmaps
# HCC
pdf(file = file.path(save_fig_folder, "8A_Heatmap_ClinicalFeatureAssociationHCC.pdf"), width = 10, height = 7)
draw(hmap_HCC)
dev.off()
# iCCA
pdf(file = file.path(save_fig_folder, "9A_Heatmap_ClinicalFeatureAssociationiCCA.pdf"), width = 9, height = 5)
draw(hmap_iCCA)
dev.off()
# CRLM
pdf(file = file.path(save_fig_folder, "9B_Heatmap_ClinicalFeatureAssociationCRLM.pdf"), width = 8, height = 4.5)
draw(hmap_CRLM)
dev.off()

#* Extended Data Figure 8B-D: Selected Kaplan-Meier plots of species presence in iCCA and CRLM patients----
# Import survival analysis results
survRes_iCCA_list <- readRDS(here("data","results","survival_resObjects_iCCA.rds"))
survRes_CRLM_list <- readRDS(here("data","results","survival_resObjects_CRLM.rds"))

# Define the species for which Kaplan-Meier plots should be generated
to_plot <- c("s__Leuconostoc citreum","g__Gemella;s__Unknown species267","g__Corynebacterium;s__Unknown species1091")

# Define the survival result objects in the order of the species
survRes_Obj_list <- list(survRes_iCCA_list,survRes_CRLM_list,survRes_CRLM_list)
dataset_name <- c("iCCA","CRLM","CRLM")
# Define figure panel numbers
panel_number <- c("B","C","D")

# select survival result object of the given taxon and plot the Kaplan-Meier curve
p <- 1
for(species in to_plot){
  
  message(species)
  survRes_Obj <- survRes_Obj_list[[p]]

  # get index of species in list
  idx <- which(str_detect(names(survRes_Obj), species))
  stopifnot(length(idx) == 1)
  species <- f_create_label(species)
  pt <- f_kaplan_meier(survRes_Obj[[idx]]$survFit,species)
  
  # save the plot
  out_name <- file.path(save_fig_folder,paste0("8",panel_number[p],"_KaplanMeier_",species,"_",dataset_name[p],".pdf"))
  pdf(out_name, width = 4, height = 5)
  print(pt, newpage = FALSE)
  dev.off()
  p <- p+1
}