#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Figure 4
# by Fabian Springer
########

# Import packages
library(tidyverse)
library(here)
require(yaml)
library(ComplexHeatmap)
library(circlize)
require(msigdbr)
library(cluster)
require(ggExtra)

source(here("src","plotting","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","Figure4")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import the raw data
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))

# Import host gene expression data ----
host_celltype_proportion_mat <- readRDS(here("data","raw","celltypeProportions_allTumors.rds"))

# Add bacterial counts per million sequencing reads, shannon diversity and genus richness to the metadata  ----
meta_test_df <- meta_combined_df %>%
  filter(Sample_ID %in% colnames(host_celltype_proportion_mat)) %>%
  dplyr::select(Sample_ID, Dataset, Batch, FastQ_libsize) %>% 
  mutate(Dataset = ifelse(str_detect(Dataset,"DKFZ"),"DKFZ",Dataset))

# Get total bacterial load (bacterial counts per million sequencing reads)
genus_counts_raw_mat <- readRDS(here("data","raw","rawScores_bulkRNAseq_genus.rds"))
# Get total bacterial scores and compute bacterial counts per million sequencing reads
total_bact_counts_df <- genus_counts_raw_mat[, meta_test_df$Sample_ID] %>%
  colSums() %>%
  enframe(name = "Sample_ID", value = "Bacteria_total") %>%
  left_join(., meta_test_df %>% dplyr::select(Sample_ID, FastQ_libsize)) %>%
  mutate(`Total bacteria (CPM)` = Bacteria_total / FastQ_libsize * 1e6)


# Import test resluts of celltype abundances vs genus abundances
test_results_celltypes_genera_df <- readRDS(here("data","results","celltypeProportion_vs_bacteria.rds")) %>% 
  mutate(genus = ifelse(genus == "Bacteria_CPM", "Total bacteria [CPM]", genus))
test_result_HCC_MetaAnalysis_df <- readRDS(here("data", "results", "HCC-meta-analysis_result_df.rds")) %>% filter(comparison == "Tumor_vs_Adj",Dataset == "Meta-analysis")



#* A: Heatmap of associations between bacterial features (genus abundances, total bacteria, Shannon div. and richness) and celltype proportions ----
# Define thresholds of p-value and adj. p-values for which significance indications should be plotted inside the hetamap tiles
p_threshold <- 0.05
fdr_threshold <- 0.2

# Generate matrices of effect sizes, pvalues and p-adj values
mat_effect_size <- test_results_celltypes_genera_df %>%
  dplyr::select(genus, celltype, effect.size) %>% 
  pivot_wider(names_from = genus, values_from = effect.size) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

mat_p <- test_results_celltypes_genera_df %>%
  dplyr::select(genus, celltype, p.val) %>% 
  pivot_wider(names_from = genus, values_from = p.val) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

mat_p_adj <- test_results_celltypes_genera_df %>%
  dplyr::select(genus, celltype, p.val_adj) %>% 
  pivot_wider(names_from = genus, values_from = p.val_adj) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

# keep only bacteria with at elast 1 p-value significant association with any celltype
column_levels <- colnames(mat_p)[colSums(mat_p<p_threshold) > 0]
row_levels <- rownames(mat_p) # keep all celltypes
# generate two groups for the columns: bacterial genera and "global" bacterial features (shannon, richness, CPM)
bac_diversity_features <- c("Total bacteria [CPM]", "Shannon diversity", "Genus richness")

# Group the bacterial features for splitting the matrix columns
column_split_vec <- ifelse(column_levels %in% bac_diversity_features, "A", "B")

# Reorder matrices
mat_effect_size <- mat_effect_size[row_levels,column_levels]
mat_p <- mat_p[row_levels,column_levels]
mat_p_adj <- mat_p_adj[row_levels,column_levels]

stopifnot(colnames(mat_effect_size) == colnames(mat_p) & colnames(mat_effect_size) == colnames(mat_p_adj))

# Normalize the score matrix for better visualization
mat_effect_size_rescaled <- f_rescale_effect_sizes(mat_effect_size)

# Define color mapping for the heatmap
range(mat_effect_size_rescaled)
max_val <- signif(max(abs(range(na.omit(mat_effect_size_rescaled)))), 1)
col_range <- c(-max_val, 0, max_val)
color_mapping <- circlize::colorRamp2(col_range, c("dodgerblue", "white", "tomato"))


# Generate column annotation:
#1) Indicate phylum to which the genus belongs
#2) Indicate the enrichment of the genus in HCC tumor tissue or adj. non-tumor

# Phylum annotation
# fetch phylum information for the genera

phylum_genus_df <- f_fetch_ncbi_taxonomy(tax_names = unique(column_levels)) %>%
  dplyr::rename(genus = query_tax)
  
# The analysed genera are from the following phyla:
phylum_genus_df$phylum %>% unique()
# rename the phylum to more common names
phylum_genus_df <-
  phylum_genus_df %>%
  mutate(phylum = case_when(
    phylum == "Actinomycetota" ~ "Actinobacteria",
    phylum == "Pseudomonadota" ~ "Proteobacteria",
    phylum == "Bacillota" ~ "Firmicutes",
    phylum == "Bacteroidota" ~ "Bacteroidetes"
  )) %>% 
  dplyr::select(genus, phylum)
phylum_colors_vec <- c("#009F4D", "#E40046","#307FE2", "#FFA300", "#8246AF", "#FFCD00", "#A8C700", "#BE5400", "#A8A99E")
names <- c(
  "Proteobacteria", "Actinobacteria", "Firmicutes", "Bacteroidetes",
  "Cyanobacteria", "Fusobacteria", "Acidobacteria", "Verrucomicrobia", "other"
)
names(phylum_colors_vec) <- names
phylum_colors_df <- phylum_genus_df %>%
  distinct() %>% 
  left_join(., phylum_colors_vec %>% enframe(name = "phylum", value = "HEX")) %>% 
  filter(!is.na(HEX))

# generate the annotation vector
anno_vec_phylum <- phylum_genus_df$phylum[match(column_levels, phylum_genus_df$genus)]

# Annotation of HCC association
HCC_assoc_df <- test_result_HCC_MetaAnalysis_df %>%
  mutate(HCC_association = case_when(
    p.val_adj < 0.2 & effect.size > 0 ~ "HCC",
    p.val_adj < 0.2 & effect.size < 0 ~ "Adj. non-tumor",
    TRUE ~ NA_character_
  ))
anno_vec_HCCassoc <- HCC_assoc_df$HCC_association[match(column_levels, HCC_assoc_df$bacteria)]

col_anno <- columnAnnotation(
    "HCC association" = anno_vec_HCCassoc,
    "Phylum" = anno_vec_phylum,
    #"Gild" = anno_vec_gild,
    gp = gpar(col = "black"),
    simple_anno_size = unit(2, "mm"),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 8),    
    col = list(
      "HCC association" = c("HCC" = "#62008c", "Adj. non-tumor" = "#CBA3D8"),
      "Phylum" = phylum_colors_df %>% select(phylum,HEX) %>% deframe()#,
     # "Gild" = gild_colors
    ),
    na_col = "white",
    annotation_legend_param = list(
      legend_gp = gpar(col = "black", lwd = 1),
      border = "black", direction = "vertical"
    )
  )

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
colLabs <- as.expression(lapply(column_levels, function(a) bquote(bolditalic(.(a)))))
rowLabs <- as.expression(lapply(row_levels, function(a) bquote(bold(.(a)))))

# plot the heatmap
hm <- ComplexHeatmap::Heatmap(mat_effect_size_rescaled,
  bottom_annotation = col_anno,
  
  cluster_rows = T,
  cluster_columns = T,

  # Column style
  column_names_rot = 45, column_names_side = "bottom", column_names_centered = F,
  column_names_gp = gpar(fontsize = 8),
  column_labels = colLabs,
  column_split = column_split_vec, # indicate whether to which group the feature belongs

  # Row style
  row_dend_side = "left",
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
  column_title=NULL,


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

pt_A <- draw(hm,annotation_legend_list = leg_combined,heatmap_legend_side = "right", annotation_legend_side = "right")
pdf(file = file.path(save_fig_folder, "A_Heatmap_CelltypeProportionsGenusAbundance.pdf"), width = 12, height = 4)
pt_A
dev.off()


#* B: selected scatterplots of celltype proportions vs bacterial rel. abundances ----


f_scatter_bac_vs_celltypes <- function(
    meta_df,
    bacterial_feature,
    celltype,
    genus_RelAbundance_corrected_mat = genus_RelAbundance_corrected_mat,
    host_celltype_proportion_mat = host_celltype_proportion_mat) {
  #* Takes rel. abundance matrix of genera and celltypes. Plots selected combination of genus vs celltype abundance. ----
  # Computes linear mixed model and returns p-value and regression coefficient

  # generate plot_df by combining bacterial relative abundaces and celltype proportions
  plot_df <- inner_join(
    genus_RelAbundance_corrected_mat[bacterial_feature, ] %>%
      enframe(name = "Sample_ID", value = "bac_rel") %>%
      inner_join(., meta_df),
    host_celltype_proportion_mat[celltype, ] %>%
      enframe(name = "Sample_ID", value = "ct_rel") %>%
      inner_join(., meta_df)
  ) %>%
    mutate(
      l10_bacteria = log10(bac_rel + 1e-5),
      l10_celltype = log10(ct_rel + 1e-5)
    ) %>%
    group_by(Dataset) %>%
    mutate(N = length(unique(Sample_ID))) %>%
    ungroup() %>%
    mutate(Dataset_label = paste0(Dataset, " (N=", N, ")"))

  y_axis_range <- seq(-5, 0, 1)
  x_axis_range <- seq(-5, 0, 1)
  x_name <- paste0(bacterial_feature, " rel. ab. [log10]")
  if (min(plot_df$l10_bacteria,na.rm = T) < min(x_axis_range) | max(plot_df$l10_bacteria,na.rm = T) > max(x_axis_range)) {
    # modify x-axis range if the range of the bacterial feature is outside the default range
    x_axis_range <- seq(0, ceiling(max(plot_df$l10_bacteria,na.rm = T)), length.out = 5) # relevant for total bacterial load
    x_name <- paste0(bacterial_feature, " [log10]")
  }

  # compute a linear-mixed model to test for association between the bacterial feature and the celltype proportion
  model <- lmerTest::lmer(l10_celltype ~ l10_bacteria + (1 | Batch), plot_df)
  # r_squared_conditional <- performance::r2_nakagawa(model)$R2_conditional
  # predict log-celltype proportions from the model
  pred_df <- data.frame(l10_bacteria = seq(range(x_axis_range)[1], range(x_axis_range)[2], length.out = 100))
  pred_df$l10_celltype <- predict(model, newdata = pred_df, re.form = NA) # use re.form = NA to only use fixed effects

  # Get p-value and regression coefficient of the linear mixed model
  p_value <- coef(summary(model))[2, 5]
  effect_size <- coef(summary(model))[2, 1]

  # Create labels
  p_label <- paste0("p-value = ", round(p_value, 3)) # round to 3 decimal places
  e_label <- paste0("regr. coeff. = ", round(effect_size, 3)) # round to 3 decimal places
  # r_sq_label <- paste0("R^2 = ",round(r_squared_conditional,2))
  annot_label <- paste0(p_label, "\n", e_label, "\n") # , r_sq_label)

  # Define name of x, y axis and plot title
  
  y_name <- paste0(celltype, " rel. ab. [log10]")
  title <- paste0(celltype, " vs\n", bacterial_feature)


  pt_tmp <- plot_df %>%
    ggplot(aes(x = l10_bacteria, y = l10_celltype, fill = Dataset_label)) +
    geom_point(pch = 21, size = 2.5, alpha = 0.75) +
    geom_line(data = pred_df, aes(x = l10_bacteria, y = l10_celltype), color = "black", inherit.aes = F, lwd = 1) +
    scale_y_continuous(limits = range(y_axis_range) * 1.05, breaks = y_axis_range, name = y_name) +
    scale_x_continuous(expand = c(0, 0), limits = range(x_axis_range) * 1.05, breaks = x_axis_range, name = x_name) +
    theme_paper +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    # geom_abline(slope=1,intercept=0)+
    ggtitle(title) +
    geom_text(
      aes(
        label = annot_label,
        x = min(x_axis_range) + 0.1, y = Inf, hjust = 0, vjust = 1.15
      ),
      color = "black", size = 4, check_overlap = TRUE
    )

  # Add marginal histograms
  pt_combined <- ggExtra::ggMarginal(pt_tmp + theme(legend.position = "none"), groupFill = FALSE, col = "black", fill = "grey", size = 7)
  pt_combined_legend <- ggExtra::ggMarginal(pt_tmp + theme(legend.position = "right"), groupFill = FALSE, col = "black", fill = "grey", size = 7)
  res_list <- list(
    pt = pt_combined,
    pt_legend = pt_combined_legend
  )

  return(res_list)
}

# Define bacteira and celltypes to be plotted
genus_selection_vec <- c("Cloacibacterium", "Tepidimonas", "Flavobacterium")
ct1 <- sort(unique(test_results_celltypes_genera_df$celltype))[15] # Pro-Inf macrophages
ct2 <- sort(unique(test_results_celltypes_genera_df$celltype))[2] # Bcells
ct3 <- sort(unique(test_results_celltypes_genera_df$celltype))[15] # Pro-Inf macrophages
celltype_selection_vec <- c(ct1, ct2, ct3)
i <- 1
meta_df <- meta_test_df %>% dplyr::select(Sample_ID,Dataset,Batch)
for(i in seq(1,length(genus_selection_vec))){
  bac <- genus_selection_vec[i]
  ct <- celltype_selection_vec[i]
  c_list <- f_scatter_bac_vs_celltypes(meta_df = meta_df, bacterial_feature = bac, celltype = ct,genus_RelAbundance_corrected_mat = genus_RelAbundance_corrected_mat, host_celltype_proportion_mat = host_celltype_proportion_mat)

  out_name <- file.path(save_fig_folder,paste0("B_Scatterplot_",bac,"_vs_",ct,".pdf"))
  ggsave(c_list$pt, file = out_name, width = 4.5, height = 4.5)
  if(i == 1){
    out_name <- file.path(save_fig_folder,paste0("B_z Only for Legend.pdf"))
    ggsave(c_list$pt_legend, file = out_name, width = 4.5, height = 4.5)
  }
}

# Plot CD8 TCells vs total bacterial load

# generate matrix with total bacterial load
bac_cpm_mat <- total_bact_counts_df %>%
  dplyr::select(Sample_ID, `Total bacteria (CPM)`) %>%
  mutate(rowN = "Total bacteria [CPM]") %>%
  pivot_wider(names_from = Sample_ID, values_from = `Total bacteria (CPM)`) %>%
  column_to_rownames("rowN") %>%
  as.matrix(drop = F)

# select CD8 TCells
cd8t <- sort(unique(test_results_celltypes_genera_df$celltype))[4]
# generate the plot
c_list <- f_scatter_bac_vs_celltypes(
  meta_df = meta_df, bacterial_feature = "Total bacteria [CPM]", celltype = cd8t,
  genus_RelAbundance_corrected_mat = bac_cpm_mat, # select for the bacterial matrix now the total bacterial load
  host_celltype_proportion_mat = host_celltype_proportion_mat
)
out_name <- file.path(save_fig_folder,paste0("B_Scatterplot_TotalBacteria_vs_",cd8t,".pdf"))
ggsave(c_list$pt, file = out_name, width = 4.5, height = 4.5)

