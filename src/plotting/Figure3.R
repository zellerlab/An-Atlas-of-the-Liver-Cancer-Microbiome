#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Figure 3
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

source(here("src","plotting","functions_plotting.R"))

# Load colors for the individual groups
parameters <- yaml::yaml.load_file(here("src", "parameters.yml"))
group_colors <- unlist(parameters$plotting)

# Define folder to store plots
save_fig_folder <- here("figures","Figure3")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import the testing result of total bacterial loads and HCC inflammation status vs gene expression and inferred pathway enrichments

test_result_genes_vs_sampleMetrics_df <- readRDS(here("data", "results", "genes_vs_sampleMetrics_df.rds"))
test_result_pathwayActivities_sampleMetrics_df <- readRDS(here("data", "results", "pathwayActivities_bySampleMetrics_df.rds"))
test_result_pathwayActivities_genera_df <- readRDS(here("data", "results", "pathwayActivities_byGenus_df.rds"))
test_result_HCC_MetaAnalysis_df <- readRDS(here("data", "results", "HCC-meta-analysis_result_df.rds")) %>% filter(comparison == "Tumor_vs_Adj",Dataset == "Meta-analysis")


#* A Forest plots of pathway activities by totla bacterial load and HCC inflammation status ----
shape_vector <- c(
    "All" = 15,
    "ISMMS-IDIBAPS" = 1,
    "INSERM" = 5,
    "DKFZ" = 4,
    "TCGA" = 2
  )
levels(shape_vector) <- names(shape_vector)

sig_color_vec <- c("sigPos" = "#E40046", "sigNeg" = "#307FE2")
sig_color_label <- c("sigPos" = "Up-regulated", "sigNeg" = "Down-regulated")

# Prepare the data for plotting
plot_df <- test_result_pathwayActivities_sampleMetrics_df %>%
  ungroup() %>%
  transmute(
    pathway = str_remove(source, "HALLMARK_"),
    pathway = str_to_title(str_replace_all(pathway, "\\_", " ")),
    feature = case_when(
      feature == "total_Bacteria" ~ "Total bacteria [CPM]",
      feature == "Inflammed_vs_nonInflamed" ~ "HCC inflammation",
    ), dataset, score, p.val_adj
  ) %>% 
  mutate(feature = factor(feature, levels = c("Total bacteria [CPM]", "HCC inflammation")))


# Define the order of the pathways based on pathway association with total bacterial load:
# Keep only pathways with FDR-significant association with total bacterial load in the meta-analysis setting
pathway_levels <-
  plot_df %>%
  filter(dataset == "All", str_detect(feature, "bacteria")) %>%
  filter(p.val_adj < 0.05) %>%
  arrange(score) %>%
  pull(pathway)

# Add levels for pathways and indicate significant directions
plot_df <- plot_df %>%
  filter(pathway %in% pathway_levels) %>% #remove non-significant associations
  mutate(pathway = factor(pathway, levels = pathway_levels)) %>% 
  mutate(sig_level = case_when(
    p.val_adj < 0.05 & score > 0 ~ "sigPos",
    p.val_adj < 0.05 & score < 0 ~ "sigNeg",
    TRUE ~ "n.s."
  ))
plot_df$dataset %>% unique()
# generate the plot
pt_A <- plot_df %>%
  ggplot(aes(x = score, y = pathway, shape = dataset, color = sig_level)) +
  geom_vline(xintercept = 0) +
  theme_paper +
  geom_point(data = plot_df %>% filter(dataset == "DKFZ"), show.legend = F) +
  geom_point(data = plot_df %>% filter(dataset == "INSERM"), show.legend = F) +
  geom_point(data = plot_df %>% filter(dataset == "TCGA"), show.legend = F) +
  geom_point(data = plot_df %>% filter(dataset == "ISMMS-IDIBAPS"), show.legend = F) +
  geom_point(data = plot_df %>% filter(dataset == "All"), show.legend = T, size = 2.5) +
  scale_shape_manual(values = shape_vector, labels = sig_color_label) +
  scale_color_manual(values = sig_color_vec, labels = sig_color_label) +
  xlab("Enrichment score") +
  labs(color = "q<0.05", shape = "Cohort") +
  facet_wrap(~feature, nrow = 1)+
  theme(axis.title.y = element_blank(),legend.position = "left")

ggsave(pt_A,filename = here(save_fig_folder,"A_ForestPlot_PathwayActivities.pdf"), width = 9, height = 8)

#* B: Forest plot of differential gene expression based on total bacterial load ----
col_vec <- c(
  "inflammatory" = "#fa7c0c",
  "metabolic" = "#1b9e77",
  "both" = "#C7A2FF"
)

col_labs <- c("inflammatory" = "in inflammatory pathways", "metabolic" = "in metabolic pathways", "both" = "in both pathways")
levels(col_vec) <- names(col_vec)
levels(col_labs) <- col_labs

# Manually classify the pathways from the forest plot into inflammatory and metabolic
inflammatory_pathways <-
  c(
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  )

test_result_pathwayActivities_sampleMetrics_df$source %>% unique() %>% sort()

metabolic_pathways <-
  c("HALLMARK_HEME_METABOLISM",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_ADIPOGENESIS",
    "HALLMARK_PEROXISOME",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_BILE_ACID_METABOLISM",
    "HALLMARK_XENOBIOTIC_METABOLISM"
  )

pathway_info_df <- bind_rows(
  inflammatory_pathways %>% enframe() %>% transmute(souurce = value, type = "inflammatory"),
  metabolic_pathways %>% enframe() %>% transmute(souurce = value, type = "metabolic")
)

# assign which genes belong to the inflammatory or metabolic hallmark pathways
hallmarks_msigdb <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
inflammatory_pathways_genes <- hallmarks_msigdb %>% filter(gs_name %in% inflammatory_pathways) %>% pull(gene_symbol) %>% unique()
metabolic_pathways_genes <- hallmarks_msigdb %>% filter(gs_name %in% metabolic_pathways) %>% pull(gene_symbol) %>% unique()
intersect(inflammatory_pathways_genes,metabolic_pathways_genes)
# Generate the plot df: Only total bacterial load vs gene expression of the meta-analysis
test_result_genes_vs_sampleMetrics_df$comparison %>% unique()

plot_df <- test_result_genes_vs_sampleMetrics_df %>%
  filter(str_detect(comparison, "Bacteria"), Dataset == "All") %>%
  transmute(gene, effect.size, p.val_adj, N_samples) %>%
  mutate(
    pathway_association = case_when(
      (gene %in% inflammatory_pathways_genes) & (gene %in% metabolic_pathways_genes) ~ "both",
      gene %in% inflammatory_pathways_genes ~ "inflammatory",
      gene %in% metabolic_pathways_genes ~ "metabolic",
      TRUE ~ "other"
    )
  ) %>% 
  mutate(pathway_association = factor(pathway_association, levels = c("inflammatory", "metabolic", "both", "other")))
pt_B <- plot_df %>%
  ggplot(aes(x = effect.size, y = -log10(p.val_adj))) +
  geom_point(data = plot_df %>% filter(pathway_association == "other"), fill = "lightgrey", color = "black", pch = 21, alpha = 0.5, size = 4) +
  geom_point(data = plot_df %>% filter(pathway_association != "other"), aes(fill = pathway_association), color = "black", pch = 21, alpha = 0.7, size = 4) +
  scale_fill_manual(values = col_vec, labels = col_labs) +
  geom_hline(yintercept = -log10(0.05)) +
  # label highly significant genes of other pathways
  ggrepel::geom_text_repel(data = plot_df %>% filter(pathway_association == "other", -log10(p.val_adj) > 8 | (effect.size > 0 & -log10(p.val_adj) > 6) | effect.size > 0.4), aes(label = gene), color = "darkgrey", size = 3)+

  # label highly significant genes of inflammatory pathways and metabolic pathways
  ggrepel::geom_label_repel(data = plot_df %>% filter(pathway_association != "other", -log10(p.val_adj) > 2), aes(label = gene, fill = pathway_association), color = "black", alpha = 0.85, show.legend = F, size = 3) +
  xlim(c(-0.5, 0.5)) +
  ylim(c(0, 15)) +
  theme_paper +
  ggtitle("Gene expression differnces depending on total bacterial load") +
  labs(fill = "") +
  ylab("-log10 q-value") +
  xlab("Enrichment effect size")

ggsave(pt_B, filename = file.path(save_fig_folder, "B_Volcano_TotalBacteriaVSgeneExpression.pdf"), width = 7, height = 8)


#* C: Heatmap of individual genera vs pathway associations ----
# Define thresholds of p-value and adj. p-values for which significance indications should be plotted inside the hetamap tiles
p_threshold <- 0.05
fdr_threshold <- 0.001

# Generate plot DF: Filter hallmark pathways with at least 10 genus associations of p.val_adj < 0.001
plot_df <- test_result_pathwayActivities_genera_df %>%
  ungroup() %>%  
  transmute(
    pathway = str_remove(source, "HALLMARK_"),
    pathway = str_to_title(str_replace_all(pathway, "\\_", " ")),
    genus = feature,
    score,
    p.val, p.val_adj
  ) %>%
  group_by(pathway) %>%
  mutate(N_sig = sum(p.val_adj < fdr_threshold)) %>%
  ungroup() %>% 
  filter(N_sig >= 10)

# Prepare the matrices of scores and p-values for plotting in the heatmap
mat_scores <- plot_df %>%
  dplyr::select(pathway, genus, score) %>% 
  pivot_wider(names_from = genus, values_from = score) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

mat_p <- plot_df %>%
  dplyr::select(pathway, genus, p.val) %>% 
  pivot_wider(names_from = genus, values_from = p.val) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

mat_p_adj <- plot_df %>%
  dplyr::select(pathway, genus, p.val_adj) %>% 
  pivot_wider(names_from = genus, values_from = p.val_adj) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Perform a clustering of the genus:pathway associations based on the scores
row_clust <- diana(mat_scores, metric = "ward")
row_levels <- rownames(mat_scores)
column_clust <- diana(t(mat_scores), metric = "ward")
column_levels <- colnames(mat_scores)

# Reorder all matrices based on the column and row levels
mat_scores <- mat_scores[row_levels, column_levels]
mat_p <- mat_p[row_levels, column_levels]
mat_p_adj <- mat_p_adj[row_levels, column_levels]

stopifnot(colnames(mat_scores) == colnames(mat_p) & colnames(mat_scores) == colnames(mat_p_adj))

# Normalize the score matrix for better visualization
# Important to note that the rescaling only affects the colors in the heatmap, not the clustering
mat_scores_rescaled <- f_rescale_effect_sizes(mat_scores)

# Define color mapping for the heatmap
range(mat_scores_rescaled)
max_val <- signif(max(abs(range(na.omit(mat_scores_rescaled)))), 1)
col_range <- c(-max_val, 0, max_val)
color_mapping <- circlize::colorRamp2(col_range, c("dodgerblue", "white", "tomato"))

# Generate the heatmap annotations


# Generate column annotation:
#1) Indicate phylum to which the genus belongs
#2) Indicate the enrichment of the genus in HCC tumor tissue or adj. non-tumor

phylum_colors_df_path <- here("data","processed","phylum_colors_df.tsv")
# check if the phylum colors df exists, if not create it
if (file.exists(phylum_colors_df_path)) {
  phylum_colors_df <- read_tsv(phylum_colors_df_path)
} else {
  # generate the phylum colors df by fetching the genus from NCBI and assigning the phylum colors
  phylum_genus_df <- f_fetch_ncbi_taxonomy(tax_names = unique(test_result_df$tax)) %>%
    dplyr::rename(genus = query_tax)

  # The analysed genera are from the following phyla:
  unique(phylum_genus_df$phylum)

  # rename the phylum to more common names
  phylum_genus_df <-
    phylum_genus_df %>%
    mutate(phylum = case_when(
      phylum == "Actinomycetota" ~ "Actinobacteria",
      phylum == "Pseudomonadota" ~ "Proteobacteria",
      phylum == "Bacillota" ~ "Firmicutes",
      phylum == "Bacteroidota" ~ "Bacteroidetes"
    ))

  phylum_colors_vec <- c("#009F4D", "#E40046", "#307FE2", "#FFA300", "#8246AF", "#FFCD00", "#A8C700", "#BE5400", "#A8A99E")
  names <- c(
    "Proteobacteria", "Actinobacteria", "Firmicutes", "Bacteroidetes",
    "Cyanobacteria", "Fusobacteria", "Acidobacteria", "Verrucomicrobia", "other"
  )
  names(phylum_colors_vec) <- names

  # Add the colors to the phylum_genus_df
  phylum_colors_df <- phylum_genus_df %>%
    dplyr::select(phylum, genus) %>%
    distinct() %>%
    left_join(., phylum_colors_vec %>% enframe(name = "phylum", value = "HEX"))
  write_tsv(phylum_colors_df, here("data", "processed", "phylum_colors_df.tsv"))
}

# generate the annotation vector
anno_vec_phylum <- phylum_colors_df$phylum[match(column_levels, phylum_colors_df$genus)]

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
leg_title <- "Association strength\n(normalized)"
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
hm <- ComplexHeatmap::Heatmap(mat_scores_rescaled,
  bottom_annotation = col_anno,
  
  cluster_rows = row_clust,
  cluster_columns = column_clust,

  # Column style
  column_names_rot = 45, column_names_side = "bottom", column_names_centered = F,
  column_names_gp = gpar(fontsize = 8),
  column_labels = colLabs,

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
pt_C <- draw(hm, annotation_legend_list = leg_combined)
pdf(file = file.path(save_fig_folder, "C_Heatmap_PathwayActivityGenusAbundance.pdf"), width = 17, height = 10)
draw(pt_C)
dev.off()
