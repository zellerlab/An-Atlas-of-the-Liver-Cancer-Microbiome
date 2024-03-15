#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Figure 2
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
save_fig_folder <- here("figures","Figure2")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

# Import Testing results
test_result_df <- readRDS(here("data", "results", "HCC-meta-analysis_result_df.rds")) %>% dplyr::rename(p.val_lm = p.val,tax = bacteria,p.val_adj_lm = p.val_adj)

# import raw data: 
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))

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

#*2B: Scatterplot of mean relative abudnances in bulk RNA-seq datasets and 5R16S seq datasets ----
threshold_for_prev_RNAseq <- 1e-3 #relative abundance in the RNA-Seq cohort for which a genus is considered prevalent
threshold_for_prev_5R16S <- 0 #same in the 5R16S cohort - In the preprocessing small values (<1e-4) were floored, so the threshold can be 0

tumor_abundance_prevalence_df <-
  genus_RelAbundance_corrected_mat %>%
  as_tibble(rownames = "genus") %>%
  gather(-genus, key = "Sample_ID", value = "rel") %>%
  left_join(., meta_combined_df %>% dplyr::select(Sample_ID, Tissue_type, Dataset)) %>%
  filter(Tissue_type == "Tumor") %>%
  mutate(thresh_for_prev = ifelse(str_detect(Dataset, "16S"), threshold_for_prev_5R16S, threshold_for_prev_RNAseq)) %>%  
  mutate(Assay = ifelse(str_detect(Dataset,"16S"),"5R16S","RNAseq")) %>%
  group_by(Assay) %>%
  mutate(N_samples = length(unique(Sample_ID))) %>%
  group_by(Assay, genus) %>%
  mutate(
    prevalence = (sum(rel > thresh_for_prev)) / N_samples,
    mean_abundance = mean(rel)
    #mean_abundance = (10^(mean(log10(rel+1e-5))))-1e-5
  ) %>%
  ungroup() %>%
  dplyr::select(genus, Assay, N_samples, prevalence, mean_abundance) %>%
  distinct() %>%
  arrange(genus)

meanAB_df <- tumor_abundance_prevalence_df %>%
  dplyr::select(genus, Assay, mean_abundance) %>%
  pivot_wider(names_from = Assay, values_from = mean_abundance) %>% 
  left_join(.,phylum_colors_df %>% dplyr::select(genus,phylum)) %>% 
  mutate(l10_RNAseq = log10(RNAseq+1e-5),
         l10_5R16S = log10(`5R16S`+1e-5))

pearson_R <- cor(meanAB_df$l10_RNAseq,meanAB_df$l10_5R16S,method = "spearman")
spearman_R <- cor(meanAB_df$l10_RNAseq,meanAB_df$l10_5R16S,method = "spearman")

pt_B <- meanAB_df %>%
  ggplot(aes(y = l10_5R16S, x = l10_RNAseq)) +
  geom_point(pch = 21,size = 2,alpha = 0.75,aes(fill = phylum)) +
  geom_abline(slope = 1, intercept = 0) +
  tune::coord_obs_pred() +
  theme_paper +
  scale_fill_manual(values = phylum_colors_df %>% dplyr::select(phylum,HEX) %>% deframe())+
  annotate("text",
    x = -5, y = 0, hjust = 0, vjust = 0.75, size = 3,
    # annotate("text", x = -2, y = -3.7, hjust = 0, vjust = 0.75,size = 3,
    label = paste("Pearson's R = ", signif(pearson_R, 2), "\nSpearman's R = ", signif(spearman_R, 2))
  ) +
  ylab("16S-Seq (N=110)\nMean rel. abundance [log10]") +
  xlab("RNA-Seq (N=750)\nMean rel. abundance [log10]") +
  #ggrepel::geom_text_repel(data = meanAB_df %>% filter(l10_RNAseq > -3 & l10_5R16S < -3.5),aes(label = genus), size = 2, nudge_x = 0.1, nudge_y = 0.1)+
  labs(fill = "Phylum")+
  ggtitle("Comparison of genus abundances\nbetween RNA-Seq and 16S-Seq")

ggsave(pt_B,filename = file.path(save_fig_folder,"B_Scatterplot_bulkRNAseq_vs_5R16S.pdf"), width = 5,height = 5)




#* d-f: Volcano plots ----
# Define size and ranges of the points in the volcano plots
size_definition <- scale_size_continuous(
        name = "Prevalence",
        range = c(1.5, 4),
        limits = c(0, 1),
        breaks = c(0.25, 0.5, 0.75, 1)
)
# Define y-range for all volcano plots
man_y_breaks <- c(-log10(c(0.05,0.01,0.001,1e-4,1e-5)))

#D: Viral vs non-viral HCCs
plot_df <- test_result_df %>% filter(comparison == "Viral_vs_nonViral_HCC", Dataset == "Meta-analysis")
range(plot_df$effect.size)
xBreaks <- c(-0.2,0.1,0,0.1,0.2,0.3)
xLims <- c(-0.23,0.35)

pt_D <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks,clean_tax_names = FALSE,
) +
  scale_fill_manual(values = group_colors) +
  size_definition


#E ALD vs MAFLD
plot_df <- test_result_df %>% filter(comparison == "ALD_vs_MAFLD", Dataset == "Meta-analysis")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.4, 0.2, 0.1), 1)
xLims <- c(-0.4,0.25)

pt_E <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks,clean_tax_names = FALSE,
) +
  scale_fill_manual(values = group_colors) +
  size_definition

#F HCV vs HB
plot_df <- test_result_df %>% filter(comparison == "HBV_vs_HCV", Dataset == "Meta-analysis")
range(plot_df$effect.size)
xBreaks <- round(seq(-0.3, 0.5, 0.1), 1)
xLims <- c(-0.3,0.5)

pt_F <- f_plot_volcano(
  plot_df = plot_df, xBreaks = xBreaks, xLims = xLims,
  man_y_breaks = man_y_breaks,clean_tax_names = FALSE,
) +
  scale_fill_manual(values = group_colors) +
  size_definition


# Plot volcanos
w <- 4.5
h <- 5  

ggsave(pt_D, filename = file.path(save_fig_folder,"D_Volcano_ViralVSnonViralHCC.pdf"), width = w,height = h)
ggsave(pt_E, filename = file.path(save_fig_folder,"E_Volcano_ALDvsMAFLDHCC.pdf"), width = w,height = h)
ggsave(pt_F, filename = file.path(save_fig_folder,"F_Volcano_HBVvsHCVHCC.pdf"), width = w,height = h)

#* C: Forest plot of tumor vs normal enrichments and heatmaps of prevalence, abundance and contamination likelihood ----
table(test_result_df$comparison,test_result_df$Dataset)

# Select tumor vs adj. non-tumor testing results of the meta-analysis and the datasets individually. Ignore DKFZ RNA-seq for individual plotting since there is only one non-tumor sample
shape_vector <- c(
    "Meta-analysis" = 22,
    "INSERM" = 5,
    "TCGA" = 2,  
    "ISMMS-IDIBAPS" = 1,    
    "DKFZ_5R16S" = 8
  )
levels(shape_vector) <- names(shape_vector)
shape_labels <-  c("All datasets [+/- 95 CI]",
                                  "INSERM",
                                  "TCGA",
                                  "ISMMS-IDIBAPS",                                  
                                  "DKFZ_5R16S")
names(shape_labels) <- names(shape_vector)

sig_color_vec <- c("sigPos" = "#E40046", "sigNeg" = "#307FE2", "n.s." = "#707372")
levels(shape_vector) <- names(shape_vector)
fdr_threshold <- 0.2

plot_df <- test_result_df %>%
  filter(comparison == "Tumor_vs_Adj", Dataset != "DKFZ_RNAseq") %>%
  mutate(sig_dir = case_when(
    p.val_adj_lm < fdr_threshold & effect.size > 0 ~ "sigPos",
    p.val_adj_lm < fdr_threshold & effect.size < 0 ~ "sigNeg",
    TRUE ~ "n.s."
  )) %>% 
  left_join(.,sig_color_vec %>% enframe(name = "sig_dir", value = "HEX")) %>% 
  mutate(sig_dir = factor(sig_dir, levels = c("sigPos", "n.s.", "sigNeg")))  %>% 
  mutate(Dataset = factor(Dataset, levels = c("Meta-analysis", "INSERM","TCGA", "ISMMS-IDIBAPS", "DKFZ_5R16S")))

# Sort bacteria according to meta-analysis significance and effect size
tax_levels <- plot_df %>%
  filter(Dataset == "Meta-analysis") %>%
  arrange(sig_dir, -effect.size) %>% 
  pull(tax)

plot_df <- plot_df %>% 
  mutate(tax = factor(tax, levels = tax_levels))


# Define max value for y-axis
max <- max(abs(plot_df$effect.size))
xRange <- c(-max, max)

pt_C1 <- plot_df %>%
  ggplot(aes(x = tax, y = effect.size, shape = Dataset, fill = sig_dir)) +
  geom_hline(yintercept = 0) +
  theme_paper +
  geom_point(data = plot_df %>% filter(Dataset == "Meta-analysis"), show.legend = T, size = 2.5, color = "black") +
  geom_linerange(data = plot_df %>% filter(Dataset == "Meta-analysis"), aes(ymin = lower95CI, ymax = upper95CI), color = "black") +
  geom_point(data = plot_df %>% filter(Dataset == "DKFZ_5R16S"), show.legend = F, fill = NA, color = "black", fill = "#707372") +
  geom_point(data = plot_df %>% filter(Dataset == "INSERM"), show.legend = F, fill = NA, color = "black", fill = "#707372") +
  geom_point(data = plot_df %>% filter(Dataset == "TCGA"), show.legend = F, fill = NA, color = "black", fill = "#707372") +
  geom_point(data = plot_df %>% filter(Dataset == "ISMMS-IDIBAPS"), show.legend = F, fill = NA, color = "black", fill = "#707372") +
  scale_y_continuous(limits = xRange, position = "right", name = "Enrichment effect size") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.99, vjust = 0.5, size = 8),
    legend.position = "right",
    axis.title.x = element_blank()
  ) +
  scale_shape_manual(values = shape_vector,labels = shape_labels) +
  scale_fill_manual(values = sig_color_vec) +
  guides(fill = "none", shape = guide_legend(reverse = TRUE))

ggsave(pt_C1,filename = file.path(save_fig_folder,"C1_ForestPlot_TumorVSAdjacent.pdf"), width = 13,height = 4)


#* C2: heatmaps of prevalence, abundance and contamination likelihood ----
hmap_parent_df <- genus_RelAbundance_corrected_mat %>%
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
  mutate(Dataset_lab = paste0(Dataset," (N=",N_samples,")")) %>% 
  mutate(Dataset = factor(Dataset,levels = c("DKFZ_5R16S","ISMMS-IDIBAPS","INSERM","TCGA","DKFZ_RNAseq"))) %>% 
  arrange(Dataset) %>% 
  mutate(Dataset_lab = as_factor(Dataset_lab),genus = factor(genus,levels = tax_levels))  

# Generate Heatmap of meanABundance
meanAB_mat <- hmap_parent_df %>%
  dplyr::select(genus, Dataset_lab, mean_abundance) %>%
  arrange(genus) %>%
  pivot_wider(names_from = genus, values_from = mean_abundance) %>%
  arrange(Dataset_lab) %>%
  column_to_rownames(var = "Dataset_lab") %>% 
  as.matrix()

col_fun_abundance = colorRamp2(c(-3.5,-0.8), c("white", "#006d65"))

# Add an empty annotation block to manually add the forest plot (C1) in illustrator
colAnno_Forest = columnAnnotation(forest = anno_empty(border = T,height = unit(4, "cm")))#,
                                # volcano = anno_empty(border = T,height = unit(4, "cm")))
hm_meanAB <- Heatmap(log10(meanAB_mat),name="mat",cluster_rows = F,cluster_columns = F,col = col_fun_abundance,top_annotation = colAnno_Forest,
                height = unit(2, "cm"),

                heatmap_legend_param = list(
                  title = "Rel. abundance", #at = linspace(-abs(round(min(mat_vals),0)),(round(max(mat_vals),0)),5),
                  legend_height = unit(3, "cm"),
                  at = seq(-4,0,1),
                  #at = 10^(seq(-4,0,1)),
                  labels = c("0.0001","0.001","0.01","0.1",1),
                  border = "black"),
                
                # appearance of labels
                column_names_rot = 90, column_names_side = "bottom",column_names_gp = gpar(fontsize = 8),                
                row_names_gp = gpar(fontsize = 8),row_names_side = c("right"),
                
                                
                cluster_row_slices = F,
                cluster_column_slices = F,
                row_gap = unit(c(3), "mm"),
                column_gap = unit(c(3), "mm"),
                border = T,
                rect_gp = gpar(col = "black", lwd = 1)
  )
  
# Generate heatmap of prevalence with an annotation for contaminant likelihood
contaminant_likelihood_df <- read_tsv(here("data","raw","contamination_likelihood_df.tsv")) %>% 
  dplyr::select(genus,`Contaminant likelohood`) %>% #contamination_likelihood_df
  #dplyr::select(genus,`Contaminant likelihood`) %>% 
  mutate(genus = factor(genus,levels = tax_levels)) %>% 
  mutate(rowN = "Contaminant\nlikelihood")
cont_likelihood_mat <- contaminant_likelihood_df %>%
  arrange(genus) %>%
  #pivot_wider(names_from = genus, values_from = `Contaminant likelihood`) %>%
  pivot_wider(names_from = genus, values_from = `Contaminant likelohood`) %>% #typo in contamination_likelihood_df
  column_to_rownames("rowN") %>%
  as.matrix(drop = F)

prev_mat <- hmap_parent_df %>%
  dplyr::select(genus, Dataset_lab, prevalence) %>%
  arrange(genus) %>%
  pivot_wider(names_from = genus, values_from = prevalence) %>%
  arrange(Dataset_lab) %>%
  column_to_rownames(var = "Dataset_lab") %>% 
  as.matrix()

col_fun_prev <- circlize::colorRamp2(c(0, 0.5, 1), c("white", "#ffb006", "#783d0c"))
bacLabs <- as.expression(lapply(colnames(prev_mat), function(a) bquote(bolditalic(.(a)))))

# Define colors of genera based on tumor vs non-tumor enrichment: Red (HCC enriched), Black (n.s.) or BLue (Adj. non-tumor enriched)
bacteria_label_colors <- pt_C1$data %>% filter(Dataset == "Meta-analysis") %>% mutate(HEX = ifelse(sig_dir == "n.s.","black",HEX)) %>% select(tax,HEX) %>% arrange(tax) %>% deframe()


# Generate heatmap annotation with contamination likelihood
cont_likelihood_mat
contaminant_color_vec <- c("0" = "#30E530",
                           "1" = "#4CBF4C",
                           "2" = "#5C995C",
                           "3" = "#5B725B",
                           "4" = "#4C4C4C")

cont_Anno <- columnAnnotation(
  "Contaminant\nlikelihood" = as.vector(cont_likelihood_mat),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 8),
  col = list("Contaminant\nlikelihood" = contaminant_color_vec),
      annotation_legend_param = list(
      legend_gp = gpar(col = "black", lwd = 1),
      at = seq(0,4,1),labels = c("0","1","2","3","4"),
      border = "black", direction = "horizontal"
    ),
  gp = gpar(col = "black")
)


# Create an expression with italic and colored row names
hm_prev <- Heatmap(prev_mat,cluster_rows = F, cluster_columns = F, col = col_fun_prev, column_labels = bacLabs,
  bottom_annotation = cont_Anno,
  height = unit(2, "cm"),
  heatmap_legend_param =
    list(
      title = "Prevalence",
      color_bar = "discrete",
      at = seq(0, 1, 0.1),
      border = "black"
    ),
  # appearance of labels
  column_names_rot = 90, column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8),
  row_names_side = c("right"),
  column_names_gp = gpar(fontsize = 8,col = bacteria_label_colors),

  # other
  cluster_column_slices = F,
  row_gap = unit(c(3), "mm"),
  column_gap = unit(c(3), "mm"),
  border = T,
  rect_gp = gpar(col = "black", lwd = 1)
)

# Combine the heatmaps
ht_list = hm_meanAB %v% hm_prev

# Draw the combined heatmaps and save plot
pdf(file = file.path(save_fig_folder,"C2_Heatmap_PrevalenceAbundance.pdf"),height = 8.27,width = 11.69)
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


