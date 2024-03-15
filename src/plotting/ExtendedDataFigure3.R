#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Reproduce the plots shown in Extended Data Figure 3
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
save_fig_folder <- here("figures","ExtendedDataFigures","ExtendedDataFigure3")
if(!dir.exists(save_fig_folder)){
  dir.create(save_fig_folder, recursive = TRUE)
}

#* Import corrected rel. abundance matrix ----
meta_combined_df <- read_tsv(here("data","metadata","meta_combined_HCC-meta-analysis.tsv"))
genus_RelAbundance_corrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_BatchEffectCorrected.rds"))
genus_RelAbundance_uncorrected_mat <- readRDS(here("data","processed","relAbundance_combined_genus_HCC-meta-analysis_raw.rds"))


# Import raw PathSeq scores
genus_counts_raw_mat <- readRDS(here("data","raw","rawScores_bulkRNAseq_genus.rds"))

# Import testing results for pairwise comparisons of sequencing reads, total bacterial reads, bacterial CPM, alpha- and beta-diversity
dataset_diversity_res_df <- read_tsv(here("data", "results", "HCC-meta-analysis_diversity_comparisons_byDataset.tsv")) %>%
  # remove RNAseq from DKFZ_RNAseq
  mutate(
    Group1 = str_remove(Group1, " RNA-Seq"),
    Group2 = str_remove(Group2, " RNA-Seq")
  )
batch_diversity_res_df <- read_tsv(here("data","results","HCC-meta-analysis_diversity_comparisons_byBatch.tsv")) %>% 
  mutate(
    Group1 = str_remove(Group1, "_RNAseq"),
    Group2 = str_remove(Group2, "_RNAseq")
  )


# Select RNAseq IDs considered in the meta-analysis and compute metrics for comparison
RNAseq_IDs <- meta_combined_df %>%
  dplyr::select(Sample_ID, Dataset, Tissue_type, FastQ_libsize) %>%
  filter(!str_detect(Dataset, "16S")) %>% 
  pull(Sample_ID)

total_bact_counts_df <- genus_counts_raw_mat[,RNAseq_IDs] %>% colSums() %>% enframe(name = "Sample_ID",value = "Bacteria_total")
dim(total_bact_counts_df)

#* Compute Shannon diversity and richness based on rarefied rawCounts, combine with Number of raw reads, number of bacterial reads, bacterial CPM ----
set.seed(1)
count_rar <- as.matrix(t(vegan::rrarefy(t(round(genus_counts_raw_mat[,RNAseq_IDs],0)),1000)))
# There is 1 sample with a rrarefied total number of bacteria smaller than 1000.
# This is because the sum of the PathSeq scores in this sample becomes less than 1000 when rounded to integers. 

count_rar_rel <- prop.table(count_rar,2)
div <- vegan::diversity(t(count_rar_rel), index = "shannon") %>% enframe(name = "Sample_ID",value = "Shannon")
rich <- colSums(count_rar_rel > 0) %>% enframe(name = "Sample_ID",value = "Richness")
diversity_richness_df <- full_join(div,rich,by = "Sample_ID")

# Combine with metadata 
diversity_richness_Dataset_df <- meta_combined_df %>%
  dplyr::select(Sample_ID, Dataset, FastQ_libsize) %>% 
  filter(Sample_ID %in% RNAseq_IDs) %>% 
  full_join(.,total_bact_counts_df) %>%
  full_join(.,diversity_richness_df) %>% 
  mutate(Bacteria_CPM = Bacteria_total/FastQ_libsize*1e6) %>% 
  left_join(.,meta_combined_df %>% dplyr::select(Sample_ID,Batch)) %>% 
  mutate(Dataset = ifelse(str_detect(Dataset,"DKFZ"), "DKFZ", Dataset)) %>% 
  mutate(Batch = ifelse(str_detect(Dataset,"DKFZ"), "DKFZ", Batch)) %>% 
  dplyr::relocate(Sample_ID,Dataset,Batch)

### Generate helper function for the boxplot
f_helper_boxplot <- function(df,metric,diversity_test_res_df,metric_clean){
  
  stopifnot(is.factor(df$Group))
  levels <- levels(df$Group)
  # Generate plot DF
  plot_df <- df %>% transmute(Sample_ID,Group,y := !!as.symbol(metric))

  # Generate boxplot
  pt_box <- f_boxplot_by_group(plot_df, xlab = "", ylab = metric_clean, corral.width = 0.49) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors)

  ### Plot matrix to indicate significance of pairwise comparisons
  # Get the p-values from the diversity test result df
  p_val_vec <- diversity_test_res_df %>%
    filter(Metric == metric_clean) %>%
    transmute(pair_comp = paste0(Group1, "_vs_", Group2), p.val_adj_wilcox) %>%
    deframe()
  # generate triangular matrix with p-values
  tri_mat <- f_create_upper_triangle_matrix(p_val_vec, condition_levels = levels)
  # For visualization: Convert "-" to linebreaks
  rownames(tri_mat) <- str_replace(rownames(tri_mat), "ISMMS-IDIBAPS$", "ISMMS-\nIDIBAPS") # dont shorten rownames when comparing batches
  colnames(tri_mat) <- str_replace(colnames(tri_mat), "ISMMS-IDIBAPS$", "ISMMS-\nIDIBAPS")

  # remove _RNAseq from names
  rownames(tri_mat) <- str_remove(rownames(tri_mat), "_RNAseq")
  colnames(tri_mat) <- str_remove(colnames(tri_mat), "_RNAseq")

  # plot significance matrix
  pt_sig_mat <- f_plot_signif_matrix(tri_mat)
  

  # return boxplot and matrix individually to manually adjust y-axis of boxplot
  return(list(box = pt_box, sig_mat = pt_sig_mat))

}

# Geenrate the boxplots

# Define metrics to compare
metrics_to_compare <- c("FastQ_libsize","Bacteria_total","Bacteria_CPM","Shannon","Richness")
# Define the clean name of the metrics (represents the Metric name found in the diversity test result df)
metrics_to_compare_clean_names <- c("Number of raw reads","Bacterial counts","Bacterial CPM","Shannon diversity","Genus richness")


#* Generate plots for Dataset comparison ----
# Define dataframes by Dataset and by Batch for comparison
levels_Dataset <- c("DKFZ","TCGA","INSERM","ISMMS-IDIBAPS")
df_byDataset <- diversity_richness_Dataset_df %>% mutate(Group = factor(Dataset,levels = levels_Dataset))

levels_Batch <- df_byDataset %>% arrange(Group,Batch) %>% pull(Batch) %>% unique()
df_byBatch <- diversity_richness_Dataset_df %>% mutate(Group = factor(Batch,levels = levels_Batch))

# Define list with plot y-axis
man_y_axis_list <- list(
  "FastQ_libsize" = list(
    "breaks" = c(5e6, 1e7, 2e7, 5e7, 1e8, 2e8), # breaks
    "labels" = c(5e6, 1e7, 2e7, 5e7, 1e8, 2e8), # labels
    "limits" = c(5e6, 2e8), # limimts
    "log-scale" = TRUE # log scale
  ),
  "Bacteria_total" = list(
    "breaks" = 10^seq(3,7,1), # breaks
    "labels" = 10^seq(3,7,1), # labels
    "limits" = 10^c(3,7), # limimts
    "log-scale" = TRUE # log scale
  ),
  "Bacteria_CPM" = list(
    "breaks" = 10^seq(1,5,1), # breaks
    "labels" = 10^seq(1,5,1), # labels
    "limits" = 10^c(1,5), # limimts
    "log-scale" = TRUE # log scale
  )
)


idx <- 1
plot_idx_dataset <- c("A","B","C","D","E")
plot_idx_batch <- c("G","H","I","J","K")
for(i in seq(1,length(metrics_to_compare))){
  metric <- metrics_to_compare[i]
  metric_clean <- metrics_to_compare_clean_names[i]
  
  #* For the Dataset comparison ----
  # generate plots
  c_plots <- f_helper_boxplot(
    df = df_byDataset,
    diversity_test_res_df = dataset_diversity_res_df,
    metric = metric,
    metric_clean = metric_clean
  )
    
  # adjust y-axis if necessary
  if(metric %in% names(man_y_axis_list)){
    y_parms_list <- man_y_axis_list[[metric]]
    if(y_parms_list[["log-scale"]]){
      c_plots$box <- c_plots$box +
        scale_y_log10(
          breaks = y_parms_list[["breaks"]],
          labels = y_parms_list[["labels"]],
          limits = y_parms_list[["limits"]]
        )

    } else {
      c_plots$box <- c_plots$box +
        scale_y_continuous(
          breaks = y_parms_list[["breaks"]],
          labels = y_parms_list[["labels"]],
          limits = y_parms_list[["limits"]]
        )
    }  
  }

  # combine plot with significance matrix
  heights = c(2, -0.25, 1)
  pt_combined_dataset <-
    (
      c_plots$box +
        theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    ) /
    plot_spacer() / c_plots$sig_mat + plot_layout(heights = heights, widths = c(1, 1, 1))


  #* For the Batch comparison ----
  # generate plots
  c_plots <- f_helper_boxplot(
    df = df_byBatch,
    diversity_test_res_df = batch_diversity_res_df,
    metric = metric,
    metric_clean = metric_clean
  )  
  c_plots$box
  # adjust y-axis if necessary
  if(metric %in% names(man_y_axis_list)){
    y_parms_list <- man_y_axis_list[[metric]]  
    if(y_parms_list[["log-scale"]]){
      c_plots$box <- c_plots$box +
        scale_y_log10(
          breaks = y_parms_list[["breaks"]],
          labels = y_parms_list[["labels"]],
          limits = y_parms_list[["limits"]]
        )

    } else {
      c_plots$box <- c_plots$box +
        scale_y_continuous(
          breaks = y_parms_list[["breaks"]],
          labels = y_parms_list[["labels"]],
          limits = y_parms_list[["limits"]]
        )
    }  
  }

  # combine plot with significance matrix
  heights = c(2, -0.2, 1.5)
  pt_combined_batch <-
    (
      c_plots$box +
        theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    ) /
    plot_spacer() / c_plots$sig_mat + plot_layout(heights = heights, widths = c(1, 1, 1))


  # save plot
  ggsave(pt_combined_dataset,filename = here(save_fig_folder,paste0(plot_idx_dataset[idx],"_",metric_clean,"_byDataset.pdf")),width = 4.5,height = 4.5)
  ggsave(pt_combined_batch,filename = here(save_fig_folder,paste0(plot_idx_batch[idx],"_",metric_clean,"_byBatch.pdf")),width = 12,height = 6)
  idx <- idx + 1
}

#* Compute PCoAs ----
### Plot PCoA based on UN-normalized rel. abudnances (after rarefication)

# For the dataset comparison
pt_pcoa <- f_plot_PCoA(
    meta_df = df_byDataset,
    mat = count_rar_rel[,df_byDataset$Sample_ID], method = "bray",prevalence_threshold = FALSE,
  ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    theme(legend.position = "none")+
    theme(aspect.ratio = 3/4)

p_val_vec <- dataset_diversity_res_df %>%
    filter(Metric == "Bray-Curtis distance") %>%
    transmute(pair_comp = paste0(Group1, "_vs_", Group2), p.val_adj_permanova) %>%
    deframe()
pt_sig_mat <- f_plot_signif_matrix(f_create_upper_triangle_matrix(p_val_vec, condition_levels = levels_Dataset))
pt_PCoA_dataset <-
    (pt_pcoa +theme(legend.position = "none")) /
    plot_spacer() / pt_sig_mat + plot_layout(heights = c(2, -0.2, 1), widths = c(1, 1, 1))
ggsave(pt_PCoA_dataset,filename = here(save_fig_folder,"F_PCoA_byDataset.pdf"),width = 4.5,height = 4.5)

# For the batch comparison
pt_pcoa <- f_plot_PCoA(
    meta_df = df_byBatch,
    mat = count_rar_rel[,df_byDataset$Sample_ID], method = "bray",prevalence_threshold = FALSE,
  ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    theme(legend.position = "none")+
    theme(aspect.ratio = 3/4)

p_val_vec <- batch_diversity_res_df %>%
    filter(Metric == "Bray-Curtis distance") %>%
    transmute(pair_comp = paste0(Group1, "_vs_", Group2), p.val_adj_permanova) %>%
    deframe()
pt_sig_mat <- f_plot_signif_matrix(f_create_upper_triangle_matrix(p_val_vec, condition_levels = levels_Batch))
pt_PCoA_Batch <-
    (pt_pcoa +theme(legend.position = "none",aspect.ratio = 3/4)) /
    plot_spacer() / (pt_sig_mat+theme(axis.text.x = element_text(angle = 90,hjust = 0.99,vjust = 0.5))) +
    plot_layout(heights = c(2.2, -0.2, 1.2), widths = c(1, 1, 1))
ggsave(pt_PCoA_Batch,filename = here(save_fig_folder,"L_PCoA_byBatch.pdf"),width = 5,height = 5.8)

