#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Analysis of DNA concentrations across cancer types and controls obtained via qPCR
# by Fabian Springer
########

library(tidyverse)
library(here)

# Import the qPCR result table
qPCR_dat_df <- read_tsv(here("data","raw","qPCR_concentrations_df.tsv"))

# Add sample numbers to the group names
qPCR_dat_df <- qPCR_dat_df %>%
    left_join(., qPCR_dat_df %>%
        group_by(Group) %>%
        summarise(N = n())) %>% 
        mutate(Group_N = paste0(Group,"\n(N=",N,")")) %>% 
        #mutate(Group_N = paste0(Group," (N=",N,")")) %>% 
        arrange(Group) %>% 
        mutate(Group_N = as_factor(Group_N))

wilcox_tests_list <- list(c("CTRL", "HCC"),c("BTC", "CRLM"),c("CTRL", "CRLM"), c("CTRL", "BTC"),c("HCC", "CRLM"))


# Do ONE-SIDED wilcoxon tests to test for higher abunadnces in the cancer samples (or in CRLM compared to primary liver cancer)
p_val_vec <- c()
for(i in seq(1,length(wilcox_tests_list))){
    ref <- wilcox_tests_list[[i]][1]
    group <- wilcox_tests_list[[i]][2]
    x <- subset(qPCR_dat_df,Group == ref)[["bacterial_cells_per_40ngDNA"]]
    y <- subset(qPCR_dat_df,Group == group)[["bacterial_cells_per_40ngDNA"]]
    res <- wilcox.test(y,x,alternative = "greater")
    p_val <- res$p.value
    names(p_val) <- paste0(ref,"_",group)
    p_val_vec <- c(p_val_vec,p_val)
}

# Duplicate the Group infor and sample number information to merge with test results
tmp_df <- qPCR_dat_df %>% 
    group_by(Group) %>% 
    summarise(median_bacterial_cells_per_40ngDNA = median(bacterial_cells_per_40ngDNA)) %>% 
    left_join(.,qPCR_dat_df %>% dplyr::select(Group,N)) %>% 
    distinct() %>% 
    mutate(Group1 = Group,Group2 = Group,
    median_bacterial_cells_per_40ngDNA_Group1 = median_bacterial_cells_per_40ngDNA,
    median_bacterial_cells_per_40ngDNA_Group2 = median_bacterial_cells_per_40ngDNA,
    N_Group1 = N,N_Group2 = N)

# Generate clean dataframe with test results of pairwise comparisons
qPCR_test_result_df <-
  enframe(p_val_vec, name = "comparison", value = "p.val_wilcox") %>%
  mutate(p.val_adj_wildox = p.adjust(p.val_wilcox, method = "fdr")) %>%
  separate(comparison, into = c("Group1", "Group2"), sep = "_") %>%
  left_join(., tmp_df %>% dplyr::select(contains("Group1"))) %>%
  left_join(., tmp_df %>% dplyr::select(contains("Group2")))

write_tsv(qPCR_test_result_df,here("data","results","qPCR_results_df.tsv"))
