#########
# An Atlas of the Human Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Run survival analyses of the 5R 16S species level data for the HCC, CRLM and iCCA samples as well as the 5R 16S immunotherapy cohort.
# by Fabian Springer
########

library(tidyverse)
library(here)
library(survival)
library(survminer)

source(here("src","analysis","functions_analysis.R"))

#* Import 16S species level data and metadata ----
meta_all_df <- read_tsv(here("data","metadata","meta_5R16S_allSamples.tsv"))
species_relAbundance_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_allSamples.rds"))
log10_species_relAbundance_mat <- log10(species_relAbundance_mat + 1e-5)

#* Prepare metadata for survival analysis
meta_all_df$Etiology %>% table()
HCC_types <- c("HBV_HCC","HCV_HCC","ALD/ASH_HCC","MAFLD/MASH_HCC","other_HCC")

meta_surival_all_df <-
  meta_all_df %>%
  dplyr::select(Sample_ID, Etiology, OverallSurvival_time_months, OverallSurvival_event_death) %>%
  mutate(cancer_type = case_when(
    Etiology %in% HCC_types ~ "HCC",
    Etiology == "CRLM" ~ "CRLM",
    Etiology == "iCCA" ~ "iCCA",
    TRUE ~ NA
  )) %>%
  filter(!is.na(cancer_type), !is.array(OverallSurvival_time_months), !is.na(OverallSurvival_event_death)) %>%
  # Set upper limit for survival time to 60 months: If event occurs after 60 months, set to 0 and 60 months as max
  mutate(
    OverallSurvival_event_death =
      case_when(
        OverallSurvival_time_months > 60 & OverallSurvival_event_death == 1 ~ 0,
        TRUE ~ OverallSurvival_event_death
      ),
    OverallSurvival_time_months = ifelse(OverallSurvival_time_months > 60, 60, OverallSurvival_time_months)
  ) %>% 
  dplyr::select(-Etiology)

table(meta_surival_all_df$cancer_type)

# Generate one metadata table for HCC, CRLM and iCCA cancers
meta_survival_HCC_df <- meta_surival_all_df %>% filter(cancer_type == "HCC")
meta_survival_CRLM_df <- meta_surival_all_df %>% filter(cancer_type == "CRLM")
meta_survival_iCCA_df <- meta_surival_all_df %>% filter(cancer_type == "iCCA")

f_survival_analysis <- function(tax, rel_mat, meta_df,prevalence_threshold = FALSE) {
  # Function to perform survival analysis for a given taxon (tax). 
  # Subsets the relative abundance matrix with the provided taxon, generates a dataframe with rel. abundance, survival time and event.
  # Computes prevalence and performs log-rank test (based on taxon presence/absence)and cox-regression (log10-transformed rel. abundances).
  # Returns list with 1) Dataframe with test results, 2) survFit object, 3) survObj object

  survTime_col <- "OverallSurvival_time_months"
  event_col <- "OverallSurvival_event_death"

  stopifnot(tax %in% rownames(rel_mat) | all(c(survTime_col, event_col) %in% colnames(meta_df)))

  #* Prepare data for survival analysis: Generate dataframe with relative abundance, prevalence and metadata
  surv_df <- suppressMessages(rel_mat[tax, ] %>%
    enframe(name = "Sample_ID", value = "relAbundance") %>%
    mutate(
      isPrev = ifelse(relAbundance > 0, "present", "absent"),
      isPrev = factor(isPrev, levels = c("present", "absent")),
      l10_rel = log10(relAbundance + 1e-5)
    ) %>%
    dplyr::select(Sample_ID, isPrev, l10_rel) %>%
    inner_join(., meta_df) %>%
    mutate(prevalence = sum(isPrev == "present") / n()))

  if (prevalence_threshold != FALSE) {
    if (unique(surv_df$prevalence) < prevalence_threshold) {
      #message("Prevalence is below threshold, skipping survival analysis")
      return(NULL)
    }
  }

  #* Run survival analysis
  # Log-Rank test for Kaplan Meyer curves based on prevalence (bacterium present/absent)
  survObj <- survival::Surv(
    time = surv_df[[survTime_col]],
    event = surv_df[[event_col]]
  )
  survFit <- survminer::surv_fit(survObj ~ isPrev, data = surv_df)

  # Cox-regression based on bacretial relative abundance (log10-transformed)
  cox_res <- coxph(survObj ~ l10_rel, data = surv_df)
  cox_res_summary <- summary(cox_res)

  # Extract hazard ratio and p-value
  hazard_ratio <- cox_res_summary$coefficients[1, "exp(coef)"]
  parm <- "l10_rel"
  # Compute upper and lower 95% confidence intervals
  HR_upper_95_CI <- exp(confint(cox_res, parm = parm, level = 0.95))[2]
  HR_lower_95_CI <- exp(confint(cox_res, parm = parm, level = 0.95))[1]
  pval_cox <- as.numeric(cox_res_summary$logtest["pvalue"])

  #* Generate output dataframe with results of the individual tests
  res_df <- survminer::surv_pvalue(fit = survFit) %>%
    transmute(p.val_LR = pval) %>%
    as_tibble() %>%
    mutate(
      p.val_Cox = pval_cox,
      HR = hazard_ratio,
      HR_upper_95_CI = HR_upper_95_CI,
      HR_lower_95_CI = HR_lower_95_CI,
      N_Samples = nrow(surv_df),
      N_Events = sum(surv_df[[event_col]] == 1),
      Prevalence = unique(surv_df$prevalence)
    )
  res_list <- list(
    "fit_stats_df" = res_df,
    "survFit" = survFit,
    "survObj" = survObj
  )

  return(res_list)
}

# Run for HCC cohort
message("Running survival analysis for HCC cohort")
survival_res_all_df <- tibble()
survival_ojects_list <- list()
pb <- progress_bar$new(total = nrow(species_relAbundance_mat))
c <- 1
for(tax in rownames(species_relAbundance_mat)){  
  res_list <- f_survival_analysis(tax, rel_mat = species_relAbundance_mat, meta_df = meta_survival_HCC_df, prevalence_threshold = 0.05)
  if(!is.null(res_list)){
    survival_res_all_df <- bind_rows(survival_res_all_df, res_list$fit_stats_df %>% mutate(tax = tax))
    survival_ojects_list[[c]] <- list('survFit' = res_list$survFit,'survObj' = res_list$survObj)
    names(survival_ojects_list)[c] <- tax
    c <- c+1
  }
  pb$tick()
}
survival_res_all_df <- survival_res_all_df %>%
  dplyr::relocate(tax) %>%
  mutate(
    p.val_adj_LR = p.adjust(p.val_LR, method = "fdr"),
    p.val_adj_Cox = p.adjust(p.val_Cox, method = "fdr")
  ) %>%
  dplyr::relocate(tax,p.val_LR, p.val_adj_LR, p.val_Cox, p.val_adj_Cox) %>% 
  arrange(p.val_LR)
write_tsv(survival_res_all_df,here("data","results","survival_results_HCC_df.tsv"))
saveRDS(survival_ojects_list,here("data","results","survival_resObjects_HCC.rds"))

# Run for CRLM cohort
message("Running survival analysis for CRLM cohort")
survival_res_all_df <- tibble()
survival_ojects_list <- list()
pb <- progress_bar$new(total = nrow(species_relAbundance_mat))
c <- 1
for(tax in rownames(species_relAbundance_mat)){  
  res_list <- f_survival_analysis(tax, rel_mat = species_relAbundance_mat, meta_df = meta_survival_CRLM_df, prevalence_threshold = 0.05)
  if(!is.null(res_list)){
    survival_res_all_df <- bind_rows(survival_res_all_df, res_list$fit_stats_df %>% mutate(tax = tax))
    survival_ojects_list[[c]] <- list('survFit' = res_list$survFit,'survObj' = res_list$survObj)
    names(survival_ojects_list)[c] <- tax
    c <- c+1
  }
  pb$tick()
}
survival_res_all_df <- survival_res_all_df %>%
  dplyr::relocate(tax) %>%
  mutate(
    p.val_adj_LR = p.adjust(p.val_LR, method = "fdr"),
    p.val_adj_Cox = p.adjust(p.val_Cox, method = "fdr")
  ) %>%
  dplyr::relocate(tax,p.val_LR, p.val_adj_LR, p.val_Cox, p.val_adj_Cox) %>% 
  arrange(p.val_LR)
write_tsv(survival_res_all_df,here("data","results","survival_results_CRLM_df.tsv"))
saveRDS(survival_ojects_list,here("data","results","survival_resObjects_CRLM.rds"))

# Run for iCCA cohort
message("Running survival analysis for iCCA cohort")
survival_res_all_df <- tibble()
survival_ojects_list <- list()
pb <- progress_bar$new(total = nrow(species_relAbundance_mat))
c <- 1
for(tax in rownames(species_relAbundance_mat)){  
  res_list <- f_survival_analysis(tax, rel_mat = species_relAbundance_mat, meta_df = meta_survival_iCCA_df, prevalence_threshold = 0.05)
  if(!is.null(res_list)){
    survival_res_all_df <- bind_rows(survival_res_all_df, res_list$fit_stats_df %>% mutate(tax = tax))
    survival_ojects_list[[c]] <- list('survFit' = res_list$survFit,'survObj' = res_list$survObj)
    names(survival_ojects_list)[c] <- tax
    c <- c+1
  }
  pb$tick()
}
survival_res_all_df <- survival_res_all_df %>%
  dplyr::relocate(tax) %>%
  mutate(
    p.val_adj_LR = p.adjust(p.val_LR, method = "fdr"),
    p.val_adj_Cox = p.adjust(p.val_Cox, method = "fdr")
  ) %>%
  dplyr::relocate(tax,p.val_LR, p.val_adj_LR, p.val_Cox, p.val_adj_Cox) %>% 
  arrange(p.val_LR)
write_tsv(survival_res_all_df,here("data","results","survival_results_iCCA_df.tsv"))
saveRDS(survival_ojects_list,here("data","results","survival_resObjects_iCCA.rds"))


#* Run the survival analysis for the HCC immunotherapy cohort ----
meta_immunotherapy_df <- read_tsv(here("data","metadata","meta_5R16S_ImmunotherapyCohort.tsv"))
species_relAbundance_immunotherapy_mat <- readRDS(here("data","raw","relAbundance_5R16S_species_ImmunotherapyCohort.rds"))

meta_surival_ICI_df <-
  meta_immunotherapy_df %>%
  dplyr::select(Sample_ID,Response_to_immunotherapy,OverallSurvival_time_months, OverallSurvival_event_death) %>%  
  filter(!is.array(OverallSurvival_time_months), !is.na(OverallSurvival_event_death)) %>%
  # Set upper limit for survival time to 60 months: If event occurs after 60 months, set to 0 and 60 months as max
  mutate(
    OverallSurvival_event_death =
      case_when(
        OverallSurvival_time_months > 60 & OverallSurvival_event_death == 1 ~ 0,
        TRUE ~ OverallSurvival_event_death
      ),
    OverallSurvival_time_months = ifelse(OverallSurvival_time_months > 60, 60, OverallSurvival_time_months)
  )

# Step1: Compare immune-responders vs. non-responders
message("Running survival analysis for ICI cohort")
tmp_mat <- meta_surival_ICI_df %>%
  dplyr::select(Sample_ID, Response_to_immunotherapy) %>%
  mutate(ICI = ifelse(Response_to_immunotherapy == "responder", 1, 0)) %>%
  dplyr::select(-Response_to_immunotherapy) %>% 
  pivot_wider(names_from = Sample_ID, values_from = ICI) %>% 
  mutate(rowN = "ICI") %>% 
  column_to_rownames("rowN") %>% 
  as.matrix(drop = F)
res_list <- f_survival_analysis("ICI", rel_mat = tmp_mat, meta_df = meta_surival_ICI_df, prevalence_threshold = 0.05)
ICI_res_list <- list(list('survFit' = res_list$survFit,'survObj' = res_list$survObj))
names(ICI_res_list) <- "RespVsNonResp"
#saveRDS(ICI_res_list,here("data","results","survival_resObjects_I.rds"))

# Perform species level surival analysis for the ICI cohort
message("Running survival analysis for ICI cohort")
survival_res_all_df <- tibble()
survival_ojects_list <- list()
pb <- progress_bar$new(total = nrow(species_relAbundance_immunotherapy_mat))
c <- 1
for(tax in rownames(species_relAbundance_immunotherapy_mat)){  
  res_list <- f_survival_analysis(tax, rel_mat = species_relAbundance_immunotherapy_mat, meta_df = meta_surival_ICI_df, prevalence_threshold = 0.05)
  if(!is.null(res_list)){
    survival_res_all_df <- bind_rows(survival_res_all_df, res_list$fit_stats_df %>% mutate(tax = tax))
    survival_ojects_list[[c]] <- list('survFit' = res_list$survFit,'survObj' = res_list$survObj)
    names(survival_ojects_list)[c] <- tax
    c <- c+1
  }
  pb$tick()
}
survival_res_all_df <- survival_res_all_df %>%
  dplyr::relocate(tax) %>%
  mutate(
    p.val_adj_LR = p.adjust(p.val_LR, method = "fdr"),
    p.val_adj_Cox = p.adjust(p.val_Cox, method = "fdr")
  ) %>%
  dplyr::relocate(tax,p.val_LR, p.val_adj_LR, p.val_Cox, p.val_adj_Cox) %>% 
  arrange(p.val_LR)

# Add results for responders vs non-responders
survival_ojects_list <- c(ICI_res_list,survival_ojects_list)
names(survival_ojects_list)[1:10]

# save everything
write_tsv(survival_res_all_df,here("data","results","survival_results_ImmunotherapyCohort_df.tsv"))
saveRDS(survival_ojects_list,here("data","results","survival_resObjects_ImmunotherapyCohort.rds"))

