#########
# An Atlas of the Liver Cancer Microbiome
# Rahbari, Springer, Zwang et al. (2024)
# Functions for the plotting scripts
# by Fabian Springer
########

require(ggrepel)
require(gghalves)
require(ggbeeswarm)
require(ggsignif)

# General plotting theme building upon the custom ggplot2 theme called "ggembl"
# ggembl can be downloaded here: https://git.embl.de/grp-zeller/ggembl
theme_paper <- ggembl::theme_presentation() +
    theme(
        axis.title = element_text(face = "bold", size = 12), # Bold axis titles
        panel.border = element_rect(fill=NA, colour='black', size=1.5),
        #axis.line = element_line(size = 1.5), # Change line thickness of axes

        axis.text = element_text(face = "bold", size = 12),
    )

f_create_label  <- function(name_vec){
    # Takes taxon strings as input and returns a cleaned label (e.g. s_Acinetobacter unknown species 1 -> Acinetobacter sp. [US1])
    require(dplyr)
    require(stringr)
    require(tidyr)
    df <- as_tibble(name_vec) %>%
        mutate(
            species = str_extract(value, "s__([^;]*)"), # Extract species
            genus = str_extract(value, "g__([^;]*)"), # Extract genus
            family = str_extract(value, "f__([^;]*)"), # Extract family
            order = str_extract(value, "o__([^;]*)"), # Extract order
            class = str_extract(value, "c__([^;]*)"), # Extract class
            phylum = str_extract(value, "p__([^;]*)"), # Extract phylum
            unknown_species = str_extract(species, "Unknown species(\\d+)"), # Extract the identifier for unknown species
            unknown_genus = str_extract(genus, "Unknown genus(\\d+)"),# Extract the identifier for unknown species
            lowest_known = case_when(
                !is.na(species) & !str_detect(species, "Unknown") ~ species,
                !is.na(genus) & !str_detect(genus, "Unknown") ~ genus,
                !is.na(family) & !str_detect(family, "Unknown") ~ family,
                !is.na(order) & !str_detect(order, "Unknown") ~ order,
                !is.na(class) & !str_detect(class, "Unknown") ~ class,
                TRUE ~ phylum
            ),
            lowest_assigned = case_when(
              !is.na(species) ~ species,
              !is.na(genus) ~ genus,
              !is.na(family) ~ family,
              !is.na(order) ~ order,
              !is.na(class) ~ class,
              TRUE ~ phylum
            ),
            # create the new label
            label = case_when(
                !is.na(species) & !str_detect(species, "Unknown species") ~ str_remove(species, "s__"), # Known species
                !is.na(species) & str_detect(species, "Unknown species") & !is.na(genus) & !str_detect(genus, "Unknown genus") ~ paste0(str_remove(genus, "g__"), " sp.", " [US", str_remove(unknown_species, "Unknown species"), "]"), # Unknown species, known genus
                !is.na(species) & str_detect(species, "Unknown species") ~ paste0(str_remove(lowest_known, "[a-z]__"), " [US", str_remove(unknown_species, "Unknown species"), "]"), # Unknown species, known genus
                #is.na(species) & !is.na(genus) ~ str_remove(genus, "g__"), # Species is NA but genus is known
                
                #for genus
                is.na(species) & !is.na(genus) & str_detect(genus, "Unknown genus") & !is.na(family) & !str_detect(family, "Unknown family") ~ paste0(str_remove(family, "f__"), " [UG", str_remove(unknown_genus, "Unknown genus"), "]"),
                
                TRUE ~ paste0(str_remove(lowest_known, "[a-z]__")) #Ignore everything, just return lowest known level
                #TRUE ~ paste0(str_remove(lowest_known, "[a-z]__"), " [US", str_remove(unknown_species, "Unknown species"), "]") # Known genus, unknown species or Unknown genus
            )
        )
    #df %>% glimpse()
    names_clean_vec <- df %>% pull(label)
    return(names_clean_vec)
}

f_convert_qval_pval = function(qvalues,fdr_threshold) {
    # This function is used to indicate in the volcano plots the p-value that would correspond to a given FDR threshold (e.g. 0.2)    
  
    #Function to approximate p-value from a distribution of q-values
    #copy-pasted from here: https://stats.stackexchange.com/questions/51070/how-can-i-convert-a-q-value-distribution-to-a-p-value-distribution

    # you need to know the estimate of pi0 used to create the q-value
    # that's the maximum q-value (or very, very close to it)
    qvalues_extended <- c(qvalues,fdr_threshold)
    pi0 = max(qvalues_extended)
    # compute m0, the estimated number of true nulls
    m0 = length(qvalues_extended) * pi0
    # then you multiply each q-value by the proportion of true nulls
    # expected to be under it (the inverse of how you get there from
    # the p-value):
    p_values_reconstructed <- qvalues_extended * rank(qvalues_extended) / m0
    #return the last element of the reconstructed p-values (corresponds to the fdr_threshold - p-value)
    return(p_values_reconstructed[length(p_values_reconstructed)])
}

p_to_symbol <- function(p) {
  # Convert p-values to asterisk symbols: takes a vector of p-values and returns a vector of symbols
  res <- character(length(p))  
  
  res[p < 0.001] <- "***"
  res[p >= 0.001 & p < 0.01] <- "**"
  res[p >= 0.01 & p < 0.1] <- "*"
  res[p >= 0.1] <- "n.s."
  
  # Handle any other cases (e.g., NAs in the original p-values)
  res[is.na(p)] <- NA
  
  return(res)
}

f_plot_volcano <- function(plot_df,xBreaks,xLims,fdr_thresh=0.2,man_y_breaks=NULL,clean_tax_names=TRUE,add_to_y_axis = 0.25){
  # Takes a dataframe with testing results and generates a volcano plot

  stopifnot(all(c("p.val_lm","effect.size","tax") %in% colnames(plot_df)))

  leg_text_size <- 8
  yName <- "P-value"
  pMethod <- "p.val_lm"

  # add column in case to a later point other p-values are used
  plot_df <- plot_df %>% mutate(p.val = p.val_lm,p.val_adj = p.val_adj_lm)

  # Get p-value equivalent of selected fdr-threshold
  q_value_threshold <- f_convert_qval_pval(qvalues = plot_df$p.val_adj_lm, fdr_threshold = fdr_thresh)


  # Define y-axis breaks if not supplied
  if(is.null(man_y_breaks)) {
    # For all except fibrosis
    man_y_breaks <- -log10(c(0.05,0.01,0.001,0.0001,0.00001))#,0.000001))    
  }  
  man_y_labs <- signif(10^-man_y_breaks,1)
  man_y_lims <- c(0,round(max(man_y_breaks+add_to_y_axis),1))

  # Check if supplied y-axis and x-axis limits are valid
  if (min(plot_df[[pMethod]]) < min(10^-(man_y_breaks))) {
    stop("smallest p-value is smaller than y-axis")
  }
  if (min(xLims) > min(plot_df$effect.size) | max(xLims) < max(plot_df$effect.size)) {
    stop("xlimits smaller than effect size")
  }

  # Clean taxon names if requested
  if (isTRUE(clean_tax_names)) {
      plot_df <- plot_df %>% mutate(lab = f_create_label(tax))
  } else {
      plot_df <- plot_df %>% mutate(lab = tax)
  }

  # extract Group names, Group numbers
  group1 <- unique(plot_df$Group1)
  group2 <- unique(plot_df$Group2)
  N_group1 <- unique(plot_df$N_Group1)
  N_group2 <- unique(plot_df$N_Group2)

  # Step1: Do some checks for significances etc
  plot_df <- plot_df %>%
    mutate(
      enriched_in =
        case_when(
          p.val < 0.05 & effect.size > 0 ~ group2,
          p.val < 0.05 & effect.size < 0 ~ group1,
          TRUE ~ "n.s."
        ),
      # Define whether a taxon is FDR significant or not
      fdr_sig = ifelse(p.val_adj < fdr_thresh, TRUE, FALSE),
    ) %>%
    mutate(
      enriched_in = factor(enriched_in, levels = c(group1, group2, "n.s.")),
      font = ifelse(fdr_sig, "bold.italic", "italic"), # for labeling
      group_prev = ifelse(effect.size > 0, Prev_Group2, Prev_Group1), # size based on prevalence in group
      lab := case_when(!!as.symbol(pMethod) < 0.05 | p.val_adj < fdr_thresh ~ lab, TRUE ~ "")
    )  
  # Step2: Do the actual plotting
  pt <- plot_df %>%
    ggplot(aes(x = effect.size, y = -log10(!!as.symbol(pMethod)))) +
    theme_paper +
    #geom_hline(yintercept = man_y_breaks[2:length(man_y_breaks) - 1], color = "grey", lty = "solid", lwd = 0.2) +
    geom_hline(yintercept = man_y_breaks[2:length(man_y_breaks)], color = "grey", lty = "solid", lwd = 0.2) +
    geom_hline(yintercept = -log10(0.05), color = "black", lty = "solid", lwd = 0.5) +
    geom_hline(yintercept = -log10(q_value_threshold), color = "darkgrey", lty = "dashed", lwd = 0.75) +
    geom_point(
      aes(
        size = group_prev,
        fill = enriched_in
      ),
      shape = 21, color = "black", alpha = 0.75
    ) +
    scale_size(range = c(2, 5), guide = guide_legend(reverse = TRUE)) + # Adjust this range as needed
    ggrepel::geom_text_repel(aes(label = lab, fontface = font),
      color = "black",
      segment.color = "black",
      size = 3,
      seed = 420
    )

  # Refine plot with legends and axis etc
  pt_res <- pt +
    theme(
      legend.box = "horizontal",
      legend.spacing.y = unit(0.1, "cm"),
      legend.position = c(.01, .99),
      legend.key.size = unit(0.75, "lines"),
      legend.justification = c(0, 1)
    ) +
    labs(fill = paste0("P < 0.05"), size = "Prevalence") +
    guides(
      color = "none",
      size = guide_legend(order = 1, reverse = T),
      fill = guide_legend(order = 2, override.aes = list(size = 3.5))
    ) +
    scale_x_continuous(breaks = xBreaks, labels = xBreaks, limits = xLims) +
    scale_y_continuous(breaks = man_y_breaks, labels = man_y_labs, name = yName, position = "left", limits = man_y_lims) +
    xlab("Enrichment effect size") + 
    annotate("text", x = xLims[1], y = -Inf, label = paste0(group1," (N=",N_group1,")"), hjust = -0.01, vjust = -0.75, fontface = "bold") +
    annotate("text", x = xLims[2], y = -Inf, label = paste0(group2," (N=",N_group2,")"), hjust = 0.99, vjust = -0.75, fontface = "bold")
  
  return(pt_res)
}

f_boxplot_by_group <- function(plot_df,xlab,ylab,corral.width = 0.49){
    # Function to create a boxplot with half-violin and beeswarm plot
    stopifnot(all(c("Group","y") %in% colnames(plot_df)))
    stopifnot(is.factor(plot_df$Group))
    #plot the basic boxplot
    pt <-
      plot_df %>%
      ggplot(aes(x = Group, y = y, fill = Group)) +
      theme_paper +
      gghalves::geom_half_boxplot(outlier.size = -1, side = "l", errorbar.length = 0.25, notch = 0, nudge = 0.025) +
      gghalves::geom_half_violin(side = "r", nudge = 0.025, alpha = 0.1, show.legend = F) + # aes(color=Group)) +
      ggbeeswarm::geom_beeswarm(alpha = 0.25, corral = "random",corral.width = corral.width,side = 1, aes(color = Group, x = as.numeric(Group) + 0.05), pch = 16, cex = 0.75,show.legend = F) + 
      ylab(ylab) +
      xlab(xlab)      

    return(pt)
}

f_kaplan_meier <- function(survFit,title){
  # Function to create a Kaplan-Meier plots from a survFit object
  pt <- survminer::ggsurvplot(survFit,
      title = title,
      conf.int = TRUE, pval = TRUE, pval.size = 3,
      risk.table = TRUE, fontsize = 3,
      legend.labs = c("present","absent"),
      legend.title = "",
      xlab = "Time post resection [months]",
      ylab = "Overall survival probability",
      conf.int.style = "step",
      # surv.median.line = "hv",
      palette = c("#E40046","grey"),
      risk.table.height = .25,
      ggtheme = theme_paper+theme(plot.title = element_text(size = 12,hjust = 0.5, face = "bold.italic")),
      risk.table.col = "strata",
      risk.table.title="",
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = FALSE
    )
    return(pt)
}

f_plot_PCoA <- function(meta_df, mat, method = "bray", threshold_for_prevalence = 0, prevalence_threshold = 0.05, point_alpha = 0.75) {
  # Takes a metadata dataframe, a species relative abundance matrix and computes a PCoA and returns a ggplot object.
  # Calls the f_compute PCoA function for computation

  stopifnot(all(c("Sample_ID", "Group") %in% colnames(meta_df)))
  pcoa_list <- f_compute_PCoA(meta_df, mat, method, threshold_for_prevalence, prevalence_threshold)
  x_lab <- paste0("PCoA1 [", pcoa_list$var_expl_df[1, 1], "%]")
  y_lab <- paste0("PCoA2 [", pcoa_list$var_expl_df[1, 2], "%]")
  z_lab <- paste0("PCoA3 [", pcoa_list$var_expl_df[1, 3], "%]")
  asp_ratio <- pcoa_list$var_expl_df$`PCo2[%]` / pcoa_list$var_expl_df$`PCo1[%]`

  pt <- pcoa_list$PCO_res %>%
    ggplot(aes(x = PCoA1, y = PCoA2, fill = Group)) +
    geom_point(size = 3, alpha = point_alpha, color = "black", pch = 21) +
    xlab(x_lab) +
    ylab(y_lab) +
    stat_ellipse(show.legend = F, aes(color = Group)) +
    theme_paper +
    theme(aspect.ratio = asp_ratio)

  return(pt)
}

f_compute_PCoA <- function(meta_df, mat, method = "bray", threshold_for_prevalence = 0, prevalence_threshold = 0.05) {
  # Takes a metadata dataframe, a species relative abundance matrix and computes a PCoA
  # Filters rel. abundance matrix to keep only feature with at least prevalence_threshold prevalence
  # A feature is considered "prevalent" in a sample if its relative abundance is above threshold_for_prevalence

  stopifnot(all(c("Sample_ID", "Group") %in% colnames(meta_df)))
  mat_clean <- mat[, meta_df$Sample_ID]

  # filter matrix to keep only features with at least 5% prevalence
  prev <- rowSums(mat_clean > threshold_for_prevalence) / ncol(mat_clean)
  relAB_sel <- mat_clean[prev >= prevalence_threshold, ]
  message(nrow(relAB_sel), " out of ", nrow(mat_clean), " features passed prevalence cutoff")

  # remove empty columns
  if (method == "bray") {
    relAB_sel <- relAB_sel[, colSums(relAB_sel) > 0]
  }
  # Compute Bray-Curtis distances
  dist <- vegan::vegdist(x = t(relAB_sel), method = method)
  pco <- labdsv::pco(dis = dist, k = 2)
  res <- as_tibble(pco$points, rownames = "Sample_ID")
  colnames(res)[2:ncol(res)] <- paste0("PCoA", 2:ncol(res) - 1)
  
  res <- res %>% left_join(., meta_df %>% dplyr::select(Group, Sample_ID))
  #* Calclualte proportion of variance explained
  var_expl <- data.frame(matrix(ncol = 3, nrow = 1))
  for (i in seq(1, 3)) {
    var_expl[1, i] <- round((pco$eig[i] / sum(pco$eig)) * 100, 1)
  }
  colnames(var_expl) <- paste0("PCo", seq(1:ncol(var_expl)), "[%]")

  # set.seed(1)
  # permanova_res <- ecole::permanova_pairwise(
  #   dist,
  #   as.character(meta_df$Group),
  #   permutations = 999,
  #   method = "bray",
  #   padj = "fdr"
  # )


  return(list(PCO_res = res, var_expl_df = var_expl))
}

f_fetch_ncbi_taxonomy <- function(tax_names) {
  #takes a vector of taxonomic names (e.g. c("Escherichia","Fusobacterium")) and queries NCBI to get full taxonomic ranks

  require(rentrez)
  require(xml2)
  require(progress)

  message("Fetching ",length(tax_names), " IDs from NCBI")

  # Initialize result_df and progress bar
  result_df <- tibble()
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = length(tax_names))
  
  max_retries <- 3    
  for (tax in tax_names) {  
    #message("\n",tax,"\n")
    success <- FALSE
    retries <- 0
    
    while (!success && retries < max_retries) {
      tryCatch({
        tax_record <- entrez_search(db = "taxonomy", term = paste0(tax, " AND (Bacteria[Organism] OR Archaea[Organism])"))
        
        if (tax_record$count > 0) {
          tax_id <- tax_record$ids[1]
          tax_fetch <- entrez_fetch(db = "taxonomy", id = tax_id, rettype = "xml")
          
          xml_doc <- read_xml(tax_fetch)
          lineage_info <- xml_find_all(xml_doc, "//LineageEx/Taxon")
          
          tax_levels <- tibble(rank=character(),scientific_name=character())          
          for (i in lineage_info) {
            scientific_name <- xml_text(xml_find_first(i, "./ScientificName"))
            rank <- xml_text(xml_find_first(i, "./Rank"))
            tax_levels <- add_row(tax_levels, rank = rank, scientific_name = scientific_name)
          }
          
          tax_levels <- add_row(tax_levels, rank = "query_tax", scientific_name = tax)
          result_df <- bind_rows(result_df, pivot_wider(tax_levels %>% filter(str_detect(rank,"phylum|class|order|family|genus|species|query_tax")), names_from = rank, values_from = scientific_name))
        } else {
          result_df <- bind_rows(result_df, tibble(query_tax = tax))
        }
        
        success <- TRUE
      }, error = function(e) {
        retries <- retries + 1
        Sys.sleep(1)
      })
    }
    
    if (!success) {
      result_df <- bind_rows(result_df, tibble(query_tax = tax))
    }
    
    Sys.sleep(0.5)
    pb$tick()
  }
  
  return(result_df)
}

f_rescale_effect_sizes <- function(mat) {
  # For heatmap plotting -> rescales effect size estimates  
  # Rescale based on 0.95 quantile (to be less sensitive to outliers)
  scale_factor <- quantile(abs(na.omit(mat)), 0.95)

  # Ensure values lie within [-1, 1]
  rescaled <- pmin(pmax(mat / scale_factor, -1), 1)
  return(rescaled)
}

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
        p.val_fisher = as.numeric(fisher_res),
        n_category_perc = 0
      ) %>%
        mutate(clin_feat_lab = paste0(clin_feat, "\n(N=", N, ")"))
    )
  }
  
  return(bar_df %>% mutate(clin_feat_lab = as_factor(clin_feat_lab)))
}

f_create_upper_triangle_matrix <- function(named_vector,condition_levels) {
  # takes a nemd vector with comparisons and returns a matrix with the comparisons in the upper triangle
  # vector must be named and comparisons separated by "_vs_"
  
  # Extracting group names
  groups <- unique(unlist(strsplit(names(named_vector), "_vs_")))
  # Creating an empty matrix
  p_matrix <- matrix(NA, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
  # Populating the matrix
  for (name in names(named_vector)) {
    parts <- unlist(strsplit(name, "_vs_"))
    p_matrix[parts[1], parts[2]] <- named_vector[name]
    p_matrix[parts[2], parts[1]] <- named_vector[name] # for symmetry
  }
  p_matrix <- p_matrix[condition_levels, condition_levels]
  p_matrix[lower.tri(p_matrix, diag = FALSE)] <- NA
  return(p_matrix)
}

f_plot_signif_matrix <- function(upper_tri_matrix, condition_levels = NULL) {
  # takes a matrix with correlation values and returns a heatmap colored by significane level

  if (is.null(condition_levels)) {
    condition_levels <- rownames(upper_tri_matrix)
  }
  sig_df <-
    upper_tri_matrix %>%
    as_tibble(rownames = "condition1") %>%
    gather(-condition1, key = "condition2", value = value) %>%
    mutate(
      condition1 = factor(condition1, levels = (condition_levels)),
      condition2 = factor(condition2, levels = rev(condition_levels))
    ) %>%
    mutate(p_sym = p_to_symbol(value))
  
  # define p-value legend
  #! Make sure this matches the assigned p-values in the p_to_symbol function
  p_val_legend <- c(
    "***" = "q<0.001",
    "**" = "q<0.01",
    "*" = "q<0.1",
    "n.s." = "q>=0.1"
  )
  p_val_colors <- c(
    "***" = "black",
    "**" = "#666666", # dark gray
    "*" = "#AAAAAA", # light gray
    "n.s." = "white"
  )  
  pt_sig_matrix <-
    sig_df %>%
    ggplot(aes(x = condition1, y = condition2)) +
    theme_paper +
    geom_tile(data = subset(sig_df, condition1 == condition2), color = NA, fill = NA) +
    geom_tile(data = subset(sig_df, !(is.na(p_sym))), aes(fill = p_sym), color = "black", size = 0.5) +
    scale_fill_manual(
      values = p_val_colors, breaks = names(p_val_legend),
      labels = p_val_legend,
      limits = names(p_val_legend), na.value = "white"
    )+
    scale_x_discrete(position = "bottom", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_abline(slope = -1, intercept = length(condition_levels)+1, size = 0.75) +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5, size = 8),
      axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(0.9, 1), legend.justification = c(1, 1),      
      legend.text = element_text(size = 6), # Adjust legend text size
      legend.title = element_blank(), # Adjust legend title size
      legend.key.size = unit(0.3, "cm")
    )     
  return(pt_sig_matrix)
}







