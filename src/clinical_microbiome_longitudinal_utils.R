# ==============================================================================
# Clinical & Microbiome Longitudinal Analysis Toolkit
# Infrastructure Functions Module
# Year: 2026
# ==============================================================================

library(tidyverse)
library(rstatix)
library(ggpubr)

# ------------------------------------------------------------------------------
# 1. Clinical Data Statistics Engine
# ------------------------------------------------------------------------------
calc_longitudinal_diff <- function(df, vars, time_col = "Time", id_col = "Clinic_ID", 
                                   comparisons, ref_time = "BL", group_col = NULL) {
  
  clinic.cate.order <-  c("K Score", "K Sub Score", "Hormone", "Glucose", "Lipid", "Liver_Renal", "Coagulation_Inflammation", "CBC", "Ultrasound")
  
  clinic.var.order <- c(
    "K_Score","FSFI","SAS","SDS","Sleep_Score",
    "Hot_flashes_and_sweating",  "Paresthesia","Insomnia","Melancholia","Nervousness" ,
    "Vertigo","Fatigue" ,"Arthralgia","Headache", "Palpitations","Formication","Sexual_dysfunction","Urinary_symptoms" ,
    
    "FSH","LH","E2", "Testosterone","TSH", 
    
    "Glu","INS", 
    "HDL_C","LDL_C" ,"TC","TG",
    "ALT","AST","BUN","Cr",
    "APTT","D_Dimer" , "Fibrinogen","INR",  "PT", "hsCRP",
    "HGB" , "Neut_Percent" ,"PLT", "RBC","WBC",
    "Left_Ovary_Diameter_1","Left_Ovary_Diameter_2","Left_Ovary_Diameter_3",                    
    "Right_Ovary_Diameter_1","Right_Ovary_Diameter_2","Right_Ovary_Diameter_3",                       
    "Uterus_Size_Diameter_1","Uterus_Size_Diameter_2","Uterus_Size_Diameter_3"  
  )
  
  if (is.null(group_col)) {
    df$Subgroup <- "Global"
    group_col <- "Subgroup"
  } else {
    df$Subgroup <- df[[group_col]]
    group_col <- "Subgroup"
  }
  
  sum_stats <- list()
  for (v in vars) {
    for (comp in comparisons) {
      temp_df <- df %>%
        select(all_of(c(id_col, time_col, group_col, v))) %>%
        rename(Time = !!sym(time_col), Group = !!sym(group_col), Value = !!sym(v)) %>%
        filter(Time %in% comp) %>%
        drop_na(Value)
      
      groups <- unique(temp_df$Group)
      for (g in groups) {
        sub_df <- temp_df %>% filter(Group == g)
        t1 <- comp[1]; t2 <- comp[2]
        
        valid_ids <- sub_df %>%
          group_by(!!sym(id_col)) %>%
          summarise(count_t1 = sum(Time == t1, na.rm = TRUE), count_t2 = sum(Time == t2, na.rm = TRUE), .groups = "drop") %>%
          filter(count_t1 == 1 & count_t2 == 1) %>% 
          pull(!!sym(id_col))
        
        if (length(valid_ids) < 3) next 
        n_paired <- length(valid_ids)
        
        plot_data <- sub_df %>%
          filter(!!sym(id_col) %in% valid_ids) %>%
          mutate(Time = factor(Time, levels = comp)) %>% 
          arrange(Time,!!sym(id_col)) 
        
        stat_res <- plot_data %>% wilcox_test(Value ~ Time, paired = TRUE,detailed = TRUE)
        wilcox_est <- if ("estimate" %in% colnames(stat_res)) stat_res$estimate else NA_real_
        
        mean_res <- plot_data %>%
          group_by(Time) %>%
          summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
          pivot_wider(names_from = Time, values_from = Mean_Value)
        
        mean1 <- mean_res[[t1]]; mean2 <- mean_res[[t2]]
        mean_diff <- mean2 - mean1
        
        res_row <- tibble(
          Factor = v, Subgroup = g, Comparison = paste(comp, collapse = "/"),
          N = n_paired,     
          Mean1 = mean1, Mean2 = mean2, F2C = log2(mean2 / mean1),
          Mean.diff.log = sign(mean_diff) * log1p(abs(mean_diff)), 
          Estimate = wilcox_est, p = stat_res$p
        )
        sum_stats[[length(sum_stats) + 1]] <- res_row
      }
    }
  }
  if (length(sum_stats) == 0) stop("No paired subjects resolved. Verify keys.")
  
  final_res <- bind_rows(sum_stats) %>%
    group_by(Subgroup) %>%                          
    mutate(p.adj = p.adjust(p, method = "fdr")) %>%  
    ungroup() %>%                                 
    mutate(
      lab.p.adj = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", p.adj < 0.05 ~ "*", TRUE ~ "ns"),
      lab.p     = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
      label     = paste(lab.p, lab.p.adj, sep = "/")
    )
  return(final_res)
}

plot_longitudinal_heatmap <- function(stat_df, phen_cate) {
  draw <- stat_df
  # draw$label <- case_when(
  #   draw$Mean.diff.log > 0 & draw$label != "ns/ns" ~ paste0(draw$label, "\n+"),
  #   draw$Mean.diff.log < 0 & draw$label != "ns/ns" ~ paste0(draw$label, "\n-"),
  #   TRUE ~ as.character(draw$label)
  # )
  draw$label <- str_replace(draw$label, "ns/ns", "")
  ## order 
  cate_levels      <- clinic.cate.order
  target_var_order <- clinic.var.order
  draw$Cate <- phen_cate$Category[match(draw$Factor, phen_cate$Phen)]
  
  draw <- draw %>%
    mutate(Cate   = factor(Cate, levels = intersect(cate_levels, unique(Cate))) ) %>%
    # mutate(Factor   = factor(Factor, levels = intersect(target_var_order, unique(Factor))) ) %>%
    arrange(Cate, match(Factor, target_var_order)  ) %>%
    mutate(Factor = fct_inorder(as.character(Factor)) )
  
  if (length(unique(draw$Subgroup)) > 1) {
    draw$Y_Axis <- paste(draw$Subgroup, draw$Comparison, sep = " - ")
  } else {
    draw$Y_Axis <- draw$Comparison
  }
  draw$Y_Axis <- factor(draw$Y_Axis, levels = unique(draw$Y_Axis))
  
  p <- ggplot(draw, aes(x = Factor, y = Y_Axis, fill = Mean.diff.log)) + 
    geom_tile(color = "lightgrey") +
    geom_text(aes(label = label), angle = 0, size = 2, color = "black") +  
    facet_grid(. ~ Cate, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#FF7F0E", midpoint = 0, name = "Log-difference") +
    scale_x_discrete(position = "top") +  scale_y_discrete(limits = rev) + 
    theme_minimal() +   
    theme(
      axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 11),
      axis.text.y     = element_text(angle = 0, size = 10), panel.grid = element_blank(),
      strip.text.x    = element_text(angle = 0, size = 9, face = "bold"), 
      strip.placement = "outside", panel.spacing.x = unit(0.3, "lines")     
    ) + labs(x = "", y = "")
  return(p)
}

run_clinic_longitudinal_pipeline <- function(df, vars_k, vars_non_k, comps_k, comps_non_k, 
                                             group_col, prefix, out_dir, phen_cate) {
  res_non_k <- calc_longitudinal_diff(df = df, vars = vars_non_k, comparisons = comps_non_k, group_col = group_col)
  if (!is.null(group_col)) res_non_k <- res_non_k %>% arrange(desc(Subgroup))
  write.csv(res_non_k, file.path(out_dir, paste0(prefix, "_phen_nonK.wilcox.csv")), row.names = FALSE)
  
  p_non_k  <- plot_longitudinal_heatmap(res_non_k, phen_cate)
  ht_non_k <- if(is.null(group_col)) 2.8 else max(3, length(unique(res_non_k$Subgroup)) * 1.5)
  ggsave(file.path(out_dir, paste0(prefix, "_wilcox.phen_nonK.heatmap.pdf")), plot = p_non_k, width = 15, height = ht_non_k)
  
  res_k <- calc_longitudinal_diff(df = df, vars = vars_k, comparisons = comps_k, group_col = group_col)
  if (!is.null(group_col)) res_k <- res_k %>% arrange(desc(Subgroup))
  write.csv(res_k, file.path(out_dir, paste0(prefix, "_phen_K.wilcox.csv")), row.names = FALSE)
  
  p_k  <- plot_longitudinal_heatmap(res_k, phen_cate)
  ht_k <- if(is.null(group_col)) 4 else max(6, length(unique(res_k$Comparison)) * length(unique(res_k$Subgroup)) * 0.5)
  ggsave(file.path(out_dir, paste0(prefix, "_wilcox.phen_K.heatmap.pdf")), plot = p_k, width = 8, height = ht_k)
}

# ------------------------------------------------------------------------------
# 2. Microbiome Profiles Data Statistics Engine
# ------------------------------------------------------------------------------
calc_microbiome_diff <- function(prof_list, meta_df, threshold_map, 
                                 comparisons, group_col = NULL, valid_4_if = TRUE) {

  if (is.null(group_col)) {
    meta_df$Subgroup <- "Global"
    group_col <- "Subgroup"
  } else {
    meta_df$Subgroup <- meta_df[[group_col]]
    group_col <- "Subgroup"
  }
  
  sum_stats <- list(); raw_data_for_plot_list <- list() 
  
  for (site in names(threshold_map)) {
    rel_ab.threshold <- threshold_map[[site]]
    
    if (!site %in% names(prof_list)) next
    abun_matrix      <- prof_list[[site]] %>% data.frame()
    if (is.null(abun_matrix)) next
    
    taxa_to_keep <- sapply(abun_matrix, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(abun_matrix))
    taxa_vars    <- names(taxa_to_keep)[taxa_to_keep]

    if (length(taxa_vars) == 0) next
    
    # abun_matrix  <- abun_matrix[, taxa_vars, drop = FALSE]
    abun_matrix_filtered  <- abun_matrix[, taxa_vars, drop = FALSE]
    
    pseu_value      <- min(abun_matrix_filtered[abun_matrix_filtered != 0]) * 0.1
    abun_matrix_log <- log10(abun_matrix_filtered + pseu_value)
    zero_log_val    <- log10(pseu_value)
    
    df_merged <- abun_matrix_log %>% as.data.frame() %>% rownames_to_column("SeqID") %>%
      inner_join(meta_df %>% dplyr::select(SeqID, Clinic_ID, Time, Subgroup), by = "SeqID") %>%
      rename(Group = Subgroup) %>% 
      dplyr::filter(!is.na(Group))
    
    groups <- unique(df_merged$Group)
    for (g in groups) {
      sub_df <- df_merged %>% dplyr::filter(Group == g)
      if (valid_4_if){
        valid_4_ids <- sub_df %>% group_by(Clinic_ID) %>% summarise(n_time = n_distinct(Time), .groups = "drop") %>% filter(n_time == 4) %>% pull(Clinic_ID)
        sub_df <- sub_df %>% filter(Clinic_ID %in% valid_4_ids)
      }
      if (nrow(sub_df) == 0) next
      
      for (comp in comparisons) {
        t1 <- comp[1]; t2 <- comp[2]
        comp_df <- sub_df %>% dplyr::filter(Time %in% comp)
        
        valid_pairs <- comp_df %>% group_by(Clinic_ID) %>% summarise(c1 = sum(Time == t1), c2 = sum(Time == t2), .groups = "drop") %>% dplyr::filter(c1 == 1 & c2 == 1) %>% pull(Clinic_ID)
        if (length(valid_pairs) < 3) next
        
        # plot_data <- comp_df %>% dplyr::filter(Clinic_ID %in% valid_pairs) %>% mutate(Time = factor(Time, levels = comp)) %>% arrange(Clinic_ID, Time) %>%
        #   pivot_longer(cols = all_of(taxa_vars), names_to = "Factor", values_to = "Value")
        plot_data <- comp_df %>% dplyr::filter(Clinic_ID %in% valid_pairs) %>% mutate(Time = factor(Time, levels = comp)) %>% arrange(Clinic_ID, Time) %>%
          pivot_longer(cols = all_of(taxa_vars), names_to = "Factor", values_to = "Value") %>%
          group_by(Factor, Clinic_ID) %>%
          dplyr::filter(!all(abs(Value - zero_log_val) < 1e-9)) %>%
          ungroup()
        
        if (nrow(plot_data) == 0) next
        
        valid_taxa_counts <- plot_data %>% 
          group_by(Factor) %>% 
          summarise(n_paired = n_distinct(Clinic_ID), .groups = "drop") %>% 
          filter(n_paired >= 3) %>% 
          pull(Factor)
        if (length(valid_taxa_counts) == 0) next
        
        # 只保留满足配对数 >= 3 的物种进行后续统计
        plot_data <- plot_data %>% filter(Factor %in% valid_taxa_counts)
        
        n_map <- plot_data %>% 
          group_by(Factor) %>% 
          summarise(N = n_distinct(Clinic_ID), .groups = "drop")
        
        stat_res <- plot_data %>% group_by(Factor) %>% wilcox_test(Value ~ Time, paired = TRUE, detailed = TRUE)
        if (!"estimate" %in% colnames(stat_res)) {
          stat_res$estimate <- NA_real_
        }
        
        mean_res <- plot_data %>% group_by(Factor, Time) %>% summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = 'drop') %>% pivot_wider(names_from = Time, values_from = Mean_Value)
        
        res_merged <- stat_res %>% left_join(mean_res, by = "Factor") %>% left_join(n_map, by = "Factor") %>%
          mutate(
            Site = site, 
            Subgroup = g, 
            Comparison = paste(comp, collapse = "/"), 
            Mean1 = !!sym(t1), 
            Mean2 = !!sym(t2), 
            Mean.diff = Mean2 - Mean1, 
            Mean.diff.log = sign(Mean.diff) * log1p(abs(Mean.diff))
          ) %>%
          dplyr::select(Site, Factor, Subgroup, Comparison, N, Mean1, Mean2, Mean.diff, Mean.diff.log, estimate, p)
        
        sum_stats[[length(sum_stats) + 1]] <- res_merged
        raw_data_cleaned <- plot_data %>% mutate(Site = site, Subgroup = g, Comparison = paste(comp, collapse = "/")) %>% 
          dplyr::select(Site, Subgroup, Comparison, Clinic_ID, SeqID, Time, Factor, Value)
        raw_data_for_plot_list[[paste(site, g, paste(comp, collapse = "vs"), sep = "_")]] <- raw_data_cleaned
      }
    }
  }
  if (length(sum_stats) == 0) return(NULL)

  final_stats <- bind_rows(sum_stats) %>% group_by(Site, Comparison, Subgroup) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    ungroup() %>%
    mutate(
      lab.p.adj = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", p.adj < 0.05 ~ "*", TRUE ~ "ns"),
      lab.p     = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
      label     = paste(lab.p, lab.p.adj, sep = "/")
    )
  
  all_raw_data <- bind_rows(raw_data_for_plot_list)
  sig_combinations <- final_stats %>% dplyr::filter(p < 0.05) %>% dplyr::select(Site, Factor, Comparison) %>% distinct()
  final_plot_raw_data <- if (nrow(sig_combinations) > 0) all_raw_data %>% inner_join(sig_combinations, by = c("Site", "Factor", "Comparison")) else data.frame()
  
  attr(final_stats, "plot_raw_data") <- final_plot_raw_data
  return(final_stats)
}

plot_microbiome_heatmap <- function(stat_df, group_levels = NULL) {
  sig_taxa <- stat_df %>% filter(p.adj < 0.05) %>% select(Site, Factor) %>% distinct()
  draw     <- stat_df %>% inner_join(sig_taxa, by = c("Site", "Factor"))
  if (nrow(draw) == 0) return(NULL)
  
  # draw$label <- case_when(
  #   draw$Mean.diff.log > 0 & draw$label != "ns/ns" ~ paste0(draw$label, "\n+"),
  #   draw$Mean.diff.log < 0 & draw$label != "ns/ns" ~ paste0(draw$label, "\n-"),
  #   TRUE ~ as.character(draw$label)
  # )
  draw$label <- str_replace(draw$label, "ns/ns", "")
  
  site_levels <- c("VA", "UR", "GUT", "TO")
  draw$Site   <- factor(draw$Site, levels = intersect(site_levels, unique(draw$Site)))
  has_subgroups <- "Subgroup" %in% colnames(draw) && length(unique(draw$Subgroup)) > 1 && !all(draw$Subgroup == "Global")
  
  if (has_subgroups) {
    draw$Subgroup   <- if (!is.null(group_levels)) factor(draw$Subgroup, levels = group_levels) else factor(draw$Subgroup, levels = unique(draw$Subgroup))
    draw$Y_Axis     <- factor(draw$Comparison, levels = unique(draw$Comparison))
    facet_layer     <- facet_grid(Subgroup ~ Site, scales = "free", space = "free", switch = "x")
  } else {
    draw$Y_Axis     <- factor(draw$Comparison, levels = unique(draw$Comparison))
    facet_layer     <- facet_grid(. ~ Site, scales = "free_x", space = "free_x", switch = "x")
  }
  
  p <- ggplot(draw, aes(x = Factor, y = Y_Axis, fill = Mean.diff.log)) + 
    geom_tile(color = "lightgrey") + geom_text(aes(label = label), angle = 0, size = 2, color = "black") +  
    facet_layer + scale_fill_gradient2(low = "#5e3c99", mid = "white", high = "#b35806", midpoint = 0, name = "Log-difference") +
    scale_x_discrete(position = "top") + scale_y_discrete(limits = rev) + theme_minimal() +   
    theme(
      axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10, face = "bold.italic"),
      axis.text.y     = element_text(angle = 0, size = 10), panel.grid = element_blank(),
      strip.text.x    = element_text(angle = 0, size = 10, face = "bold"), strip.text.y = element_text(angle = -90, size = 10, face = "bold"),
      strip.placement = "outside", panel.spacing = unit(0.5, "lines")     
    ) + labs(x = "", y = "")
  return(p)
}

plot_microbiome_boxplot_sig <- function(boxplot_raw_data, out_dir, prefix) {
  if (is.null(boxplot_raw_data) || nrow(boxplot_raw_data) == 0) return(invisible(NULL))
  slices <- boxplot_raw_data %>% dplyr::select(Site, Comparison) %>% distinct()
  
  for (i in 1:nrow(slices)) {
    cur_site <- as.character(slices$Site[i]); cur_comp <- as.character(slices$Comparison[i])
    data_sub <- boxplot_raw_data %>% dplyr::filter(Site == cur_site, Comparison == cur_comp)
    if (nrow(data_sub) == 0) next
    
    sample_counts <- data_sub %>% group_by(Subgroup) %>% summarise(n_patients = n_distinct(Clinic_ID), .groups = "drop")
    n_nr          <- sample_counts %>% filter(Subgroup == "NR") %>% pull(n_patients) %>% {if(length(.) == 0) 0 else .}
    n_r           <- sample_counts %>% filter(Subgroup == "R") %>% pull(n_patients) %>% {if(length(.) == 0) 0 else .}
    
    p_box <- ggplot(data_sub, aes(x = Time, y = Value, fill = Time)) +
      # geom_line(aes(group = Clinic_ID), color = "grey85", alpha = 0.5, linewidth = 0.4) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.4, color = "black") +
      geom_jitter(width = 0.12, size = 1, alpha = 0.7, aes(color = Time)) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8, fill = "#D95F02", color = "black", stroke = 0.8, show.legend = FALSE) +
      facet_grid(Subgroup ~ gsub("^GUT_", "", Factor), scales = "free_y", axes = "all") +
      scale_fill_manual(values = c("#9ecae1", "#3182bd")) + scale_color_manual(values = c("#9ecae1", "#3182bd")) +
      labs(x = "Timeline", y = "Microbial Abundance (Log10)", title = paste("Significant Longitudinal Shifts -", cur_site, paste0("(", cur_comp, ")")), caption = paste0("Sample size: NR (n = ", n_nr, "), R (n = ", n_r, ")")) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5), panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 8.5, face = "italic", color = "grey30", hjust = 1, margin = margin(t = 10)),
        strip.background = element_rect(fill = "white", color = "black"), strip.text.x = element_text(face = "bold.italic", size = 9), strip.text.y = element_text(face = "bold", size = 10)
      )
    ggsave(file.path(out_dir, paste0(prefix, ".Boxplot.", cur_site, "_", gsub("/", "_vs_", cur_comp), ".pdf")), plot = p_box, width = max(5, n_distinct(data_sub$Factor) * 2 + 1.5), height = if(length(unique(data_sub$Subgroup)) > 1) 5.5 else 3.5, limitsize = FALSE)
  }
}




run_microbiome_pipeline = function(prof_list, meta_df, threshold_map, 
                                   comparisons, group_col, prefix, out_dir, valid_4_if = TRUE) {
  res <- calc_microbiome_diff(prof_list, meta_df, threshold_map, comparisons, group_col, valid_4_if = valid_4_if)
  if (is.null(res)) return(invisible(NULL))
  
  grp_levels <- if (!is.null(group_col) && is.factor(meta_df[[group_col]])) levels(meta_df[[group_col]]) else NULL
  res$Subgroup <- if (!is.null(group_col) && !is.null(grp_levels)) factor(res$Subgroup, levels = grp_levels) else res$Subgroup
  res <- res %>% arrange(if(!is.null(group_col)) Subgroup else Site, Site, Factor)
  write.csv(res, file.path(out_dir, paste0(prefix, ".TimeDiff.wilcox.csv")), row.names = FALSE)
  
  p_heatmap <- plot_microbiome_heatmap(res, group_levels = grp_levels)
  if (!is.null(p_heatmap)) ggsave(file.path(out_dir, paste0(prefix, ".TimeDiff.wilcox.heatmap.pdf")), plot = p_heatmap, width = max(7, n_distinct(res %>% filter(p.adj < 0.05) %>% pull(Factor)) * 0.25 + 2), height = if(is.null(group_col)) 4 else max(5, length(unique(res$Subgroup)) * length(comparisons) * 0.5), limitsize = FALSE)
  # if (!is.null(attr(res, "plot_raw_data")) && nrow(attr(res, "plot_raw_data")) > 0) plot_microbiome_boxplot_sig(attr(res, "plot_raw_data"), out_dir, prefix)
}