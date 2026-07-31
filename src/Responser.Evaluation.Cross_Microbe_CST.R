# ==============================================================================
# Pipeline: Cross-Sectional Microbiome Features & CST Stratified Analysis
# Modules: Abundance Wilcoxon, CST Proportions (Fisher/Chi-sq) & Logistic Models
# Year: 2026
# ==============================================================================

library(tidyverse)
library(rstatix)
library(ggpubr)
library(broom)
library(conflicted)
library(vegan)
library(colorspace)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

# Shared meta-subgroup mapping configuration matrices
subgroups_config <- list(
  "Hormone_E2" = phen.gr.E2.all %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "mKI"        = phen.gr.K      %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))))

# filter_strategy    <- "core_BL_24"
filter_strategy    <- "all"

all_sites          <- c("VA", "UR", "GUT", "TO")
time_points_vec    <- c("BL", "T04", "T12", "T24")
cst_ordered_factor <- c("GUT-P.cop", "GUT-P.v", "TO-N.s", "TO-P.m", "UROG-Div", "UROG-G.v", "UROG-L.i", "UROG-L.c")
threshold_map      <- list("GUT" = 0.001, "VA" = 0.0001, "TO" = 0.001, "UR" = 0.0001)

# Custom color palette configuration
cst_colors <- c(
  "GUT-P.v"   = "#66c2a4", "GUT-P.cop" = "#b2e2e2",
  "TO-P.m"    = "#fc8d59", "TO-N.s"   = "#fdcc8a",
  "UROG-L.c"  = "#117733", "UROG-L.i" = "#74C476", 
  "UROG-G.v"  = "#FF7F0E", "UROG-Div" = "#2171B5"
)

get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"   
  )}
# ==============================================================================
# Helper Functions Infrastructure Module
# ==============================================================================
calc_outcome_diff_by_time <- function(prof_list, meta_df, threshold_map, time_points, outcome_col) {
  sum_stats <- list(); raw_data_for_plot_list <- list()
  
  patient_time_status <- meta_df %>%
    dplyr::filter(!is.na(Time)) %>%
    group_by(Clinic_ID) %>%
    summarise(available_times = list(unique(Time)), .groups = "drop") %>%
    mutate(keep_sample = sapply(available_times, function(times) {
      keep <- FALSE
      if (filter_strategy == "strict_4") {
        if (all(c("BL", "T04", "T12", "T24") %in% times)) keep <- TRUE
      } else if (filter_strategy == "core_BL_24") {
        if (all(c("BL", "T24") %in% times)) keep <- TRUE
      } else if (filter_strategy == "all") {
        keep <- TRUE
      } else {
        stop(paste("Unknown filter strategy:", filter_strategy))
      }
      return(keep)
    })) %>%
  dplyr::filter(keep_sample == TRUE) %>%
  dplyr::select(Clinic_ID)  
  
  for (site in names(threshold_map)) {
    rel_ab.threshold <- threshold_map[[site]]
    abun_matrix      <- prof_list[[site]] %>% data.frame()
    if (is.null(abun_matrix)) next
    
    taxa_to_keep <- sapply(abun_matrix, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(abun_matrix))
    taxa_vars    <- names(taxa_to_keep)[taxa_to_keep]
    abun_matrix  <- abun_matrix[, taxa_vars, drop = FALSE]
    if (length(taxa_vars) == 0) next
    
    pseu_value      <- min(abun_matrix[abun_matrix != 0]) * 0.1
    abun_matrix_log <- log10(abun_matrix + pseu_value)
    
    df <- abun_matrix_log %>% as.data.frame() %>% rownames_to_column("SeqID") %>%  
      inner_join(meta_df %>% dplyr::select(SeqID, Clinic_ID, Time, !!sym(outcome_col)), by = "SeqID") %>%
      rename(Outcome = !!sym(outcome_col)) %>%
      inner_join(patient_time_status, by = "Clinic_ID") %>%
      dplyr::filter(!is.na(Outcome), Time %in% time_points)
    
    for (tp in time_points) {
      tp_df <- df %>% dplyr::filter(Time == tp)
      if (length(unique(tp_df$Outcome)) < 2) next
      
      group_counts <- tp_df %>% group_by(Outcome) %>% summarise(n = n(), .groups = "drop")
      if (any(group_counts$n < 3)) next
      
      plot_data <- tp_df %>% mutate(Outcome = factor(Outcome, levels = levels(meta_df[[outcome_col]]))) %>%
        pivot_longer(cols = all_of(taxa_vars), names_to = "Factor", values_to = "Value")
      
      stat_res <- plot_data %>% group_by(Factor) %>% wilcox_test(Value ~ Outcome, paired = FALSE, detailed = TRUE) %>% dplyr::select(Factor, p)
      mean_levels <- levels(plot_data$Outcome)
      
      mean_res <- plot_data %>% group_by(Factor, Outcome) %>% summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
        pivot_wider(names_from = Outcome, values_from = Mean_Value)
      
      for(lvl in mean_levels) { if (!lvl %in% colnames(mean_res)) mean_res[[lvl]] <- NA }
      
      # 动态提取对应两个比对组的真实样本量大小
      n_g1 <- group_counts$n[group_counts$Outcome == mean_levels[1]]
      n_g2 <- group_counts$n[group_counts$Outcome == mean_levels[2]]
      
      res_merged <- stat_res %>% left_join(mean_res, by = "Factor") %>%
        mutate(
          Groups        = paste0("Group.", site),
          Group1        = mean_levels[1],
          Group2        = mean_levels[2],
          Sites         = site,
          Species       = Factor,
          Timepoint     = tp,
          SampleNum1    = n_g1,
          SampleNum2    = n_g2,
          `Mean1 (Log10)` = !!sym(mean_levels[1]),
          `Mean2 (Log10)` = !!sym(mean_levels[2]),
          `P value` = p, 
          
          # 2. 完全保留原有底层老列名，作为核心底盘无损输送给下游可视化函数
          p             = p,
          Site          = site, 
          TimePoint     = tp, 
          Factor        = Factor,
          Mean_R        = !!sym(mean_levels[1]), 
          Mean_NR       = !!sym(mean_levels[2]),
          Mean.diff     = Mean_R - Mean_NR, 
          Mean.diff.log = sign(Mean.diff) * log1p(abs(Mean.diff))
        )
      
      sum_stats[[length(sum_stats) + 1]] <- res_merged
      raw_data_for_plot_list[[paste(site, tp, sep = "_")]] <- plot_data %>% mutate(Site = site, TimePoint = tp) %>%
        dplyr::select(Site, TimePoint, Clinic_ID, SeqID, Outcome, Factor, Value)
    }
  }
  if (length(sum_stats) == 0) return(NULL)
  final_stats <- bind_rows(sum_stats) %>%
    group_by(Site) %>% 
    mutate(Adjusted_P_value = p.adjust(p, method = "BH")) %>% 
    ungroup() %>%
    mutate(
      p_star       = get_significance_stars(p),
      q_star       = get_significance_stars(Adjusted_P_value),
      raw_mark = paste(p_star, q_star, sep = "/"),
      Significance = if_else(p > 0.05, "", raw_mark),
      
      # 兼容旧代码热图逻辑所需的映射标签
      p.adj        = Adjusted_P_value,
      p_mark       = Significance
    ) %>%
    # 3. 按照新截图格式强行控制第一视角前置列序，作图私有变量优雅顺延至尾部
    dplyr::select(
      Groups, Group1, Group2, Sites, Species, Timepoint, SampleNum1, SampleNum2, 
      `Mean1 (Log10)`, `Mean2 (Log10)`, `P value`, `Adjusted P value` = Adjusted_P_value, Significance,
      everything()
    )
  
  all_raw_data <- bind_rows(raw_data_for_plot_list)
  sig_combinations <- final_stats %>% dplyr::filter(`P value` < 0.05) %>% dplyr::select(Site, TimePoint, Factor) %>% distinct()
  attr(final_stats, "plot_raw_data") <- if (nrow(sig_combinations) > 0) all_raw_data %>% inner_join(sig_combinations, by = c("Site", "TimePoint", "Factor")) else data.frame()
  
  return(final_stats)
}

plot_outcome_heatmap <- function(stat_df) {
  sig_taxa <- stat_df %>% dplyr::filter(p < 0.05) %>% dplyr::select(Site, Factor) %>% distinct()
  draw     <- stat_df %>% inner_join(sig_taxa, by = c("Site", "Factor"))
  if (nrow(draw) == 0) return(NULL)
  
  draw$Site         <- factor(draw$Site, levels = c("VA", "UR", "TO", "GUT"))
  draw$TimePoint    <- factor(draw$TimePoint, levels = c("BL", "T04", "T12", "T24"))
  draw$clean_factor <- gsub("^GUT_|^VA_|^UR_|^TO_", "", draw$Factor)
  
  p <- ggplot(draw, aes(x = clean_factor, y = TimePoint, fill = Mean.diff.log)) + 
    geom_tile(color = "lightgrey") +
    geom_text(aes(label = p_mark), size = 1.8, color = "black", fontface = "bold", lineheight = 0.85) +  
    facet_grid(. ~ Site, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_gradient2(low = "#5e3c99", mid = "white", high = "#b35806", midpoint = 0, name = "Log Diff\n(R - NR)") +
scale_x_discrete(position = "top") + scale_y_discrete(limits = rev) + theme_minimal() +   
  theme(
    axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 9, face = "bold.italic"),
    axis.text.y     = element_text(angle = 0, size = 10, face = "bold"), panel.grid = element_blank(),
    strip.text.x    = element_text(angle = 0, size = 11, face = "bold"), strip.placement = "outside", panel.spacing = unit(0.5, "lines")     
  ) + labs(x = "", y = "Timeline")
return(p)
}

plot_outcome_boxplot_sig <- function(boxplot_raw_data, out_dir, prefix, box_palette) {
  if (is.null(boxplot_raw_data) || nrow(boxplot_raw_data) == 0) return(invisible(NULL))
  slices <- boxplot_raw_data %>% dplyr::select(Site, TimePoint) %>% distinct()
  
  for (i in 1:nrow(slices)) {
    cur_site <- as.character(slices$Site[i]); cur_tp <- as.character(slices$TimePoint[i])
    data_sub <- boxplot_raw_data %>% dplyr::filter(Site == cur_site, TimePoint == cur_tp)
    if (nrow(data_sub) == 0) next
    
    sample_counts         <- data_sub %>% group_by(Outcome) %>% summarise(n_patients = n_distinct(Clinic_ID), .groups = "drop")
    mean_levels           <- levels(data_sub$Outcome)
    n_nr                  <- sample_counts %>% filter(Outcome == mean_levels[2]) %>% pull(n_patients) %>% {if(length(.) == 0) 0 else .}
    n_r                   <- sample_counts %>% filter(Outcome == mean_levels[1]) %>% pull(n_patients) %>% {if(length(.) == 0) 0 else .}
    data_sub$clean_factor <- gsub("^GUT_|^VA_|^UR_|^TO_", "", data_sub$Factor)
    
    p_box <- ggplot(data_sub, aes(x = Outcome, y = Value, fill = Outcome)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.4, color = "black") +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, aes(color = Outcome)) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8, fill = "#D95F02", color = "black", stroke = 0.8, show.legend = FALSE) +
      facet_wrap(~ clean_factor, scales = "free_y", ncol = 4) +
      theme_bw() +
      scale_fill_manual(values = box_palette) + scale_color_manual(values = box_palette) +
      labs(x = "Clinical Response Status", y = "Abundance (Log10)", 
           title = paste("Cross-Sectional Abundance Shifts -", cur_site, paste0("(", cur_tp, ")")), 
           caption = sprintf("Sample size: %s (n = %d), %s (n = %d)", mean_levels[2], n_nr, mean_levels[1], n_r)) +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5), panel.grid.minor = element_blank(), legend.position = "none",
        plot.caption = element_text(size = 8.5, face = "italic", color = "grey30", hjust = 1, margin = margin(t = 10)),
        strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(face = "bold.italic", size = 9.5)
      )
    ggsave(file.path(out_dir, paste0(prefix, ".Boxplot.", cur_site, "_", cur_tp, ".pdf")), plot = p_box, 
           width = max(4.5, n_distinct(data_sub$Factor) * 2.0 + 1.2), height = 3.8, limitsize = FALSE)
  }
}

# ==============================================================================
# 3. Main Operational Pipeline Iteration Execution
# ==============================================================================
for (sg_name in names(subgroups_config)) {
  
  # Establish consolidated file output tracking directories
  target_dir_cst <- file.path(base_dir, "Results/Responser", sg_name, "Cross_Microbe_CST")
  if(!dir.exists(target_dir_cst)) dir.create(target_dir_cst, recursive = TRUE)
  setwd(target_dir_cst)
  
  current_gr_df <- subgroups_config[[sg_name]]
  target_levels <- levels(current_gr_df$group)
  box_palette   <- setNames(c("#C7A6FA", "#E36C5B")[1:length(target_levels)], target_levels)
  
  # Align primary microbiome metadata objects
  mg_meta_outcome <- Microbe.phen.prof %>%
    dplyr::select(-any_of(c("group"))) %>% 
    left_join(current_gr_df, by = "Clinic_ID") %>%
    filter(!is.na(group))
  
  # ----------------------------------------------------------------------------
  # PART I: Cross-Sectional Taxa Abundance Wilcoxon Profiling Pipeline #####################
  # ----------------------------------------------------------------------------
  res_taxa <- calc_outcome_diff_by_time(
    prof_list = prof_filtered, meta_df = mg_meta_outcome, 
    threshold_map = threshold_map, time_points = time_points_vec, outcome_col = "group"
  )
  
  if (!is.null(res_taxa)) {
    write.csv(res_taxa %>% arrange(Site, TimePoint, p), sprintf("%s_OutcomeDiff.wilcox.%s.csv", sg_name, filter_strategy), row.names = FALSE)
    
    p_heat <- plot_outcome_heatmap(res_taxa)
    if (!is.null(p_heat)) {
      ggsave(sprintf("%s_OutcomeDiff.wilcox.heatmap.%s.pdf", sg_name, filter_strategy), plot = p_heat, 
             width = max(7, n_distinct(res_taxa %>% filter(p < 0.05) %>% pull(Factor)) * 0.25 + 2), 
             height = length(time_points_vec) * 0.8 + 1.5, limitsize = FALSE)
    }
    
    # plot_outcome_boxplot_sig(attr(res_taxa, "plot_raw_data"), target_dir_cst, sg_name, box_palette)
  }
  
  # ----------------------------------------------------------------------------
  # PART II: CST Distribution Cross-Proportion Framework (Fisher/Chi-Sq) #######################
  # ----------------------------------------------------------------------------
  df_cst_pre <- mg_meta_outcome %>%
    left_join(prof_diversity %>% select(SeqID, CST), by = "SeqID") %>%
    dplyr::select(Clinic_ID, Time, Site, CST, group) %>%
    drop_na(Time, Site, CST, group)
  
  cst_whitelist <- df_cst_pre %>%
    group_by(Site, Clinic_ID) %>%
    summarise(available_times = list(unique(Time)), .groups = "drop") %>%
    mutate(keep_sample = sapply(available_times, function(times) {
      keep <- FALSE
      if (filter_strategy == "strict_4") {
        if (all(c("BL", "T04", "T12", "T24") %in% times)) keep <- TRUE
      } else if (filter_strategy == "core_BL_24") {
        if (all(c("BL", "T24") %in% times)) keep <- TRUE
      } else if (filter_strategy == "all") {
        keep <- TRUE
      } else {
        stop(paste("Unknown filter strategy:", filter_strategy))
      }
      return(keep)
    })) %>%
    dplyr::filter(keep_sample == TRUE) %>%
    dplyr::select(Site, Clinic_ID) 
  
  df_cst <- df_cst_pre %>%
    inner_join(cst_whitelist, by = c("Site", "Clinic_ID")) %>%
    mutate(
      Time = factor(Time, levels = c("BL", "T04", "T12", "T24")),
      Site = factor(Site, levels = all_sites),
      CST  = factor(CST, levels = cst_ordered_factor)
    ) %>% 
    drop_na(Time)
  
  cst_stats <- df_cst %>%
    group_by(Site, Time) %>%
    group_modify(~ {
      current_df <- .x %>% droplevels()
      tab <- table(current_df$group, current_df$CST)
      if (nrow(tab) < 2 || ncol(tab) < 2 || any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
        return(data.frame(p = NA, method = "None"))
      }
      exp_freq <- suppressWarnings(chisq.test(tab)$expected)
      
      if (any(exp_freq < 5)) {
        p_val <- fisher.test(tab)$p.value
        return(data.frame(p = p_val, method = "Fisher"))
      } else {
        p_val <- suppressWarnings(chisq.test(tab)$p.value)
        return(data.frame(p = p_val, method = "Chi-square"))
      }
    }) %>%
    ungroup() %>%
    group_by(Site) %>%
    mutate(p.adj = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    mutate(
      p_star   = get_significance_stars(p),
      q_star   = get_significance_stars(p.adj),
      raw_mark = paste(p_star, q_star, sep = "/"),
      p_mark   = if_else(p > 0.05, "", raw_mark),
      y_pos    = 1.03
    )
  
  df_plot_cst <- df_cst %>%
    group_by(Site, Time, group, CST) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(Total = sum(count),Percentage = count / sum(count)) %>%
    ungroup()
  
  p_cst_ratio <- ggplot(df_plot_cst, aes(x = group, y = Percentage, fill = CST)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
    facet_grid(Site ~ Time, scales = "free_y", axes = "all") +
    scale_fill_manual(values = cst_colors, name = "CST Type", drop = TRUE) +
    scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Relative Proportion of CST", title = sprintf("CST Compositional Proportions [%s]", sg_name)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5), strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 11, face = "bold"), axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 9, color = "black"), axis.ticks.y = element_line(color = "black"), panel.grid = element_blank()
    ) +
    geom_text(data = cst_stats, aes(x = 1.5, y = y_pos, label = p_mark), inherit.aes = FALSE, vjust = -0.1, size = 2.2, fontface = "bold", color = "black", lineheight = 0.85)
  
  ggsave(sprintf("%s_CST_Cross_Proportion_Filtered_%s.pdf", sg_name,filter_strategy), plot = p_cst_ratio, width = 9, height = 8.5)
  write.csv(df_plot_cst, sprintf("%s_CST_Cross_Outcome_Prop_%s.csv", sg_name,filter_strategy), row.names = FALSE)
  write.csv(cst_stats, sprintf("%s_CST_Cross_Outcome_Prop_chiq_%s.csv", sg_name, filter_strategy), row.names = FALSE)
  # ----------------------------------------------------------------------------
  # PART III: Cross-Sectional CST Binary Logistic Regression Modeling ##########################
  # ----------------------------------------------------------------------------
  logistic_whitelist <- df_cst_pre %>%
    group_by(Site, Clinic_ID) %>%
    summarise(available_times = list(unique(Time)), .groups = "drop") %>%
    mutate(keep_sample = sapply(available_times, function(times) {
      keep <- FALSE
      if (filter_strategy == "strict_4") {
        if (all(c("BL", "T04", "T12", "T24") %in% times)) keep <- TRUE
      } else if (filter_strategy == "core_BL_24") {
        if (all(c("BL", "T24") %in% times)) keep <- TRUE
      } else if (filter_strategy == "all") {
        keep <- TRUE
      } else {
        stop(paste("Unknown filter strategy:", filter_strategy))
      }
      return(keep)
    })) %>%
    dplyr::filter(keep_sample == TRUE) %>%
    dplyr::select(Site, Clinic_ID)
  
  df_logistic_prep <- df_cst_pre %>%
    inner_join(logistic_whitelist, by = c("Site", "Clinic_ID")) %>%
    mutate(
      Time = factor(Time, levels = c("BL", "T04", "T12", "T24")),
      Site = factor(Site, levels = all_sites), 
      CST  = as.character(CST),

      Outcome_Binary = if_else(group == target_levels[1], 1, 0)
    ) %>% 
    drop_na(Time)
  
  
  logistic_results_list <- list()
  
  for (st in levels(df_logistic_prep$Site)) {
    for (tp in levels(df_logistic_prep$Time)) {
      slice_df <- df_logistic_prep %>% dplyr::filter(Site == st, Time == tp)
      if (nrow(slice_df) == 0) next
      
      unique_csts <- unique(slice_df$CST)
      if (length(unique_csts) < 2) next 
      
      # 📌 修复1：使用原生安全计数，抓取当前切片大盘总响应(1)与不响应(0)人数
      n_total_outcome_1 <- sum(slice_df$Outcome_Binary == 1)
      n_total_outcome_0 <- sum(slice_df$Outcome_Binary == 0)
      
      # 如果大盘切片本身就缺乏某一个分类（比如全都是响应者），则无法建立逻辑回归，跳过
      if (n_total_outcome_1 == 0 || n_total_outcome_0 == 0) next
      
      for (target_cst in unique_csts) {
        run_df <- slice_df %>% mutate(CST_Binary = if_else(CST == target_cst, 1, 0))
        
        # 📌 修复2：用标准逻辑条件求和，彻底打掉不稳定的 table 矩阵字符索引漏洞
        # TP: 携带该 CST 且 临床响应
        # FP: 携带该 CST 但 临床不响应
        TP <- sum(run_df$CST_Binary == 1 & run_df$Outcome_Binary == 1)
        FP <- sum(run_df$CST_Binary == 1 & run_df$Outcome_Binary == 0)
        
        # 🚨 安全过滤：如果在目标 CST 携带者中，有任意一组人数为 0（Zero cell 导致 OR 无法计算或无穷大），则跳过
        if (TP == 0 || FP == 0) next
        
        # 📌 核心新增：基于四格表联动推导诊断效能指标
        FN <- n_total_outcome_1 - TP  # 未携带该 CST 但 临床响应
        TN <- n_total_outcome_0 - FP  # 未携带该 CST 且 临床不响应
        
        sensitivity <- TP / (TP + FN)  # 敏感性
        specificity <- TN / (TN + FP)  # 特异性
        ppv         <- TP / (TP + FP)  # 阳性预测值（即该 CST 携带者的响应率）
        npv         <- TN / (TN + FN)  # 阴性预测值
        
        # 建立一元二分类 Logistic 模型
        model_stats <- tidy(glm(Outcome_Binary ~ CST_Binary, data = run_df, family = binomial(link = "logit")), 
                            exponentiate = TRUE, conf.int = TRUE) %>%
          dplyr::filter(term == "CST_Binary") 
        
        if (nrow(model_stats) == 0) next
        
        logistic_results_list[[paste(st, tp, target_cst, sep = "_")]] <- data.frame(
          Site = st, 
          TimePoint = tp, 
          Target_CST = target_cst, 
          
          # 样本量追踪字段
          N_Target_Outcome_0 = FP,  # 对应原表名：目标 CST 中 Outcome 0 人数
          N_Target_Outcome_1 = TP,  # 对应原表名：目标 CST 中 Outcome 1 人数
          N_Total_Outcome_0  = n_total_outcome_0,
          N_Total_Outcome_1  = n_total_outcome_1,
          
          # 📌 核心新增：预测效能/诊断效能指标输出
          Sensitivity = round(sensitivity, 3),
          Specificity = round(specificity, 3),
          PPV         = round(ppv, 3),
          NPV         = round(npv, 3),
          
          Odds_Ratio = model_stats$estimate,
          CI_Lower   = model_stats$conf.low, 
          CI_Upper   = model_stats$conf.high, 
          p_value    = model_stats$p.value, 
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(logistic_results_list) > 0) {
    all_logistic_report <- bind_rows(logistic_results_list) %>%
      group_by(Site) %>%
      mutate(p.adj = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>% 
      arrange(p_value) %>%
      # 📌 调整表格列次序：让样本量和审稿人最关心的【预测效能四驾马车】紧跟在分类信息之后
      dplyr::select(Site, TimePoint, Target_CST, 
                    N_Target_Outcome_0, N_Target_Outcome_1, 
                    N_Total_Outcome_0, N_Total_Outcome_1, 
                    Sensitivity, Specificity, PPV, NPV,
                    Odds_Ratio, CI_Lower, CI_Upper, p_value, p.adj,
                    everything())
    
    write.csv(all_logistic_report, sprintf("%s_CST_Cross_Outcome_Logistic_Filtered_%s.csv",sg_name, filter_strategy), row.names = FALSE)
  }
    
    plot_forest_df <- all_logistic_report %>% 
      mutate(
        Plot_Label    = factor(Target_CST, levels = cst_ordered_factor),
        TimePoint     = factor(TimePoint, levels = c("BL", "T04", "T12", "T24")),
        Site          = factor(Site, levels = all_sites),
        # p_fmt         = if_else(p_value < 0.001, sprintf("%.1e", p_value), sprintf("%.3f", p_value)),
        # padj_fmt      = if_else(p.adj < 0.001, sprintf("%.1e", p.adj), sprintf("%.3f", p.adj)),
        # p_mark        = if_else(p_value < 0.05, sprintf(":%s/:%s", p_fmt, padj_fmt), "")
        p_star       = get_significance_stars(p_value),
        q_star       = get_significance_stars(p.adj),
        raw_mark = paste(p_star, q_star, sep = "/"),
        p_mark = if_else(p_value > 0.05, "", raw_mark),
      )
    
    p_forest <- ggplot(plot_forest_df, aes(x = Odds_Ratio, y = Plot_Label)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "#2B8CBE", linewidth = 0.7) +
      geom_point(color = "#E41A1C", size = 2.5) +
      geom_text(aes(label = p_mark), color = "black", size = 2.2, fontface = "italic", vjust = -0.8, hjust = -0.05) +
      facet_grid(Site ~ TimePoint, scales = "free_y", space = "free", axes = "all") +
      scale_x_log10() + theme_bw() +
      labs(x = "Odds Ratio (Log Scale, 95% CI)", y = "CST Types", title = sprintf("CST and Treatment Response Forest Matrix [%s]", sg_name)) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5), strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 11, face = "bold"), axis.text = element_text(color = "black", size = 9), panel.grid.minor = element_blank()
      )
    ggsave(sprintf("%s_CST_Cross_Logistic_Forest_Filtered_%s.pdf",sg_name, filter_strategy), plot = p_forest, width = 11, height = 8)
    
         
  }                                                          
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                       