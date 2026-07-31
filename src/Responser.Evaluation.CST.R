# ==============================================================================
# Pipeline: Menopause Hormone Therapy (MHT) Multi-Omic Subgroup Evaluation
# Modules: Stratified CST Longitudinal Ratio & Logistic Regression Forest Plots
# Year: 2026
# ==============================================================================

library(tidyverse)
library(ggpubr)
library(rstatix)
library(conflicted)
library(broom)
library(openxlsx)
library(vegan)
library(colorspace)
library(RColorBrewer)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

# Define optimized shared meta-subgroup mapping configuration matrices
subgroups_config <- list(
  "Hormone_E2" = phen.gr.E2.all %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "mKI"        = phen.gr.K      %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "Gut_Type"   = phen.gr.gut    %>% select(Clinic_ID, group = CST)   %>% mutate(group = factor(group, levels = c("GUT-P.v", "GUT-P.cop")))
)

filter_strategy    <- "None"

all_sites          <- c("VA", "UR", "GUT", "TO")
follow_ups         <- c("T04", "T12", "T24")
cst_ordered_factor <- c("GUT-P.cop", "GUT-P.v", "TO-N.s", "TO-P.m", "UROG-Div", "UROG-G.v", "UROG-L.i", "UROG-L.c")

# Base metadata harmonization
mg_meta_outcome_base <- Microbe.phen.prof %>%
  left_join(prof_diversity %>% select(SeqID, CST), by = "SeqID")

# ==============================================================================
# CST -  SUBGROUP ITERATION WORKFLOW  ###################################
# ==============================================================================
for (sg_name in names(subgroups_config)) {
  
  # Dynamic sub-directory mapping configuration
  target_dir_cst <- file.path(base_dir, "Results/Responser", sg_name, "CST")
  if(!dir.exists(target_dir_cst)) dir.create(target_dir_cst, recursive = TRUE)
  setwd(target_dir_cst)
  
  current_subgroup_df <- subgroups_config[[sg_name]]
  
  # Inject stratified subgroup mapping factor back into metadata frame
  mg_meta_outcome <- mg_meta_outcome_base %>%
    select(-any_of("group")) %>%
    left_join(current_subgroup_df, by = "Clinic_ID") %>%
    filter(!is.na(group))
  
  df_cst_base <- mg_meta_outcome %>%
    dplyr::select(Clinic_ID, Time, Site, CST, group) %>%
    drop_na(Time, Site, CST, group)
  
  # Initialize processing storage containers per subgroup
  longitudinal_stats     <- list()
  logistic_results_list  <- list()
  valid_rows_list        <- list()
  
  # ============================================================================
  # 1. Automated Stratified Longitudinal Statistics (Fisher & Logistic)
  # ============================================================================
  for (st in all_sites) {
    for (fu in follow_ups) {
      
      comp_df_all <- df_cst_base %>%
        dplyr::filter(Site == st, Time %in% c("BL", fu))
      
      # Perform calculations independently inside each subgroup factor level
      for (g_level in levels(df_cst_base$group)) {
        comp_df <- comp_df_all %>% filter(group == g_level)
        if (nrow(comp_df) == 0) next
        
        # 1.1 Pairwise longitudinal tracking alignment
        paired_ids <- comp_df %>%
          group_by(Clinic_ID) %>%
          summarise(has_bl = "BL" %in% Time, has_fu = fu %in% Time, .groups = "drop") %>%
          dplyr::filter(has_bl & has_fu) %>%
          pull(Clinic_ID)
        
        if (filter_strategy == "strict_4") {
          strict_4_ids <- df_cst_base %>%
            dplyr::filter(Site == st, group == g_level) %>%
            group_by(Clinic_ID) %>%
            summarise(n_time = n_distinct(Time), .groups = "drop") %>%
            dplyr::filter(n_time == 4) %>% 
            pull(Clinic_ID)
          paired_ids <- intersect(paired_ids, strict_4_ids)
        }
        
        current_n <- length(paired_ids)
        if (current_n < 3) next 
        
        pair_df <- comp_df %>% dplyr::filter(Clinic_ID %in% paired_ids)
        valid_rows_list[[paste(st, fu, g_level, sep = "_")]] <- pair_df
        
        # 1.2 Contingency Table & Fisher's Exact Test
        tab <- table(pair_df$CST, pair_df$Time[, drop = TRUE])
        tab <- tab[, c("BL", fu), drop = FALSE] 
        tab <- tab[rowSums(tab) > 0, , drop = FALSE]
        
        if (nrow(tab) >= 2) {
          test_res <- fisher.test(tab, simulate.p.value = TRUE, B = 2000)
          longitudinal_stats[[paste(st, fu, g_level, sep = "_")]] <- data.frame(
            Site = st, Subgroup = g_level, Comparison = paste0(fu, "_vs_BL"), 
            p = test_res$p.value, N_Samples = current_n, stringsAsFactors = FALSE
          )
        }

        # # 1.3 Conditional Logistic Regression Engine (clogit Paired Design)
        # ==============================================================================
        library(survival)
        
        # 构建条件回归基盘：时间转换为干预二分类（BL=0, 随访点=1）
        slice_df    <- pair_df %>% mutate(Time_Binary = if_else(Time == fu, 1, 0))
        unique_csts <- unique(slice_df$CST)
        if (length(unique_csts) < 2) next
        
        for (target_cst in unique_csts) {
          run_df <- slice_df %>% mutate(CST_Binary = if_else(CST == target_cst, 1, 0))
          contingency_table <- table(run_df$CST_Binary, run_df$Time_Binary)
          if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) next
          tryCatch({
            model <- clogit(CST_Binary ~ Time_Binary + strata(Clinic_ID), data = run_df)
            model_stats <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% 
              dplyr::filter(term == "Time_Binary")
            if (nrow(model_stats) == 0 || is.na(model_stats$estimate)) next
            logistic_results_list[[paste(st, fu, g_level, target_cst, sep = "_")]] <- data.frame(
              Site          = st, 
              Subgroup      = g_level, 
              FollowUpPoint = fu, 
              Target_CST    = target_cst, 
              Odds_Ratio    = model_stats$estimate,   
              CI_Lower      = model_stats$conf.low, 
              CI_Upper      = model_stats$conf.high, 
              p_value       = model_stats$p.value, 
              N_Patients    = current_n, 
              stringsAsFactors = FALSE
            )
            
          }, error = function(e) {
            cat(sprintf("      [INFO] CST [%s] 在 [%s_vs_BL] 条件下无足够变异，已跳过。\n", target_cst, fu))
          })
        }
        
        
      }
    }
  }
  
  if (length(longitudinal_stats) == 0) {
    message(sprintf(">>> Skip Subgroup [%s]: Insufficient sample sizes resolved.", sg_name))
    next
  }
  
  # ============================================================================
  # 2. Composition Visualization (Stacked Proportions & Forest Plots)
  # ============================================================================
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
    )
  }
  # 2.1 Format Stratified Fisher Outputs
  cst_stats_long <- bind_rows(longitudinal_stats) %>%
    group_by(Site) %>%
    mutate(p.adj = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    mutate(
      p_fmt    = if_else(p < 0.0001, sprintf("%.1e", p), sprintf("%.3f", p)),
      padj_fmt = if_else(p.adj < 0.0001, sprintf("%.1e", p.adj), sprintf("%.3f", p.adj)),
      p_label  = sprintf("%s/%s\nn=%d", p_fmt, padj_fmt, N_Samples),
      y_pos    = 1.03, 
      Time     = factor(gsub("_vs_BL", "", Comparison), levels = c("BL", "T04", "T12", "T24")),
      Subgroup = factor(Subgroup, levels = levels(df_cst_base$group)),
      group    = Subgroup
    )
  write.csv(cst_stats_long, sprintf("CST_Longitudinal_Proportion_%s_%s.csv", sg_name, filter_strategy), row.names = FALSE)
  
  df_cst <- bind_rows(valid_rows_list) %>%
    distinct() %>%
    mutate(
      Time = factor(Time, levels = c("BL", "T04", "T12", "T24")),
      Site = factor(Site, levels = all_sites),
      CST  = factor(CST, levels = cst_ordered_factor),
      group = factor(group, levels = levels(df_cst_base$group))
    )
  
  # 2.2 Stratified Stacked Proportion Plot (Facet Grid: Subgroup x Site)
  df_plot_cst <- df_cst %>%
    group_by(Site, group, Time, CST) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(Percentage = count / sum(count)) %>%
    ungroup()
  
  p_cst_long_ratio <- ggplot(df_plot_cst, aes(x = Time, y = Percentage, fill = CST)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
    facet_grid( Site ~group, scales = "free_y", axes = "all") + 
    scale_fill_manual(values = cst_colors, name = "CST Type", drop = TRUE) +
    scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Relative Proportion of CST", 
         title = sprintf("Longitudinal Shifts of CST Proportion [%s]", sg_name)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 9, color = "black"), 
      axis.ticks.y = element_line(color = "black"),
      panel.grid = element_blank()
    ) +
    geom_text(
      data        = cst_stats_long, 
      aes(x = Time, y = y_pos, label = p_label), 
      inherit.aes = FALSE, vjust = 0.5, size = 2.2, color = "black", lineheight = 0.85
    )
  ggsave(sprintf("CST_Longitudinal_Proportions_%s_%s.pdf", sg_name, filter_strategy), plot = p_cst_long_ratio, width = 3 + length(levels(df_cst_base$group))*2 , height =12 )
  
  # 2.3 Format & Export Stratified Logistic Report
  if (length(logistic_results_list) > 0) {
    all_logistic_report <- bind_rows(logistic_results_list) %>%
      group_by(Site) %>%
      mutate(p.adj = p.adjust(p_value, method = "BH"))%>%
      ungroup() %>%
      mutate(
        p_star = get_significance_stars(p_value),
        q_star = get_significance_stars(p.adj),
        raw_mark = paste(p_star, q_star, sep = "/"),
        p_mark = if_else(p_value > 0.05, "", raw_mark)
      ) %>%
      dplyr::select(-p_star, -q_star,-raw_mark) %>%
      arrange(p_value)
    write.csv(all_logistic_report, sprintf("CST_Longitudinal_Logistic_%s_%s.csv", sg_name, filter_strategy), row.names = FALSE)
    
    # 2.4 Stratified Forest Plot Matrix Layout (Site-Internal BH Correction)
    plot_forest_df <- all_logistic_report %>% 
      mutate(
        Plot_Label    = factor(Target_CST, levels = cst_ordered_factor),
        FollowUpPoint = factor(FollowUpPoint, levels = follow_ups),
        Site          = factor(Site, levels = all_sites),
        Subgroup      = factor(Subgroup, levels = levels(df_cst_base$group))
      )%>%
      filter(!is.na(CI_Upper) & !is.na(CI_Lower),
              !is.infinite(CI_Upper) & !is.infinite(CI_Upper))
    
    # 森林图核心图层渲染
    p_forest <- ggplot(plot_forest_df, aes(x = Odds_Ratio, y = Plot_Label, color = Target_CST)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, linewidth = 0.7) +
      geom_point(size = 2.5) +
      
      # 动态文本映射：将复合校正指标以精细斜体平铺在几率比置信区间（Confidence Interval）右上方
      geom_text(aes(label = p_mark), size = 3, fontface = "italic", vjust = -0.6, hjust = -0.05, show.legend = FALSE) +
      
      # 矩阵分面布局 (行: 部位 + 亚组, 列: 随访节点)
      facet_grid(Site + Subgroup ~ FollowUpPoint, scales = "free_y", space = "free_y", axes = "all") +
      
      scale_x_log10() + 
      scale_color_manual(values = cst_colors, name = "CST Type", drop = TRUE) +
      labs(
        x     = "Odds Ratio vs Baseline (Log Scale, 95% CI)", 
        y     = "CST Types", 
        title = sprintf("Longitudinal Logistic Odds Ratio of CST Types [%s]", sg_name)
      ) +
      theme_bw() +
      theme(
        plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text       = element_text(size = 11, face = "bold"),
        axis.text.x      = element_text(color = "black", size = 9),
        axis.text.y      = element_text(color = "black", size = 9),
        axis.ticks       = element_line(color = "black", linewidth = 0.4), 
        panel.grid.minor = element_blank(), 
        legend.position  = "bottom"
      )
    
    # 高清图像文件输出
    ggsave(sprintf("CST_Longitudinal_Logistic_Forest_%s_%s.pdf", sg_name, filter_strategy), 
           plot = p_forest, width = 10, height = 4 + length(levels(df_cst_base$group)) * 3.5)
    }
}



# ==============================================================================
# Diversity -  SUBGROUP ITERATION WORKFLOW  ###############################
# ==============================================================================

cols <- c("BL/T4" = "#D9534F", "BL/T12" = "#F0AD4E", "BL/T24" = "#337AB7", "T4/T12" = "#5CB85C", "T4/T24" = "#9467BD", "T12/T24" = "#8C564B")
time_order_bray <- c("BL", "T4", "T12", "T24") 

for (sg_name in names(subgroups_config)) {
  
  # 1. 动态子目录切换
  target_dir_cst <- file.path(base_dir, "Results/Responser", sg_name, "CST")
  if(!dir.exists(target_dir_cst)) dir.create(target_dir_cst, recursive = TRUE)
  setwd(target_dir_cst)
  
  # 2. 抽取并洗净当前亚组对应的临床响应数据
  current_subgroup_df <- subgroups_config[[sg_name]]
  subgroup_statuses   <- levels(factor(current_subgroup_df$group))
  
  mg_meta_outcome <- mg_meta_outcome_base %>%
    select(-any_of("group")) %>%
    left_join(current_subgroup_df, by = "Clinic_ID") %>%
    filter(!is.na(group))
  
  # ----------------------------------------------------------------------------
  # MODULE 5. PERMANOVA Sub-Stratified Pairwise Execution (Longitudinal strict_4_ids) ####################
  # ----------------------------------------------------------------------------
  permanova_stat_list <- list()
  
  for (site_name in names(prof_filtered)) {
    
    mat_raw <- data.frame(prof_filtered[[site_name]])
    
    # 动态咬合当前部位的微生物与亚组表型
    mat_with_meta_site <- mat_raw %>%
      mutate(
        Group     = mg_meta_outcome$Time[match(rownames(mat_raw), mg_meta_outcome$SeqID)],
        Clinic_ID = mg_meta_outcome$Clinic_ID[match(rownames(mat_raw), mg_meta_outcome$SeqID)],
        Resp_Grp  = mg_meta_outcome$group[match(rownames(mat_raw), mg_meta_outcome$SeqID)]
      ) %>%
      mutate(Group = recode(Group, "T04" = "T4")) %>%
      filter(Group %in% time_order_bray & !is.na(Resp_Grp))
    
    if (nrow(mat_with_meta_site) == 0) next

    if (filter_strategy == "strict_4") {
      strict_4_ids <- mat_with_meta_site %>%
        group_by(Clinic_ID) %>%
        summarise(n_time = n_distinct(Group), .groups = "drop") %>%
        filter(n_time == 4) %>% 
        pull(Clinic_ID)
      
      mat_with_meta_site <- mat_with_meta_site %>% filter(Clinic_ID %in% strict_4_ids)
    }
    
    pure_abundance_all <- mat_raw[rownames(mat_with_meta_site), , drop = FALSE] # 随之更新丰度表矩阵切片
    site_pairs_list    <- list()
    
    # 拆分响应状态层级进行体内纵向两两比对
    for (r_status in subgroup_statuses) {
      mat_with_meta <- mat_with_meta_site %>% filter(Resp_Grp == r_status)
      if (nrow(mat_with_meta) == 0) next
      
      for (pair_name in names(cols)) {
        tps        <- strsplit(pair_name, "/")[[1]]
        tp1        <- tps[1]
        tp2        <- tps[2]
        comp_label <- paste0(tp1, "_vs_", tp2)
        
        # 严格交集控线：必须在当前响应状态下，两个时间点都具备采样的患者
        ids_in_tp1 <- mat_with_meta %>% filter(Group == tp1) %>% pull(Clinic_ID) %>% unique()
        ids_in_tp2 <- mat_with_meta %>% filter(Group == tp2) %>% pull(Clinic_ID) %>% unique()
        paired_clinic_ids <- intersect(ids_in_tp1, ids_in_tp2)
        if (length(paired_clinic_ids) < 3) next
        
        sub_meta <- mat_with_meta %>% 
          filter(Group %in% c(tp1, tp2) & Clinic_ID %in% paired_clinic_ids)
        
        n_g1 <- sum(sub_meta$Group == tp1)
        n_g2 <- sum(sub_meta$Group == tp2)
        
        sub_abundance <- pure_abundance_all[rownames(sub_meta), , drop = FALSE]
        sub_abundance <- sub_abundance[, colSums(sub_abundance) > 0, drop = FALSE] 
        
        # 距离矩阵分盘及属性测度
        sub_bray    <- vegan::vegdist(sub_abundance, method = "bray")
        dist_vector <- as.vector(sub_bray)
        med_val     <- round(median(dist_vector), 3)
        mean_val    <- round(mean(dist_vector), 3)
        sd_val      <- round(sd(dist_vector), 3)
        mean_sd_str <- paste0(mean_val, "_", sd_val)
        
        tryCatch({
          ad_res  <- vegan::adonis2(sub_bray ~ Group, data = sub_meta, permutations = 999)
          p_val   <- ad_res$`Pr(>F)`[1]
          r2_val  <- ad_res$R2[1]
          f_model <- ad_res$F[1]
          
          uniq_p_key <- paste(r_status, pair_name, sep = "___")
          site_pairs_list[[uniq_p_key]] <- data.frame(
            Sites            = site_name,
            Subgroup_Status  = r_status,
            Comparisons      = paste(tp1, "vs", tp2, sep = "_"),
            time_pair        = pair_name,
            N_Group1         = n_g1, 
            N_Group2         = n_g2,
            Median           = med_val,
            `Mean±SD`        = mean_sd_str,
            R2               = r2_val,
            F.model          = f_model,
            p_value          = p_val,
            stringsAsFactors = FALSE
          )
        }, error = function(e) {})
      }
    }
    
    if (length(site_pairs_list) == 0) next
    
    # 局部多重校正：基于解剖部位在各自响应组内独立实施 FDR 校正
    site_permanova_df <- bind_rows(site_pairs_list) %>%
      group_by(Subgroup_Status) %>%
      mutate(Adjusted_p = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>%
      mutate(
        p_star   = get_significance_stars(p_value),
        q_star   = get_significance_stars(Adjusted_p),
        Significance = paste(p_star, q_star, sep = "/"),
        Method   = "PERMANOVA"
      ) %>%
      dplyr::select(-p_star, -q_star)
    
    permanova_stat_list[[site_name]] <- site_permanova_df
  }
  
  if (length(permanova_stat_list) > 0) {
    final_permanova_all <- bind_rows(permanova_stat_list) %>%
      mutate(
        Sites     = factor(Sites, levels = names(prof_filtered)),
        time_pair = factor(time_pair, levels = names(cols))
      ) %>%
      arrange(Sites, Subgroup_Status, time_pair)
    
    write.csv(final_permanova_all, sprintf("Pairwise_PERMANOVA_Stats_Subgroup_%s_Four_Sites.csv", sg_name), row.names = FALSE)
  }
  
  # ----------------------------------------------------------------------------
  # MODULE 6. Alpha Diversity Longitudinal Paired Wilcoxon Test (Supp Table 13) ##################
  # ----------------------------------------------------------------------------
  alpha_wilcox_list <- list()
  
  prof_diversity_sub <- prof_diversity %>%
    inner_join(current_subgroup_df, by = "Clinic_ID") %>%
    mutate(Time = recode(Time, "T04" = "T4")) %>%
    filter(Time %in% time_order_bray) %>%
    drop_na(Shannon, Site, Time, group)

  for (site_name in unique(prof_diversity_sub$Site)) {
    site_alpha_data <- prof_diversity_sub %>% filter(Site == site_name)
    
    
    if (filter_strategy == "strict_4") {
      strict_4_alpha_ids <- site_alpha_data %>%
        group_by(Clinic_ID) %>%
        summarise(n_time = n_distinct(Time), .groups = "drop") %>%
        filter(n_time == 4) %>% 
        pull(Clinic_ID)
      
      site_alpha_data <- site_alpha_data %>% filter(Clinic_ID %in% strict_4_alpha_ids)
    }
    
    for (r_status in subgroup_statuses) {
      status_alpha_data <- site_alpha_data %>% filter(group == r_status)
      if (nrow(status_alpha_data) == 0) next
      
      for (pair_name in names(cols)) {
        tps <- strsplit(pair_name, "/")[[1]]
        tp1 <- tps[1]
        tp2 <- tps[2]
        
        # 严格交集控线：必须在当前状态下，两个时间点都具备 Shannon 值的同一批受试者
        ids_tp1    <- status_alpha_data %>% filter(Time == tp1) %>% pull(Clinic_ID)
        ids_tp2    <- status_alpha_data %>% filter(Time == tp2) %>% pull(Clinic_ID)
        paired_ids <- intersect(ids_tp1, ids_tp2)
        if (length(paired_ids) < 3) next
        
        calc_df <- status_alpha_data %>% 
          filter(Time %in% c(tp1, tp2) & Clinic_ID %in% paired_ids) %>%
          mutate(Time = factor(Time, levels = c(tp1, tp2))) %>%
          arrange(Clinic_ID, Time) 
        
        shannon_tp1 <- calc_df %>% filter(Time == tp1) %>% pull(Shannon)
        shannon_tp2 <- calc_df %>% filter(Time == tp2) %>% pull(Shannon)
        
        mean1_val <- round(mean(shannon_tp1), 6)
        mean2_val <- round(mean(shannon_tp2), 6)
        est_val   <- round(median(shannon_tp1 - shannon_tp2), 6)
        
        tryCatch({
          w_test <- calc_df %>% wilcox_test(Shannon ~ Time, paired = TRUE, detailed = FALSE)
          p_val  <- w_test$p
          
          uniq_w_key <- paste(site_name, r_status, pair_name, sep = "___")
          alpha_wilcox_list[[uniq_w_key]] <- data.frame(
            Groups           = sg_name,         
            Subgroups        = r_status,        
            Sites            = site_name,       
            Comparisons      = paste0(tp1, "/", tp2),
            SampleNum1       = length(shannon_tp1),
            SampleNum2       = length(shannon_tp2),
            Mean1            = mean1_val,
            Mean2            = mean2_val,
            Estimate         = est_val,
            p_value          = p_val,
            stringsAsFactors = FALSE
          )
        }, error = function(e) {})
      }
    }
  }
  
  if (length(alpha_wilcox_list) > 0) {
    final_alpha_wilcox_report <- bind_rows(alpha_wilcox_list) %>%
      group_by(Sites, Subgroups) %>%
      mutate(Adjusted_P_value = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>%
      mutate(
        p_star       = get_significance_stars(p_value),
        q_star       = get_significance_stars(Adjusted_P_value),
        Significance = paste(p_star, q_star, sep = "/")
      ) %>%
      dplyr::select(-p_star, -q_star) %>%
      arrange(factor(Sites, levels = c("VA", "UR", "TO", "GUT")), Subgroups, factor(Comparisons, levels = names(cols)))
    
    write.csv(final_alpha_wilcox_report, sprintf("Supplementary_Table_13_Subgroup_%s_Alpha_Wilcox_Paired.csv", sg_name), row.names = FALSE)
  }
}

