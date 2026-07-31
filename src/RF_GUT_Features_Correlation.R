# ==============================================================================
# Pipeline: Baseline Gut Microbiome vs Hormones & KI Clinical Scores
# Modules: Wilcoxon Group Analysis, Spearman Heatmaps & Segmented Scatter Plots
# Year: 2026
# ==============================================================================

library(tidyverse)
library(ggpubr)
library(rstatix)
library(psych)
library(pheatmap)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("margin", "ggplot2")

# ==============================================================================
# 1. 全局数据加载与预洗
# ==============================================================================
cat("\n[INIT] 加载全局基础数据...\n")
base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))
phen.Cate.BL     <- phen.Cate %>% subset(Phen_Time == "BL")

# 加载特征和机器学习 Master 盘
base_in_dir     <- paste0(base_dir,"/Results/RF_Model_Blood_3Model/mKI/CLR/Wilcoxon/Gut/Target_Response_K/")
gut_pred_feature <- read.csv(paste0(base_in_dir, "4_RF_Importance_Features_Microbe_Clinical.csv"))
master_data      <- read.csv(paste0(base_in_dir, "0_Master_Data_Response_K.csv"))

target_dir <- paste0(base_dir, "/Results/RF_GUT_Feature/")
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)
# 提取基线目标菌群丰度
target_features      <- gut_pred_feature$Feature[grep("GUT_",gut_pred_feature$Feature)]
matched_cols         <- intersect(c("Response", target_features), colnames(master_data))
merged_data_baseline <- master_data %>% select(Clinic_ID, all_of(matched_cols))
microbe_cols         <- colnames(merged_data_baseline)[3:ncol(merged_data_baseline)]

get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"   
  )}
# ==============================================================================
# 2. 所有基线特征组间差异分析 (Wilcoxon) 与箱线图 ##################################
# ==============================================================================
cat("\n[RUN] 正在计算菌株基线 Response 组间差异并刻画富集方向...\n")
wilcox_results_list <- list()
boxplot_list        <- list()

for (mb in c(gut_pred_feature$Feature)) {
  df_sub <- master_data %>%
    dplyr::select(Response, all_of(mb)) %>%
    drop_na(Response, all_of(mb)) %>%
    mutate(Response = factor(Response, levels = c("NR", "R")))
  
  if (length(unique(df_sub$Response)) < 2) next
  
  # ----------------------------------------------------------------------------
  # 核心新增步骤：提取两组描述性均值，硬核刻画群体层面的生物学富集方向
  # ----------------------------------------------------------------------------
  mean_stats <- df_sub %>%
    group_by(Response) %>%
    summarise(mean_val = mean(.data[[mb]], na.rm = TRUE), .groups = "drop")
  
  mean_nr <- mean_stats$mean_val[mean_stats$Response == "NR"]
  mean_r  <- mean_stats$mean_val[mean_stats$Response == "R"]
  
  # 根据均值判定富集朝向
  enrich_direction <- case_when(
    mean_r > mean_nr ~ "Enriched in R",
    mean_r < mean_nr ~ "Enriched in NR",
    TRUE             ~ "No Difference"
  )
  
  # 防御机制：安全包裹变量名以兼容含特殊字符的列名
  safe_mb_name <- sprintf("`%s`", mb)
  formula_str  ~ Response
  
  stat_res <- df_sub %>%
    wilcox_test(as.formula(paste(safe_mb_name, "~ Response"))) %>%
    mutate(
      Feature = mb,
      Mean_NR = mean_nr,
      Mean_R  = mean_r,
      Enrichment_Direction = enrich_direction 
    )
  
  wilcox_results_list[[mb]] <- stat_res
  
  # 仅为 P < 0.2 趋势或显著的特征建立散点箱线图
  if (stat_res$p < 0.2) {
    clean_title <- gsub("^GUT_", "", mb)
    
    p <- ggplot(df_sub, aes(x = Response, y = .data[[mb]], fill = Response)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
      geom_jitter(aes(color = Response), width = 0.15, alpha = 0.6, size = 1.2) +
      theme_bw() +
      scale_fill_manual(values = c("NR" = "#C7A6FA", "R" = "#E36C5B")) +
      scale_color_manual(values = c("NR" = "#C7A6FA", "R" = "#E36C5B")) +
      labs(title = clean_title, x = "Response Group", y = "CLR Abundance") +
      stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.35, size = 3.5) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 9), axis.text = element_text(color = "black", size = 8),
        panel.grid.minor = element_blank(), legend.position = "none"
      )
    boxplot_list[[mb]] <- p
  }
}

# ----------------------------------------------------------------------------
# 3. 结果合并、多重校正与大归档盘出库
# ----------------------------------------------------------------------------
all_stats <- bind_rows(wilcox_results_list) %>% 
  mutate(Adjusted_P_value = p.adjust(p, method = "BH")) %>% 
  dplyr::select(Feature, group1, group2, Mean_NR, Mean_R, Enrichment_Direction, 
                statistic, P_value = p, Adjusted_P_value, everything()) %>%
  mutate(
    p_star           = get_significance_stars(P_value),
    q_star           = get_significance_stars(Adjusted_P_value),
    Significance     = paste(p_star, q_star, sep = "/") 
  )%>%
  select(-p_star,-q_star)

write.csv(all_stats, "Baseline_Response_Wilcoxon.csv", row.names = FALSE)

# if (length(boxplot_list) > 0) {
#   combined_microbe_box <- ggarrange(plotlist = boxplot_list, ncol = 4, nrow = 3, align = "hv")
#   ggexport(combined_microbe_box, filename = "Combined_Baseline_Response_Boxplots.pdf", width = 15, height = 11)
# }


# RF 特征量重新作图 - 双向 ########################
df_gut_plot <- all_stats %>%
  left_join(gut_pred_feature%>%select(Feature,MeanDecreaseAccuracy), by = "Feature") %>%
  arrange(-MeanDecreaseAccuracy) %>%
  slice_head(n = 30) %>%
  mutate(
    Clean_Feature = gsub("^GUT_", "", Feature),
    
    Plot_Accuracy = if_else(Enrichment_Direction == "Enriched in NR", 
                            -MeanDecreaseAccuracy, 
                            MeanDecreaseAccuracy)
  )%>%
  mutate(
    Sign.col = case_when(Adjusted_P_value<0.05~ "Adjusted_P < 0.05",
                           P_value       <0.05 ~"P_value < 0.05",
                           TRUE ~ "ns")
  ) 
df_gut_plot$Significance[df_gut_plot$Significance == "ns/ns"] <- ""
df_gut_plot$Significance[df_gut_plot$Significance == "ns/ns"] <- ""


write.csv(df_gut_plot%>%select(-Plot_Accuracy), "Baseline_Response_Wilcoxon_Top30withImp.csv", row.names = FALSE)

df_gut_plot$Clean_Feature <- factor(df_gut_plot$Clean_Feature, 
                                    levels = df_gut_plot$Clean_Feature[order(abs(df_gut_plot$Plot_Accuracy))])

p_dual_lollipop <- ggplot(df_gut_plot, aes(x = Plot_Accuracy, y = Clean_Feature, color = Sign.col)) +

  geom_vline(xintercept = 0, color = "darkgray", linetype = "dashed", size = 0.5) +
  geom_segment(aes(x = 0, xend = Plot_Accuracy, y = Clean_Feature, yend = Clean_Feature), 
               size = 0.8, alpha = 0.6) +

  geom_point(aes(size = MeanDecreaseAccuracy), alpha = 0.9) +
  scale_color_manual(values = c("Adjusted_P < 0.05" = "#1F78B4", "P_value < 0.05" = "#A6CEE3"  ,"ns" = "grey" )) +

  scale_x_continuous(
    labels = function(x) abs(x),
    limits = c(-max(df_gut_plot$MeanDecreaseAccuracy)*1.1, max(df_gut_plot$MeanDecreaseAccuracy)*1.1)
  ) +

    scale_size_continuous(range = c(3, 6)) +

  theme_pubclean() +
  labs(
    title = "Random Forest Feature Importance & Enrichment Direction",
    x = "Importance (Mean Decrease Accuracy)",
    y = "Top 30 Predictive Features (Clinical & Microbiome)",
    color = "Enrichment"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, face = "italic", color = "gray30", hjust = 0.5),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"), # 纵轴特征名加粗
    axis.text.x = element_text(size = 9, color = "black"),
    # axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25,  color = "gray80"),
    legend.position = "bottom",             
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  guides(size = "none")

print(p_dual_lollipop)

ggsave("Dual_Direction_Plot.pdf", plot = p_dual_lollipop, width = 6, height = 7)



# ==============================================================================
# 3. 基线肠道菌 vs 激素, K score 关联分析 ####################################
# ==============================================================================
# 核心通用函数一：相关性计算、FDR校正 =========
run_microbe_vs_clinical_pipeline <- function(merged_df, microbe_cols, clinical_cols, 
                                             output_csv = "Stats_Table.csv", 
                                             output_pdf = "Heatmap.pdf", 
                                             heatmap_title = "Correlation Heatmap") {
  
  # 1. 提取外部已经标准化好的临床矩阵与对应的微生物矩阵
  clinical_z_mat  <- merged_df %>% dplyr::select(all_of(clinical_cols))
  microbiome_mat  <- merged_df %>% dplyr::select(all_of(microbe_cols))
  
  # 2. 核心相关性计算 (Unadjusted P values)
  cor_res <- psych::corr.test(x = microbiome_mat, y = clinical_z_mat, method = "spearman", adjust = "none")
  
  # 3. 展平高维矩阵为 SCI 二维标准表
  cor_r_df <- as.data.frame(as.table(cor_res$r)) %>% rename(Microbiome = Var1, Clinical_Factor = Var2, R_value = Freq)
  cor_p_df <- as.data.frame(as.table(cor_res$p)) %>% rename(Microbiome = Var1, Clinical_Factor = Var2, P_value = Freq)
  
  # 4. 融汇并执行全局独立 BH (FDR) 假阳性控制惩罚
  final_stats_table <- cor_r_df %>%
    inner_join(cor_p_df, by = c("Microbiome", "Clinical_Factor")) %>%
    mutate(
      Adjusted_P_value = p.adjust(P_value, method = "BH"),
      p_star           = get_significance_stars(P_value),
      q_star           = get_significance_stars(Adjusted_P_value),
      Significance     = paste(p_star, q_star, sep = "/") 
    ) %>%
    arrange(P_value)
  
  # 导出包含完整多重校正指标的标准 CSV 附表（供审稿人查阅）
  write.csv(final_stats_table %>% dplyr::select(-p_star, -q_star), output_csv, row.names = FALSE)
  
  # 5. 逆向还原高保真符号矩阵
  sig_matrix_full <- final_stats_table %>%
    dplyr::select(Microbiome, Clinical_Factor, p_star) %>%
    pivot_wider(names_from = Clinical_Factor, values_from = p_star) %>%
    column_to_rownames("Microbiome") %>%
    as.matrix()
  
  # 确保顺序与原始计算矩阵一致
  sig_matrix_full <- sig_matrix_full[rownames(cor_res$r), colnames(cor_res$r), drop = FALSE]
  
  # 6. 【核心新增核心拦截逻辑】: 动态识别并提取“至少有一个 P < 0.05”的微生物和临床指标
  # 在原始 P 值矩阵中，找出至少有一个值小于 0.05 的行和列
  p_mat <- cor_res$p
  
  keep_rows <- apply(p_mat, 1, function(row) any(row < 0.05))
  keep_cols <- apply(p_mat, 2, function(col) any(col < 0.05))
  
  # 提取合格的特征名字
  valid_microbes <- names(keep_rows)[keep_rows]
  valid_clinicals <- names(keep_cols)[keep_cols]
  
  if (length(valid_microbes) == 0 || length(valid_clinicals) == 0) {
    warning("no P < 0.05, skip")
    return(final_stats_table)
  }
  
  # 7. 矩阵瘦身：切片过滤相关系数矩阵和符号星标矩阵
  r_matrix_filtered   <- cor_res$r[valid_microbes, valid_clinicals, drop = FALSE]
  sig_matrix_filtered <- sig_matrix_full[valid_microbes, valid_clinicals, drop = FALSE]
  
  # 视觉升级：彻底拿掉碍眼的 "ns/ns" 或 "ns"，让不显著的单元格保持高雅留白
  sig_matrix_filtered[sig_matrix_filtered == "ns"] <- ""
  sig_matrix_filtered[is.na(sig_matrix_filtered)]  <- ""
  
  # 8. 绘图颜色与高级细网格线定义
  my_color <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)
  
  # 9. 渲染过滤后的清爽学术热图
  pheatmap(r_matrix_filtered, 
           display_numbers = sig_matrix_filtered, 
           fontsize_number = 12,                
           color = my_color,
           cluster_rows = TRUE, 
           cluster_cols = FALSE, 
           angle_col = 45,
           border_color = "gray92",               
           breaks = seq(-1, 1, length.out = 101), 
           main = heatmap_title,
           filename = output_pdf, 
           width = max(6, length(valid_clinicals) * 0.8), 
           height = max(4, length(valid_microbes) * 0.4))  
  
  return(final_stats_table %>% dplyr::select(-p_star, -q_star))
}

# 核心通用函数二：散点图 =========
generate_pipeline_scatter_plots <- function(stats_table, plot_data_df, 
                                            regex_bl, regex_t24, regex_delta, 
                                            pdf_prefix = "Combined_Scatter_") {
  
  sig_pairs <- stats_table %>%
    dplyr::filter(P_value < 0.05) %>%
    rename(p_value = P_value, r_value = R_value)
  
  scatter_vault <- list("BL" = list(), "T24" = list(), "Delta" = list())
  
  if (nrow(sig_pairs) > 0) {
    for (i in 1:nrow(sig_pairs)) {
      mb_name   <- as.character(sig_pairs$Microbiome[i])
      fac_name  <- as.character(sig_pairs$Clinical_Factor[i])
      pair_key  <- paste(mb_name, fac_name, sep = "_vs_")
      
      current_stage <- case_when(
        grepl(regex_bl, fac_name)    ~ "BL",
        grepl(regex_t24, fac_name)   ~ "T24",
        grepl(regex_delta, fac_name) ~ "Delta",
        TRUE                         ~ "Unknown"
      )
      
      if (current_stage == "Unknown" || length(scatter_vault[[current_stage]]) >= 24) next
      
      clean_fac_name <- gsub(paste0("^(", regex_bl, "|", regex_t24, "|", regex_delta, ")"), "", fac_name)
      clean_fac_name <- gsub(regex_bl, "", clean_fac_name) 
      clean_fac_name <- gsub(regex_t24, "", clean_fac_name)
      
      y_lab_parsed <- if(current_stage == "Delta") bquote(Delta * " " * .(clean_fac_name) * " (T24 - BL)") else clean_fac_name
      
      p_scatter <- ggplot(plot_data_df, aes(x = .data[[mb_name]], y = .data[[fac_name]])) +
        geom_point(color = "#377EB8", alpha = 0.7, size = 1.5) +
        geom_smooth(method = "lm", color = "#E41A1C", fill = "#E41A1C", alpha = 0.15) +
        theme_bw() +
        labs(title = paste(gsub("^GUT_", "", mb_name), "vs", clean_fac_name), x = "CLR Abundance", y = y_lab_parsed) +
        stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3.5) +
        theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
      
      scatter_vault[[current_stage]][[pair_key]] <- p_scatter
    }
    
    for (stage_name in names(scatter_vault)) {
      if (length(scatter_vault[[stage_name]]) > 0) {
        grid_plot <- ggarrange(plotlist = scatter_vault[[stage_name]], ncol = 4, nrow = 3, align = "hv")
        grid_plot <- grid_plot + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
        ggexport(grid_plot, filename = paste0(pdf_prefix, stage_name, ".pdf"), width = 16, height = 12)
      }
    }
  }
}


# hormone  #################################################################################
cat("\n[RUN] 正在处理激素与微生物联合分析（外部标准化版）...\n")

# 1. 原始清洗与计算 Delta
hormones_target <- Phen.Seq %>% 
  subset(Time %in% c("BL", "T24")) %>%
  select(Clinic_ID, Time, E2, FSH, LH) %>%
  drop_na(E2, FSH, LH) %>% 
  group_by(Clinic_ID) %>%
  filter(all(c("BL", "T24") %in% Time)) %>%
  ungroup() %>%
  pivot_wider(names_from = Time, values_from = c(E2, FSH, LH)) %>%
  mutate(
    Delta_E2  = E2_T24 - E2_BL,
    Delta_FSH = FSH_T24 - FSH_BL,
    Delta_LH  = LH_T24 - LH_BL
  )

# 定义需要分析的列名
hormone_cols <- c("E2_BL", "FSH_BL", "LH_BL", "E2_T24", "FSH_T24", "LH_T24", "Delta_E2", "Delta_FSH", "Delta_LH")

# 2. 【核心修改：外部无损标准化】
clinical_z_df <- hormones_target %>% 
  select(all_of(hormone_cols)) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(Clinic_ID = hormones_target$Clinic_ID)

# 3. 合并
merged_hormone <- inner_join(clinical_z_df, merged_data_baseline, by = "Clinic_ID")
merged_hormone_raw <- inner_join(hormones_target, merged_data_baseline, by = "Clinic_ID")

sub_microbe_cols <- intersect(target_features, colnames(master_data)) 

# 4. 计算与画图
stats_hormone <- run_microbe_vs_clinical_pipeline(
  merged_df = merged_hormone, 
  microbe_cols = sub_microbe_cols, 
  clinical_cols = hormone_cols,
  output_csv = "STx_Microbe_vs_Hormones_Correlation_Stats.csv",
  output_pdf = "Spearman_Global_Heatmap_Hormone.pdf",
  heatmap_title = "Spearman Correlation: Microbes vs Hormones"
)

generate_pipeline_scatter_plots(
  stats_table = stats_hormone, 
  plot_data_df = merged_hormone_raw, 
  regex_bl = "_BL$", regex_t24 = "_T24$", regex_delta = "^Delta_",
  pdf_prefix = "Combined_Spearman_Scatter_Stage_"
)

# ki_scores ######################################################################
cat("\n[RUN] 正在处理 KI 临床评分与微生物联合分析（外部标准化版）...\n")

ki_scores <- phen.Cate.BL$Phen[phen.Cate.BL$Category %in% c("K Score")]

scores_target <- Phen.Seq %>% 
  dplyr::filter(Time %in% c("BL", "T24")) %>%
  dplyr::select(Clinic_ID, Time, all_of(ki_scores)) %>%
  pivot_longer(cols = all_of(ki_scores), names_to = "Score_Type", values_to = "value") %>%
  drop_na(value) %>%
  group_by(Clinic_ID, Score_Type) %>%
  dplyr::filter(all(c("BL", "T24") %in% Time)) %>%
  ungroup() %>%
  pivot_wider(names_from = Time, values_from = value) %>%
  mutate(Delta = T24 - BL)

scores_matrix_wide  <- scores_target %>% dplyr::select(Clinic_ID, Score_Type, BL, T24) %>% pivot_wider(names_from = Score_Type, values_from = c(BL, T24), names_sep = "_")
scores_matrix_delta <- scores_target %>% dplyr::select(Clinic_ID, Score_Type, Delta) %>% pivot_wider(names_from = Score_Type, values_from = Delta, names_prefix = "Delta_")
ki_scores_final      <- inner_join(scores_matrix_wide, scores_matrix_delta, by = "Clinic_ID")

ki_cols <- ki_scores_final %>% dplyr::select(starts_with("BL_"), starts_with("T24_"), starts_with("Delta_")) %>% colnames()

# 【核心步骤：外部标准化】
ki_z_df <- ki_scores_final %>% 
  select(all_of(ki_cols)) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(Clinic_ID = ki_scores_final$Clinic_ID)

merged_ki     <- inner_join(ki_z_df, merged_data_baseline, by = "Clinic_ID")
merged_ki_raw <- inner_join(ki_scores_final, merged_data_baseline, by = "Clinic_ID")

# 4. 计算与画图
stats_ki <- run_microbe_vs_clinical_pipeline(
  merged_df = merged_ki, 
  microbe_cols = sub_microbe_cols, 
  clinical_cols = ki_cols,
  output_csv = "STx_Microbe_vs_KI_Scores_Correlation_Stats.csv",
  output_pdf = "Spearman_Global_Heatmap_KI_Scores.pdf",
  heatmap_title = "Spearman Correlation: Microbes vs KI Clinical Scores"
)

generate_pipeline_scatter_plots(
  stats_table = stats_ki, 
  plot_data_df = merged_ki_raw, # 📌 传入未标准化的原始量表分值大表，让散点图展示真实的临床改良K评分数值区间
  regex_bl = "^BL_", regex_t24 = "^T24_", regex_delta = "^Delta_",
  pdf_prefix = "Combined_Spearman_Scatter_KI_Stage_"
)