library(tidyverse)
library(rstatix)
library(ggpubr)



# ==============================================================================
# 1. 初始化变量域与参数配置定义
# ==============================================================================
base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))
setwd(file.path(base_dir, "Results/Responser"))

non_k_vars <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & !phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()
k_vars     <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()

subgroups_config <- list(
  "Hormone_E2" = phen.gr.E2.all %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "mKI"        = phen.gr.K      %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR")))
)

all_target_features <- unique(c(non_k_vars, k_vars))

phen_core_data <- Phen.Seq %>% 
  drop_na(Time, K_Score) %>%
  dplyr::select(Clinic_ID, Time, all_of(all_target_features)) %>%
  mutate(Time = as.character(Time))

# filter_strategy <- "core_BL_24" # 可选值: "strict_4", "core_BL_24", "all"
filter_strategy <- "all" # 可选值: "strict_4", "core_BL_24", "all"

get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"   
  )}

# ==============================================================================
# 2. 核心分析引擎：执行统计与描述性均值计算（按时间点独立校正）
# ==============================================================================
cat("\n[RUN] 正在开始跨时间点独立校正的横截面非配对统计规整计算...\n")

global_stats_vault <- list()
global_data_vault  <- list() 

patient_whitelist <- phen_core_data %>%
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

for (sg_name in names(subgroups_config)) {
  current_subgroup_df <- subgroups_config[[sg_name]]
  target_levels        <- levels(current_subgroup_df$group)
  
  df_sg_merged <- phen_core_data %>%
    inner_join(patient_whitelist, by = "Clinic_ID") %>%  
    inner_join(current_subgroup_df, by = "Clinic_ID") %>%
    filter(!is.na(group))
  
  df_sg_long <- df_sg_merged %>%
    pivot_longer(cols = all_of(all_target_features), names_to = "Variable", values_to = "Value") %>%
    drop_na(Value, group)
  
  df_sg_cleaned <- df_sg_long %>%
    filter(
      (Variable %in% non_k_vars & Time %in% c("BL", "T24")) |
        (Variable %in% k_vars     & Time %in% c("BL", "T04", "T12", "T24"))
    )
  
  if (nrow(df_sg_cleaned) == 0) next
  
  global_data_vault[[sg_name]] <- df_sg_cleaned %>% mutate(Subgroup_Type = sg_name)
  
  # =========================================================================
  # 2.1 基础非配对 Wilcoxon 检验计算（自带 n1, n2 列）
  # =========================================================================
  sg_stat_res <- df_sg_cleaned %>%
    group_by(Variable, Time) %>%
    filter(length(unique(group)) == 2 & n() >= 4) %>%
    wilcox_test(Value ~ group, paired = FALSE) %>%
    ungroup()
  
  # =========================================================================
  # 2.2 描述性均值提取（顺手计算真实的独立观测样本量）
  # =========================================================================
  mean_res <- df_sg_cleaned %>%
    group_by(Variable, Time, group) %>%
    # 📌 不仅算均值，还把去除缺失值后的样本量计算出来
    summarise(
      Mean_Val = mean(Value, na.rm = TRUE),
      N_Count  = sum(!is.na(Value)),
      .groups  = "drop"
    ) %>%
    pivot_wider(names_from = group, values_from = c(Mean_Val, N_Count), names_sep = "_")
  
  # =========================================================================
  # 2.3 按时间点(Time)独立校正并拼装样本量列
  # =========================================================================
  sg_stat_final <- sg_stat_res %>%
    left_join(mean_res, by = c("Variable", "Time")) %>%
    mutate(
      # 1) 动态映射两组的均值
      Mean1         = !!sym(paste0("Mean_Val_", target_levels[1])),
      Mean2         = !!sym(paste0("Mean_Val_", target_levels[2])),
      
      # 2) 📌 核心新增：提取并整合总样本量、组1样本量、组2样本量
      N_Group1      = n1, # 继承自 wilcox_test 的组1观测数
      N_Group2      = n2, # 继承自 wilcox_test 的组2观测数
      Total_N       = n1 + n2, # 自动咬合两组总样本量
      
      Subgroup_Type = sg_name
    ) %>%
    group_by(Time) %>%  
    mutate(
      p.adj        = p.adjust(p, method = "BH"),
      p_star       = get_significance_stars(p),
      q_star       = get_significance_stars(p.adj),
      raw_mark     = paste(p_star, q_star, sep = "/"),
      Significance = if_else(p > 0.05, "", raw_mark)
    ) %>%
    ungroup() %>%
    # 3) 📌 高级规范：重组列序，把审稿人最关心的【样本量三驾马车】和【均值】推向视觉中心
    dplyr::select(Variable, Time, Subgroup_Type, 
                  Total_N, N_Group1, N_Group2, 
                  Mean1, Mean2, 
                  statistic, p, p.adj, Significance, everything())
  
  global_stats_vault[[sg_name]] <- sg_stat_final
}

# 整合大表并导出根目录备份
final_global_table <- bind_rows(global_stats_vault)
final_global_data  <- bind_rows(global_data_vault)

if (nrow(final_global_table) > 0) {
  final_global_table <- final_global_table %>%
    dplyr::select(
      Subgroup_Type, Variable, Time, 
      Group_1 = group1, Group_2 = group2, 
      N_Group1, N_Group2, 
      Mean1, Mean2,
      Statistic = statistic, 
      P_value = p, 
      FDR_p.adj = p.adj, 
      Significance
    ) %>%
    arrange(Subgroup_Type, factor(Time, levels = c("BL", "T04", "T12", "T24")), Variable)
  
  # 在主运行空间留下全局总归档
  write.csv(final_global_table, sprintf("Responser_Cross_Clinic_Wilcox_Stats_Global_%s.csv", filter_strategy), row.names = FALSE)
  cat("\n[SUCCESS] 全局统一大归档数据盘导出完成。\n")
}

# ==============================================================================
# 3. 外置独立的基线显著性筛选作图函数（完美融合亚组独占目录路由机制）
# ==============================================================================
#' Plot Baseline Significance Boxplots inside Subgroup Directories
#' @param stats_df 全局统计结果大表 (data.frame 或路径)
#' @param raw_data_df 清洗后的原始表达/得分丰度大表
#' @param root_path 传导主系统的工作根路径（即你的 base_dir）
#' @param p_cutoff 筛选作图的原始 P 值阈值，默认 0.05
plot_baseline_sig_scores <- function(stats_df, raw_data_df, root_path, p_cutoff = 0.05) {
  cat("\n[PLOT] 激活外置作图路由引擎，开始执行亚组隔离分流导出...\n")
  
  if (is.character(stats_df)) stats_df <- read.csv(stats_df, check.names = FALSE)
  
  # 1. 锁死基线(BL)且满足显著性阈值的记录切片
  bl_sig_meta <- stats_df %>%
    filter(Time == "BL" & P_value < p_cutoff)
  
  # 2. 动态扫描大表中所有被激活的亚组分类群
  unique_subgroups <- unique(stats_df$Subgroup_Type)
  
  for (sg in unique_subgroups) {
    
    # --------------------------------------------------------------------------
    # 【动态目录管理核心】：强行切入当前亚组专属的物理磁盘空间
    # --------------------------------------------------------------------------
    target_dir_hormone <- file.path(root_path, "Results/Responser", sg, "Cross_Clinic")
    if(!dir.exists(target_dir_hormone)) dir.create(target_dir_hormone, recursive = TRUE)
    setwd(target_dir_hormone)
    
    # 3. 剥离属于当前亚组的单维度全时段 Wilcox 统计文件并就地存盘归档
    sg_full_stats <- stats_df %>% filter(Subgroup_Type == sg)
    write.csv(sg_full_stats, sprintf("Cross_Sectional_Wilcox_Stats_%s.csv", filter_strategy), row.names = FALSE)
    cat(sprintf("[FILE] 亚组 [%s] 的全随访差异表格已输出至专属 Cross_Clinic 目录。\n", sg))
    
    # 4. 提取当前亚组下在基线表现出真显著的变量群
    sg_sig_vars <- bl_sig_meta %>% filter(Subgroup_Type == sg) %>% pull(Variable) %>% unique()
    
    # 作图防线拦截：若当前亚组在基线无任何指标显著，直接静默放行去往下一个亚组，防止空图抛错
    if (length(sg_sig_vars) == 0) {
      cat(sprintf("[INFO] 亚组 [%s] 基线无显著特征指标，跳过图形绘制步骤。\n", sg))
      next
    }
    
    # 5. 裁剪属于当前亚组空间的表达子集
    df_plot <- raw_data_df %>% 
      filter(Subgroup_Type == sg & Time == "BL" & Variable %in% sg_sig_vars)
    
    stat_plot <- bl_sig_meta %>% 
      filter(Subgroup_Type == sg & Variable %in% sg_sig_vars) %>%
      rename(group1 = Group_1, group2 = Group_2, p = P_value)
    
    df_plot <- df_plot %>% mutate(group = factor(group))
    
    if (exists("clinic.var.order")) {
      df_plot   <- df_plot   %>% mutate(Variable = factor(Variable, levels = clinic.var.order))
      stat_plot <- stat_plot %>% mutate(Variable = factor(Variable, levels = clinic.var.order))
    }

    target_levels <- levels(df_plot$group)
    
    # 6. 动态推算各分面专属的 Y 轴弹性上限位置，防止括号跟星标跟数据打架
    y_positions <- df_plot %>%
      group_by(Variable) %>%
      summarise(y_pos = max(Value, na.rm = TRUE) * 0.98, .groups = "drop")
    
    # 公式记忆补齐（显式注入 Value ~ group 挽救丢失的 formula 属性属性）
    stat_plot <- stat_plot %>%
      left_join(y_positions, by = "Variable") %>%
      rstatix::add_xy_position(formula = Value ~ group, data = df_plot, 
                               x = "group", scales = "free_y", step.increase = 0) %>%
      mutate(y.position = y_pos)
    
    # 学术配色方案绑定
    box_palette_fill   <- setNames(c("#C7A6FA", "#E36C5B", "#9ecae1", "#a1d99b")[1:length(target_levels)], target_levels)
    box_palette_border <- c(box_palette_fill, "Between_Group" = "grey30")
    
    # 7. 构建 ggplot 艺术版面
    p_box <- ggplot(df_plot, aes(x = group, y = Value)) +
      geom_boxplot(aes(fill = group), alpha = 0.6, outlier.shape = NA, width = 0.45) +
      geom_jitter(aes(color = group), width = 0.18, alpha = 0.5, size = 1.2) +
      facet_wrap(~ Variable, scales = "free_y", axes = "all", nrow = ifelse(length(sg_sig_vars) > 6, 2, 1)) +
      theme_bw() +
      scale_fill_manual(values = box_palette_fill) +
      scale_color_manual(values = box_palette_border) +
      labs(
        x = "Clinical Stratification Domain", 
        y = "Baseline Value",
        title = sprintf("Baseline Cross-Sectional Shifts [%s]", sg)
      ) +
      theme(
        plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text       = element_text(size = 9.5, face = "bold"),
        legend.position  = "right"
      ) +
      stat_pvalue_manual(
        data        = stat_plot,
        label       = "Significance",
        xmin        = "group1",
        xmax        = "group2",
        y.position  = "y.position",
        tip.length  = 0.015,
        bracket.size = 0.45,
        size        = 2,
        color       = "black",
        inherit.aes = FALSE
      )
    
    calc_width <- max(4.2, min(14, length(sg_sig_vars) * 2.5))
    file_name  = sprintf("Baseline_Only_Scores_Boxplot_%s.pdf", filter_strategy)
    
    ggsave(file_name, plot = p_box, width = calc_width, height = 5.0)
    cat(sprintf("[PLOT] 亚组 [%s] 的显著性 PDF 图盘渲染完毕并导出成功。\n", filter_strategy))
  }
}

# ==============================================================================
# 4. 运行外置作图引擎示例
# ==============================================================================
# 将主路径 base_dir 传导进去，它会自动解析并在各个专属文件夹中就地吐出对应的分析结果
plot_baseline_sig_scores(stats_df = final_global_table, raw_data_df = final_global_data, root_path = base_dir, p_cutoff = 0.05)
