library(tidyverse)
library(ggpubr)
library(rstatix)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ==============================================================================
# 1. 环境初始化与全局配置加载
# ==============================================================================
base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

# 创建并锁定主输出根目录
main_output_dir <- file.path(base_dir, "Results/Responser")
if(!dir.exists(main_output_dir)) dir.create(main_output_dir, recursive = TRUE)
setwd(main_output_dir)

# 动态提取 continuous 且不属于 K Score 的变量 (只具备 BL 和 T24)
non_k_vars <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & !phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()
# 动态提取 continuous 且属于 K Score 的变量 (具备 BL, T04, T12, T24)
k_vars     <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()

# 统合全量目标变量底盘
all_target_features <- unique(c(non_k_vars, k_vars))

subgroups_config <- list(
  "Hormone_E2" = phen.gr.E2.all %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "mKI"        = phen.gr.K      %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "Gut_Type"   = phen.gr.gut    %>% select(Clinic_ID, group = CST)   %>% mutate(group = factor(group, levels = c("GUT-P.v", "GUT-P.cop")))
)

# 提取元数据核心数据列
phen_core_data <- Phen.Seq %>% 
  dplyr::select(Clinic_ID, Time, all_of(all_target_features)) %>%
  mutate(Time = as.character(Time))

# 全量随访窗口状态字典配置
time_matrix_config <- list(
  "T04" = c("BL", "T04"),
  "T12" = c("BL", "T12"),
  "T24" = c("BL", "T24")
)

# ==============================================================================
# 2. 核心分析引擎：跨亚组、跨随访时段的 Absolute Delta 自动化循环
# ==============================================================================

all_stats_list <- list()

for (sg_name in names(subgroups_config)) {
  
  current_subgroup_df <- subgroups_config[[sg_name]]
  target_levels        <- levels(current_subgroup_df$group)
  box_palette_fill    <- setNames(c("#C7A6FA", "#E36C5B", "#9ecae1", "#a1d99b")[1:length(target_levels)], target_levels)
  
  # 融合当前亚组表型标签
  df_sg_merged <- phen_core_data %>%
    inner_join(current_subgroup_df, by = "Clinic_ID") %>%
    filter(!is.na(group))
  
  # 遍历不同的随访终点窗口 (T04, T12, T24)
  for (time_label in names(time_matrix_config)) {
    
    target_times <- time_matrix_config[[time_label]]
    t_base   <- target_times[1]
    t_follow <- target_times[2]
    
    # 动态确定在当前时间视窗下，哪些变量具备有效性
    current_vars_pool <- k_vars
    if (time_label == "T24") {
      current_vars_pool <- unique(c(k_vars, non_k_vars))
    }
    
    # 裁剪当前时段下的原始数据长表
    df_long_delta <- df_sg_merged %>%
      dplyr::filter(Time %in% target_times) %>%
      dplyr::select(Clinic_ID, Time, group, all_of(current_vars_pool)) %>%
      pivot_longer(cols = all_of(current_vars_pool), names_to = "Variable", values_to = "Value") %>%
      drop_na(Value, group) %>%
      # 严格受试者卡控：必须基线和当前随访点数据全齐
      group_by(Clinic_ID, Variable) %>%
      filter(all(target_times %in% Time)) %>%
      ungroup() %>%
      # 展开为宽表计算绝对差值 Delta
      pivot_wider(names_from = Time, values_from = Value) %>%
      mutate(
        Delta    = !!sym(t_follow) - !!sym(t_base),
        Variable = factor(Variable, levels = current_vars_pool),
        group    = factor(group, levels = target_levels)
      )
    
    if (nrow(df_long_delta) == 0) next
    
    # 📌 新增 2.0：计算每组 Delta 的均值 (Mean)
    # 1. 动态获取当前的两个分组名称，以便精确命名列名
    active_groups <- levels(df_long_delta$group)
    g1_name <- active_groups[1]
    g2_name <- active_groups[2]
    
    # 2. 统计每个指标在两组中的均值并转化为宽表
    mean_stats <- df_long_delta %>%
      group_by(Variable, group) %>%
      summarise(Mean_Val = mean(Delta, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = Mean_Val, names_prefix = "Mean_")
    
    # 2.1 核心统计：针对变化量 Delta 进行组间 Wilcoxon 非配对检验
    stat_delta_res <- df_long_delta %>%
      group_by(Variable) %>%
      filter(length(unique(group)) == 2 & n() >= 4) %>%
      wilcox_test(Delta ~ group, paired = FALSE) %>%
      ungroup()
    
    if (nrow(stat_delta_res) == 0) next
    
    stat_delta_res <- stat_delta_res %>%
      left_join(mean_stats, by = "Variable") %>%
      mutate(
        p.adj        = p.adjust(p, method = "BH"),
        p_star       = get_significance_stars(p),
        q_star       = get_significance_stars(p.adj),
        Significance = paste(p_star, q_star, sep = "/")
      ) %>%
      dplyr::select(Variable, n1, n2, starts_with("Mean_"), p, p.adj, Significance, everything())
    
    # 2.3 严格图盘卡控：提取当前随访终点下，Delta 差异显著 (P < 0.05) 的指标
    sig_delta_vars <- stat_delta_res %>% 
      filter(p < 0.05) %>% 
      pull(Variable) %>% 
      unique()

    if (length(sig_delta_vars) > 0) {
      
      df_plot   <- df_long_delta %>% filter(Variable %in% sig_delta_vars)
      stat_plot <- stat_delta_res %>% filter(Variable %in% sig_delta_vars)
      
      y_positions <- df_plot %>%
        group_by(Variable) %>%
        summarise(y_pos = max(Delta, na.rm = TRUE) * 1, .groups = "drop")

      stat_plot <- stat_plot %>%
        left_join(y_positions, by = "Variable") %>%
        rstatix::add_xy_position(formula = Delta ~ group, data = df_plot, x = "group", scales = "free_y", step.increase = 0) %>%
        mutate(y.position = y_pos)
      
      if (exists("clinic.var.order")) {
        df_plot   <- df_plot   %>% mutate(Variable = factor(Variable, levels = clinic.var.order))
        stat_plot <- stat_plot %>% mutate(Variable = factor(Variable, levels = clinic.var.order))
      }

      p_delta <- ggplot(df_plot, aes(x = group, y = Delta)) +
        geom_boxplot(aes(fill = group), alpha = 0.6, outlier.shape = NA, width = 0.45) +
        geom_jitter(aes(color = group), width = 0.18, alpha = 0.5, size = 1.2) +
        facet_wrap(~ Variable, scales = "free_y", nrow = ifelse(length(sig_delta_vars) > 8, 2, 1)) +
        theme_bw() +
        scale_fill_manual(values = box_palette_fill) +
        scale_color_manual(values = box_palette_fill) +
        labs(
          x = "Stratified Grouping Vector", 
          y = sprintf("Delta Deviation Shift Value (\\Delta %s - %s)", t_follow, t_base),
          title = sprintf("Global Scores Variance Delta Profiles [%s - %s]", sg_name, time_label)
        ) +
        theme(
          plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text       = element_text(size = 9.5, face = "bold"),
          legend.position  = "right"
        ) +
        stat_pvalue_manual(
          data = stat_plot, label = "Significance", xmin = "group1", xmax = "group2",
          y.position = "y.position", tip.length = 0.015, bracket.size = 0.45, size = 3.3, inherit.aes = FALSE
        )

      calc_width <- max(4.5, min(14, length(sig_delta_vars) * 2.6))

      target_sub_dir <- file.path(main_output_dir, sg_name, "Clinic")
      if(!dir.exists(target_sub_dir)) dir.create(target_sub_dir, recursive = TRUE)
      
      ggsave(file.path(target_sub_dir, sprintf("%s_All_Features_Delta_%s.pdf", sg_name, time_label)), 
             plot = p_delta, width = calc_width, height = 5.0)
    }
    
    target_sub_dir <- file.path(main_output_dir, sg_name, "Clinic")
    if(!dir.exists(target_sub_dir)) dir.create(target_sub_dir, recursive = TRUE)
    
    
    stat_delta_res_tagged <- stat_delta_res %>%
      mutate(
        Subgroup  = sg_name,     
        Timepoint = time_label,  
        .before = 1              
      )

    all_stats_list[[paste(sg_name, time_label, sep = "_")]] <- stat_delta_res_tagged
    write.csv(stat_delta_res, file.path(target_sub_dir, sprintf("%s_All_Features_Delta_Stats_%s.csv", sg_name, time_label)), row.names = FALSE)
  }
}

## combination
master_stats_table <- bind_rows(all_stats_list)

write.csv(master_stats_table, 
          file.path(main_output_dir, "Responser_Clinic_All_Features_Delta_Stats.csv"), 
          row.names = FALSE)
