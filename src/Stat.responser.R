

library(ggplot2)
library(tidyverse)
library(ggpubr) # 用于最后的拼图

# ==============================================================================
# 0. 自动化百分比堆叠图核心绘图函数 (复制即可用)
# ==============================================================================
plot_stacked_percentage <- function(data, y_var, fill_var, y_label, fill_legend_title, colors) {
  
  # 1. 动态重命名变量，方便统一处理
  plot_df <- data %>%
    select(Y = all_of(y_var), Fill = all_of(fill_var)) %>%
    na.omit() %>%
    mutate(Y = as.factor(Y), Fill = as.factor(Fill))
  
  total_n <- nrow(plot_df)
  
  # 2. 自动化执行卡方检验获取学术 P 值
  chisq_res <- chisq.test(table(plot_df$Y, plot_df$Fill))
  p_val <- chisq_res$p.value
  p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  test_label <- paste0("Chi-squared test, ", p_label, " (Total N = ", total_n, ")")
  
  # 3. 计算频数和百分比，用于在柱子内部写字 "72% (136)"
  anno_df <- plot_df %>%
    group_by(Y, Fill) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(Y) %>%
    mutate(
      total = sum(count),
      percentage = (count / total) * 100,
      # 计算文本标签
      label_text = paste0(sprintf("%.1f%%", percentage), " (", count, ")")
    )
  
  # 4. 开始使用 ggplot2 构图
  p <- ggplot(anno_df, aes(x = count, y = Y, fill = Fill)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) +
    geom_text(aes(label = label_text), 
              position = position_fill(vjust = 0.5), 
              color = "white", size = 3, fontface = "bold") +
    # 格式化 X 轴为百分比标尺
    scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    # 高雅的学术配色
    scale_fill_manual(values = colors, name = fill_legend_title) +
    # 纯白学术主题与细节微调
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_blank(),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      plot.margin = margin(t = 10, r = 18, b = 10, l = 10, unit = "pt")
    ) +
    labs(
      x = "Proportion",
      y = y_label,
      caption = test_label 
    )
  
  return(p)
}

# ==============================================================================
# 1. 模拟/准备你的真实数据集 (请让同学用包含这三列的真实大表替换)
# ==============================================================================

# Load datasets (Ensure the paths are correct)
base_dir   <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

setwd("D:/WorkProjects/Demo-MHT 2026/Results/Responser")
final_clinical_data <- merge(
  merge(phen.gr.E2.all %>% mutate(Hormone_E2 = recode_factor(Group,"NR" = "H_NR","R"  = "H_R")) %>% select(-Group),
        phen.gr.K      %>% mutate(mKI = recode_factor(Group,"NR" = "mKI_NR","R"  = "mKI_R"))%>% select(-Group),all = T),
        phen.gr.gut    %>% mutate(Gut_Type = recode_factor(CST, "GUT-P.cop"="GUT-P.cop","GUT-P.v"="GUT-P.v"))%>% select(-CST),all = T) 

dim(final_clinical_data)

# 定义两组经典的临床对比色
colors_set1 <- c("H_R" = "#2ecc71", "H_NR" = "#d35400")          # 绿 vs 橙（对应你图中的颜色）
colors_set2 <- c("mKI_R" = "#1a5276", "mKI_NR" = "#bdc3c7") # 蓝 vs 灰（标准临床速度对比色）
colors_set3 <- c( "GUT-P.v"   = "#66c2a4", "GUT-P.cop" = "#b2e2e2")

# ==============================================================================
# 2. 一键批量生成 3 张学术切片图
# ==============================================================================
# 图 A: 肠道基线 vs 激素应答 (还原你的 image_887339.png)
fig_A <- plot_stacked_percentage(
  data = final_clinical_data, 
  y_var = "Gut_Type", fill_var = "Hormone_E2", 
  y_label = "Baseline Gut Type", fill_legend_title = "Hormonal Response", 
  colors = colors_set1
)
fig_A

# 图 B: 肠道基线 vs 临床见效速度
fig_B <- plot_stacked_percentage(
  data = final_clinical_data, 
  y_var = "Gut_Type", fill_var = "mKI", 
  y_label = "Baseline Gut Type", fill_legend_title = "Clinical Response", 
  colors = colors_set2
)
fig_B


final_merged_figure <- ggarrange(
  fig_A, fig_B,
  ncol = 1, nrow = 2,align = "hv"
)

ggsave("MHT_Response_Group_Tri-Matrix.pdf", plot = final_merged_figure, 
       width = 6, height = 14)


# Change.dot Hormone_E2 ######################
plot_trajectory_1p <- function(data, var_prefix, group_col, title_suffix, endpoint = "T24") {
  bl_col <- paste0(var_prefix, "_BL")
  # 动态生成终点列名（例如：var_prefix_T24 或 var_prefix_T04）
  ep_col <- paste0(var_prefix, "_", endpoint)
  
  # 剔除缺失值
  plot_df <- data %>%
    select(Clinic_ID, BL = !!sym(bl_col), Endpoint = !!sym(ep_col), Group = !!sym(group_col)) %>%
    filter(!is.na(BL) & !is.na(Endpoint) & !is.na(Group))
  
  # 【新增】自动计算各分组的人数，生成 caption 文本
  group_counts <- plot_df %>% count(Group)
  caption_text <- paste("Sample Size:", paste(group_counts$Group, "=", group_counts$n, collapse = " | "))
  
  # 对每个分组按基线(BL)降序排列，生成个体序号，呈现瀑布流水般的视觉效果
  plot_df <- plot_df %>%
    group_by(Group) %>%
    arrange(desc(BL)) %>%
    mutate(Patient_Idx = row_number()) %>%
    ungroup()
  
  # 定义临床变量是否"越高越好"以决定参考线的样式或解释
  is_positive_indicator <- var_prefix == "FSFI"
  
  p <- ggplot(plot_df) + 
    geom_segment(aes(x = Patient_Idx, xend = Patient_Idx, y = BL, yend = Endpoint), 
                 color = "grey70", alpha = 0.5, linewidth = 0.2,
                 arrow = arrow(length = unit(0.15, "cm"), ends = "last", type = "closed")) +
    geom_point(aes(x = Patient_Idx, y = BL), color = "#4477AA", size = 1.2, alpha = 0.8) +
    geom_point(aes(x = Patient_Idx, y = Endpoint), color = "#FF7F0E", size = 1.2, alpha = 0.8) +
    facet_grid(~ Group, scales = "free_x", space = "free") +
    labs(
      title = paste(var_prefix, "-", title_suffix),
      x = "Individual Patients (Ordered by Baseline)", 
      y = paste(var_prefix, ""),
      caption = paste(caption_text, "(Blue dot: Baseline | Orange dot:", endpoint,")")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      plot.caption = element_text(hjust = 1, face = "italic", color = "black", size = 11),  
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "white", color = "white")
    )
  return(p)
}


final_clinical_data <- merge(
  merge(phen.gr.E2.all %>% mutate(Hormone_E2 = recode_factor(Group,"NR" = "H_NR","R"  = "H_R")) %>% select(-Group),
        phen.gr.K      %>% mutate(mKI = recode_factor(Group,"NR" = "mKI_NR","R"  = "mKI_R"))%>% select(-Group),all = T),
  phen.gr.gut    %>% mutate(Gut_Type = recode_factor(CST, "GUT-P.cop"="GUT-P.cop","GUT-P.v"="GUT-P.v"))%>% select(-CST),all = T)

phen_merged <- phen.enroll %>% select(Clinic_ID,K_Score_BL,K_Score_T04,E2_BL,E2_T24,FSH_BL,FSH_T24)%>%
  left_join(final_clinical_data)
  

# plot_trajectory_1p(phen_merged, "K_Score", "mKI", "Clinical Response",endpoint = "T04")
# ggsave("Change.dots.mKI.pdf",width = 12,height = 5)

plot_trajectory_1p(phen_merged, "E2", "Hormone_E2", "Hormonal Response")
ggsave("Change.dots.Hormone_E2.pdf",width = 12,height = 5)


#  Change.dot mKI with Hormone_E2 and mKI #######################################################
library(tidyverse)
 
var_prefix <- "K_Score"
group_col  <- "mKI"
title_suffix <- "Clinical Response"
endpoint   <- "T04"
 
final_clinical_data <- merge(
  merge(
    phen.gr.E2.all %>% 
      mutate(Hormone_E2 = recode_factor(Group, "NR" = "H_NR", "R" = "H_R")) %>% 
      select(-Group),
    phen.gr.K %>% 
      mutate(mKI = recode_factor(Group, "NR" = "mKI_NR", "R" = "mKI_R")) %>% 
      select(-Group), 
    all = TRUE
  ),
  phen.gr.gut %>% 
    mutate(Gut_Type = recode_factor(CST, "GUT-P.cop" = "GUT-P.cop", "GUT-P.v" = "GUT-P.v")) %>% 
    select(-CST), 
  all = TRUE
)

phen_merged <- phen.enroll %>% 
  select(Clinic_ID, K_Score_BL, K_Score_T04, E2_BL, E2_T24, FSH_BL, FSH_T24) %>%
  left_join(final_clinical_data, by = "Clinic_ID")
 
bl_col <- paste0(var_prefix, "_BL")
ep_col <- paste0(var_prefix, "_", endpoint)

plot_df <- phen_merged %>%
  select(
    Clinic_ID, 
    BL         = !!sym(bl_col), 
    Endpoint   = !!sym(ep_col), 
    Group      = !!sym(group_col), 
    Hormone_E2 = Hormone_E2,
    Gut_Type   = Gut_Type
  ) %>%
  filter(!is.na(BL) & !is.na(Endpoint) & !is.na(Group)) %>% 
  mutate(
    Hormone_E2_str = ifelse(is.na(Hormone_E2), "NA", as.character(Hormone_E2)),
    Gut_Type_str   = ifelse(is.na(Gut_Type), "NA", as.character(Gut_Type))
  )

# 自动计算各分组人数
group_counts <- plot_df %>% count(Group)
caption_text <- paste("Sample Size:", paste(group_counts$Group, "=", group_counts$n, collapse = " | "))

# 组内按基线降序排列生成 Patient_Idx
plot_df <- plot_df %>%
  group_by(Group) %>%
  arrange(desc(BL)) %>%
  mutate(Patient_Idx = row_number()) %>%
  ungroup()

 
#  定义颜色映射表（加入 "NA" 对应灰色）
e2_colors  <- c("H_R" = "#2ecc71", "H_NR" = "#d35400")
gut_colors <- c("GUT-P.v" = "#66c2a4", "GUT-P.cop" = "#b2e2e2")
all_fills  <- c(e2_colors, gut_colors, "NA" = "grey30")

p <- ggplot(plot_df) + 
  # 1. 轨迹线段 (映射 Hormone_E2 颜色，NA 时默认为灰色)
  geom_segment(
    aes(
      x = Patient_Idx, xend = Patient_Idx, 
      y = BL, yend = Endpoint
    ), 
    color = "grey",
    alpha = 0.6, 
    linewidth = 0.3,
    arrow = arrow(length = unit(0.12, "cm"), ends = "last", type = "closed")
  ) +
  geom_point(aes(x = Patient_Idx, y = BL), color = "#4477AA", size = 1.2, alpha = 0.8) +
  geom_point(aes(x = Patient_Idx, y = Endpoint), color = "#FF7F0E", size = 1.2, alpha = 0.8) +
  
  geom_tile(
    data = plot_df,
    aes(x = Patient_Idx, y = -1, fill = Hormone_E2_str),
    height = 0.8,       # 控制色块高度（可根据 y 轴刻度微调）
    color = "white",    # 色块边框颜色，设为 white 会有网格隔断感，若要连续填色可设为 NA
    linewidth = 0.1     # 边框线宽
  ) +
  geom_tile(
    data = plot_df,
    aes(x = Patient_Idx, y = -2, fill = Gut_Type_str),
    height = 0.8,       # 控制色块高度
    color = "white",    # 色块边框颜色
    linewidth = 0.1
  ) +
  
  facet_grid(~ Group, scales = "free_x", space = "free") +
  
  # 线段颜色控制
  scale_color_manual(
    name = "Hormone E2 (Trajectory)", 
    values = e2_colors, 
    na.value = "grey70"
  ) +
  
  # 🌟 填充颜色映射控制
  scale_fill_manual(
    name = "Annotation Tracks", 
    values = all_fills,
    breaks = c("H_R", "H_NR", "GUT-P.v", "GUT-P.cop", "NA"),
    labels = c("H_R", "H_NR", "GUT-P.v", "GUT-P.cop", "Missing (NA)")
  ) +
  labs(
    title = paste(var_prefix, "-", title_suffix),
    x = "Individual Patients (Ordered by Baseline)", 
    y = var_prefix,
    caption = paste(caption_text, "(Light dot: Baseline | Dark dot:", endpoint, ")")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(hjust = 0.5, color = "grey70"),
    plot.caption       = element_text(hjust = 1, face = "italic", color = "black", size = 11),  
    axis.text.x        = element_blank(),  
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold"),
    strip.background   = element_rect(fill = "white", color = "white"),
    legend.position    = "bottom",
    legend.box         = "vertical"
  )

print(p)
ggsave("Change.dots.mKI.new.pdf",width = 15,height = 6)

