library(tidyverse)
library(reshape2)
library(scales)
library(ggbeeswarm)
library(ggpubr)

setwd("D:/WorkProjects/Demo-MHT 2026/Results/")

# ==============================================================================
# 1. 数据清洗与转移矩阵频数/比例精细计算
# ==============================================================================
# 显式重组转换长表至宽表，卡死成对存在基线与 T24 的受试者
df_matrix_wide <- prof_diversity %>%
  dplyr::filter(Site == "VA" & Time %in% c("BL", "T24")) %>%
  dplyr::select(Clinic_ID, Time, CST) %>%
  reshape2::dcast(Clinic_ID ~ Time, value.var = "CST") %>%
  na.omit()

# 统一卡死学术分型因子级别（保证行列次序绝对镜像对齐）
cst_ordered_levels <- c("UROG-L.c", "UROG-L.i", "UROG-G.v", "UROG-Div")

# 建立原始观察频数数据框
df_counts <- xtabs(~ BL + T24, data = df_matrix_wide) %>% 
  as.data.frame() %>% 
  rename(Count = Freq)

# 建立按行（基线底盘）归一化的转移概率百分比数据框
df_ratio <- xtabs(~ BL + T24, data = df_matrix_wide) %>% 
  prop.table(margin = 1) %>% 
  as.data.frame() %>% 
  rename(Ratio = Freq)

# 强力咬合频数与比例，构建完美绘图矩阵底盘
df_plot_final <- df_counts %>%
  inner_join(df_ratio, by = c("BL", "T24")) %>%
  mutate(
    # X轴正序排布
    T24 = factor(T24, levels = rev(cst_ordered_levels)),
    # 【高阶技巧】：Y轴（基线）倒序因子化，迫使高阶 CST 在热图上方聚集，符合从左上到右下的转移阅读习惯
    BL  = factor(BL,  levels = rev(cst_ordered_levels)),
    # 生成带频数和比例的双重硬核标签：例如 "72.5% (29)"
    Label_Text = if_else(Ratio > 0, paste0(sprintf("%.1f%%", Ratio * 100), "\n(n=", Count, ")"), "")
  )

# ==============================================================================
# 2. 气泡转移矩阵热图（Bubble Shift Matrix Heatmap）
# ==============================================================================
p_cst_shift <- ggplot(df_plot_final, aes(x = T24, y = BL)) +
  geom_tile(color = "grey90", fill = "grey98", linewidth = 0.5) +
  geom_point(aes(size = Ratio, color = Ratio), alpha = 0.85) +
  scale_size_continuous(range = c(2, 16), guide = "none") +
  scale_color_distiller(palette = "Blues", direction = 1, labels = scales::percent_format(), name = "Shift Ratio") +
  geom_text(aes(label = Label_Text), size = 2.8, fontface = "bold", color = "grey10", vjust = 0.5) +
  theme_bw() +
  labs(
    x = "Follow-up (T24)",
    y = "Baseline"
  ) +
  theme(
    plot.title        = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title        = element_text(face = "bold", size = 9.5),
    axis.text         = element_text(size = 9, color = "black", face = "bold"),
    axis.text.x       = element_text(angle = 30, vjust = 1, hjust = 1), # 倾斜30度防爆重叠
    legend.title      = element_text(face = "bold", size = 8.5),
    legend.text       = element_text(size = 8),
    legend.position   = "right"
  )
p_cst_shift

ggsave("CST_Microbiome_Transition_Matrix.pdf", plot = p_cst_shift, width = 6.2, height = 5.2)



# FSFI & L.c############################
my_mean <- function(x) mean(x, na.rm = TRUE)

my_custom_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text        = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.line        = element_line(linewidth = 0.5),
      legend.position   = "bottom",
      plot.title       = element_text(face = "bold", size = 15),
      plot.subtitle    = element_text(color = "grey30"),
      strip.text       = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.spacing    = unit(1.5, "lines")
    )
}

df.diversity <- prof_diversity %>%
  mutate(
    Time = Microbe.phen.prof$Time[match(SeqID, Microbe.phen.prof$SeqID)],
    Site = Microbe.phen.prof$Site[match(SeqID, Microbe.phen.prof$SeqID)],
    vt   = CST
  ) %>%
  filter(Site %in% c("UR", "VA"))

group_data_t24 <- df.diversity %>%
  filter(Time == "T24", Site == "VA") %>%
  drop_na(vt) %>%
  mutate(
    Group = ifelse(vt == "UROG-L.c", "UROG-L.c at T24", "Non-UROG-L.c at T24"),
    Group = factor(Group, levels = c("UROG-L.c at T24", "Non-UROG-L.c at T24"))
  )

# ==============================================================================
# Part 1: FSFI VS T24 UROG-L.c 纵向分析
# ==============================================================================
df_fsfi <- Phen.Seq %>%
  select(Clinic_ID, Time, FSFI) %>%
  mutate(
    Group = group_data_t24$Group[match(Clinic_ID, group_data_t24$Clinic_ID)],
    value = FSFI
  ) %>%
  drop_na(Group, value) %>%
  mutate(
    Time  = factor(Time, levels = unique(as.character(Time))),
    Group = factor(Group, levels = c("UROG-L.c at T24", "Non-UROG-L.c at T24"))
  )

# 四个时间点完整的患者
df_fsfi_complete <- df_fsfi %>%
  group_by(Clinic_ID) %>%
  filter(n() == 4) %>%
  ungroup()

# Friedman 检验与标签制作
friedman_fsfi <- df_fsfi_complete %>%
  group_by(Group) %>%
  summarise(
    p_val = friedman.test(value ~ Time | Clinic_ID)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    sig   = case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ " "),
    label = paste0(Group, ", P = ", round(p_val, 3), " ", sig)
  )

# 绘图 1
ggplot(df_fsfi, aes(x = Time, y = value, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.3, linewidth = 0.5) +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.5, size = 0.5) +
  stat_summary(
    fun = my_mean, geom = "line", aes(group = Group), 
    position = position_dodge(width = 0.5), alpha = 0.8, linewidth = 1.5
  ) +
  stat_summary(
    fun = my_mean, geom = "point", shape = 23, size = 2, fill = "white",
    position = position_dodge(0.8)
  ) +
  stat_compare_means(
    method = "wilcox.test", aes(label = ..p.signif..),
    step.increase = 0.08, label.y = max(df_fsfi$value) * 1.05,
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  ) +
  scale_fill_manual(values = c("#CC6677", "#332288")) +
  scale_color_manual(values = c("#CC6677", "#332288")) +
  labs(caption = paste0(friedman_fsfi$label, collapse = "\n"), x = "Time Point", y = "FSFI") +
  my_custom_theme()

ggsave("FSFI_L.c_group_T24.pdf", width = 6, height = 5)

# ==============================================================================
# Part 2: Lactobacillus crispatus 丰度纵向分析 (多部位 Facet)
# ==============================================================================
# 合并 VA 和 UR 丰度数据，减少重复代码
df_taxa_long <- bind_rows(
  prof_filtered[["VA"]] %>% data.frame() %>% rownames_to_column("SeqID") %>% mutate(Site = "VA"),
  prof_filtered[["UR"]] %>% data.frame() %>% rownames_to_column("SeqID") %>% mutate(Site = "UR")
) %>%
  select(SeqID, Site, Lactobacillus_crispatus, Gardnerella_vaginalis) %>%
  mutate(
    Clinic_ID = Microbe.phen.prof$Clinic_ID[match(SeqID, Microbe.phen.prof$SeqID)],
    Time      = Microbe.phen.prof$Time[match(SeqID, Microbe.phen.prof$SeqID)],
    Group     = group_data_t24$Group[match(Clinic_ID, group_data_t24$Clinic_ID)],
    value     = Lactobacillus_crispatus
  ) %>%
  filter(!is.na(Group)) %>%
  drop_na(value) %>%
  mutate(
    Time  = factor(Time, levels = c("BL", "T04", "T12", "T24")),
    Group = factor(Group, levels = c("UROG-L.c at T24", "Non-UROG-L.c at T24")),
    Site  = factor(Site, levels = unique(as.character(Site)))
  )

# 四点共存控制
df_taxa_complete <- df_taxa_long %>%
  group_by(Site, Clinic_ID) %>%
  filter(n() == 4) %>%
  ungroup()

# Friedman 检验
friedman_taxa <- df_taxa_complete %>%
  group_by(Site, Group) %>%
  summarise(
    p_val = friedman.test(value ~ Time | Clinic_ID)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    sig   = case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ " "),
    label = paste0(Group, ", P = ", round(p_val, 3), " ", sig)
  )

# 绘图 2
ggplot(df_taxa_complete, aes(x = Time, y = value, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.3, linewidth = 0.5) +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.5, size = 0.5) +
  stat_summary(
    fun = my_mean, geom = "line", aes(group = interaction(Group, Site)), 
    position = position_dodge(width = 0.5), alpha = 0.8, linewidth = 1.5
  ) +
  facet_grid(~Site) +
  stat_summary(
    fun = my_mean, geom = "point", shape = 23, size = 2, fill = "white",
    position = position_dodge(0.8)
  ) +
  stat_compare_means(
    method = "wilcox.test", aes(label = ..p.signif..),
    step.increase = 0.08, label.y = max(df_taxa_long$value) * 1.05,
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  ) +
  scale_fill_manual(values = c("#CC6677", "#332288")) +
  scale_color_manual(values = c("#CC6677", "#332288")) +
  labs(caption = paste0(friedman_taxa$label, collapse = "\n"), x = "Time Point", y = "Lactobacillus crispatus") +
  my_custom_theme()

ggsave("UR_VA_L.c_group_T24.pdf", width = 8, height = 6)

# ==============================================================================
# Part 3: Baseline Vagina Type 与 FSFI 横向关联分析
# ==============================================================================
# 获取基线 CST 亚组分类
group_data_bl <- df.diversity %>%
  filter(Time == "BL", Site == "VA") %>%
  drop_na(vt) %>%
  mutate(
    Group = factor(vt, levels = c("UROG-L.c", "UROG-L.i", "UROG-G.v", "UROG-Div"))
  )

df_baseline <- Phen.Seq %>%
  select(Clinic_ID, Time, FSFI) %>%
  mutate(
    Group = group_data_bl$Group[match(Clinic_ID, group_data_bl$Clinic_ID)],
    value = FSFI
  ) %>%
  drop_na(Group, value) %>%
  mutate(Time = factor(Time, levels = unique(as.character(Time))))

# 两两比较组合矩阵
my_com <- combn(levels(df_baseline$Group), 2, simplify = FALSE)

# 绘图 3
ggplot(df_baseline, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.3, linewidth = 0.5) +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.5, size = 0.5) +
  facet_grid(~Time) +
  stat_compare_means(
    method = "wilcox.test", aes(label = ..p.signif..),
    comparisons = my_com, step.increase = 0.05, tip.length = 0.01, size = 3,
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  ) +
  scale_fill_manual(values = c("#CC6677", "#807DBA", "#9467BD", "#3F007D")) +
  scale_color_manual(values = c("#CC6677", "#807DBA", "#9467BD", "#3F007D")) +
  labs(x = "Time Point", y = "FSFI") +
  my_custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 仅针对此图旋转 x 轴文本

ggsave("FSFI_va_Baseline_group.pdf", width = 10, height = 6)