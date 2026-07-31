
# ==============================================================================
# Pipeline: Vaginal (VA) Cohort & External Reference Micro-Ecology Analysis
# Year: 2026
# ==============================================================================

library(dplyr)
library(tibble)
library(tidyr)
library(rstatix)
library(ggplot2)
library(stringr)
library(vegan)
library(ape)
library(purrr)
library(openxlsx)
library(ggpubr)
library(reshape2)
library(forcats)
library(RColorBrewer)
library(patchwork)

# ------------------------------------------------------------------------------
# 0. Global Parameters & Directory Configuration
# ------------------------------------------------------------------------------
base_dir <- "D:/WorkProjects/Demo-MHT 2026"
out_dir  <- "D:/WorkProjects/Demo-MHT 2026/Results/OverallEva/VA_REF_Microbe"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

# 加载全局对象映像
load(file.path(base_dir, "data/MHT.demo.RData"))

prefix           <- "VA_Ref_Only"
rel_ab.threshold <- 0.0001

# ==============================================================================
# 1. Dataset I Preparation: VA Longitudinal Cohort
# ==============================================================================
message(">>> 正在提取并过滤 VA 纵向随访队列数据集...")

data_va <- prof_filtered[["VA"]] %>% data.frame()

# Taxa Filter: 流行度卡控（保留在 >10% 的样本中超过阈值的物种）
taxa_keep_va <- sapply(data_va, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(data_va))
taxa.vars_va <- names(taxa_keep_va)[taxa_keep_va]
data_va      <- data_va[, taxa.vars_va, drop = FALSE]

df_va <- data_va %>%
  mutate(
    Group       = Microbe.phen.prof$Time[match(rownames(data_va), Microbe.phen.prof$SeqID)],
    Clinic_ID   = Microbe.phen.prof$Clinic_ID[match(rownames(data_va), Microbe.phen.prof$SeqID)],
    Site.source = "VA" 
  )

# 严格配对：筛选出具有 4 个随访时间点完整记录的受试者
Clinic_IDs_4_va <- df_va %>% 
  group_by(Clinic_ID, Group) %>% 
  group_by(Clinic_ID) %>% 
  summarise(pairs = n(), .groups = "drop") %>% 
  filter(pairs == 4) %>% 
  pull(Clinic_ID)

df_va <- df_va %>% subset(Clinic_ID %in% Clinic_IDs_4_va)

# ==============================================================================
# 2. Dataset II Preparation: External Reference Cohort
# ==============================================================================
message(">>> 正在提取并过滤外部 Reference 对照组数据集...")

phen.ref <- phen.va.ref
prof.ref <- prof.va.ref %>% column_to_rownames("X")

# 共有的样本对齐
intersect_samples <- intersect(phen.ref$SeqID, rownames(prof.ref))
prof.ref          <- prof.ref[intersect_samples, , drop = FALSE]
prof.ref          <- prof.ref[rowSums(prof.ref) != 0, colSums(prof.ref) != 0, drop = FALSE]

data_ref      <- as.data.frame(prof.ref)
taxa_keep_ref <- sapply(data_ref, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(data_ref))
taxa.vars_ref <- names(taxa_keep_ref)[taxa_keep_ref]
data_ref      <- data_ref[, taxa.vars_ref, drop = FALSE]

df_ref <- data_ref %>%
  mutate(
    Group       = phen.ref$Group[match(rownames(data_ref), phen.ref$SeqID)],
    Clinic_ID   = phen.ref$Cohort[match(rownames(data_ref), phen.ref$SeqID)],
    Site.source = "Ref" 
  )

# ==============================================================================
# 3. 核心双表扁平化合并 (取物种并集 ➕ 缺失值 NA 自动补 0)
# ==============================================================================
message(">>> 正在执行双表联合（取并集物种并补0）...")

# # 纵向合并数据，缺失的格子由 bind_rows 自动生成 NA
cols_to_keep <- names(df_va)
cat("df_ref delete:", setdiff(names(df_ref), cols_to_keep), "\n")

df.merged <- dplyr::bind_rows(
  df_va,
  df_ref[intersect(names(df_ref), cols_to_keep)]
)

# 统一时间/对照组标签（T04 -> T4, Menopause -> H_M, Reproductive -> H_R）
df.merged <- df.merged %>%
  mutate(Group = recode(Group, 
                        "T04" = "T4", "T12" = "T12", "T24" = "T24",
                        "Menopause" = "H_M", "Reproductive" = "H_R"))

# 将双表合并产生的因并集空缺导致的 NA 批量安全填充为 0
df_numeric_cols <- sapply(df.merged, is.numeric)
df.merged[, df_numeric_cols][is.na(df.merged[, df_numeric_cols])] <- 0

# 建立全局因子顺序
df.merged$Group <- factor(df.merged$Group, levels = c("BL", "T4", "T12", "T24", "H_R", "H_M"))
top_factors     <- levels(df.merged$Group)

# 重新校正映射分类：Ref 在生物学本质上归属于阴道样本 (VA)
df.merged$Site <- gsub("Ref", "VA", df.merged$Site.source) %>% factor(levels = c("VA"))
df.merged$CST  <- prof_diversity$CST[match(rownames(df.merged), prof_diversity$SeqID)] %>%
  factor(levels = rev(c("UROG-Div", "UROG-G.v", "UROG-L.i", "UROG-L.c")))

print("合并完成后的各观测组样本分布大盘（已移除UR）：")
print(xtabs(~ Site.source + Group, df.merged))


# ==============================================================================
# 6. Pipeline II: Wilcoxon Differential Analysis & Evaluation Heatmap
# ==============================================================================
message(">>> Starting Wilcoxon statistical comparisons for VA and Reference...")

# 6.0 丰度矩阵 Log10 转换对齐 ---------------------------------------------------
# 提取纯微生物群落相对丰度列
abundance_cols <- setdiff(colnames(df.merged), c("Group", "Clinic_ID", "Site.source", "Site", "CST"))

df_meta_labels <- df.merged %>% select(Group, Clinic_ID)
df_numeric_log <- df.merged[, abundance_cols]

# 加伪计数进行稳定对数转化
pseu_value <- min(df_numeric_log[df_numeric_log != 0]) * 0.1
df_numeric_log <- log10(df_numeric_log + pseu_value)

# 重新组装用于统计的平铺数据框
df_stat_ready <- cbind(df_meta_labels, df_numeric_log)
top_factors   <- levels(df_stat_ready$Group) # 沿用全局因子顺序

# 6.1 定义标准的对比对子 (完美对齐新版归口名称：T4, H_M, H_R) -------------------
my_comparisons <- list(
  c("BL", "T4"),  c("BL", "T12"), c("BL", "T24"), 
  c("BL", "H_M"), c("BL", "H_R"), 
  c("T24", "H_M"), c("T24", "H_R")
)

stats_results     <- list()
plot_data_results <- list()

for (i in 1:length(my_comparisons)) {
  t1 <- my_comparisons[[i]][1]
  t2 <- my_comparisons[[i]][2]
  
  # 索引 1:3 为纵向随访配对设计，4:7 涉及与外部对照组的非配对设计
  if (i %in% 1:3) {
    # 严格锁死配对样本
    pairs.samp <- df_stat_ready %>%
      filter(Group %in% my_comparisons[[i]]) %>% 
      group_by(Clinic_ID) %>%
      summarise(pairs = n(), .groups = "drop") %>%
      filter(pairs == 2) %>% 
      pull(Clinic_ID)
    
    plot_data <- df_stat_ready %>%
      filter(Clinic_ID %in% pairs.samp & Group %in% my_comparisons[[i]]) %>%
      arrange(Clinic_ID, Group) %>%
      reshape2::melt(id.vars = c("Group", "Clinic_ID"), variable.name = "Factor", value.name = "Value") %>%
      mutate(Group = fct_relevel(Group, top_factors))
    
    stat_data <- plot_data %>%
      group_by(Factor) %>%
      wilcox_test(Value ~ Group, paired = TRUE, detailed = TRUE) %>%  
      adjust_pvalue(method = "fdr")
    
  } else {
    # 外部横向独立样本检验
    plot_data <- df_stat_ready %>%
      filter(Group %in% my_comparisons[[i]]) %>%
      arrange(Group) %>%
      reshape2::melt(id.vars = c("Group", "Clinic_ID"), variable.name = "Factor", value.name = "Value") %>%
      mutate(Group = fct_relevel(Group, top_factors))
    
    stat_data <- plot_data %>%
      group_by(Factor) %>%
      wilcox_test(Value ~ Group, paired = FALSE, detailed = TRUE) %>%  
      adjust_pvalue(method = "fdr")
  }
  
  # 精确提取各组均值，规避硬编码索引带来的计算颠倒风险
  stat_plot <- plot_data %>%
    group_by(Factor, Group) %>%
    summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    reshape2::dcast(Factor ~ Group, value.var = "Mean_Value") %>%
    mutate(
      Mean1         = !!sym(t1), 
      Mean2         = !!sym(t2),
      Mean.diff     = Mean2 - Mean1,
      Mean.diff.log = sign(Mean.diff) * log1p(abs(Mean.diff))
    ) %>%
    dplyr::select(Factor, Mean1, Mean2, Mean.diff.log)
  
  stat_data <- merge(stat_data, stat_plot, by = "Factor", all = TRUE)
  comp_name <- paste(my_comparisons[[i]], collapse = "/")
  
  stats_results[[i]]     <- stat_data %>% mutate(comparisons = comp_name)
  plot_data_results[[i]] <- plot_data %>% mutate(comparisons = comp_name)
}

# 6.2 整合结果并追加多重检验显著性星号 -----------------------------------------
stats_results_df     <- bind_rows(stats_results)
plot_data_results_df <- bind_rows(plot_data_results)

stats_results_df <- stats_results_df %>%
  mutate(
    lab.p.adj = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", p.adj < 0.05 ~ "*", TRUE ~ "ns"),
    lab.p     = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
    label     = paste(lab.p, lab.p.adj, sep = "/")
  )

sum_stats_results     <- stats_results_df %>% mutate(Sites = "VA") %>% na.omit()
sum_plot_data_results <- plot_data_results_df %>% mutate(Sites = "VA") %>% na.omit()

# 导出原始统计数据 CSV
write.csv(sum_stats_results, paste0(prefix, ".TimeDiff.wilcox.csv"), row.names = FALSE)

# ==============================================================================
# 7. Visualization: Heatmap Grid Rendering (SCI Layout)
# ==============================================================================
message(">>> Filtering significant records and plotting validation heatmap...")

# 挑选出在至少一个对比组中 FDR 调整后显著 (p.adj < 0.05) 的所有关键菌株
sig_taxa <- sum_stats_results %>% filter(p.adj < 0.05) %>% select(Factor) %>% distinct()

if (nrow(sig_taxa) > 0) {
  draw <- sum_stats_results %>% 
    inner_join(sig_taxa, by = "Factor") %>% 
    mutate(groups = comparisons)
  
  # 剔除无意义的完全不显著标签符号显示
  draw$label  <- gsub("ns/ns", "", draw$label)
  
  # 规范设置横轴对比组因子的标准 SCI 排版顺序 (已经替换为了新版的 H_M, H_R)
  target_comparisons <- c("BL/T4", "BL/T12", "BL/T24", "BL/H_M", "BL/H_R", "T24/H_M", "T24/H_R")
  draw$groups        <- factor(draw$groups, levels = target_comparisons)
  
  # 开始绘制学术矩阵热图
  p_heat <- ggplot(draw, aes(x = groups, y = Factor, fill = Mean.diff.log)) +
    geom_tile(color = "lightgrey") +
    geom_text(aes(label = label), size = 2, color = "black") +  
    geom_vline(xintercept = c(3.5, 5.5), linetype = "dashed", color = "grey60") +
    facet_grid(Sites ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#5e3c99", mid = "grey98", high = "#b35806", midpoint = 0) +
    labs(x = NULL, y = NULL, fill = "Log-difference") +
    theme_minimal() +   
    theme(
      axis.text.y  = element_text(face = "bold.italic", size = 9),
      axis.text.x  = element_text(size = 10, angle = 30, hjust = 1, face = "bold"),
      panel.grid   = element_blank(),
      strip.text.y = element_text(angle = -90, face = "bold", size = 11),  
      panel.spacing.y = unit(0.5, "lines")  
    )
  
  print(p_heat)
  ggsave(paste0(prefix, ".TimeDiff.wilcox.heatmap.pdf"), plot = p_heat, width = 7, height = 6)
  write.csv(draw, paste0(prefix, ".TimeDiff.wilcox.heatmap.csv"), row.names = FALSE)
  message(">>> Statistical verification and Heatmap generated completely.")
} else {
  message("提示: 未发现任何经 FDR 校正后 p.adj < 0.05 的显著物种，自动跳过热图绘制。")
}

