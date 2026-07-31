# ==============================================================================
# 0. 加载所有必要的学术绘图与分析包
# ==============================================================================
library(dplyr)
library(tidyverse)
library(lme4)
library(mediation)
library(ggeffects)
library(ggplot2)
# if (!requireNamespace("compositions", quietly = TRUE)) install.packages("compositions")
library(compositions) 
library(conflicted)

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))
Phen.Seq.BL     <- Phen.Seq%>% subset(Time == "BL")

target_dir <- paste0(base_dir, "/Results/Mediation/")
if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)
 
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(mediation::mediate)

# ==============================================================================
# 1. 原始丰度矩阵的 CLR 预处理 (自动化数据准备)
# ==============================================================================
# 提取带有 SeqID 的原始数据
rel_ab.threshold <- 0.0001
data_va <- prof_filtered$VA%>%data.frame(.)

taxa_keep_va <- sapply(data_va, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(data_va))
taxa.vars_va <- names(taxa_keep_va)[taxa_keep_va]
taxa_matrix  <- data_va[, taxa.vars_va, drop = FALSE]

# 自动计算零置换伪值 (用最小非零值的 10% 补零)
pseu_value <- min(taxa_matrix[taxa_matrix > 0], na.rm = TRUE) * 0.1
data_imputed <- taxa_matrix + pseu_value

# 执行中心对数比值变换 (CLR Transformation)
data_trans <- as.data.frame(unclass(compositions::clr(data_imputed)))

# 将加工好的 CLR 丰度矩阵重新融合回带 SeqID 的标准作图格式 dat.VA
dat.VA <- data_trans %>% rownames_to_column("SeqID")


# ==============================================================================
# 2. 自定义你想批量循环分析的核心阴道菌列表
# ==============================================================================
microbe_list <- c("Lactobacillus_crispatus",  "Lactobacillus_iners", 
                  "Gardnerella_vaginalis",    "Fannyhessea_vaginae",
                  "Aerococcus_christensenii", "Dialister_micraerophilus")

# ==============================================================================
# 3. 【核心自动化】多菌循环分析引擎 (构建存储器)
# ==============================================================================
# 创建用于收集所有菌中介结果的作图空表格
mediation_forest_table <- data.frame()

# 开始自动化多菌大循环
for (target_microbe in microbe_list) {
  
  cat("\n====================================================================\n")
  cat("👉 正在执行核心自动化中介分析，目标靶点菌:", target_microbe, "\n")
  cat("====================================================================\n")
  
  # 检查数据中是否存在该菌，防止拼写错误中断循环
  if (!target_microbe %in% colnames(dat.VA)) {
    cat("⚠️ 警告: 菌名", target_microbe, "不在丰度表中，已自动跳过。\n")
    next
  }
  
  # 3.1 动态组装特定菌的长数据 metadata
  long_data <- dat.VA %>% 
    select(SeqID, Target_Microbe = all_of(target_microbe)) %>%
    left_join(Microbe.phen.prof %>% select(SeqID, Clinic_ID, Time), by = "SeqID") %>%
    left_join(Phen.Seq.BL %>% select(Clinic_ID, Age, BMI), by = "Clinic_ID") %>%
    left_join(Phen.Seq %>% select(Clinic_ID, Time, FSFI), by = c("Clinic_ID", "Time")) %>%
    filter(Time %in% c("BL", "T04")) %>%
    na.omit()
  
  # 规范化因子的时点基准线
  long_data$Time <- factor(long_data$Time, levels = c("BL", "T04"))
  
  set.seed(123)
  # 3.2 运行双混合效应模型 (LMM)
  model.M_mixed <- lmer(Target_Microbe ~ Time + Age + BMI + (1 | Clinic_ID), data = long_data)
  model.Y_mixed <- lmer(FSFI ~ Time + Target_Microbe + Age + BMI + (1 | Clinic_ID), data = long_data)
  
  # 3.3 运行配对中介分析
  results_mixed <- mediate(model.M_mixed, model.Y_mixed, 
                           treat = "Time", mediator = "Target_Microbe", sims = 1000)
  
  # 实时展示当前菌的简要统计看板
  cat("\n[统计汇总 - ", target_microbe, "]\n")
  print(summary(results_mixed))
  
  # 3.4 【自动化信息抓取】提取森林图所需的 4 大金刚核心指标
  med_summary <- summary(results_mixed)
  
  current_res <- data.frame(
    Microbe     = target_microbe,
    Estimate    = med_summary$d0,                # ACME 估计值
    CI_Lower    = med_summary$d0.ci[1],           # 95% 置信区间下限
    CI_Upper    = med_summary$d0.ci[2],           # 95% 置信区间上限
    P_Value     = med_summary$d0.p,                  # ACME p值
    Prop_Med    = med_summary$n0,                # 中介占比估计值
    Prop_P      = med_summary$n0.p,                   # 中介占比p值
    
    ADE.Estimate    = med_summary$z0,                # ADE 估计值
    ADE.CI_Lower    = med_summary$z0.ci[1],           # 95% 置信区间下限
    ADE.CI_Upper    = med_summary$z0.ci[2],           # 95% 置信区间上限
    ADE.P_Value     = med_summary$z0.p
  )
  
  # 将当前菌的战果合并到总表
  mediation_forest_table <- rbind(mediation_forest_table, current_res)
  
  # 3.5 【选做保存】自动为当前菌绘制并导出对齐校正后的偏残差图
  eff_Path_b <- ggpredict(model.Y_mixed, terms = "Target_Microbe")
  plot_data_line <- as.data.frame(eff_Path_b)
  
  beta_microbe <- fixef(model.Y_mixed)["Target_Microbe"]
  intercept    <- fixef(model.Y_mixed)["(Intercept)"]
  beta_age     <- fixef(model.Y_mixed)["Age"]
  beta_bmi     <- fixef(model.Y_mixed)["BMI"]
  mean_age     <- mean(long_data$Age, na.rm = TRUE)
  mean_bmi     <- mean(long_data$BMI, na.rm = TRUE)
  
  y_matched_residuals <- (beta_microbe * long_data$Target_Microbe) + 
    intercept + (beta_age * mean_age) + (beta_bmi * mean_bmi) + 
    residuals(model.Y_mixed, type = "response")
  
  plot_data_points <- data.frame(x = long_data$Target_Microbe, y = y_matched_residuals)
  
  p_partial <- ggplot() +
    geom_ribbon(data = plot_data_line, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#bdc3c7", alpha = 0.4) +
    geom_point(data = plot_data_points, aes(x = x, y = y, color = y), size = 2, alpha = 0.6) +
    geom_line(data = plot_data_line, aes(x = x, y = predicted), color = "#1a5276", size = 1.2) +
    scale_color_viridis_c(option = "D", name = "Adjusted FSFI") +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(title = paste("Path b Partial Residual Plot:", target_microbe),
         x = paste("CLR Abundance of", target_microbe), y = "Partial Residuals of FSFI")
  
  # 自动导出图片文件到本地磁盘，省时省力
  ggsave(filename = paste0("Partial_Residual_Plot_", target_microbe, ".pdf"), plot = p_partial, width = 7, height = 6)
}


# ==============================================================================
# 4. 输出并查看最终绘制森林图的整合表格 (作图表格输出完成)
# ==============================================================================
cat("\n\n🎉 循环全部结束！已为您生成终极作图数据汇总表格：\n")


mediation_forest_table <- mediation_forest_table %>%
  mutate(FDR = p.adjust(P_Value, method = "BH"),
         ADE.FDR = p.adjust(ADE.P_Value, method = "BH"),
         Prop_P.FDR = p.adjust(Prop_P , method = "BH")) %>%
  
  mutate(P_Label = ifelse(FDR < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", FDR)))%>%
  mutate(ADE.P_Label = ifelse(ADE.FDR < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", ADE.FDR)))


print(mediation_forest_table)
write.csv(mediation_forest_table, "Mediation_Results_Total.csv") 


# ==============================================================================
# 5. 统一实现高颜值多菌联合学术森林图 (Forest Plot)
# ==============================================================================
# 固化其物理排列顺序（防止默认字母乱序排列）
mediation_forest_table$Microbe <- factor(mediation_forest_table$Microbe, levels = rev(microbe_list))

# 自动依据正负中介效应标记不同分类上色（正向促进 vs 负向压制）
mediation_forest_table <- mediation_forest_table %>%
  mutate(Effect_Type = ifelse(Estimate > 0 & P_Value < 0.05, "Positive Mediation",
                              ifelse(Estimate < 0 & P_Value < 0.05, "Negative Mediation", "Not Significant")))


# 查看包含 FDR 的新表格
print(mediation_forest_table)
ggplot(mediation_forest_table, aes(x = Estimate, y = Microbe, color = Effect_Type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray20", size = 0.8) +
  geom_pointrange(aes(xmin = CI_Lower, xmax = CI_Upper), size = 0.8, fatten = 4) +
  geom_text(aes(x = Estimate, label = paste0(P_Label)), 
            hjust = -0.1,vjust = -0.5, color = "black", fontface = "italic", size = 2.5) +
  scale_color_manual(values = c("Positive Mediation" = "#1a5276", 
                                "Negative Mediation" = "#b03a2e", 
                                "Not Significant"    = "#7f8c8d")) + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.y = element_text(face = "italic"),                   
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Mediation Effect of Vaginal Microbes",
    x = "Indirect Mediation Effect Size",
    y = "",
    color = "Statistical Inference"
  ) +
  xlim(min(mediation_forest_table$CI_Lower) - 0.05, max(mediation_forest_table$CI_Upper) + 0.2)

ggsave(filename = paste0("FOREST_ACME_", "FSFI_Vaginal" , ".pdf"),  width = 7, height = 6)

