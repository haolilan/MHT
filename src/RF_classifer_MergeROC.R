library(pROC)
library(dplyr)
library(ggplot2)


setwd("D:/WorkProjects/Demo-MHT 2026/Results/RF_Model_Blood_3Model/mKI/CLR/Wilcoxon/Gut/Target_Response_K/")

# 1. 根据截图定义这三个文件的名称，并给它们指定简短的模型标签
file_list <- c(
  "Clinical_Only"    = "SourceData_Test_Predictions_Features_Clinical_Only.csv",
  "Microbe_Only"     = "SourceData_Test_Predictions_Features_Microbe_only.csv",
  "Microbe_Clinical" = "SourceData_Test_Predictions_Features_Microbe_Clinical.csv"
)

# 初始化空列表，用于存储提取的数据和 ROC 对象
roc_data_list <- list()
roc_objects   <- list() # 顺便保留原始 roc_obj，方便以后做 DeLong 检验

# 2. 循环读取文件并提取 ROC 坐标
for (model_name in names(file_list)) {
  
  # 读取对应的预测结果
  pred_data <- read.csv(file_list[model_name])
  
  # 计算 ROC
  roc_obj <- roc(
    response = pred_data$Actual_Response, 
    predictor = pred_data$Pred_Prob_R, 
    levels = c("NR", "R"), 
    direction = "<", 
    quiet = TRUE
  )
  
  roc_objects[[model_name]] <- roc_obj
  
  # 提取当前模型的图源数据，新增一列 "Model" 以区分数据来源
  roc_source_data <- data.frame(
    Model               = model_name,
    Threshold           = roc_obj$thresholds,
    Sensitivity         = roc_obj$sensitivities,
    Specificity         = roc_obj$specificities,
    `1_Minus_Specificity` = 1 - roc_obj$specificities,
    check.names         = FALSE # 防止 R 语言自动把以数字开头的列名加上 "X" 前缀
  )
  
  roc_data_list[[model_name]] <- roc_source_data
}

# 3. 将三个模型的 ROC 坐标数据合并为一张总表
all_roc_source_data <- bind_rows(roc_data_list) %>%
  # 【关键修复 1】：强制按 X 轴和 Y 轴从小到大排序，彻底消除锯齿和线条折返
  arrange(Model, `1_Minus_Specificity`, Sensitivity)

# 导出整合后的源数据
write.csv(all_roc_source_data, "Combined_ROC_Curve_SourceData.csv", row.names = FALSE, quote = FALSE)

# ==============================================================================
# 4. 直接使用合并后的数据画出三条对比曲线
# ==============================================================================
p_roc <- ggplot(all_roc_source_data, aes(x = `1_Minus_Specificity`, y = Sensitivity, color = Model, fill = Model)) +
  geom_area(alpha = 0.1, position = "identity", color = NA) +
  
  geom_path(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") + 

  scale_color_manual(
    values = c(
      "Clinical_Only" = "grey20", 
      "Microbe_Only" = "#2ca02c", 
      "Microbe_Clinical" = "#ff7f0e"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Clinical_Only" = "grey20", 
      "Microbe_Only" = "#2ca02c", 
      "Microbe_Clinical" = "#ff7f0e"
    )
  ) +
  labs(
    title = "Test ROC Curves Comparison",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = c(0.75, 0.25), 
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by = 0.2))

print(p_roc)
ggsave("Combined_ROC_Curves_Fixed.pdf", plot = p_roc, width = 6, height = 6)
