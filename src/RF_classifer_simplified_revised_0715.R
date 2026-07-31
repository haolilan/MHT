# ==============================================================================
# PHASE 0: 环境与依赖配置 & 全局参数设定
# ==============================================================================
library(tidyverse)
library(randomForest)
library(caret)
library(pROC)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(ggrepel)
library(data.table)
library(glmnet)
library(compositions) # 用于 CLR
library(broom)        # 用于 GLM tidy 输出
library(conflicted)

conflicts_prefer(compositions::cor)
conflicts_prefer(psych::alpha)
conflicts_prefer(dplyr::filter)


base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

out_dir  <- "D:/WorkProjects/Demo-MHT 2026/Results/RF_Model_Blood_3Model/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)



Phen.Seq.BL <- Phen.Seq %>% subset(Time == "BL")%>%
  left_join(phen.enroll%>%select(Clinic_ID,Hospital=Hospital_BL))
phen.Cate.BL <- phen.Cate %>% subset(Phen_Time == "BL")

# 定义全局共享的静态常量参数（所有亚组一致）
IF_70_30                  <- TRUE  
SAMPLING_METHOD            <- "none"
TRANSFORM_METHODS          <- c("CLR")
FEATURE_SELECTION_METHODS  <- "Wilcoxon"
site_definitions           <- list("Gut" = "GUT") 
# target_list = "Response_K" 


vars.cate_filter <- c("K Score",  "Hormone","Glucose", "Lipid","K Sub Score")
vars.cate_filter_H <- c("K Score", "Hormone","Glucose", "Lipid", "K Sub Score")

subgroups_master_config <- list(
  "mKI" = list(
    cate_filter =vars.cate_filter,
    data_df     = phen.gr.K %>% 
      select(Clinic_ID, Response_K = Group) %>% 
      mutate(Response_K = factor(Response_K, levels = c("R", "NR"))),
    target_list = "Response_K" 
    
  ) ,
  "Hormone_E2" = list(
    cate_filter = vars.cate_filter_H,
    data_df     = phen.gr.E2.all %>%
      select(Clinic_ID, Response_E2 = Group) %>%
      mutate(Response_E2 = factor(Response_E2, levels = c("R", "NR"))),
    target_list = "Response_E2"
  )
  
)

# ==============================================================================
# 3. 自动化流水线总控大循环
# ==============================================================================
for (sg_name in names(subgroups_master_config)) {
  
  cat(sprintf("\n>>> [PIPELINE] Commencing automated execution grid for: [%s] ...\n", sg_name))
  
  # 3.1 动态解包提取当前亚组的特定参数
  current_cfg <- subgroups_master_config[[sg_name]]
  sub_data_df <- current_cfg$data_df
  sub_filter  <- current_cfg$cate_filter
  
  target_list <- current_cfg$target_list
  
  # 3.2 动态切换并创建与当前亚组完全对齐的保存目录
  base_out_dir <- file.path(out_dir, sg_name)
  if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE)
  setwd(base_out_dir)
  
  # 3.3 动态过滤临床特征子集：精准适配不同亚组的科学假说，并剔除胰岛素(INS),因为有效数据太少
  clinic.KH.vars <- phen.Cate.BL %>%
    filter(Category %in% sub_filter) %>%
    mutate(Category = factor(Category, levels = sub_filter)) %>%
    arrange(Category) %>%
    pull(Phen) %>%
    setdiff("INS") %>%
    unique()
  
  # 3.4 清洗样本目标集，排除不完整缺失
  phen_targets <- Phen.Seq.BL %>% 
    dplyr::select(Clinic_ID) %>%  
    left_join(sub_data_df, by = "Clinic_ID") %>%
    na.omit()
  
  # 3.5 质量卡控：打印当前的就绪简报，防止底层建模因为零频数闪退
  cat(sprintf("     Path Set: %s\n", base_out_dir))
  cat(sprintf("     Resolved Targets: N = %d subjects\n", nrow(phen_targets)))
  cat(sprintf("     Injected Clinical Features: N = %d vars [%s...]\n", length(clinic.KH.vars), clinic.KH.vars[1]))
  print(table(phen_targets$group))
  
  ## 以下为开始：######################################################################################
  threshold_map <- list("GUT" = 0.001)
  prevalence_cutoff <- 0.1 
  
  # 初始化全局结果记录器
  global_summary_list <- list()
  
  # ==============================================================================
  # PHASE 2 & 3 
  # ==============================================================================
  for(trans_method in TRANSFORM_METHODS) {
    cat(paste0("\n", rep("=", 80), collapse = ""), "\n")
    cat(sprintf("🌀 [TRANSFORMATION] 正在处理数据转换模式: %s\n", trans_method))
    
    # --------------------------------------------------------------------------
    # 2.1 针对当前转换方法生成 Profile (避免在内部循环重复计算)
    # --------------------------------------------------------------------------
    prof_BL_highpreval <- list()
    for (site_key in names(threshold_map)) {
      rel_ab.threshold <- threshold_map[[site_key]]
      data <- prof_filtered[[site_key]] %>%data.frame()
      
      Seq_id <- Microbe.phen.prof$SeqID[Microbe.phen.prof$Site == site_key & Microbe.phen.prof$Time == "BL"]
      data <- data %>% subset(rownames(data) %in% Seq_id)
      
      taxa_to_keep <- sapply(data, function(x) { sum(x > rel_ab.threshold) >= prevalence_cutoff * nrow(data) })
      data <- data[, taxa_to_keep, drop = FALSE]
      
      pseu_value <- min(data[data > 0], na.rm = TRUE) * 0.1
      data_imputed <- data + pseu_value
      
      if(trans_method == "CLR") {
        data_trans <- as.data.frame(unclass(clr(data_imputed)))
      } else {
        data_trans <- log10(data_imputed)
      }
      
      colnames(data_trans) <- paste(site_key, colnames(data_trans), sep = "_")
      data_trans$Clinic_ID <- Microbe.phen.prof$Clinic_ID[match(rownames(data_trans), Microbe.phen.prof$SeqID)]
      prof_BL_highpreval[[site_key]] <- data_trans
    }
    
    # 【核心修改 1】：合并多部位数据，但保留原生的 NA，不再填补任何人为制造的全局最小值
    merged_prof_BL_highpreval <- Reduce(function(x, y) merge(x, y, by = "Clinic_ID", all = TRUE), prof_BL_highpreval)
    
    taxa.species.vars_master <- setdiff(colnames(merged_prof_BL_highpreval), "Clinic_ID")
    
    # --------------------------------------------------------------------------
    # 2.2 进入特征筛选循环
    # --------------------------------------------------------------------------
    for(fs_method in FEATURE_SELECTION_METHODS) {
      cat(sprintf("  -> 🧪 [SELECTION] 当前特征筛选方法: %s\n", fs_method))
      
      # ------------------------------------------------------------------------
      # 2.3 进入身体部位循环
      # ------------------------------------------------------------------------
      for(site_name in names(site_definitions)) {
        current_site_codes <- site_definitions[[site_name]]
        
        exclude_sites <- setdiff(c("GUT", "VA", "TO", "UR"), current_site_codes)
        regex_pattern <- if(length(exclude_sites) > 0) paste0("^", exclude_sites, "_", collapse = "|") else "^X_X_X"
        
        taxa_site_cands <- taxa.species.vars_master[!grepl(regex_pattern, taxa.species.vars_master)]
        prof_BL_input_site <- merged_prof_BL_highpreval[, !grepl(regex_pattern, colnames(merged_prof_BL_highpreval))]
        
        # 【核心修改 2】：针对当前 site 筛选出的有效特征，强力剔除该范围内全是 NA 的无用样本
        # 排除 Clinic_ID 列后，计算当前样本在激活物种列上的非 NA 数量，只要为 0 则意味着该部位完全缺失
        valid_sample_idx <- rowSums(!is.na(prof_BL_input_site %>% select(all_of(taxa_site_cands)))) > 0
        prof_BL_input_site <- prof_BL_input_site[valid_sample_idx, , drop = FALSE]
        
        # ----------------------------------------------------------------------
        # 2.4 进入预测靶标循环
        # ----------------------------------------------------------------------
        for(current_target in target_list) {
          
          # 建立并进入嵌套子目录：base/TRANS/FS/SITE/TARGET/
          target_dir <- file.path(base_out_dir, trans_method, fs_method, site_name, paste0("Target_", current_target))
          if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
          setwd(target_dir)
          
          cat(sprintf("     >>> Target: %s | Site: %s | Pipeline: %s + %s <<<\n", 
                      current_target, site_name, trans_method, fs_method))
          
          phen.gr.current <- phen_targets %>% 
            select(Clinic_ID, Response = !!sym(current_target)) %>%
            filter(!is.na(Response)) %>%
            mutate(Response = factor(Response, levels = c("NR", "R")))
          
          if(nrow(phen.gr.current) < 20 || length(unique(phen.gr.current$Response)) < 2) next
          
          # 构建 Master Data
          basic.vars <- NULL 
          
          if (grepl("K", current_target, ignore.case = TRUE)) {
            basic.vars <- c(basic.vars, "K_Score","SAS","SDS","FSH","LH","Sleep_Score","Palpitations","Sexual_dysfunction") %>% intersect(clinic.KH.vars)
            cat(c("     [Info] ",paste0(basic.vars,collapse = " "),"纳入 basic.vars\n"))
          } else if (grepl("FSFI", current_target, ignore.case = TRUE)) {
            basic.vars <- c(basic.vars, "FSFI") %>% intersect(clinic.KH.vars)
            cat(c("     [Info] ",basic.vars,"纳入 basic.vars。\n"))
          } else if (grepl("FSH", current_target, ignore.case = TRUE)) {
            basic.vars <- c(basic.vars, "FSH") %>% intersect(clinic.KH.vars)
            cat(c("     [Info] ",basic.vars,"纳入 basic.vars。\n"))
          } else if (grepl("E2", current_target, ignore.case = TRUE)) {
            basic.vars <- c(basic.vars, "E2") %>% intersect(clinic.KH.vars)
            cat(c("     [Info] ",basic.vars,"纳入 basic.vars。\n"))
          }
          
          # clinical_features_master <- unique(c(basic.vars, clinic.KH.vars))
          clinical_features_master <- unique(c(clinic.KH.vars, basic.vars))
          
          # 数据集划分
          if(IF_70_30) {
            master_data <- Phen.Seq.BL %>% 
              select(Clinic_ID, all_of(clinical_features_master)) %>%
              left_join(prof_BL_input_site, by = "Clinic_ID") %>%
              left_join(phen.gr.current, by = "Clinic_ID") %>%
              select(Clinic_ID, Response, all_of(clinical_features_master), all_of(taxa_site_cands)) %>%
              filter(!is.na(Response)) %>% 
              na.omit() %>%
              mutate(across(all_of(clinical_features_master), ~ as.numeric(scale(.))))
            
            write.csv(master_data, paste0("0_Master_Data_", current_target, ".csv"), row.names = FALSE)
            
            set.seed(123)
            train_index <- createDataPartition(master_data$Response, p = 0.7, list = FALSE)
          } 
          
          train_data_full <- master_data[train_index, ]
          test_data_full  <- master_data[-train_index, ] 
          
          # 1. 计算整体的 Response 分布
          total_dist <- table(master_data$Response)
          train_dist <- table(train_data_full$Response)
          test_dist  <- table(test_data_full$Response)
          distribution_summary <- data.frame(
            Class = names(total_dist),
            Total_Count = as.numeric(total_dist),
            Train_Count = as.numeric(train_dist),
            Test_Count  = as.numeric(test_dist),
            Train_Ratio = round(as.numeric(train_dist) / as.numeric(total_dist), 3) 
          )
          cat("[Data Distribution] 训练集(Train)分布:\n")
          print(train_dist)
          cat("[Data Distribution] 测试集(Test/非Train)分布:\n")
          print(test_dist)
          write.csv(distribution_summary, 
                    paste0("0_Response_Distribution_", current_target, ".csv"), 
                    row.names = FALSE)
          
          # --------------------------------------------------------------------
          # 步骤 3：特征预选 (仅在 Train Set)
          # --------------------------------------------------------------------
          taxa_selected <- taxa_site_cands # 初始化
          
          if(fs_method == "Wilcoxon") {
            # --- Wilcoxon 筛选逻辑 ---
            long_data <- train_data_full %>% select(Response, all_of(taxa_site_cands)) %>%
              pivot_longer(cols = -Response, names_to = "Taxa", values_to = "Abundance")
            
            stats_summary <- long_data %>% group_by(Taxa, Response) %>%
              summarise(n = n(), Mean = mean(Abundance), Median = median(Abundance), SD = sd(Abundance), .groups = 'drop') %>%
              pivot_wider(names_from = Response, values_from = c(n, Mean, Median, SD), names_glue = "{.value}_{Response}")
            
            p_values <- sapply(taxa_site_cands, function(t) {
              ab <- train_data_full[[t]]; gp <- train_data_full$Response
              if(length(unique(ab))<=1 || length(unique(gp))<2) return(NA)
              wilcox.test(ab ~ gp, exact = FALSE)$p.value
            })
            
            taxa_stats <- data.frame(Taxa = names(p_values), P_value = p_values) %>%
              filter(!is.na(P_value)) %>% left_join(stats_summary, by = "Taxa") %>%
              mutate(FDR = p.adjust(P_value, method = "BH")) %>% arrange(P_value)
            
            write.csv(taxa_stats, paste0("1_Wilcoxon_Stats_", current_target, ".csv"), row.names = FALSE)
            
            sig_taxa <- taxa_stats %>% filter(P_value < 0.2)
            if(nrow(sig_taxa) < 15) taxa_selected <- taxa_stats$Taxa[1:min(15, nrow(taxa_stats))]
            else if(nrow(sig_taxa) > 50) taxa_selected <- sig_taxa$Taxa[1:50]
            else taxa_selected <- sig_taxa$Taxa
            
          } 
          
          # --------------------------------------------------------------------
          # 步骤 4：构建 RF 训练集并建模 (加入了纯临床对照模型)
          # --------------------------------------------------------------------
          RF.data.list <- list(
            "Features_Clinical_Only" = master_data %>% select(Response,  all_of(clinical_features_master)),
            "Features_Microbe_only" = master_data %>% select(Response, all_of(taxa_selected)),
            "Features_Microbe_Clinical" = master_data %>% select(Response,all_of(clinical_features_master), all_of(taxa_selected))
          )
          

          results_rf <- list()
          
          # 开启并行计算加速调参
          cl <- makeCluster(max(1, min(10,detectCores() - 2)))
          registerDoParallel(cl)
          
          for (model_name in names(RF.data.list)) {
            cat(sprintf("      -> Training RF: %s\n", model_name))
            
            data_sub <- RF.data.list[[model_name]]
            train_data  <- data_sub[train_index, ]
            test_data   <- data_sub[-train_index, ]
            
            # 模型调参控制 #################################################################
            # 0. 确保安装了需要的底层依赖包
            if(!require(smotefamily)) install.packages("smotefamily")
            if(!require(pROC)) install.packages("pROC")
            library(pROC)
            library(randomForest)
            
            # ====================================================================
            # 1. 升级自定义 RFE 函数集：让它完整支持 ROC/AUC 的计算
            # ====================================================================
            custom_rfFuncs <- rfFuncs
            
            # 彻底重写 pred 函数：确保同时预测类别和概率 (完美兼容并行)
            custom_rfFuncs$pred <- function(modelFit, testX) {
              res <- predict(modelFit, testX)
              prob <- predict(modelFit, testX, type = "prob")
              
              out <- data.frame(pred = res)
              if(!is.null(prob)) out <- cbind(out, prob)
              out
            }
            
            # 将 RFE 内部的评估机制绑定为 twoClassSummary
            custom_rfFuncs$summary <- twoClassSummary
            
            # ====================================================================
            # 2. 动态计算特征范围（必须提前计算，因为生成种子需要用到长度）
            # ====================================================================
            total_feats <- ncol(train_data) - 1
            
            if (total_feats > 80) {
              sizes_to_test <- unique(c(1:80, seq(85, min(total_feats, 100), by = 5), total_feats))
            } else {
              sizes_to_test <- 1:total_feats
            }
            
            # ====================================================================
            # 为并行 RFE 严格构建全套随机种子矩阵
            # ====================================================================
            set.seed(224) # 基础随机源   ## 224 is the best 0.771 
            
            # 5折CV对应5个常规Fold + 1个最终完整重训，共需要 6 个 List 元素
            rfe_seeds <- vector(mode = "list", length = 6) 
            
            # 前5个Fold，每个Fold内部需要为每一种 sizes_to_test 的取值分配一个独立种子
            for(i in 1:5) {
              rfe_seeds[[i]] <- sample.int(n = 10000, size = length(sizes_to_test))
            }
            # 第6个元素是给最后锁定最佳特征量后，在全训练集上重训最终模型使用的种子
            rfe_seeds[[6]] <- sample.int(n = 10000, size = 1)
            saveRDS(rfe_seeds, file = paste0("RFE_Gold_Seeds_",  model_name, ".rds"))
            ##
            names(rfe_seeds) <- c(paste0("Fold_", 1:5), "Final_Re-train")
            seeds_df <- tibble::enframe(rfe_seeds, name = "Fold_Index", value = "Seed_Value") %>%
              tidyr::unnest_longer(Seed_Value)
            write.csv(seeds_df, paste0("RFE_Gold_Seeds_", model_name,".csv"), row.names = FALSE)
            
            # ====================================================================
            # 4. 配置针对 ROC 优化且【锁死种子】的 rfeControl
            # ====================================================================
            if (SAMPLING_METHOD == "none") {
              rfe_ctrl <- rfeControl(functions = custom_rfFuncs, 
                                     method = "cv", 
                                     number = 5, 
                                     returnResamp = "all",      
                                     saveDetails = TRUE,   
                                     allowParallel = TRUE,       
                                     verbose = FALSE,
                                     seeds = rfe_seeds)          
            } else {
              rfe_ctrl <- rfeControl(functions = custom_rfFuncs, 
                                     method = "cv", 
                                     number = 5, 
                                     returnResamp = "all", 
                                     saveDetails = TRUE,
                                     allowParallel = TRUE, 
                                     verbose = FALSE,
                                     sampling = SAMPLING_METHOD, 
                                     seeds = rfe_seeds)         
            }
            
            # 3. 运行带有 5 折交叉验证的递归特征筛选 (RFECV)
            cat(sprintf("      -> 🚀 正在执行 5-Fold CV 递归特征筛选 (基于 ROC 评估)... 候选特征量: [%s]\n", 
                        paste(sizes_to_test, collapse = ", ")))
            
            # 分离特征和标签
            X_train <- train_data %>% select(-Response)
            Y_train <- train_data$Response
            
            # 运行筛选
            rf_rfe_profile <- rfe(x = X_train, 
                                  y = Y_train,
                                  sizes = sizes_to_test, 
                                  rfeControl = rfe_ctrl,
                                  metric = "ROC",           
                                  ntree = 500)
            
            # ====================================================================
            # 4. 挑选出的“最佳特征数量”和“最佳特征列表”
            # ====================================================================
            best_feature_count <- rf_rfe_profile$bestSubset
            taxa_selected_cv <- predictors(rf_rfe_profile)
            write_rds(rf_rfe_profile, paste0("1_RFECV_Selected_Features_", model_name, ".RDS"))
            
            cat(sprintf("         [📊 对比] 绝对最高 ROC 对应特征量：%d 个\n", best_feature_count))
            
            ##绝对最高 ROC 分布 #################################
            # 1. 提取 RFE 内部所有特征数量候选集的 CV 结果
            rfe_results <- rf_rfe_profile$results
            
            # 2. 锁定最佳特征量所在的行，用于在图上做高亮红点标记
            best_perf <- rfe_results %>% filter(Variables == best_feature_count)
            
            # ====================================================================
            # 3. 使用 ggplot2 绘制学术级趋势折线图
            # ====================================================================
            p_rfe_trend <- ggplot(rfe_results, aes(x = Variables, y = ROC)) +
              geom_ribbon(aes(ymin = ROC - ROCSD, ymax = ROC + ROCSD),
                          fill = "#1f4e79", alpha = 0.1) +
              geom_line(color = "#134074", size = 1) +
              geom_point(color = "#134074", size = 1.5, alpha = 0.7) +
              geom_point(data = best_perf, aes(x = Variables, y = ROC),
                         color = "#d84315", size = 3.5, shape = 16) +
              geom_vline(xintercept = best_feature_count, linetype = "dashed", color = "#d84315", size = 0.5) +
              geom_text(data = best_perf,
                        aes(x = Variables, y = ROC,
                            label = sprintf(" Best Size = %d\n (CV ROC = %.3f)", Variables, ROC)),
                        hjust = -0.1, vjust = -1.1, size = 3, fontface = "bold", color = "#d84315") +
              scale_x_continuous(breaks = seq(0, max(rfe_results$Variables), by = 10)) +
              scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, by = 0.1)) +
              labs(x = "Number of Microbe Features (Variables)",
                   y = "5-Fold CV Mean ROC (AUC)",
                   title = paste("Feature Selection Curve -", model_name),
                   caption = "Shaded region represents ±1 SD across 5 cross-validation folds") +
              theme_bw() +
              theme(
                plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0.5),
                axis.title = element_text(face = "bold", size = 10),
                panel.grid.major =  element_blank(),
                panel.grid.minor = element_blank()
              )
            ggsave(filename = paste0("1_Feature_Count_Curve_", model_name, ".pdf"),
                   plot = p_rfe_trend, width = 6.5, height = 5)
            
            # --------------------------------------------------------------------
            # 展示 Train Set 5个 Fold 分别的 ROC 曲线
            # --------------------------------------------------------------------
            library(pROC)
            library(ggplot2)
            library(dplyr)
            library(tidyr)
            
            # 1. 提取出最佳特征数量下的所有 CV 预测结果
            cv_predictions <- rf_rfe_profile$pred %>% 
              filter(Variables == best_feature_count) # 只看最佳特征组合下的表现
            
            # ====================================================================
            # 计算全训练集总体的 ROC 以及 95% 灰色置信区间 (CI)
            # ====================================================================
            total_roc <- roc(cv_predictions$obs, cv_predictions$R, levels = c("NR", "R"), direction = "<", quiet = TRUE)
            total_auc <- as.numeric(auc(total_roc))
            
            # ====================================================================
            # 【提取多折曲线】准备 5 个 Fold 的各自折线数据
            # ====================================================================
            cv_folds <- unique(cv_predictions$Resample)
            df_folds_list <- list()
            fold_auc_list <- c()
            
            for (i in seq_along(cv_folds)) {
              current_fold <- cv_folds[i]
              fold_data <- cv_predictions %>% filter(Resample == current_fold)
              
              fold_roc <- roc(fold_data$obs, fold_data$R, levels = c("NR", "R"), direction = "<", quiet = TRUE)
              fold_auc_list <- c(fold_auc_list, as.numeric(auc(fold_roc)))
              
              df_folds_list[[current_fold]] <- data.frame(
                Specificity = fold_roc$specificities,
                Sensitivity = fold_roc$sensitivities,
                Fold        = sprintf("%s (AUC = %.3f)", current_fold, as.numeric(auc(fold_roc)))
              )
            }
            df_folds_all <- do.call(rbind, df_folds_list)
            
            # ====================================================================
            # 【ggplot2 精英绘图】应用低饱和度、非同色系
            # ====================================================================
            morandi_palette <- c("#7A8B99", "#A77E71", "#70877F", "#96859A", "#C4A46F")
            
            mean_label <- sprintf("Overall Mean (AUC = %.3f)", mean(fold_auc_list))
            
            p_roc <- ggplot() +
              geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), linetype = "dashed", color = "gray60", size = 0.5) +
              geom_path(data = df_folds_all, aes(x = Specificity, y = Sensitivity, color = Fold), size = 0.9) +
              scale_x_reverse(expand = c(0.01, 0.01), limits = c(1, 0)) +
              scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
              scale_color_manual(values = morandi_palette) +
              labs(x = "Specificity (True Negative Rate)", 
                   y = "Sensitivity (True Positive Rate)", 
                   title = paste("5-Fold CV ROC Curves -", model_name),
                   subtitle = paste("Best Features Number =",best_feature_count,"|",mean_label),
                   color = "Cross-Validation Folds") +
              theme_bw() + 
              theme(
                plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0.5),
                axis.title = element_text(face = "bold", size = 10),
                legend.title = element_text(face = "bold", size = 9),
                legend.text = element_text(size = 8),
                legend.position = c(0.72, 0.22), 
                # legend.background = element_rect(fill = alpha("white", 0.7), color = "gray90"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
              )
            p_roc
            ggsave(filename = paste0("1_Train_BestFeatureCount_5Fold_CV_ROC_", model_name, ".pdf"), 
                   plot = p_roc, width = 6, height = 6)
            
            cat(sprintf("         [🎉 ROC输出成功] 平均 5-Fold CV AUC: %.4f\n", mean(fold_auc_list)))
            
            # --------------------------------------------------------------------
            # 5. 使用挑选出的最佳特征子集，重新进行最终模型的网格调参 (mtry 优化)
            # --------------------------------------------------------------------
            train_data_selected <- train_data %>% select(Response, all_of(taxa_selected_cv))
            test_data_selected  <- test_data %>% select(Response, all_of(taxa_selected_cv))
            
            # 重新配置标准最终训练的控制参数
            if (SAMPLING_METHOD == "none") {
              ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, allowParallel = TRUE, 
                                   summaryFunction = twoClassSummary)
            } else {
              ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, allowParallel = TRUE, 
                                   summaryFunction = twoClassSummary, sampling = SAMPLING_METHOD)
            }
            
            # 动态动态设置 mtry 网格（防止特征太少时报错）
            tune_grid <- expand.grid(mtry = seq(1, max(1, min(10, length(taxa_selected_cv) - 1)), by = 1))
            
            set.seed(123)
            rf_caret <- train(Response ~ ., 
                              data = train_data_selected, 
                              method = "rf", 
                              metric = "ROC", 
                              ntree = 1000,         
                              importance = TRUE, 
                              trControl = ctrl, 
                              tuneGrid = tune_grid)
            
            # 6. 提取内置的最终最优模型
            final_model <- rf_caret$finalModel
            write_rds(final_model, paste0("2_Train_Final_Model_", model_name, ".RDS"))
            
            # --------------------------------------------------------------------
            # 7. 在 Test 独立测试集上评估性能并输出真正的验证 ROC
            # --------------------------------------------------------------------
            pred_prob <- predict(rf_caret, test_data_selected, type = "prob")[, "R"]
            roc_obj <- roc(test_data_selected$Response, pred_prob, levels = c("NR", "R"), direction = "<", quiet = TRUE)
            
            pred_class <- predict(rf_caret, test_data_selected)
            cm <- confusionMatrix(data = pred_class, reference = test_data_selected$Response, positive = "R")
            
            # 提取最终筛选后特征的 MeanDecreaseAccuracy 重要性
            imp_data <- importance(final_model) %>% 
              as.data.frame() %>% 
              rownames_to_column("Feature") %>% 
              arrange(desc(MeanDecreaseAccuracy))
            
            # 保存该模型的核心结果
            results_rf[[model_name]] <- list(
              AUC = as.numeric(auc(roc_obj)), 
              Accuracy = as.numeric(cm$overall["Accuracy"]),
              Sensitivity = as.numeric(cm$byClass["Sensitivity"]), 
              Specificity = as.numeric(cm$byClass["Specificity"]),
              Importance = imp_data
            )
            cat(sprintf("         [🎉 ROC输出成功] TEST AUC: %.4f\n", results_rf[[model_name]]$AUC))
            
            # --------------------------------------------------------------------
            # 8. 绘制并输出高分辨率测试集独立验证 ROC 曲线 PDF
            # --------------------------------------------------------------------
            pdf(file = paste0("3_Test_ROC_Curve_", model_name, ".pdf"), width = 6, height = 6)
            plot(roc_obj, 
                 main = paste("Test ROC -", model_name),
                 col = "#e74c3c", lwd = 3, legacy.axes = TRUE, print.auc = TRUE,
                 auc.polygon = TRUE, auc.polygon.col = "#fdedec")
            dev.off()
            
            # 9. 重要性特征输出与你的棒棒糖图落盘 (自动适应新筛选出的特征数量)
            plot_data <- results_rf[[model_name]]$Importance
            write.csv(plot_data, paste0("4_RF_Importance_", model_name, ".csv"), row.names = FALSE)
            
            p_rf_imp <- ggplot(plot_data %>% head(min(30, nrow(plot_data))), aes(x = MeanDecreaseAccuracy, y = reorder(Feature, MeanDecreaseAccuracy))) +
              geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = Feature, yend = Feature), color = "gray70", linetype = "dashed") +
              geom_point(aes(color = MeanDecreaseAccuracy), size = 4) + 
              scale_color_gradient(low = "#85c1e9", high = "#e74c3c") +
              labs(x = "Mean Decrease Accuracy (Importance)", y = "Feature", color = "Accuracy", 
                   title = paste("Features:", model_name), 
                   subtitle = sprintf("Test set: AUC=%.3f\nBest Feature Count = %d", 
                                      results_rf[[model_name]]$AUC, best_feature_count),
                   caption = sprintf("Test: Acc=%.3f | Sens=%.3f | Spec=%.3f", 
                                     results_rf[[model_name]]$Accuracy,
                                     results_rf[[model_name]]$Sensitivity, results_rf[[model_name]]$Specificity)) +
              theme_minimal() + 
              theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10, face = "bold"))
            
            ggsave(paste0("4_RF_Importance_Plot_", model_name, ".pdf"), plot = p_rf_imp, width = 8, height = 8)
            
            # 将汇总结果压入全局内存记录器 (用于最终的上帝视角图)
            global_summary_list[[length(global_summary_list) + 1]] <- data.frame(
              Trans_Method = trans_method,
              FS_Method = fs_method,
              Folder_Name = site_name,
              Target_Grouping = current_target,
              Model_Type = model_name,
              Total_Features = ncol(train_data_selected) - 1,
              AUC = results_rf[[model_name]]$AUC,
              Accuracy = results_rf[[model_name]]$Accuracy,
              Sensitivity = results_rf[[model_name]]$Sensitivity,
              Specificity = results_rf[[model_name]]$Specificity
            )
          }
          
          try(stopCluster(cl), silent = TRUE)
        }
      }
    }
  }
  # ==============================================================================
  # PHASE 4
  # ==============================================================================
  cat("\n>>> 🎉 所有计算完成，正在生成终极上帝视角报告图表...\n")
  setwd(base_out_dir)
  
  final_master_summary <- bind_rows(global_summary_list) %>% arrange(desc(AUC))
  
  write.csv(final_master_summary, "0_Global_All_Models_AUC_Summary.csv", row.names = FALSE)
  
  
  # 绘图整理数据
  colors_line <- c(
    "Features_Clinical_Only" = "grey20",
    "Features_Microbe_only"          = "#2ca02c", 
    "Features_Microbe_Clinical"      = "#ff7f0e"
  )
  final_master_summary$Model_Type <- factor(final_master_summary$Model_Type, levels = names(colors_line))
  final_master_summary$Folder_Name <- factor(final_master_summary$Folder_Name, levels = names(site_definitions))
  final_master_summary$Target_Grouping <- factor(final_master_summary$Target_Grouping, levels = target_list)
  
  
  # 创建终极多维分面柱状图 (Grid: Trans + FS ~ Site)
  p_gods_eye <- ggplot(final_master_summary, aes(x = Model_Type, y = AUC, fill = Model_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "white") +
    geom_text(aes(label = sprintf("%.3f", AUC)), position = position_dodge(width = 0.8), vjust = -0.5, size = 2.8, fontface = "bold") +
    geom_text(aes(y = 0.35, label = paste0("F.N=", Total_Features)), position = position_dodge(width = 0.8), size = 2, color = "white", fontface="bold") +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey60", linewidth = 0.5, alpha = 0.5) +
    coord_cartesian(ylim = c(0.3, max(final_master_summary$AUC, na.rm=TRUE) + 0.1)) +
    scale_fill_manual(values = colors_line) +
    labs(title = "Global Pipeline Performance Matrix", x = "", y = "Test Set AUC", fill = "") +
    facet_grid(~Trans_Method)+
    theme_bw(base_size = 14) + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = "right", panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "#2C3E50", color = "white"),
      strip.text = element_text(face = "bold", size = 11, color = "white")
    ) 
  p_gods_eye
  
  ggsave("0_Global_AUC_Pipeline_Matrix.pdf", plot = p_gods_eye, width = 6, height = 5)
  
  cat(sprintf("🚀 任务全部结束！请前往 %s 查看全局摘要 (0_Global_All_Models_AUC_Summary.csv) 和最终矩阵图表！\n", base_out_dir))
  
}

