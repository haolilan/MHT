
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
library(compositions) 
library(broom)    
library(conflicted)
conflicts_prefer(compositions::cor)
conflicts_prefer(psych::alpha)
conflicts_prefer(dplyr::filter)

###
base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))
Phen.Seq.BL <- Phen.Seq %>% subset(Time == "BL")%>%
  left_join(phen.enroll%>%select(Clinic_ID,Hospital=Hospital_BL))
phen.Cate.BL <- phen.Cate %>% subset(Phen_Time == "BL")

###
out_dir  <- "D:/WorkProjects/Demo-MHT 2026/Results/RF_Model_Blood_3Model/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

# Global configuration: cohort partitioning and preprocessing settings
IF_70_30                  <- TRUE   # TRUE: 70/30 stratified random split; FALSE: hospital-based split
SAMPLING_METHOD           <- "none" # Options: smote, down, up, none
TRANSFORM_METHODS          <- c("CLR")
FEATURE_SELECTION_METHODS  <- "Wilcoxon" 
site_definitions           <- list("Gut" = "GUT") 

vars.cate_filter <- c("K Score",  "Hormone","Glucose", "Lipid","K Sub Score")
subgroups_master_config <- list(
  "mKI" = list(
    cate_filter =vars.cate_filter,
    data_df     = phen.gr.K %>% 
      select(Clinic_ID, Response_K = Group) %>% 
      mutate(Response_K = factor(Response_K, levels = c("R", "NR"))),
    target_list = "Response_K" 
    
  ),
  "Hormone_E2" = list(
    cate_filter = vars.cate_filter,
    data_df     = phen.gr.E2.all %>%
      select(Clinic_ID, Response_E2 = Group) %>%
      mutate(Response_E2 = factor(Response_E2, levels = c("R", "NR"))),
    target_list = "Response_E2"
  )
)

for (sg_name in names(subgroups_master_config)) {
  
  cat(sprintf("\n>>> [PIPELINE] Commencing automated execution grid for: [%s] ...\n", sg_name))

  current_cfg <- subgroups_master_config[[sg_name]]
  sub_data_df <- current_cfg$data_df
  sub_filter  <- current_cfg$cate_filter
  target_list <- current_cfg$target_list

  base_out_dir <- file.path(out_dir, sg_name)
  if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE)
  setwd(base_out_dir)
  
  # Select baseline clinical features matching target categories; exclude INS due to low coverage
  clinic.KH.vars <- phen.Cate.BL %>%
    filter(Category %in% sub_filter) %>%
    mutate(Category = factor(Category, levels = sub_filter)) %>%
    arrange(Category) %>%
    pull(Phen) %>%
    setdiff("INS") %>%
    unique() 
  
  phen_targets <- Phen.Seq.BL %>% 
    dplyr::select(Clinic_ID) %>%  
    left_join(sub_data_df, by = "Clinic_ID") %>%
    na.omit()
 
  cat(sprintf("     Path Set: %s\n", base_out_dir))
  cat(sprintf("     Resolved Targets: N = %d subjects\n", nrow(phen_targets)))
  cat(sprintf("     Injected Clinical Features: N = %d vars [%s...]\n", length(clinic.KH.vars), clinic.KH.vars[1]))
  print(table(phen_targets$group))
 
  threshold_map <- list("GUT" = 0.001)
  prevalence_cutoff <- 0.1 
 
  global_summary_list <- list()
  
  # ==============================================================================
  # PHASE 2 & 3 
  # ==============================================================================
  for(trans_method in TRANSFORM_METHODS) {
    cat(paste0("\n", rep("=", 80), collapse = ""), "\n")
    
    # Precompute abundance transformations per anatomical site (prevalence > 10%)
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
    
    # Merge profiles across sites while preserving native NA values
    merged_prof_BL_highpreval <- Reduce(function(x, y) merge(x, y, by = "Clinic_ID", all = TRUE), prof_BL_highpreval)
    
    taxa.species.vars_master <- setdiff(colnames(merged_prof_BL_highpreval), "Clinic_ID")
 
    for(fs_method in FEATURE_SELECTION_METHODS) {
      cat(sprintf("  -> 🧪 [SELECTION] Method: %s\n", fs_method))
      
      for(site_name in names(site_definitions)) {
        current_site_codes <- site_definitions[[site_name]]
        
        exclude_sites <- setdiff(c("GUT", "VA", "TO", "UR"), current_site_codes)
        regex_pattern <- if(length(exclude_sites) > 0) paste0("^", exclude_sites, "_", collapse = "|") else "^X_X_X"
        
        taxa_site_cands <- taxa.species.vars_master[!grepl(regex_pattern, taxa.species.vars_master)]
        prof_BL_input_site <- merged_prof_BL_highpreval[, !grepl(regex_pattern, colnames(merged_prof_BL_highpreval))]
        
        # Exclude samples lacking measurements for the target anatomical site
        valid_sample_idx <- rowSums(!is.na(prof_BL_input_site %>% select(all_of(taxa_site_cands)))) > 0
        prof_BL_input_site <- prof_BL_input_site[valid_sample_idx, , drop = FALSE]
        
        for(current_target in target_list) {
          
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
          
          # clinical_features_master
          basic.vars <- NULL 
          
          if (grepl("K", current_target, ignore.case = TRUE)) {
            basic.vars <- c(basic.vars, "K_Score","SAS","SDS","FSH","LH","Sleep_Score","Palpitations","Sexual_dysfunction") %>% intersect(clinic.KH.vars)
            cat(c("     [Info] ",paste0(basic.vars,collapse = " "),"纳入 basic.vars\n"))
          } 
          
          clinical_features_master <- unique(c(clinic.KH.vars, basic.vars))
          
          # Cohort partitioning: 70/30 stratified split or geographic multicenter split
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
          cat("[Data Distribution] Train:\n")
          print(train_dist)
          cat("[Data Distribution] Test:\n")
          print(test_dist)
          write.csv(distribution_summary, 
                    paste0("0_Response_Distribution_", current_target, ".csv"), 
                    row.names = FALSE)
          
          # Feature pre-filtering using Wilcoxon rank-sum test on the training set
          taxa_selected <- taxa_site_cands  
          
          if(fs_method == "Wilcoxon") {
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
          
          # Construct feature sets: Clinical-only, Microbe-only, and Combined models
          RF.data.list <- list(
            "Features_Clinical_Only" = master_data %>% select(Response,  all_of(clinical_features_master)),
            "Features_Microbe_only" = master_data %>% select(Response, all_of(taxa_selected)),
            "Features_Microbe_Clinical" = master_data %>% select(Response,all_of(clinical_features_master), all_of(taxa_selected))
          )
          
          results_rf <- list()
           
          cl <- makeCluster(max(1, min(10,detectCores() - 2)))
          registerDoParallel(cl)
          
          for (model_name in names(RF.data.list)) {
            cat(sprintf("      -> Training RF: %s\n", model_name))
            
            data_sub <- RF.data.list[[model_name]]
            train_data  <- data_sub[train_index, ]
            test_data   <- data_sub[-train_index, ]
 
            if(!require(smotefamily)) install.packages("smotefamily")
            if(!require(pROC)) install.packages("pROC")
            library(pROC)
            library(randomForest)
            
            # Configure custom RFE prediction functions to retain class probabilities
            custom_rfFuncs <- rfFuncs
            custom_rfFuncs$pred <- function(modelFit, testX) {
              res <- predict(modelFit, testX)
              prob <- predict(modelFit, testX, type = "prob")
              
              out <- data.frame(pred = res)
              if(!is.null(prob)) out <- cbind(out, prob)
              out
            }
            
            custom_rfFuncs$summary <- twoClassSummary
            total_feats <- ncol(train_data) - 1
            
            if (total_feats > 80) {
              sizes_to_test <- unique(c(1:80, seq(85, min(total_feats, 100), by = 5), total_feats))
            } else {
              sizes_to_test <- 1:total_feats
            }
            
            # Generate deterministic random seeds for 5-fold CV RFE
            set.seed(224)  
            
            rfe_seeds <- vector(mode = "list", length = 6) 
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
            
            # Configure RFE control settings
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
            
            # Execute 5-fold CV recursive feature elimination based on ROC metric
            cat(sprintf("      -> 5-Fold CV features counts: [%s]\n", 
                        paste(sizes_to_test, collapse = ", ")))
             
            X_train <- train_data %>% select(-Response)
            Y_train <- train_data$Response
            
            rf_rfe_profile <- rfe(x = X_train, 
                                  y = Y_train,
                                  sizes = sizes_to_test, 
                                  rfeControl = rfe_ctrl,
                                  metric = "ROC",             
                                  ntree = 500)
 
            best_feature_count <- rf_rfe_profile$bestSubset
            taxa_selected_cv <- predictors(rf_rfe_profile)
            write_rds(rf_rfe_profile, paste0("1_RFECV_Selected_Features_", model_name, ".RDS"))
            
            cat(sprintf("         [📊 对比] 绝对最高 ROC 对应特征量：%d 个\n", best_feature_count))
 
            rfe_results <- rf_rfe_profile$results
            best_perf <- rfe_results %>% filter(Variables == best_feature_count)
            
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
            
            # Plot cross-validation ROC curves across individual folds
            library(pROC)
            library(ggplot2)
            library(dplyr)
            library(tidyr)
             
            cv_predictions <- rf_rfe_profile$pred %>% 
              filter(Variables == best_feature_count) 
            
            total_roc <- roc(cv_predictions$obs, cv_predictions$R, levels = c("NR", "R"), direction = "<", quiet = TRUE)
            total_auc <- as.numeric(auc(total_roc))
            
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
            
            # Retrain final model with optimal feature subset and tune mtry parameter
            train_data_selected <- train_data %>% select(Response, all_of(taxa_selected_cv))
            test_data_selected  <- test_data %>% select(Response, all_of(taxa_selected_cv))

            if (SAMPLING_METHOD == "none") {
              ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, allowParallel = TRUE, 
                                   summaryFunction = twoClassSummary)
            } else {
              ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, allowParallel = TRUE, 
                                   summaryFunction = twoClassSummary, sampling = SAMPLING_METHOD)
            }
            
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
            
            final_model <- rf_caret$finalModel
            write_rds(final_model, paste0("2_Train_Final_Model_", model_name, ".RDS"))
            
            # Evaluate model performance on the independent test set
            pred_prob <- predict(rf_caret, test_data_selected, type = "prob")[, "R"]
            roc_obj <- roc(test_data_selected$Response, pred_prob, levels = c("NR", "R"), direction = "<", quiet = TRUE)
            
            pred_class <- predict(rf_caret, test_data_selected)
            cm <- confusionMatrix(data = pred_class, reference = test_data_selected$Response, positive = "R")
        
            imp_data <- importance(final_model) %>% 
              as.data.frame() %>% 
              rownames_to_column("Feature") %>% 
              arrange(desc(MeanDecreaseAccuracy))

            results_rf[[model_name]] <- list(
              AUC = as.numeric(auc(roc_obj)), 
              Accuracy = as.numeric(cm$overall["Accuracy"]),
              Sensitivity = as.numeric(cm$byClass["Sensitivity"]), 
              Specificity = as.numeric(cm$byClass["Specificity"]),
              Importance = imp_data
            )
            cat(sprintf("         [ROC] TEST AUC: %.4f\n", results_rf[[model_name]]$AUC))

            saveRDS(rf_caret, paste0("Model_Caret_Full_", model_name, ".RDS"))
            write.csv(imp_data, paste0("SourceData_Feature_Importance_", model_name, ".csv"), row.names = FALSE, quote = FALSE)

            test_predictions <- data.frame(
              Actual_Response = test_data_selected$Response,
              Pred_Prob_R = pred_prob,
              Pred_Class = pred_class
            )
            write.csv(test_predictions, paste0("SourceData_Test_Predictions_", model_name, ".csv"), row.names = TRUE, quote = FALSE)

            roc_source_data <- data.frame(
              Threshold = roc_obj$thresholds,
              Sensitivity = roc_obj$sensitivities,     
              Specificity = roc_obj$specificities,
              `1_Minus_Specificity` = 1 - roc_obj$specificities  
            )
            write.csv(roc_source_data, paste0("SourceData_ROC_Curve_", model_name, ".csv"), row.names = FALSE, quote = FALSE)
            
            # metrics
            metrics_df <- data.frame(
              Model = model_name,
              AUC = as.numeric(auc(roc_obj)), 
              Accuracy = as.numeric(cm$overall["Accuracy"]),
              Sensitivity = as.numeric(cm$byClass["Sensitivity"]), 
              Specificity = as.numeric(cm$byClass["Specificity"])
            )
            write.csv(metrics_df, paste0("Table_Metrics_", model_name, ".csv"), row.names = FALSE, quote = FALSE)

            # ROC & PDF
            pdf(file = paste0("3_Test_ROC_Curve_", model_name, ".pdf"), width = 6, height = 6)
            plot(roc_obj, 
                 main = paste("Test ROC -", model_name),
                 col = "#e74c3c", lwd = 3, legacy.axes = TRUE, print.auc = TRUE,
                 auc.polygon = TRUE, auc.polygon.col = "#fdedec")
            dev.off()
            
            # Importance
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
  # PHASE 4 OUTPUT
  # ==============================================================================
  setwd(base_out_dir)
  
  final_master_summary <- bind_rows(global_summary_list) %>% arrange(desc(AUC))
  
  write.csv(final_master_summary, "0_Global_All_Models_AUC_Summary.csv", row.names = FALSE)
  
  colors_line <- c(
    "Features_Clinical_Only" = "grey20",
    "Features_Microbe_only"          = "#2ca02c", 
    "Features_Microbe_Clinical"      = "#ff7f0e",
    "Features_Microbe_BL"            = "#0077B6"
  )
  final_master_summary$Model_Type <- factor(final_master_summary$Model_Type, levels = names(colors_line))
  # Grid: Trans + FS ~ Site
  p_gods_eye <- ggplot(final_master_summary, aes(x = Model_Type, y = AUC, fill = Model_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "white") +
    geom_text(aes(label = sprintf("%.3f", AUC)), position = position_dodge(width = 0.8), vjust = -0.5, size = 2.8, fontface = "bold") +
    geom_text(aes(y = 0.35, label = paste0("F.N=", Total_Features)), position = position_dodge(width = 0.8), size = 2, color = "white", fontface="bold") +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey60", linewidth = 0.5, alpha = 0.5) +
    coord_cartesian(ylim = c(0.3, max(final_master_summary$AUC, na.rm=TRUE) + 0.1)) +
    scale_fill_manual(values = colors_line) +
    labs(title = "Global Pipeline Performance Matrix", x = "", y = "Test Set AUC", fill = "") +
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
  
}

