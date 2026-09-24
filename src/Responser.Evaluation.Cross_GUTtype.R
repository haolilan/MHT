library(tidyverse)
library(broom)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

subgroups_config <- list(
  "Hormone_E2" = phen.gr.E2.all %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "mKI"        = phen.gr.K      %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR")))
)

time_points_vec <- c("BL", "T04", "T12", "T24")

format_p <- function(p) {
  ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p))
}

df_cst_pre <- Microbe.phen.prof %>%
  left_join(prof_diversity %>% select(SeqID, CST), by = "SeqID") %>%
  dplyr::select(Clinic_ID, Time, Site, CST)

# ==============================================================================
# Longitudinal Gut CST Logistic Regression (Ref: GUT-P.cop)
# ==============================================================================
for (sg_name in names(subgroups_config)) {
  
  target_dir <- file.path(base_dir, "Results/Responser", sg_name, "Logistic_GUTtype")
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  
  current_gr_df <- subgroups_config[[sg_name]]
  target_levels <- levels(current_gr_df$group)
  
  df_gut <- df_cst_pre %>%
    filter(Site == "GUT", Time %in% time_points_vec, CST %in% c("GUT-P.cop", "GUT-P.v")) %>%
    inner_join(current_gr_df, by = "Clinic_ID") %>%
    filter(!is.na(group), !is.na(CST)) %>%
    mutate(
      Outcome_Binary = if_else(group == target_levels[1], 1, 0),
      CST = factor(CST, levels = c("GUT-P.cop", "GUT-P.v"))
    )
  
  logistic_results <- list()
  
  for (tp in time_points_vec) {
    slice_df <- df_gut %>% filter(Time == tp)
    if (nrow(slice_df) == 0) next
 
    if (n_distinct(slice_df$CST) < 2 || n_distinct(slice_df$Outcome_Binary) < 2) next
    
    n_total_outcome_1 <- sum(slice_df$Outcome_Binary == 1)
    n_total_outcome_0 <- sum(slice_df$Outcome_Binary == 0)
     
    TP <- sum(slice_df$CST == "GUT-P.v" & slice_df$Outcome_Binary == 1)
    FP <- sum(slice_df$CST == "GUT-P.v" & slice_df$Outcome_Binary == 0)
    
    if (TP == 0 || FP == 0) next
    
    FN <- n_total_outcome_1 - TP
    TN <- n_total_outcome_0 - FP
 
    fit <- tryCatch(
      glm(Outcome_Binary ~ CST, data = slice_df, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    model_stats <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "CSTGUT-P.v")
    
    if (nrow(model_stats) == 0 || is.na(model_stats$estimate)) next
    
    logistic_results[[tp]] <- tibble(
      Subgroup           = sg_name,
      Site               = "GUT",
      TimePoint          = tp,
      Comparison         = "GUT-P.v (vs. GUT-P.cop)",
      N_Target_Outcome_0 = FP,
      N_Target_Outcome_1 = TP,
      N_Total_Outcome_0  = n_total_outcome_0,
      N_Total_Outcome_1  = n_total_outcome_1,
      Sensitivity        = round(TP / (TP + FN), 3),
      Specificity        = round(TN / (TN + FP), 3),
      PPV                = round(TP / (TP + FP), 3),
      NPV                = round(TN / (TN + FN), 3),
      Odds_Ratio         = model_stats$estimate,
      CI_Lower           = model_stats$conf.low,
      CI_Upper           = model_stats$conf.high,
      p_value            = model_stats$p.value
    )
  }
  
  if (length(logistic_results) == 0) next
  
  res_table <- bind_rows(logistic_results) %>%
    mutate(
      TimePoint = factor(TimePoint, levels = time_points_vec),
      p_fmt     = format_p(p_value),
      p_label   = paste0("P: ", p_fmt)
    ) %>%
    arrange(TimePoint)
  
  write.csv(res_table, file.path(target_dir, sprintf("%s_GUT_Logistic_Results.csv", sg_name)), row.names = FALSE)
 
  p_forest <- ggplot(res_table, aes(x = Odds_Ratio, y = Comparison)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.18, color = "#2B8CBE", linewidth = 0.7) +
    geom_point(color = "#E41A1C", size = 2.8) +
    geom_text(aes(label = p_label), size = 3.0, vjust = -1.0, hjust = 0.5) +
    facet_grid(. ~ TimePoint, scales = "free_y") +
    scale_x_log10() +
    theme_bw() +
    labs(
      x = "Odds Ratio (Log Scale, 95% CI)",
      y = NULL,
      title = sprintf("Gut CST (P.v vs. P.cop) vs. Treatment Response Across Timeline [%s]", sg_name)
    ) +
    theme(
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text       = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(color = "black", size = 9.5, face = "bold"),
      axis.text.x      = element_text(color = "black", size = 8.5)
    )
  
  ggsave(file.path(target_dir, sprintf("%s_GUT_Logistic_Forest.pdf", sg_name)), plot = p_forest, width = 9.5, height = 2.6)
}