library(dplyr)
library(tidyverse)
library(lme4)
library(mediation)
library(ggeffects)
library(ggplot2)
library(compositions)
library(conflicted)

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))
Phen.Seq.BL <- Phen.Seq %>% subset(Time == "BL")

target_dir <- file.path(base_dir, "Results/Mediation")
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(mediation::mediate)

# Prevalence filtering (>10% samples) and CLR transformation for vaginal taxa
rel_ab.threshold <- 0.0001
data_va <- prof_filtered$VA %>% data.frame()

taxa_keep_va <- sapply(data_va, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(data_va))
taxa.vars_va <- names(taxa_keep_va)[taxa_keep_va]
taxa_matrix  <- data_va[, taxa.vars_va, drop = FALSE]

pseu_value   <- min(taxa_matrix[taxa_matrix > 0], na.rm = TRUE) * 0.1
data_imputed <- taxa_matrix + pseu_value
data_trans   <- as.data.frame(unclass(compositions::clr(data_imputed)))
dat.VA       <- data_trans %>% rownames_to_column("SeqID")

microbe_list <- c("Lactobacillus_crispatus",  "Lactobacillus_iners", 
                  "Gardnerella_vaginalis",    "Fannyhessea_vaginae",
                  "Aerococcus_christensenii", "Dialister_micraerophilus")

mediation_forest_table <- data.frame()

for (target_microbe in microbe_list) {
  if (!target_microbe %in% colnames(dat.VA)) next
  
  long_data <- dat.VA %>% 
    select(SeqID, Target_Microbe = all_of(target_microbe)) %>%
    left_join(Microbe.phen.prof %>% select(SeqID, Clinic_ID, Time), by = "SeqID") %>%
    left_join(Phen.Seq.BL %>% select(Clinic_ID, Age, BMI), by = "Clinic_ID") %>%
    left_join(Phen.Seq %>% select(Clinic_ID, Time, FSFI), by = c("Clinic_ID", "Time")) %>%
    filter(Time %in% c("BL", "T04")) %>%
    na.omit()
  
  long_data$Time <- factor(long_data$Time, levels = c("BL", "T04"))
  
  set.seed(123)
  # Longitudinal mediation models controlling for Age, BMI, and repeated measures per participant
  model.M_mixed <- lmer(Target_Microbe ~ Time + Age + BMI + (1 | Clinic_ID), data = long_data)
  model.Y_mixed <- lmer(FSFI ~ Time + Target_Microbe + Age + BMI + (1 | Clinic_ID), data = long_data)
  
  results_mixed <- mediate(model.M_mixed, model.Y_mixed, 
                           treat = "Time", mediator = "Target_Microbe", sims = 1000)
  
  med_summary <- summary(results_mixed)
  
  current_res <- data.frame(
    Microbe      = target_microbe,
    Estimate     = med_summary$d0,
    CI_Lower     = med_summary$d0.ci[1],
    CI_Upper     = med_summary$d0.ci[2],
    P_Value      = med_summary$d0.p,
    Prop_Med     = med_summary$n0,
    Prop_P       = med_summary$n0.p,
    ADE.Estimate = med_summary$z0,
    ADE.CI_Lower = med_summary$z0.ci[1],
    ADE.CI_Upper = med_summary$z0.ci[2],
    ADE.P_Value  = med_summary$z0.p
  )
  
  mediation_forest_table <- rbind(mediation_forest_table, current_res)
  
  # Path b partial residuals adjusted for covariates (Age, BMI)
  eff_Path_b     <- ggpredict(model.Y_mixed, terms = "Target_Microbe")
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
    geom_line(data = plot_data_line, aes(x = x, y = predicted), color = "#1a5276", linewidth = 1.2) +
    scale_color_viridis_c(option = "D", name = "Adjusted FSFI") +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    labs(
      title = paste("Path b Partial Residual Plot:", target_microbe),
      x = paste("CLR Abundance of", target_microbe), 
      y = "Partial Residuals of FSFI"
    )
  
  ggsave(filename = paste0("Partial_Residual_Plot_", target_microbe, ".pdf"), plot = p_partial, width = 7, height = 6)
}

mediation_forest_table <- mediation_forest_table %>%
  mutate(
    # FDR         = p.adjust(P_Value, method = "BH"),
    # ADE.FDR     = p.adjust(ADE.P_Value, method = "BH"),
    # Prop_P.FDR  = p.adjust(Prop_P, method = "BH"),
    P_Label     = ifelse(P_Value < 0.001, "P < 0.001", sprintf("P = %.3f", P_Value)),
    ADE.P_Label = ifelse(ADE.P_Value < 0.001, "ADE.P < 0.001", sprintf("ADE.P = %.3f", ADE.P_Value)),
    Microbe     = factor(Microbe, levels = rev(microbe_list)),
    Effect_Type = case_when(
      Estimate > 0 & P_Value < 0.05 ~ "Positive Mediation",
      Estimate < 0 & P_Value < 0.05 ~ "Negative Mediation",
      TRUE                          ~ "Not Significant"
    )
  )

write.csv(mediation_forest_table, "Mediation_Results_Total.csv", row.names = FALSE)

p_forest <- ggplot(mediation_forest_table, aes(x = Estimate, y = Microbe, color = Effect_Type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray20", linewidth = 0.8) +
  geom_pointrange(aes(xmin = CI_Lower, xmax = CI_Upper), linewidth = 0.8, fatten = 4) +
  geom_text(aes(label = P_Label), hjust = -0.1, vjust = -0.5, color = "black", fontface = "italic", size = 2.5) +
  scale_color_manual(values = c("Positive Mediation" = "#1a5276", 
                                "Negative Mediation" = "#b03a2e", 
                                "Not Significant"    = "#7f8c8d")) + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(face = "bold", size = 12),
    axis.text          = element_text(size = 11, color = "black"),
    axis.text.y        = element_text(face = "italic"),                    
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position    = "right"
  ) +
  labs(
    title = "Mediation Effect of Vaginal Microbes",
    x = "Indirect Mediation Effect Size (ACME)",
    y = NULL,
    color = "Statistical Inference"
  ) +
  xlim(min(mediation_forest_table$CI_Lower) - 0.05, max(mediation_forest_table$CI_Upper) + 0.2)

ggsave(filename = "FOREST_ACME_FSFI_Vaginal.pdf", plot = p_forest, width = 7, height = 6)

