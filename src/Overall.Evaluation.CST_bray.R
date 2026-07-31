# ==============================================================================
# Pipeline: Menopause Hormone Therapy (MHT) Multi-Omic Evaluation
# Modules: CST Longitudinal Ratio, Logistic Regression & Longitudinal Bray Dynamics
# Year: 2026
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Global Setup & Library Dependencies
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
library(conflicted)
library(broom)
library(openxlsx)
library(vegan)
library(colorspace)
library(RColorBrewer)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
out_dir  <- "D:/WorkProjects/Demo-MHT 2026/Results/OverallEva/CST_Bray"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

# Load global dataset workspace image
load(file.path(base_dir, "data/MHT.demo.RData"))

# Define sample inclusion strategy: "strict_4" or "None"
# filter_strategy <- "None"
filter_strategy <- "strict_4"

# Core Configuration Arrays
all_sites        <- c("VA", "UR", "GUT", "TO")
follow_ups       <- c("T04", "T12", "T24")
time_order_bray  <- c("BL", "T4", "T12", "T24")
cst_ordered_factor <- c("GUT-P.cop", "GUT-P.v", "TO-N.s", "TO-P.m", "UROG-Div", "UROG-G.v", "UROG-L.i", "UROG-L.c")

# Base metadata harmonization
mg_meta_outcome <- Microbe.phen.prof %>%
  left_join(prof_diversity %>% select(SeqID, CST), by = "SeqID")

df_cst_base <- mg_meta_outcome %>%
  dplyr::select(Clinic_ID, Time, Site, CST) %>%
  drop_na(Time, Site, CST)

# ==============================================================================
# 1. Automated Longitudinal Composition (Fisher's Exact Test & Logistic) #################
# ==============================================================================
longitudinal_stats     <- list()
logistic_results_list  <- list()
valid_rows_list        <- list()

for (st in all_sites) {
  for (fu in follow_ups) {
    
    # 1.1 Pairwise longitudinal cleaning & alignment
    comp_df <- df_cst_base %>%
      dplyr::filter(Site == st, Time %in% c("BL", fu))
    
    paired_ids <- comp_df %>%
      group_by(Clinic_ID) %>%
      summarise(has_bl = "BL" %in% Time, has_fu = fu %in% Time, .groups = "drop") %>%
      dplyr::filter(has_bl & has_fu) %>%
      pull(Clinic_ID)
    
    if (filter_strategy == "strict_4") {
      strict_4_ids <- df_cst_base %>%
        dplyr::filter(Site == st) %>%
        group_by(Clinic_ID) %>%
        summarise(n_time = n_distinct(Time), .groups = "drop") %>%
        dplyr::filter(n_time == 4) %>% 
        pull(Clinic_ID)
      paired_ids <- intersect(paired_ids, strict_4_ids)
    }
    
    current_n <- length(paired_ids)
    if (current_n < 3) next 
    
    pair_df <- comp_df %>% dplyr::filter(Clinic_ID %in% paired_ids)
    valid_rows_list[[paste(st, fu, sep = "_")]] <- pair_df
    
    # 1.2 Contingency Table & Fisher's Exact Test
    tab <- table(pair_df$CST, pair_df$Time[, drop = TRUE])
    tab <- tab[, c("BL", fu), drop = FALSE] 
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    
    if (nrow(tab) >= 2) {
      test_res <- fisher.test(tab, simulate.p.value = TRUE, B = 2000)
      longitudinal_stats[[paste(st, fu, sep = "_")]] <- data.frame(
        Site = st, Comparison = paste0(fu, "_vs_BL"), p = test_res$p.value, N_Samples = current_n, stringsAsFactors = FALSE
      )
    }

    # # 1.3 Conditional Logistic Regression Engine (clogit Paired Design)
    # ==============================================================================
    library(survival)
    
    # 构建条件回归基盘：时间转换为干预二分类（BL=0, 随访点=1）
    slice_df    <- pair_df %>% mutate(Time_Binary = if_else(Time == fu, 1, 0))
    unique_csts <- unique(slice_df$CST)
    if (length(unique_csts) < 2) next
    
    for (target_cst in unique_csts) {
      run_df <- slice_df %>% mutate(CST_Binary = if_else(CST == target_cst, 1, 0))
      contingency_table <- table(run_df$CST_Binary, run_df$Time_Binary)
      if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) next
      tryCatch({
        model <- clogit(CST_Binary ~ Time_Binary + strata(Clinic_ID), data = run_df)
        model_stats <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% 
          dplyr::filter(term == "Time_Binary")
        if (nrow(model_stats) == 0 || is.na(model_stats$estimate)) next
        logistic_results_list[[paste(st, fu, target_cst, sep = "_")]] <- data.frame(
          Site          = st, 
          FollowUpPoint = fu, 
          Target_CST    = target_cst, 
          Odds_Ratio    = model_stats$estimate,   
          CI_Lower      = model_stats$conf.low, 
          CI_Upper      = model_stats$conf.high, 
          p_value       = model_stats$p.value, 
          N_Patients    = current_n, 
          stringsAsFactors = FALSE
        )
        
      }, error = function(e) {
        cat(sprintf("      [INFO] CST [%s] 在 [%s_vs_BL] 条件下无足够变异，已跳过。\n", target_cst, fu))
      })
    }
    
  }
}

# ==============================================================================
# 2. Composition Reports & Visualization (Stacked Bars & Forest Plots) ####################
# ==============================================================================
cst_colors <- c(
  "GUT-P.v"   = "#66c2a4", "GUT-P.cop" = "#b2e2e2",
  "TO-P.m"    = "#fc8d59", "TO-N.s"   = "#fdcc8a",
  "UROG-L.c"  = "#117733", "UROG-L.i" = "#74C476", 
  "UROG-G.v"  = "#FF7F0E", "UROG-Div" = "#2171B5"
)
# 2.1 Format Fisher Statistics Outputs
cst_stats_long <- bind_rows(longitudinal_stats) %>%
  mutate(
    p_txt = if_else(p < 0.0001, sprintf("p=%.2e", p), sprintf("p=%.3f", p)),
    p_label = sprintf("%s\n(n=%d)", p_txt, N_Samples), y_pos = 1.02,
    Time = factor(gsub("_vs_BL", "", Comparison), levels = c("BL", "T04", "T12", "T24"))
  )

df_cst <- bind_rows(valid_rows_list) %>%
  distinct() %>%
  mutate(
    Time = factor(Time, levels = c("BL", "T04", "T12", "T24")),
    Site = factor(Site, levels = all_sites),
    CST  = factor(CST, levels = cst_ordered_factor)
  )

print("Longitudinal composition differences and sample sizes:")
print(cst_stats_long %>% dplyr::select(Site, Comparison, p_txt, N_Samples))

# 2.2 Stacked Proportion Plot
df_plot_cst <- df_cst %>%
  group_by(Site, Time, CST) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(Percentage = count / sum(count)) %>%
  ungroup()

p_cst_long_ratio <- ggplot(df_plot_cst, aes(x = Time, y = Percentage, fill = CST)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
  
  facet_wrap(~ Site, nrow = 1, scales = "free_y") + 
  scale_fill_manual(values = cst_colors, name = "CST Type", drop = TRUE) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = " ", 
    y = "Relative Proportion of CST", 
    title = paste0("Longitudinal Shifts of CST Proportion (Strategy: ", filter_strategy, ")")
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(size = 11, face = "bold"),
    axis.text.x      = element_text(size = 10, face = "bold"),
    axis.text.y      = element_text(size = 9, color = "black"), 
    axis.ticks.y     = element_line(color = "black"), # 补齐右侧子图的 Y 轴小刻度线
    
    panel.grid       = element_blank()
  ) +
  geom_text(
    data = cst_stats_long, 
    aes(x = Time, y = y_pos, label = p_label), 
    inherit.aes = FALSE, vjust = 0.5, size = 2.5, fontface = "bold", color = "black", lineheight = 0.8
  )
ggsave(sprintf("CST_Longitudinal_Proportions_%s.pdf", filter_strategy), plot = p_cst_long_ratio, width = 11, height = 5)

# 2.3 Format & Save Logistic Regression Reports
get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"   
  )
}
all_logistic_report <- bind_rows(logistic_results_list) %>%
  group_by(Site) %>%
  mutate(p.adj = p.adjust(p_value, method = "BH"))%>%
  ungroup() %>%
  mutate(
    p_star = get_significance_stars(p_value),
    q_star = get_significance_stars(p.adj),
    raw_mark = paste(p_star, q_star, sep = "/"),
    p_mark = if_else(p_value > 0.05, "", raw_mark)
  ) %>%
  dplyr::select(-p_star, -q_star,-raw_mark) %>%
  arrange(p_value)

write.csv(all_logistic_report, sprintf("CST_Longitudinal_Logistic_%s.csv", filter_strategy), row.names = FALSE)

# 2.4 Forest Plot Visualization
plot_forest_df <- all_logistic_report %>% 
  mutate(
    Plot_Label = factor(Target_CST, levels = cst_ordered_factor),
    FollowUpPoint = factor(FollowUpPoint, levels = follow_ups),
    Site = factor(Site, levels = all_sites)
  )

if (nrow(plot_forest_df) > 0) {
  p_forest <- ggplot(plot_forest_df, aes(x = Odds_Ratio, y = Plot_Label, color = Target_CST)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, linewidth = 0.7) +
    geom_point(size = 2.5) +
    geom_text(aes(label = p_mark), size = 2.5, fontface = "italic", vjust = -0.6, hjust = -0.05) +
    
    facet_grid(Site ~ FollowUpPoint, scales = "free_y", space = "free_y", axes = "all") +
    
    scale_x_log10() + 
    scale_color_manual(values = cst_colors, name = "CST Type", drop = TRUE) +
    labs(
      x = "Odds Ratio vs Baseline (Log Scale, 95% CI)", 
      y = "CST Types", 
      title = paste0("Longitudinal Logistic Odds Ratio of CST Types (Strategy: ", filter_strategy, ")")
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text       = element_text(size = 11, face = "bold"),
      axis.text.x      = element_text(color = "black", size = 9),
      axis.text.y      = element_text(color = "black", size = 9),
      axis.ticks       = element_line(color = "black", linewidth = 0.4), 
      
      panel.grid.minor = element_blank(), 
      legend.position  = "bottom"
    )
  ggsave(sprintf("CST_Longitudinal_Logistic_Forest_%s.pdf", filter_strategy), plot = p_forest, width = 11, height = 8)
}

# ==============================================================================
# 3. CST Transition Sankey Flow Graphs ##############################################
# ==============================================================================
message(">>> Starting Batch CST Transition Tracking Plots...")

for (site_name in names(prof_filtered)) {
  phen.input_raw <- prof_diversity %>%
    filter(Site == site_name & Time %in% c("BL", "T04", "T12", "T24")) %>%
    mutate(Time = recode(Time, "T04" = "T4")) %>%
    select(Clinic_ID, Time, CST) %>%
    drop_na()
  
  if (nrow(phen.input_raw) == 0) next
  
  site_present_csts <- cst_ordered_factor[cst_ordered_factor %in% unique(phen.input_raw$CST)]
  if (length(site_present_csts) == 0) next
  
  phen.input <- phen.input_raw %>%
    filter(CST %in% site_present_csts) %>%
    mutate(CST = factor(CST, levels = site_present_csts))
  
  data.point <- phen.input %>% count(Time, CST, name = "Freq") %>% mutate(Time = factor(Time, levels = time_order_bray))
  
  data.seg <- phen.input %>%
    inner_join(phen.input, by = "Clinic_ID", suffix = c("", ".V2"),, relationship = "many-to-many") %>%
    filter((Time == "BL" & Time.V2 == "T4") | (Time == "T4" & Time.V2 == "T12") | (Time == "T12" & Time.V2 == "T24")) %>%
    count(CST, Time, CST.V2, Time.V2, name = "count") %>%
    mutate(
      Time = factor(Time, levels = time_order_bray), Time.V2 = factor(Time.V2, levels = time_order_bray),
      CST = factor(CST, levels = site_present_csts), CST.V2 = factor(CST.V2, levels = site_present_csts)
    )
  
  total_labels <- data.point %>% group_by(Time) %>% summarise(SamNum = sum(Freq), .groups = "drop")
  
  p_trans <- ggplot() +
    geom_segment(data = data.seg, aes(x = Time, y = CST, xend = Time.V2, yend = CST.V2, size = count, color = count), alpha = 0.8, lineend = "round") +
    geom_point(data = data.point, aes(x = Time, y = CST, size = Freq), color = "#FFD700", alpha = 0.95) +
    geom_text(data = data.point, aes(x = Time, y = CST, label = Freq), color = "black", size = 3.5, fontface = "bold") +
    geom_text(data = total_labels, aes(x = Time, y = site_present_csts[1], label = paste0("N=", SamNum)), vjust = 4.5, size = 4, color = "grey15", inherit.aes = FALSE) +
    scale_x_discrete(limits = time_order_bray, expand = c(0.12, 0.12)) +
    scale_y_discrete(limits = site_present_csts) + 
    scale_color_gradient(low = "#e0ecf4", high = "purple4", name = "Flow Count") +
    scale_size_area(max_size = 14, guide = "none") + 
    labs(title = paste0("CST Longitudinal Transitions - ", site_name), x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      panel.grid.major.y = element_line(color = "grey93", linetype = "dotted"),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      legend.position = "right", legend.title = element_text(size = 10, face = "bold")
    )
  
  ggsave(sprintf("%s_CST_Longitudinal_Transition.pdf", site_name), plot = p_trans, width = 6.5, height = min(5.2, 3.5 + length(site_present_csts)*0.4))
  write.csv(data.seg, sprintf("%s_CST_Transition_Matrix.csv", site_name), row.names = FALSE)
}

# ==============================================================================
# 4. Individual Longitudinal Bray-Curtis Distance Dispersal #############################
# ==============================================================================
message(">>> Commencing Individual Bray-Curtis Tracking & Stat Plots...")

filter_strategy <- "strict_4" # 可选："strict_4" 或 "none"

cols               <- c("BL/T4" = "#D9534F", "BL/T12" = "#F0AD4E", "BL/T24" = "#337AB7", "T4/T12" = "#5CB85C", "T4/T24" = "#9467BD", "T12/T24" = "#8C564B")
time_border_colors <- colorspace::darken(cols, 0.3)
site_bray_list     <- list()
site_stat_list     <- list()

for (site_name in names(prof_filtered)) {
  mat_raw <- data.frame(prof_filtered[[site_name]])
  
  mat_with_meta <- mat_raw %>%
    mutate(
      Group     = Microbe.phen.prof$Time[match(rownames(mat_raw), Microbe.phen.prof$SeqID)],
      Clinic_ID = Microbe.phen.prof$Clinic_ID[match(rownames(mat_raw), Microbe.phen.prof$SeqID)]
    ) %>%
    mutate(Group = recode(Group, "T04" = "T4")) %>%
    filter(Group %in% time_order_bray) %>%
    select(Clinic_ID, Group, everything())
  
  if (filter_strategy == "strict_4") {
    strict_4_ids <- mat_with_meta %>%
      group_by(Clinic_ID) %>%
      summarise(n_time = n_distinct(Group), .groups = "drop") %>%
      filter(n_time == 4) %>% 
      pull(Clinic_ID)
    
    mat_with_meta <- mat_with_meta %>% filter(Clinic_ID %in% strict_4_ids)
  }
  
  clinic_list <- unique(mat_with_meta$Clinic_ID)
  person_bray_list <- list()
  
  for (id in clinic_list) {
    person_sub <- mat_with_meta %>% filter(Clinic_ID == id)
    if (nrow(person_sub) < 2) next 
    
    abundance_matrix <- person_sub %>% select(-Clinic_ID, -Group) %>% as.matrix()
    rownames(abundance_matrix) <- person_sub$Group
    
    bc_matrix <- vegdist(abundance_matrix, method = "bray")
    bc_table  <- as.data.frame(as.table(as.matrix(bc_matrix)))
    colnames(bc_table) <- c("time1", "time2", "bray_distance")
    
    bc_table_clean <- bc_table %>%
      filter(time1 %in% time_order_bray & time2 %in% time_order_bray) %>%
      filter(match(time1, time_order_bray) < match(time2, time_order_bray))
    
    if (nrow(bc_table_clean) == 0) next
    bc_table_clean$Clinic_ID <- id
    person_bray_list[[id]] <- bc_table_clean
  }
  
  if(length(person_bray_list) == 0) next
  
  site_bray_df <- bind_rows(person_bray_list) %>%
    mutate(time_pair = factor(paste(time1, time2, sep = "/"), levels = names(cols)), Site = site_name)
  
  stat_test <- site_bray_df %>%
    wilcox_test(bray_distance ~ time_pair) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(Site = site_name, y.position = seq(max(site_bray_df$bray_distance) * 1.05, max(site_bray_df$bray_distance) * 1.45, length.out = n()))
  
  site_bray_list[[site_name]] <- site_bray_df
  site_stat_list[[site_name]] <- stat_test
}

final_bray_all  <- bind_rows(site_bray_list) %>% mutate(Site = factor(Site, levels = names(prof_filtered)))
final_stats_all <- bind_rows(site_stat_list) %>% mutate(Site = factor(Site, levels = names(prof_filtered)))

write.csv(final_bray_all, "Four_Sites_Bray_Distance_RawData.csv", row.names = FALSE)
write.csv(final_stats_all, "Four_Sites_Bray_Distance_stat.csv", row.names = FALSE)

# 4.1 Generate Single Plots without Facet (Filtering 'ns')
for (site_name in unique(final_bray_all$Site)) {
  site_bray_data     <- final_bray_all %>% filter(Site == site_name)
  site_stat_filtered <- final_stats_all %>% filter(Site == site_name & p.adj.signif != "ns")
  
  if (nrow(site_stat_filtered) > 0) {
    max_val            <- max(site_bray_data$bray_distance, na.rm = TRUE)
    site_stat_filtered <- site_stat_filtered %>% mutate(y.position = seq(max_val * 1.05, max_val * 1.40, length.out = n()))
    y_max_limit        <- max_val * 1.48
  } else {
    y_max_limit        <- max(site_bray_data$bray_distance, na.rm = TRUE) * 1.15
  }
  
  p_single <- ggplot(site_bray_data, aes(x = time_pair, y = bray_distance)) +
    stat_boxplot(geom = "errorbar", width = 0.35, linewidth = 0.5, aes(color = time_pair)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, size = 0.4, alpha = 0.8, aes(fill = time_pair, color = time_pair)) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = time_border_colors) +
    scale_y_continuous(limits = c(0, y_max_limit), expand = c(0, 0)) +
    labs(title = sprintf("%s - Longitudinal", site_name), y = "Bray-Curtis Dissimilarity", x = NULL) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 11, color = "black"), legend.position = "none"
    )
  
  if (nrow(site_stat_filtered) > 0) {
    p_single <- p_single +
      ggpubr::stat_pvalue_manual(data = site_stat_filtered, label = "p.adj.signif", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0.015, bracket.size = 0.45, size = 4.2, color = "grey20", inherit.aes = FALSE)
  }
  
  ggsave(sprintf("%s_bray_distance_boxplot.pdf", site_name), plot = p_single, width = 4.8, height = 5.5)
  write.csv(site_stat_filtered, sprintf("%s_bray_distance_stats.csv", site_name), row.names = FALSE)
}

# 5 - 7. permanova and alpha - shannon #################################################################
 
get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"
  )
}

cols <- c("BL/T4" = "#D9534F", "BL/T12" = "#F0AD4E", "BL/T24" = "#337AB7", "T4/T12" = "#5CB85C", "T4/T24" = "#9467BD", "T12/T24" = "#8C564B")
time_order_bray <- c("BL", "T4", "T12", "T24") 

permanova_stat_list <- list()
alpha_wilcox_list    <- list()

for (site_name in names(prof_filtered)) {
  
  mat_raw <- data.frame(prof_filtered[[site_name]])
  
  # 动态咬合临床时间盘与受试者 ID 字典
  mat_with_meta <- mat_raw %>%
    mutate(
      Group     = Microbe.phen.prof$Time[match(rownames(mat_raw), Microbe.phen.prof$SeqID)],
      Clinic_ID = Microbe.phen.prof$Clinic_ID[match(rownames(mat_raw), Microbe.phen.prof$SeqID)]
    ) %>%
    mutate(Group = recode(Group, "T04" = "T4")) %>%
    filter(Group %in% time_order_bray)
  
  if (filter_strategy == "strict_4") {
    strict_4_ids <- mat_with_meta %>%
      group_by(Clinic_ID) %>%
      summarise(n_time = n_distinct(Group), .groups = "drop") %>%
      filter(n_time == 4) %>% 
      pull(Clinic_ID)
    
    mat_with_meta <- mat_with_meta %>% filter(Clinic_ID %in% strict_4_ids)
  }
  
  pure_abundance_all <- mat_raw[rownames(mat_with_meta), , drop = FALSE] # 随之更新丰度表矩阵切片
  site_pairs_list    <- list()
  
  # ============================================================================
  # 5. PERMANOVA 整体配对分析模块###########################
  # ============================================================================
  for (pair_name in names(cols)) {
    
    tps        <- strsplit(pair_name, "/")[[1]]
    tp1        <- tps[1]
    tp2        <- tps[2]
    # comp_label <- paste0(tp1, "_vs_", tp2)
    
    # 严格受试者配对卡控线：寻找双时间点均存在的交集个体
    ids_in_tp1 <- mat_with_meta %>% filter(Group == tp1) %>% pull(Clinic_ID) %>% unique()
    ids_in_tp2 <- mat_with_meta %>% filter(Group == tp2) %>% pull(Clinic_ID) %>% unique()
    paired_clinic_ids <- intersect(ids_in_tp1, ids_in_tp2)
    if (length(paired_clinic_ids) < 3) next
    
    sub_meta <- mat_with_meta %>% 
      filter(Group %in% c(tp1, tp2) & Clinic_ID %in% paired_clinic_ids)
    
    n_g1 <- sum(sub_meta$Group == tp1)
    n_g2 <- sum(sub_meta$Group == tp2)
    
    sub_abundance <- pure_abundance_all[rownames(sub_meta), , drop = FALSE]
    sub_abundance <- sub_abundance[, colSums(sub_abundance) > 0, drop = FALSE] 
    
    sub_bray    <- vegan::vegdist(sub_abundance, method = "bray")
    dist_vector <- as.vector(sub_bray)
    
    med_val     <- round(median(dist_vector), 3)
    mean_val    <- round(mean(dist_vector), 3)
    sd_val      <- round(sd(dist_vector), 3)
    mean_sd_str <- paste0(mean_val, "_", sd_val)
    
    tryCatch({
      ad_res  <- vegan::adonis2(sub_bray ~ Group, data = sub_meta, permutations = 999)
      p_val   <- ad_res$`Pr(>F)`[1]
      r2_val  <- ad_res$R2[1]
      f_model <- ad_res$F[1]
      
      site_pairs_list[[pair_name]] <- data.frame(
        Sites            = site_name,
        Comparisons      = paste(tp1, "vs", tp2, sep = "_"),
        time_pair        = pair_name,
        N_Group1         = n_g1, 
        N_Group2         = n_g2,
        Median           = med_val,
        `Mean±SD`        = mean_sd_str,
        R2               = r2_val,
        F.model          = f_model,
        p_value          = p_val,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {})
  }
  
  if (length(site_pairs_list) > 0) {
    # 局部多重校正 (按部位独立实施 FDR)
    permanova_stat_list[[site_name]] <- bind_rows(site_pairs_list) %>%
      mutate(Adjusted_p = p.adjust(p_value, method = "BH")) %>%
      mutate(
        p_star   = get_significance_stars(p_value),
        q_star   = get_significance_stars(Adjusted_p),
        Signific = paste(p_star, q_star, sep = "/"),
        Method   = "PERMANOVA"
      ) %>%
      dplyr::select(-p_star, -q_star)
  }
  
  # ============================================================================
  # 6. Alpha 多样性 (Shannon) 整体纵向严格配对 Wilcoxon 检验模块##################
  # ============================================================================

  site_alpha_all <- prof_diversity %>%
    filter(Site == site_name) %>%
    mutate(Time = recode(Time, "T04" = "T4")) %>%
    filter(Time %in% time_order_bray) %>%
           drop_na(Shannon, Time, Clinic_ID)
           
  if (filter_strategy == "strict_4") {
    strict_4_alpha_ids <- site_alpha_all %>%
      group_by(Clinic_ID) %>%
      summarise(n_time = n_distinct(Time), .groups = "drop") %>%
      filter(n_time == 4) %>% 
      pull(Clinic_ID)
    
    site_alpha_all <- site_alpha_all %>% filter(Clinic_ID %in% strict_4_alpha_ids)
  }
  
  site_alpha_pairs <- list()
  
  for (pair_name in names(cols)) {
    tps <- strsplit(pair_name, "/")[[1]]
    tp1 <- tps[1]
    tp2 <- tps[2]
    
    ids_tp1    <- site_alpha_all %>% filter(Time == tp1) %>% pull(Clinic_ID)
    ids_tp2    <- site_alpha_all %>% filter(Time == tp2) %>% pull(Clinic_ID)
    paired_ids <- intersect(ids_tp1, ids_tp2)
    if (length(paired_ids) < 3) next
    
    calc_df <- site_alpha_all %>%
      filter(Time %in% c(tp1, tp2) & Clinic_ID %in% paired_ids) %>%
      mutate(Time = factor(Time, levels = c(tp1, tp2))) %>%
      arrange(Clinic_ID, Time)
    
    shannon_tp1 <- calc_df %>% filter(Time == tp1) %>% pull(Shannon)
    shannon_tp2 <- calc_df %>% filter(Time == tp2) %>% pull(Shannon)
    
    mean1_val <- round(mean(shannon_tp1), 6)
    mean2_val <- round(mean(shannon_tp2), 6)
    est_val   <- round(median(shannon_tp1 - shannon_tp2), 6) # 配对中位数差值
    
    tryCatch({
      w_test <- calc_df %>% wilcox_test(Shannon ~ Time, paired = TRUE, detailed = FALSE)
      
      site_alpha_pairs[[pair_name]] <- data.frame(
        Groups           = "Whole",
        Subgroups        = "-",
        Sites            = site_name,
        Comparisons      = paste0(tp1, "/", tp2),
        SampleNum1       = length(shannon_tp1),
        SampleNum2       = length(shannon_tp2),
        Mean1            = mean1_val,
        Mean2            = mean2_val,
        Estimate         = est_val,
        p_value          = w_test$p,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {})
  }
  
  if (length(site_alpha_pairs) > 0) {
    alpha_wilcox_list[[site_name]] <- bind_rows(site_alpha_pairs) %>%
      mutate(Adjusted_P_value = p.adjust(p_value, method = "BH")) %>%
      mutate(
        p_star       = get_significance_stars(p_value),
        q_star       = get_significance_stars(Adjusted_P_value),
        Significance = paste(p_star, q_star, sep = "/")
      ) %>%
      dplyr::select(-p_star, -q_star)
  }
}

# ----------------------------------------------------------------------------
# 最终数据平铺聚合、物理重排与写入写出
# ----------------------------------------------------------------------------
# 1. 导出全队列两两 PERMANOVA 大盘
final_permanova_all <- bind_rows(permanova_stat_list) %>%
  filter(Sites %in% c("GUT","TO"))%>%
  mutate(
    Sites     = factor(Sites, levels = names(prof_filtered)),
    time_pair = factor(time_pair, levels = names(cols))
  ) %>%
  arrange(Sites, time_pair)
write.csv(final_permanova_all, "Pairwise_PERMANOVA_Stats_GUT_TO_Sites_Time.csv", row.names = FALSE)


# 2. 导出全队列 Alpha 配对 Wilcoxon 检验报告 (完美的 Supp Table 13 样式)
final_alpha_wilcox_all <- bind_rows(alpha_wilcox_list) %>%
  arrange(factor(Sites, levels = names(prof_filtered)), factor(Comparisons, levels = names(cols)))
write.csv(final_alpha_wilcox_all, "Supplementary_Table_13_Whole_Cohort_Alpha_Wilcox_Paired.csv", row.names = FALSE)
