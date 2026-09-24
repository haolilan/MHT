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

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
out_dir  <- "D:/WorkProjects/Demo-MHT 2026/Results/OverallEva/VA_REF_Microbe"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

load(file.path(base_dir, "data/MHT.demo.RData"))

prefix           <- "VA_Ref_Only"
rel_ab.threshold <- 0.0001

# Filter longitudinal cohort taxa (prevalence > 10% above threshold) and subjects with all 4 visits
data_va <- prof_filtered[["VA"]] %>% data.frame()
taxa_keep_va <- sapply(data_va, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(data_va))
taxa.vars_va <- names(taxa_keep_va)[taxa_keep_va]
data_va      <- data_va[, taxa.vars_va, drop = FALSE]

df_va <- data_va %>%
  mutate(
    Group       = Microbe.phen.prof$Time[match(rownames(data_va), Microbe.phen.prof$SeqID)],
    Clinic_ID   = Microbe.phen.prof$Clinic_ID[match(rownames(data_va), Microbe.phen.prof$SeqID)],
    Site.source = "VA" 
  )

Clinic_IDs_4_va <- df_va %>% 
  group_by(Clinic_ID) %>% 
  summarise(pairs = n(), .groups = "drop") %>% 
  filter(pairs == 4) %>% 
  pull(Clinic_ID)

df_va <- df_va %>% filter(Clinic_ID %in% Clinic_IDs_4_va)

# Process external reference cohort
phen.ref <- phen.va.ref
prof.ref <- prof.va.ref %>% column_to_rownames("X")

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

# Merge longitudinal cohort with external reference
cols_to_keep <- names(df_va)
df.merged <- dplyr::bind_rows(
  df_va,
  df_ref[intersect(names(df_ref), cols_to_keep)]
) %>%
  mutate(Group = recode(Group, 
                        "T04" = "T4", "T12" = "T12", "T24" = "T24",
                        "Menopause" = "H_M", "Reproductive" = "H_R"))

df_numeric_cols <- sapply(df.merged, is.numeric)
df.merged[, df_numeric_cols][is.na(df.merged[, df_numeric_cols])] <- 0

df.merged$Group <- factor(df.merged$Group, levels = c("BL", "T4", "T12", "T24", "H_R", "H_M"))
top_factors     <- levels(df.merged$Group)

df.merged$Site <- gsub("Ref", "VA", df.merged$Site.source) %>% factor(levels = c("VA"))
df.merged$CST  <- prof_diversity$CST[match(rownames(df.merged), prof_diversity$SeqID)] %>%
  factor(levels = rev(c("UROG-Div", "UROG-G.v", "UROG-L.i", "UROG-L.c")))

abundance_cols <- setdiff(colnames(df.merged), c("Group", "Clinic_ID", "Site.source", "Site", "CST"))
df_meta_labels <- df.merged %>% select(Group, Clinic_ID)
df_numeric_log <- df.merged[, abundance_cols]

pseu_value     <- min(df_numeric_log[df_numeric_log != 0]) * 0.1
df_numeric_log <- log10(df_numeric_log + pseu_value)
df_stat_ready  <- cbind(df_meta_labels, df_numeric_log)

my_comparisons <- list(
  c("BL", "T4"),  c("BL", "T12"), c("BL", "T24"), 
  c("BL", "H_M"), c("BL", "H_R"), 
  c("T24", "H_M"), c("T24", "H_R")
)

stats_results     <- list()
plot_data_results <- list()

# Paired Wilcoxon for longitudinal visits (1:3); unpaired Wilcoxon against reference groups (4:7)
for (i in seq_along(my_comparisons)) {
  t1 <- my_comparisons[[i]][1]
  t2 <- my_comparisons[[i]][2]
  
  if (i %in% 1:3) {
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

sum_stats_results <- bind_rows(stats_results) %>%
  mutate(
    lab.p.adj = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", p.adj < 0.05 ~ "*", TRUE ~ "ns"),
    lab.p     = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
    label     = paste(lab.p, lab.p.adj, sep = "/"),
    Sites     = "VA"
  ) %>%
  na.omit()

write.csv(sum_stats_results, paste0(prefix, ".TimeDiff.wilcox.csv"), row.names = FALSE)

# Heatmap of taxa with at least one significant contrast after FDR correction (p.adj < 0.05)
sig_taxa <- sum_stats_results %>% 
  filter(p.adj < 0.05) %>% 
  select(Factor) %>% 
  distinct()

if (nrow(sig_taxa) > 0) {
  draw <- sum_stats_results %>% 
    inner_join(sig_taxa, by = "Factor") %>% 
    mutate(
      groups = factor(comparisons, levels = c("BL/T4", "BL/T12", "BL/T24", "BL/H_M", "BL/H_R", "T24/H_M", "T24/H_R")),
      label  = gsub("ns/ns", "", label)
    )
  
  p_heat <- ggplot(draw, aes(x = groups, y = Factor, fill = Mean.diff.log)) +
    geom_tile(color = "lightgrey") +
    geom_text(aes(label = label), size = 2, color = "black") +  
    geom_vline(xintercept = c(3.5, 5.5), linetype = "dashed", color = "grey60") +
    facet_grid(Sites ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#5e3c99", mid = "grey98", high = "#b35806", midpoint = 0) +
    labs(x = NULL, y = NULL, fill = "Log-difference") +
    theme_minimal() +   
    theme(
      axis.text.y     = element_text(face = "bold.italic", size = 9),
      axis.text.x     = element_text(size = 10, angle = 30, hjust = 1, face = "bold"),
      panel.grid      = element_blank(),
      strip.text.y    = element_text(angle = -90, face = "bold", size = 11),  
      panel.spacing.y = unit(0.5, "lines")  
    )
  
  ggsave(paste0(prefix, ".TimeDiff.wilcox.heatmap.pdf"), plot = p_heat, width = 7, height = 6)
  write.csv(draw, paste0(prefix, ".TimeDiff.wilcox.heatmap.csv"), row.names = FALSE)
}