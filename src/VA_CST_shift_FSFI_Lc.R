library(tidyverse)
library(reshape2)
library(scales)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)

setwd("D:/WorkProjects/Demo-MHT 2026/Results/")

cst_ordered_levels <- c("UROG-L.c", "UROG-L.i", "UROG-G.v", "UROG-Div")

# ==============================================================================
# 1. Vaginal CST Transition Matrix (BL to T24) ###################################
# ==============================================================================
# Require paired complete samples at BL and T24
df_matrix_wide <- prof_diversity %>%
  dplyr::filter(Site == "VA" & Time %in% c("BL", "T24")) %>%
  dplyr::select(Clinic_ID, Time, CST) %>%
  reshape2::dcast(Clinic_ID ~ Time, value.var = "CST") %>%
  na.omit()

# Compute counts and row-normalized baseline transition probabilities
df_counts <- xtabs(~ BL + T24, data = df_matrix_wide) %>% 
  as.data.frame() %>% 
  rename(Count = Freq)

df_ratio <- xtabs(~ BL + T24, data = df_matrix_wide) %>% 
  prop.table(margin = 1) %>% 
  round(3) %>% 
  as.data.frame() %>% 
  rename(Ratio = Freq)

df_plot_final <- df_counts %>%
  inner_join(df_ratio, by = c("BL", "T24")) %>%
  mutate(
    T24        = factor(T24, levels = rev(cst_ordered_levels)),
    BL         = factor(BL, levels = rev(cst_ordered_levels)),
    Label_Text = if_else(Ratio > 0, paste0(sprintf("%.1f%%", Ratio * 100), "\n(n=", Count, ")"), "")
  )

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
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 9.5),
    axis.text        = element_text(size = 9, color = "black", face = "bold"),
    axis.text.x      = element_text(angle = 30, vjust = 1, hjust = 1),
    legend.title     = element_text(face = "bold", size = 8.5),
    legend.text      = element_text(size = 8),
    legend.position  = "right"
  )

ggsave("CST_Microbiome_Transition_Matrix.pdf", plot = p_cst_shift, width = 6.2, height = 5.2)
write.csv(df_plot_final %>% select(-Label_Text), "CST_Microbiome_Transition_Matrix.csv", row.names = FALSE)

# ==============================================================================
# Functions & Theme Configuration
# ==============================================================================
my_mean <- function(x) mean(x, na.rm = TRUE)

format_p_val <- function(p) {
  p <- as.numeric(p)
  ifelse(p < 0.001, sprintf("p = %.2e", p), sprintf("p = %.3f", p))
}

my_custom_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text        = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.line        = element_line(linewidth = 0.5),
      legend.position  = "bottom",
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

# Define stratification groups based on T24 vaginal CST status
group_data_t24 <- df.diversity %>%
  filter(Time == "T24", Site == "VA") %>%
  drop_na(vt) %>%
  mutate(
    Group = ifelse(vt == "UROG-L.c", "UROG-L.c at T24", "Non-UROG-L.c at T24"),
    Group = factor(Group, levels = c("UROG-L.c at T24", "Non-UROG-L.c at T24"))
  )

# ==============================================================================
# 2. Longitudinal FSFI Dynamics Stratified by T24 UROG-L.c #####################
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

# Retain participants with complete observations across all 4 timepoints
df_fsfi_complete <- df_fsfi %>%
  group_by(Clinic_ID) %>%
  filter(n() == 4) %>%
  ungroup()

# Non-parametric Friedman test for repeated measures within groups
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
# Pre-calculate Wilcoxon rank-sum test between groups per Timepoint
stat_fsfi_wilcox <- df_fsfi %>%
  group_by(Time) %>%
  wilcox_test(value ~ Group, paired = FALSE) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(
    p_fmt     = format_p_val(p),
    p_adj_fmt = format_p_val(p.adj)
  )

write.csv(stat_fsfi_wilcox, "FSFI_L.c_group_T24_Wilcox_Stats.csv", row.names = FALSE)

# Generate plotting coordinates for significance brackets
stat_fsfi_plot <- stat_fsfi_wilcox %>%
  add_xy_position(x = "Time", dodge = 0.8)

p_fsfi <- ggplot(df_fsfi, aes(x = Time, y = value, color = Group, fill = Group)) +
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
  stat_pvalue_manual(
    stat_fsfi_plot,
    label = "p_fmt",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3
  ) +
  scale_fill_manual(values = c("#CC6677", "#332288")) +
  scale_color_manual(values = c("#CC6677", "#332288")) +
  labs(caption = paste0(friedman_fsfi$label, collapse = "\n"), x = "Time Point", y = "FSFI") +
  my_custom_theme()

ggsave("FSFI_L.c_group_T24.pdf", plot = p_fsfi, width = 6, height = 5)
write.csv(df_fsfi %>% select(Time, FSFI, Group, value), "FSFI_L.c_group_T24.csv", row.names = FALSE)

# ==============================================================================
# 3. Longitudinal L. crispatus Abundance in VA and UR Sites ############################
# ==============================================================================
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

df_taxa_complete <- df_taxa_long %>%
  group_by(Site, Clinic_ID) %>%
  filter(n() == 4) %>%
  ungroup()

# Friedman test across timepoints within each Site and Group
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

# Pre-calculate Wilcoxon rank-sum test between groups per Site and Timepoint
stat_taxa_wilcox <- df_taxa_complete %>%
  group_by(Site, Time) %>%
  wilcox_test(value ~ Group, paired = FALSE) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(
    p_fmt     = format_p_val(p),
    p_adj_fmt = format_p_val(p.adj)
  )

write.csv(stat_taxa_wilcox, "UR_VA_L.c_group_T24_Wilcox_Stats.csv", row.names = FALSE)

# Generate plotting coordinates for significance brackets
stat_taxa_plot <- stat_taxa_wilcox %>%
  add_xy_position(x = "Time", dodge = 0.8)
p_taxa <- ggplot(df_taxa_complete, aes(x = Time, y = value, color = Group, fill = Group)) +
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
  stat_pvalue_manual(
    stat_taxa_plot,
    label = "p_fmt",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3
  ) +
  scale_fill_manual(values = c("#CC6677", "#332288")) +
  scale_color_manual(values = c("#CC6677", "#332288")) +
  labs(caption = paste0(friedman_taxa$label, collapse = "\n"), x = "Time Point", y = "Lactobacillus crispatus") +
  my_custom_theme()

ggsave("UR_VA_L.c_group_T24.pdf", plot = p_taxa, width = 10, height = 7)
write.csv(df_taxa_complete %>% select(Site, Time, Group, value), "UR_VA_L.c_group_T24.csv", row.names = FALSE)

# ==============================================================================
# 4. Cross-Sectional FSFI Profiles Stratified by Baseline Vaginal CST#####################
# ==============================================================================
group_data_bl <- df.diversity %>%
  filter(Time == "BL", Site == "VA") %>%
  drop_na(vt) %>%
  mutate(Group = factor(vt, levels = c("UROG-L.c", "UROG-L.i", "UROG-G.v", "UROG-Div")))

df_baseline <- Phen.Seq %>%
  select(Clinic_ID, Time, FSFI) %>%
  mutate(
    Group = group_data_bl$Group[match(Clinic_ID, group_data_bl$Clinic_ID)],
    value = FSFI
  ) %>%
  drop_na(Group, value) %>%
  mutate(Time = factor(Time, levels = unique(as.character(Time))))

my_com <- combn(levels(df_baseline$Group), 2, simplify = FALSE)

# Pre-calculate pairwise Wilcoxon tests across CST categories per Timepoint
stat_baseline_wilcox <- df_baseline %>%
  group_by(Time) %>%
  wilcox_test(value ~ Group, comparisons = my_com, paired = FALSE) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(
    p_fmt     = format_p_val(p),
    p_adj_fmt = format_p_val(p.adj)
  )

write.csv(stat_baseline_wilcox, "FSFI_va_Baseline_group_Wilcox_Stats.csv", row.names = FALSE)

# Filter for significant pairwise comparisons to avoid visual clutter
stat_baseline_plot <- stat_baseline_wilcox %>%
  add_xy_position(x = "Group", step.increase = 0.4) %>%
  filter(p < 0.05)

p_baseline <- ggplot(df_baseline, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.3, linewidth = 0.5) +
  geom_quasirandom(dodge.width = 0.8, alpha = 0.5, size = 0.5) +
  facet_grid(~Time) +
  stat_pvalue_manual(
    stat_baseline_plot,
    label = "p_fmt",
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3
  ) +
  scale_fill_manual(values = c("#CC6677", "#807DBA", "#9467BD", "#3F007D")) +
  scale_color_manual(values = c("#CC6677", "#807DBA", "#9467BD", "#3F007D")) +
  labs(x = "CST Classification", y = "FSFI") +
  my_custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("FSFI_va_Baseline_group.pdf", plot = p_baseline, width = 10, height = 6)
write.csv(df_baseline %>% select(FSFI, Time, Group, value), "FSFI_va_Baseline_group.csv", row.names = FALSE)
