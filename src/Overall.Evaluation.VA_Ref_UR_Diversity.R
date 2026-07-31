
# Pipeline: Combined 3 Datasets (UR + VA + Reference) Micro-Ecology Analysis 
# Pipeline: Combined 2 Datasets (UR + VA  ) Micro-Ecology Analysis 

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

# pacman::p_load( broom, survival)

# ------------------------------------------------------------------------------
# 0. Global Parameters & Directory Configuration #################################
# ------------------------------------------------------------------------------
base_dir <- "D:/WorkProjects/Demo-MHT 2026"
out_dir  <- "D:/WorkProjects/Demo-MHT 2026/Results/OverallEva/VA_REF_UR"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

# Load global dataset workspace image
load(file.path(base_dir, "data/MHT.demo.RData"))


filter_strategy <- "strict_4" 

prefix           <- "Combined_UR_VA_Ref"

# 4. Core Matrix Merge ##################################

## 1. Dataset I Preparation: UR Longitudinal Cohort ##################

message(">>> Extracting and filtering UR longitudinal table...")

data_ur <- prof_filtered[["UR"]] %>% data.frame()

df_ur <- data_ur %>%
  mutate(
    Group       = Microbe.phen.prof$Time[match(rownames(data_ur), Microbe.phen.prof$SeqID)],
    Clinic_ID   = Microbe.phen.prof$Clinic_ID[match(rownames(data_ur), Microbe.phen.prof$SeqID)],
    Site.source = "UR" 
  )

if (filter_strategy == "strict_4") {
  Clinic_IDs_4_ur <- df_ur %>% 
    group_by(Clinic_ID, Group) %>% 
    group_by(Clinic_ID) %>% 
    summarise(pairs = n(), .groups = "drop") %>% 
    filter(pairs == 4) %>% 
    pull(Clinic_ID)
  
  df_ur <- df_ur %>% subset(Clinic_ID %in% Clinic_IDs_4_ur)
}

## 2. Dataset II Preparation: VA Longitudinal Cohort ################################

message(">>> Extracting and filtering VA longitudinal table...")

data_va <- prof_filtered[["VA"]] %>% data.frame()

df_va <- data_va %>%
  mutate(
    Group       = Microbe.phen.prof$Time[match(rownames(data_va), Microbe.phen.prof$SeqID)],
    Clinic_ID   = Microbe.phen.prof$Clinic_ID[match(rownames(data_va), Microbe.phen.prof$SeqID)],
    Site.source = "VA" 
  )

if (filter_strategy == "strict_4") {
  Clinic_IDs_4_va <- df_va %>% 
    group_by(Clinic_ID, Group) %>% 
    group_by(Clinic_ID) %>% 
    summarise(pairs = n(), .groups = "drop") %>% 
    filter(pairs == 4) %>% 
    pull(Clinic_ID)
  
  df_va <- df_va %>% subset(Clinic_ID %in% Clinic_IDs_4_va)
}

## 3. Dataset III Preparation: External Reference Cohort###########################

message(">>> Extracting and filtering external reference table...")

phen.ref <- phen.va.ref
prof.ref <- prof.va.ref %>% column_to_rownames("X")

intersect_samples <- intersect(phen.ref$SeqID, rownames(prof.ref))
prof.ref          <- prof.ref[intersect_samples, , drop = FALSE]
data_ref          <- prof.ref[rowSums(prof.ref) != 0, colSums(prof.ref) != 0, drop = FALSE]

df_ref <- data_ref %>%
  mutate(
    Group       = phen.ref$Group[match(rownames(data_ref), phen.ref$SeqID)],
    Clinic_ID   = phen.ref$Cohort[match(rownames(data_ref), phen.ref$SeqID)],
    Site.source = "Ref" 
  )


## 4. Core Matrix Merging & NA Zero Imputation ##################################

message(">>> Merging datasets on species union and mapping metadata...")

df.merged <- dplyr::bind_rows(df_ur, df_va, df_ref)

df.merged <- df.merged %>%
  mutate(Group = recode(Group, 
                        "T04" = "T4", "T12" = "T12", "T24" = "T24",
                        "Menopause" = "H_M", "Reproductive" = "H_R"))

df_numeric_cols <- sapply(df.merged, is.numeric)
df.merged[, df_numeric_cols][is.na(df.merged[, df_numeric_cols])] <- 0

df.merged$Group <- factor(df.merged$Group, levels = c("BL", "T4", "T12", "T24", "H_R", "H_M"))
top_factors     <- levels(df.merged$Group)

df.merged$Site <- gsub("Ref", "VA", df.merged$Site.source) %>% factor(levels = c("VA", "UR"))
df.merged$CST  <- prof_diversity$CST[match(rownames(df.merged), prof_diversity$SeqID)] %>%
  factor(levels = rev(c("UROG-Div", "UROG-G.v", "UROG-L.i", "UROG-L.c")))

print("Distribution of combined observational groups:")
print(xtabs(~ Site.source + Group, df.merged))


# 2. Analysis I: Taxonomic Composition Barplot ####################

message(">>> Generating taxonomic composition stacked barplot...")

abundance_cols <- setdiff(colnames(df.merged), c("Group", "Clinic_ID", "Site.source", "Site", "CST"))
df             <- df.merged[, abundance_cols]
rownames(df)   <- rownames(df.merged)

group_info       <- df.merged %>% select(Group, Clinic_ID, Site)
group_info$SeqID <- rownames(df.merged)

df <- df[, sapply(df, is.numeric)]
df <- data.frame(t(df))

df$mean_abundance <- rowMeans(df)
top_taxa <- df %>%
  arrange(desc(mean_abundance)) %>%
  head(30) %>%
  rownames()

top_data <- df[top_taxa, -ncol(df)]
others   <- 1 - colSums(df[rownames(df) %in% top_taxa, -ncol(df)])
top_data <- rbind(top_data, Others = others)

group_avg <- top_data %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SeqID") %>%
  left_join(group_info, by = "SeqID") %>%
  select(-SeqID, -Clinic_ID) %>%
  group_by(Group, Site) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  melt(id.vars = c("Group", "Site"), variable.name = "Taxonomy", value.name = "Abundance")

other_colors <- c(
  "#D55E00", "#E69F00", "#CC79A7", "#332288", "#DDCC77", "#AA4499", "#F7B6D2", 
  "#393B79", "#DBDB8D", "#882255", "#661100", "#AA4466", "#999933", "#9467BD",
  "#E58606", "#ED645A", "#CC3A8E", "#A5AA99", "#6A3D9A", "#B15928", "#8DD3C7", 
  "#FDB462", "#BEBADA"
)

taxa.vars <- group_avg$Taxonomy %>% unique()
manual_colors <- c(
  colorRampPalette(brewer.pal(9, "Greens")[5:9])(length(grep("Lactobacillus", taxa.vars))),
  colorRampPalette(brewer.pal(9, "Blues")[5:9])(length(grep("Prevotella", taxa.vars))),
  "#C7C7C7", 
  other_colors
)
names(manual_colors) <- c(
  taxa.vars[grep("Lactobacillus", taxa.vars)],
  taxa.vars[grep("Prevotella", taxa.vars)],
  taxa.vars[grep("Others", taxa.vars)],
  taxa.vars[grep("Lactobacillus|Prevotella|Others", taxa.vars, invert = TRUE)]
)

ggplot(group_avg, aes(x = Group, y = Abundance, fill = Taxonomy)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = manual_colors) +
  labs(x = "Group", y = "Relative Abundance (%)", title = "Top 30 Microbial Taxa") +
  scale_y_continuous(labels = scales::percent_format(scale = 100), expand = c(0, 0)) +
  facet_grid(~ Site, scales = "free", space = "free") +
  guides(fill = guide_legend(ncol = 2)) +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 11),
    axis.title         = element_text(size = 12),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  )

ggsave("VA_UR_barplot.pdf", width = 10, height = 6)
write.csv(group_avg, "VA_UR_barplot_avg.csv", row.names = FALSE, quote = FALSE)


# 3. Analysis II: Alpha Diversity (Shannon Index Calculation & Plotting) ##########################
get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"
  )
}

message(">>> Starting Alpha diversity pipeline...")
abundance_matrix <- as.matrix(df.merged[, abundance_cols])
df.merged$Shannon <- vegan::diversity(abundance_matrix, index = "shannon")

## summary shannon
summary_table <- df.merged %>%
  drop_na(Shannon,Site, Group) %>%
  group_by(Site, Group) %>%
  summarise(
    NumberofSample     = n(),
    Median             = round(median(Shannon), 3),
    Mean_val           = round(mean(Shannon), 3),
    SD_val             = round(sd(Shannon), 3),
    .groups            = "drop"
  ) %>%
  mutate(
    Mean_SD = paste0(Mean_val, "_", SD_val)
  ) %>%
  dplyr::select(Sites = Site, Timepoint = Group, NumberofSample, Median, Mean_SD)
write.csv(summary_table, "VA_UR_alpha_summary.csv", row.names = FALSE, quote = FALSE)


paired_test_by_site_time <- function(data_input, target_site, t1, t2) {
  tmp <- data_input %>%
    filter(Site.source == target_site, Group %in% c(t1, t2)) %>%
    group_by(Clinic_ID) %>%
    filter(n_distinct(Group) == 2) %>%
    ungroup()
  
  if (nrow(tmp) == 0) return(NULL)
  
  # 物理重排与对齐：物理防区线
  t1_df <- tmp %>% filter(Group == t1) %>% arrange(Clinic_ID)
  t2_df <- tmp %>% filter(Group == t2) %>% arrange(Clinic_ID)
  
  # 确保受试者完全一一匹配
  common_ids <- intersect(t1_df$Clinic_ID, t2_df$Clinic_ID)
  t1_df <- t1_df %>% filter(Clinic_ID %in% common_ids)
  t2_df <- t2_df %>% filter(Clinic_ID %in% common_ids)
  
  n_samples <- nrow(t1_df)
  if (n_samples < 3) return(NULL)
  
  val1 <- t1_df$Shannon
  val2 <- t2_df$Shannon
  
  # 配对 Wilcoxon 检验
  test <- wilcox.test(val1, val2, paired = TRUE, exact = FALSE)
  
  # 计算指标：Mean1, Mean2 以及配对差值中位数（Estimate）
  mean1_val <- mean(val1, na.rm = TRUE)
  mean2_val <- mean(val2, na.rm = TRUE)
  est_val   <- median(val2 - val1, na.rm = TRUE) # 对应图片中 Estimate 算法
  
  data.frame(
    Groups           = "Whole",
    Subgroups        = "-",
    Sites            = target_site,
    Comparisons      = paste0(t1, "/", t2),
    SampleNum1       = n_samples,
    SampleNum2       = n_samples,
    Mean1            = round(mean1_val, 6),
    Mean2            = round(mean2_val, 6),
    Estimate         = round(est_val, 6),
    p_value          = test$p.value,
    stringsAsFactors = FALSE
  )
}

unpaired_test_by_site <- function(data_input, target_site, g1, g2) {
  tmp <- data_input %>% filter((Site.source == target_site & Group == g1) | (Site.source == "Ref" & Group == g2))
  x   <- tmp$Shannon[tmp$Site.source == target_site & tmp$Group == g1]
  y   <- tmp$Shannon[tmp$Site.source == "Ref" & tmp$Group == g2]
  
  n1  <- length(x)
  n2  <- length(y)
  
  if (n1 < 3 || n2 < 3) return(NULL)
  
  test <- wilcox.test(x, y, paired = FALSE, exact = FALSE)
  
  mean1_val <- mean(x, na.rm = TRUE)
  mean2_val <- mean(y, na.rm = TRUE)
  est_val   <- median(y, na.rm = TRUE) - median(x, na.rm = TRUE) 
  
  data.frame(
    Groups           = "Whole",
    Subgroups        = "-",
    Sites            = target_site,
    Comparisons      = paste0(g1, "/", g2),
    SampleNum1       = n1,
    SampleNum2       = n2,
    Mean1            = round(mean1_val, 6),
    Mean2            = round(mean2_val, 6),
    Estimate         = round(est_val, 6),
    p_value          = test$p.value,
    stringsAsFactors = FALSE
  )
}

## paired
time_pairs       <- list(c("BL", "T4"), c("BL", "T12"), c("BL", "T24"), c("T4", "T12"), c("T4", "T24"), c("T12", "T24"))
stat_paired_list <- list()

for (st in c("UR", "VA")) {
  for (pair in time_pairs) {
    res <- paired_test_by_site_time(df.merged, st, pair[1], pair[2])
    if (!is.null(res)) stat_paired_list[[length(stat_paired_list) + 1]] <- res
  }
}
stat_paired_df <- bind_rows(stat_paired_list)

## unpaired
follow_up          <- c("BL", "T4", "T12", "T24")
ref_groups         <- c("H_R", "H_M")
unpaired_grid      <- expand_grid(follow_up, ref_groups)
stat_unpaired_list <- list()

for (st in c("VA")) {
  for (r in 1:nrow(unpaired_grid)) {
    res <- unpaired_test_by_site(df.merged, st, unpaired_grid$follow_up[r], unpaired_grid$ref_groups[r])
    if (!is.null(res)) stat_unpaired_list[[length(stat_unpaired_list) + 1]] <- res
  }
}
stat_unpaired_df <- bind_rows(stat_unpaired_list)

## combine
stat_all_df <- bind_rows(stat_paired_df, stat_unpaired_df) %>%
  group_by(Sites) %>% 
  mutate(Adjusted_P_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    p_star       = get_significance_stars(p_value),
    q_star       = get_significance_stars(Adjusted_P_value),
    Significance = paste(p_star, q_star, sep = "/") 
  ) 

write.csv(stat_all_df%>%dplyr::select(-p_star, -q_star), paste0(prefix, ".ST13_Ref_Cohort_Alpha_Wilcox_Paired_Full.csv"), row.names = FALSE)


# 4. data for plot 
df_plot_alpha <- df.merged %>% mutate(Site = factor(Site, levels = c("VA", "UR")))

# 转换 Comparisons 回 group1 / group2 结构，供 ggpubr 绘图读取
plot_stat_df <- stat_all_df %>%
  filter(Significance != "ns/ns") %>%  
  mutate(
    group1 = sapply(strsplit(Comparisons, "/"), `[`, 1),
    group2 = sapply(strsplit(Comparisons, "/"), `[`, 2)
  ) %>%
  group_by(Sites) %>% 
  mutate(
    y.position = max(df_plot_alpha$Shannon, na.rm = TRUE) * (0.95 + 0.04 * row_number())
  ) %>%
  ungroup() %>%
  rename(Site = Sites, label = q_star)

# 5. plot

time_colors   <- c("BL" = "#E15759", "T4" = "#F1C40F", "T12" = "#56B4E9", "T24" = "#59A14F", "H_R" = "#9B59B6", "H_M" = "#7F7F7F")
border_colors <- sapply(time_colors, function(x) rgb(t(col2rgb(x) * 0.7) / 255)) 

p_alpha <- ggplot(df_plot_alpha, aes(x = Group, y = Shannon, fill = Group, color = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = time_colors) +
  scale_color_manual(values = border_colors) +
  facet_grid(. ~ Site, space = "free", scales = "free_x") + 
  theme_bw(base_size = 14) +
  labs(y = "Shannon index", x = NULL, title = "Alpha Diversity (Shannon Index) across Sites") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    axis.text.x      = element_text(face = "bold", angle = 30, hjust = 1),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none"
  )

if (nrow(plot_stat_df) > 0) {
  p_alpha <- p_alpha +
    ggpubr::stat_pvalue_manual(
      data        = plot_stat_df,        
      label       = "label",            
      xmin        = "group1",            
      xmax        = "group2",            
      y.position  = "y.position",       
      tip.length  = 0.01,              
      size        = 3.2,                
      color       = "grey30",            
      inherit.aes = FALSE                
    )
}

print(p_alpha)
ggsave(paste0(prefix, ".Shannon_Timepoint.pdf"), plot = p_alpha, width = 9, height = 6)


## supplot.Shannon comparison between sites split by CST -------------------------------
ggplot(df.merged %>% filter(Site.source != "Ref"), aes(x = CST, y = Shannon, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("UR" = "#ADD8E6", "VA" = "#B19CD9")) +  
  labs(title = "Shannon Index Comparison", x = " ", y = "Shannon Index", fill = "Site") +
  theme_light(base_size = 13) +
  theme(
    legend.position  = "right",
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  stat_compare_means(aes(group = Site), method = "wilcox.test", label = "p.signif")

ggsave(paste0(prefix, ".Shannon_CST.pdf"), width = 6, height = 4)


# 4. Analysis III: Beta Diversity Analysis (PCoA & Ellipses Hardening) ##################
message(">>> Commencing PCoA Beta-Diversity")

pure_abundance <- df.merged %>%
  select(all_of(abundance_cols)) %>%
  .[rowSums(.) > 0, colSums(.) > 0]

df.merged.pure <- df.merged %>% filter(rownames(df.merged) %in% rownames(pure_abundance))

bray_dist   <- vegdist(pure_abundance, method = "bray")
pcoa_out    <- cmdscale(bray_dist, k = 2, eig = TRUE)
axes_pct    <- round(pcoa_out$eig[1:2] / sum(pcoa_out$eig) * 100, 1)

pcoa_plot_df <- data.frame(
  PC1       = pcoa_out$points[, 1],
  PC2       = pcoa_out$points[, 2],
  Group     = df.merged.pure$Group,
  Site      = df.merged.pure$Site,
  Clinic_ID = df.merged.pure$Clinic_ID,
  CST       = df.merged.pure$CST
) %>% mutate(Site = factor(Site, levels = c("VA", "UR")))


## Colored by Timepoint --------------------------------------------------------
p_pcoa <- ggplot(pcoa_plot_df, aes(x = PC1, y = PC2, color = Group, shape = Site)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(linetype = Site), type = "t", level = 0.75, size = 0.6) +
  labs(title = "Beta Diversity - PCoA (Bray-Curtis)", x = paste0("PC1 (", axes_pct[1], "%)"), y = paste0("PC2 (", axes_pct[2], "%)")) +
  scale_color_manual(values = time_colors) +
  scale_shape_manual(values = c("UR" = 17, "VA" = 19)) +
  scale_linetype_manual(values = c("UR" = "dashed", "VA" = "solid")) +
  theme_bw() +
  theme(
    plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text        = element_text(color = "black", size = 9),
    panel.grid.minor = element_blank()
  )

print(p_pcoa)
ggsave(paste0(prefix, ".PCoA_Timepoint.pdf"), plot = p_pcoa, width = 8, height = 6)

## Colored by CST (with Reference groups integrated safely) --------------------
cst_colors <- c(
  "UROG-L.c" = "#117733", "UROG-G.v" = "#FF7F0E", 
  "UROG-L.i" = "#74C476", "UROG-Div" = "#2171B5",
  "H_R"      = "#9B59B6", "H_M"      = "#7F7F7F"
)

pcoa_plot_clean <- pcoa_plot_df %>%
  mutate(
    # Map the reference groups directly into the color grouping to prevent unaligned drops
    CST_plot = ifelse(is.na(CST) | CST == "", as.character(Group), as.character(CST)),
    CST_plot = factor(CST_plot, levels = c("H_R", "H_M","UROG-L.c", "UROG-G.v", "UROG-L.i", "UROG-Div"))
  )%>%arrange(CST_plot)

p_pcoa_cst <- ggplot(pcoa_plot_clean, aes(x = PC1, y = PC2, color = CST_plot, shape = Site)) +
  geom_point(size = 2.2, alpha = 0.6) +
  stat_ellipse(data = pcoa_plot_clean%>%filter(!(CST_plot%in%c("H_R", "H_M"))),aes(linetype = Site,color = CST_plot), 
               type = "norm", level = 0.75, size = 0.5) +
  labs(title = "Beta Diversity - PCoA (Bray-Curtis)", x = paste0("PC1 (", axes_pct[1], "%)"), y = paste0("PC2 (", axes_pct[2], "%)")) +
  scale_color_manual(values = cst_colors, name = "") +
  scale_shape_manual(values = c("UR" = 17, "VA" = 19)) +
  scale_linetype_manual(values = c("UR" = "dashed", "VA" = "solid")) +
  theme_bw() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text         = element_text(color = "black", size = 9),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank()
  )
print(p_pcoa_cst)
ggsave(paste0(prefix, ".PCoA_CST.pdf"), plot = p_pcoa_cst, width = 8, height = 6)

message(">>> Multi-omic diversity calculations finished.")

## 4.1 . PERMANOVA 整体两两分析 VA /ref/UR ###########################

permanova_stat_list <- list()

for (site_name in c("VA", "UR")) {
  
  site_pairs_list <- list()
  
  if (site_name == "UR") {
    time_order_bray <- c("BL", "T4", "T12", "T24")
  } else {
    time_order_bray <- c("BL", "T4", "T12", "T24", "H_R", "H_M")
  }
  
  # 基于当前 site 的时间点生成组合矩阵
  pair_matrix <- combn(time_order_bray, 2)
  
  # 提取当前 Site 的元数据 (以 df.merged.pure 为基准，确保行名与 pure_abundance 一致)
  site_meta <- df.merged.pure %>% 
    filter(Site == site_name)
  
  if (nrow(site_meta) == 0) next
  
  for (i in 1:ncol(pair_matrix)) {
    
    tp1 <- pair_matrix[1, i]
    tp2 <- pair_matrix[2, i]
    
    pair_name  <- paste0(tp1, "_vs_", tp2)
    
    # 🌟 移除配对限制：直接提取当前 Site 下属于 tp1 或 tp2 的所有样本
    sub_meta <- site_meta %>% 
      filter(Group %in% c(tp1, tp2))
    
    n_g1 <- sum(sub_meta$Group == tp1)
    n_g2 <- sum(sub_meta$Group == tp2)
    
    # 保证两组都有足够的样本量（至少各 3 个）
    if (n_g1 < 3 || n_g2 < 3) next
    
    # 提取对应的丰度子矩阵并剔除全零物种
    sub_abundance <- pure_abundance[rownames(sub_meta), , drop = FALSE]
    sub_abundance <- sub_abundance[, colSums(sub_abundance) > 0, drop = FALSE]
    
    if (ncol(sub_abundance) == 0) next
    
    # 计算 Bray-Curtis 距离矩阵及描述性统计量
    sub_bray    <- vegan::vegdist(sub_abundance, method = "bray")
    dist_vector <- as.vector(sub_bray)
    
    med_val     <- round(median(dist_vector), 3)
    mean_val    <- round(mean(dist_vector), 3)
    sd_val      <- round(sd(dist_vector), 3)
    mean_sd_str <- paste0(mean_val, " ± ", sd_val)
    
    # PERMANOVA 检验
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
    }, error = function(e) {
      warning(paste("PERMANOVA failed for", site_name, pair_name, ":", e$message))
    })
  }
  # 局部多重校正 (按部位独立实施 BH FDR 校正)
  if (length(site_pairs_list) > 0) {
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
}

final_permanova_all <- bind_rows(permanova_stat_list)
write.csv(final_permanova_all, "Pairwise_PERMANOVA_Stats_VA_REF_UR_Time.csv", row.names = FALSE)


# 5. Analysis IV: UR & VA CST Consistence Analysis  ####################

col_cons      <- c("UROG-L.i_Same" = "#74C476", "UROG-L.c_Same" = "#117733", "UROG-G.v_Same" = "#FF7F0E", "UROG-Div_Same" = "#2171B5", "Lacto_Both" = "#C7E9C0", "Lacto_NonLacto" = "#BEBADA", "NonLacto_Both" = "#CC79A7")
target_levels <- names(col_cons)
time_points   <- c("BL", "T04", "T12", "T24")

dat_cst_wider <- prof_diversity %>%
  filter(Site %in% c("UR", "VA") & Time %in% time_points) %>%
  pivot_wider(id_cols = c(Clinic_ID, Time), names_from = Site, values_from = CST) %>%
  drop_na(UR, VA) %>%
  mutate(
    new_var = case_when(
      UR == VA ~ paste0(UR, "_Same"),
      grepl("^UROG-L", UR) & grepl("^UROG-L", VA)   ~ "Lacto_Both",
      !grepl("^UROG-L", UR) & !grepl("^UROG-L", VA) ~ "NonLacto_Both",
      xor(grepl("^UROG-L", UR), grepl("^UROG-L", VA))~ "Lacto_NonLacto",
      TRUE ~ NA_character_
    ))
         
dat_cst_filtered <- dat_cst_wider

count_data <- dat_cst_filtered %>%
  count(Time, new_var, name = "Count") %>%
  group_by(Time) %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Time       = factor(Time, levels = time_points),
    new_var    = factor(new_var, levels = target_levels)
  ) %>%
  ungroup()

total_labels <- count_data %>% 
  group_by(Time) %>% 
  summarise(Total = sum(Count), .groups = "drop")

p_cons <- ggplot(count_data, aes(x = Time, y = Percentage, fill = new_var)) +
  geom_col(position = position_stack(), width = 0.65, alpha = 0.9, color = "white", size = 0.2) +
  geom_text(aes(label = sprintf("%d, %.1f%%", Count, Percentage)), position = position_stack(vjust = 0.5), color = "black", size = 3) +
  geom_text(data = total_labels, aes(x = Time, y = 104, label = paste0("n = ", Total)), inherit.aes = FALSE, size = 3.8, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = col_cons, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 108)) +
  labs(title = "Consistence of CSTs", x = " ", y = "Percentage (%)", fill = "Consistence") +
  theme_bw(base_size = 12) +
  theme(
    plot.title        = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle     = element_text(size = 11, color = "grey40"),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x= element_blank(), 
    legend.position   = "right",
    legend.title      = element_text(face = "bold")
  )

print(p_cons)
ggsave("Consistence_of_CSTs_4samp.pdf", plot = p_cons, width = 7.5, height = 7)
         

# 6. Analysis V: UR & VA CST - Conditional Logistic Regression (clogit) ##############

pacman::p_load(dplyr, purrr, broom, survival, ggplot2, patchwork)

time_points <- c("BL", "T04", "T12", "T24")
cst_vars    <- c("UROG-L.c", "UROG-G.v", "UROG-L.i", "UROG-Div")
cst_colors  <- c("UROG-L.c" = "#E15759", "UROG-G.v" = "#F1C40F", "UROG-L.i" = "#56B4E9", "UROG-Div" = "#59A14F")

phen_clean <- prof_diversity %>%
  filter(Time %in% time_points) %>%
  mutate(CST_label = CST) %>%
  filter(!is.na(CST_label) & Site %in% c("UR", "VA"))

analyze_tp <- function(tp) {
  data_sub <- phen_clean %>%
    filter(Time == tp) %>%
    group_by(Clinic_ID) %>%
    filter(any(Site == "UR") & any(Site == "VA")) %>%
    ungroup()
  
  cst_props <- data_sub %>%
    group_by(Site) %>%
    mutate(Total_Samples = n()) %>% 
    group_by(Site, CST_label, Total_Samples) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Pct = Count / Total_Samples * 100,
      Display_Str = sprintf("%d/%d (%.1f%%)", Count, Total_Samples, Pct)
    ) %>%
    select(Site, CST_label, Display_Str) %>%
    tidyr::pivot_wider(names_from = Site, values_from = Display_Str, names_prefix = "Prop_")
  
  res <- map_dfr(cst_vars, function(var) {
    sub_model <- data_sub %>% mutate(y = as.integer(CST_label == var))
    
    if (sum(sub_model$y) == 0 || sum(sub_model$y) == nrow(sub_model)) {
      return(data.frame(term = "SiteVA", estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_, CST = var))
    }
    
    contingency_table <- table(sub_model$Site, sub_model$y)
    if (any(rowSums(contingency_table) == 0) || any(colSums(contingency_table) == 0)) {
      return(data.frame(term = "SiteVA", estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_, CST = var))
    }
    
    tryCatch({
      fit <- clogit(y ~ Site + strata(Clinic_ID), data = sub_model)
      tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "SiteVA") %>% mutate(CST = var)
    }, error = function(e) {
      message(sprintf(">>> Warning: Fitting failed for CST [%s] at time [%s] due to numerical issues.", var, tp))
      return(data.frame(term = "SiteVA", estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_, CST = var))
    })
  }) %>%
    mutate(
      p.adjust = p.adjust(p.value, "BH"),
      label    = ifelse(!is.na(p.adjust) & p.adjust <= 0.05, sprintf("%.3f", p.adjust), ""),
      CST      = factor(CST, levels = rev(cst_vars))
    ) %>%
    left_join(cst_props, by = c("CST" = "CST_label")) %>%
    mutate(
      Prop_UR = ifelse(is.na(Prop_UR), "0%", Prop_UR),
      Prop_VA = ifelse(is.na(Prop_VA), "0%", Prop_VA)
    ) %>%
    select(CST, Prop_UR, Prop_VA, everything())
  
  res_plot_df <- res %>% filter(!is.na(estimate) & !is.na(conf.low) & !is.na(conf.high),
                                !is.infinite(estimate) & !is.infinite(conf.low) & !is.infinite(conf.high))
  
  p <- ggplot(res_plot_df, aes(x = CST, y = estimate, ymin = conf.low, ymax = conf.high, color = CST)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(width = 0.15, size = 0.8) +
    geom_point(size = 3) +
    geom_text(aes(label = label), hjust = -0.15, vjust = -0.4, size = 3.2, fontface = "italic", show.legend = FALSE) +
    coord_flip() +
    scale_color_manual(values = cst_colors, drop = FALSE) + 
    scale_y_log10() + 
    labs(title = tp, y = "Odds Ratio", x = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), panel.grid.minor = element_blank(), legend.position = "none")
  
  list(res = res, plot = p)
}

out_all     <- map(time_points, analyze_tp) %>% set_names(time_points)
final_stats <- map_dfr(out_all, "res", .id = "Timepoint")
write.csv(final_stats,"CST_clogit_facet_forest_plot.csv")

## plot 
res_plot_all <- final_stats %>% 
  filter(
    !is.na(estimate) & !is.na(conf.low) & !is.na(conf.high),
    !is.infinite(estimate) & !is.infinite(conf.low) & !is.infinite(conf.high),
    conf.low > 0.00001,
    conf.high < 100000
  ) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("BL", "T04", "T12", "T24"), 
                       labels = c("BL", "T4", "T12", "T24"))
  )

p_combined_facet <- ggplot(res_plot_all, aes(x = CST, y = estimate, ymin = conf.low, ymax = conf.high, color = CST)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
  geom_errorbar(width = 0.15, size = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = label), hjust = -0.15, vjust = -0.4, size = 3.2, fontface = "italic", show.legend = FALSE) +
  coord_flip() +
  scale_color_manual(values = cst_colors, drop = FALSE) + 

  facet_grid(. ~ Timepoint, scales = "fixed", space = "fixed") + 
  
  labs(
    title   = "Comparison of CST distribution between UR and VA at each time point", 
    y       = "Odds Ratio", 
    x       = NULL,
    caption = "Filtered out incomplete separations (Inf/0 conf intervals)"
  ) +
  theme_light(base_size = 12) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.background   = element_rect(fill = "white", color = "white"), 
    strip.text         = element_text(face = "bold", size = 11,color = "black"),         
    panel.grid.minor   = element_blank(),
    panel.grid.major = element_blank(),                               
    legend.position    = "none"
  )

# 3. 打印图表并高清输出
print(p_combined_facet)
ggsave("CST_clogit_facet_forest_plot.pdf", plot = p_combined_facet, width = 11, height = 4.5)