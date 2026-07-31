# ==============================================================================
# Pipeline: Streamlined PAM Clustering, PCoA & Top 6 Taxa Profiling
# Groups: [UR+VA Combined], [TO Independent], [GUT Independent]
# Year: 2026
# ==============================================================================
library(vegan)
library(cluster)
library(factoextra)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

out_dir <- "D:/WorkProjects/Demo-MHT 2026/Results/OverallEva/CST_PAM_Clustering"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

warm_colors <- c("#E15759", "#F1C40F", "#56B4E9", "#59A14F", "#9B59B6", "#FF9E4A", "#B07AA1")
target_datasets <- list()

# 数据准备 ###########################
# Load global dataset workspace image
load(file.path(base_dir, "data/MHT.demo.RData"))

df_ur_raw <- data.frame(prof_filtered[["UR"]])
df_va_raw <- data.frame(prof_filtered[["VA"]])
df_ur_va  <- dplyr::bind_rows(df_ur_raw, df_va_raw)
df_ur_va[is.na(df_ur_va)] <- 0

target_datasets[["UR_VA"]] <- df_ur_va
target_datasets[["TO"]]    <- data.frame(prof_filtered[["TO"]])
target_datasets[["GUT"]]   <- data.frame(prof_filtered[["GUT"]])

## clustering #######################################
compiled_cluster_records <- list()

for (group_name in names(target_datasets)) {
  data_mat <- target_datasets[[group_name]]
  dist_mat <- vegan::vegdist(data_mat, method = "bray")
  
  sil_plot <- fviz_nbclust(data_mat, cluster::pam, method = "silhouette", diss = dist_mat, k.max = 10) +
    labs(title = sprintf("Optimal Number of Clusters (PAM) - %s", group_name)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  best_k <- which.max(sil_plot$data$y)
  ggsave(sprintf("%s_best_cluster_silhouette.pdf", group_name), plot = sil_plot, width = 5, height = 3.8)
  
  pam_res <- cluster::pam(dist_mat, k = best_k)
  
  current_cluster_list <- data.frame(
    SeqID       = rownames(data_mat),
    CST_Cluster = paste0(group_name, "_Cluster_", pam_res$clustering),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Site = Microbe.phen.prof$Site[match(SeqID, Microbe.phen.prof$SeqID)]) %>%
    select(SeqID, Site, CST_Cluster)
  
  compiled_cluster_records[[group_name]] <- current_cluster_list
  
  pcoa_res  <- ape::pcoa(dist_mat)
  data_plot <- data.frame(PCoA1 = pcoa_res$vectors[, 1], PCoA2 = pcoa_res$vectors[, 2], Cluster = factor(pam_res$clustering))
  var_exp   <- pcoa_res$values$Relative_eig * 100
  
  p_pcoa <- ggplot(data_plot, aes(x = PCoA1, y = PCoA2, color = Cluster, fill = Cluster)) +
    geom_point(size = 2, alpha = 0.8) + 
    stat_ellipse(geom = "polygon", alpha = 0.2, linewidth = 0.4) +
    labs(title = sprintf("PCoA - %s", group_name), x = sprintf("PCoA1 (%.2f%%)", var_exp[1]), y = sprintf("PCoA2 (%.2f%%)", var_exp[2])) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank())
  ggsave(sprintf("%s_PCoA.pdf", group_name), plot = p_pcoa, width = 5.5, height = 4.2)
  
  data_mat$cluster  <- factor(pam_res$clustering)
  data_long_samples <- data_mat %>% rownames_to_column("SeqID") %>% reshape2::melt(id.vars = c("SeqID", "cluster"), variable.name = "Taxon", value.name = "Abundance")
  
  taxa_global_order <- data_long_samples %>% group_by(Taxon) %>% summarise(Global_Mean = mean(Abundance, na.rm = TRUE), .groups = "drop") %>% arrange(desc(Global_Mean)) %>% pull(Taxon)
  top6_species      <- head(taxa_global_order, 6)
  plot_boxplot_df   <- data_long_samples %>% filter(Taxon %in% top6_species) %>% mutate(Taxon = factor(Taxon, levels = top6_species), Percentage = Abundance * 100)
  
  p_taxa_boxplot <- ggplot(plot_boxplot_df, aes(x = Taxon, y = Percentage, fill = cluster, color = cluster)) +
    stat_boxplot(geom = "errorbar", size = 1.2, width = 0.3, position = position_dodge(0.75)) +
    geom_boxplot(width = 0.6, size = 1.2, alpha = 0.6, outlier.shape = 1, position = position_dodge(0.75)) +
    scale_fill_manual(values = warm_colors[1:best_k], name = "Cluster") +
    scale_color_manual(values = warm_colors[1:best_k], name = "Cluster") +
    labs(title = sprintf("Top 6 Species Boxplot - %s", group_name), x = NULL, y = "Relative Abundance (%)") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1, face = "bold.italic"), panel.grid.major.x = element_blank())
  ggsave(sprintf("%s_Dominant_Species_Boxplot.pdf", group_name), plot = p_taxa_boxplot, width = 8, height = 5)
}

final_integrated_cluster_list <- dplyr::bind_rows(compiled_cluster_records)
print(head(final_integrated_cluster_list, 10))
write.csv(final_integrated_cluster_list, "Consolidated_Sample_CST_Cluster_List.csv", row.names = FALSE)

# refine CST name and plot #####################################
cst_colors <- c(
  "GUT-P.v"   = "#66c2a4", "GUT-P.cop" = "#b2e2e2",
  "TO-P.m"    = "#fc8d59", "TO-N.s"   = "#fdcc8a",
  "UROG-L.c"  = "#117733", "UROG-L.i" = "#74C476", 
  "UROG-G.v"  = "#FF7F0E", "UROG-Div" = "#2171B5"
)
cst_border_colors <- colorspace::darken(cst_colors, 0.3)

final_cst_mapped_list <- final_integrated_cluster_list %>%
  mutate(
    CST = case_when(
      CST_Cluster == "UR_VA_Cluster_1" ~ "UROG-L.c",
      CST_Cluster == "UR_VA_Cluster_2" ~ "UROG-G.v",
      CST_Cluster == "UR_VA_Cluster_3" ~ "UROG-L.i",
      CST_Cluster == "UR_VA_Cluster_4" ~ "UROG-Div",
      CST_Cluster == "TO_Cluster_1"    ~ "TO-N.s",
      CST_Cluster == "TO_Cluster_2"    ~ "TO-P.m",
      CST_Cluster == "GUT_Cluster_1"   ~ "GUT-P.v",
      CST_Cluster == "GUT_Cluster_2"   ~ "GUT-P.cop",
      TRUE                             ~ NA_character_
    )
  ) %>%
  filter(!is.na(CST))

write.csv(final_cst_mapped_list, "Final_Consolidated_CST_Mapped_List.csv", row.names = FALSE)

group_mapping <- list("UR" = "UR", "VA" = "VA", "TO" = "TO", "GUT" = "GUT")

for (group_name in names(group_mapping)) {
  target_sites <- group_mapping[[group_name]]
  
  data_mat_raw <- data.frame(prof_filtered[[group_name]])
  data_mat_raw[is.na(data_mat_raw)] <- 0
  
  sub_cst_meta      <- final_cst_mapped_list %>% filter(Site %in% target_sites)
  intersect_samples <- intersect(rownames(data_mat_raw), sub_cst_meta$SeqID)
  
  if (length(intersect_samples) == 0) next
  
  data_mat_clean     <- data_mat_raw[intersect_samples, , drop = FALSE]
  data_mat_clean$CST <- sub_cst_meta$CST[match(rownames(data_mat_clean), sub_cst_meta$SeqID)]
  
  taxa_cols         <- setdiff(colnames(data_mat_clean), "CST")
  data_long_samples <- data_mat_clean %>%
    rownames_to_column("SeqID") %>%
    reshape2::melt(id.vars = c("SeqID", "CST"), variable.name = "Taxon", value.name = "Abundance") %>%
    mutate(Taxon = as.character(Taxon))
  
  taxa_global_order <- data_long_samples %>%
    group_by(Taxon) %>%
    summarise(Global_Mean = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Global_Mean)) %>%
    pull(Taxon)
  
  top6_species <- head(taxa_global_order, 6)
  
  plot_boxplot_df <- data_long_samples %>%
    filter(Taxon %in% top6_species) %>%
    mutate(
      Taxon      = factor(Taxon, levels = top6_species),
      CST        = factor(CST, levels = intersect(names(cst_colors), data_long_samples$CST)),
      Percentage = Abundance * 100
    )
  
  p_taxa_boxplot <- ggplot(plot_boxplot_df, aes(x = Taxon, y = Percentage, fill = CST, color = CST)) +
    stat_boxplot(geom = "errorbar", size = 1.2, width = 0.3, position = position_dodge(0.75)) +
    geom_boxplot(width = 0.6, size = 1.2, alpha = 0.8, outlier.shape = NA, outlier.size = 1.5, position = position_dodge(0.75)) +
    scale_fill_manual(values = cst_colors, name = "CST", drop = TRUE) +
    scale_color_manual(values = cst_border_colors, name = "CST", drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)), labels = function(x) paste0(x, "%")) +
    geom_vline(xintercept = seq(1.5, 5.5), color = "grey") +
    labs(title = sprintf("Top 6 Dominant Species - %s", group_name), x = NULL, y = "Relative Abundance (%)") +
    theme_bw() + 
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.text.x       = element_text(angle = 30, hjust = 1, face = "bold.italic", color = "black", size = 10), 
      axis.text.y       = element_text(size = 10, color = "black"),
      axis.title        = element_text(face = "bold", size = 11),
      panel.grid.minor  = element_blank(), 
      panel.grid.major  = element_blank(),
      legend.position   = "right",
      legend.title      = element_text(face = "bold")
    )
  
  print(p_taxa_boxplot)
  ggsave(sprintf("CST_%s_Dominant_Species_Boxplot.pdf", group_name), plot = p_taxa_boxplot, width = 10, height = 4)
}