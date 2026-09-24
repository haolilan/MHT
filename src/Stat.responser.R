library(ggplot2)
library(tidyverse)
library(ggpubr) 

base_dir   <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

setwd("D:/WorkProjects/Demo-MHT 2026/Results/Responser")

# MHT_Response_Group  ##################
plot_stacked_percentage <- function(data, y_var, fill_var, y_label, fill_legend_title, colors) {
  
  plot_df <- data %>%
    select(Y = all_of(y_var), Fill = all_of(fill_var)) %>%
    na.omit() %>%
    mutate(Y = as.factor(Y), Fill = as.factor(Fill))
  
  total_n <- nrow(plot_df)
  
  chisq_res <- chisq.test(table(plot_df$Y, plot_df$Fill))
  p_val <- chisq_res$p.value
  p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  test_label <- paste0("Chi-squared test, ", p_label, " (Total N = ", total_n, ")")
  
  anno_df <- plot_df %>%
    group_by(Y, Fill) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(Y) %>%
    mutate(
      total = sum(count),
      percentage = (count / total) * 100,
      label_text = paste0(sprintf("%.1f%%", percentage), " (", count, ")")
    )
  
  p <- ggplot(anno_df, aes(x = count, y = Y, fill = Fill)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) +
    geom_text(aes(label = label_text), 
              position = position_fill(vjust = 0.5), 
              color = "white", size = 3, fontface = "bold") +
    scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = colors, name = fill_legend_title) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_blank(),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(t = 10, r = 18, b = 10, l = 10, unit = "pt")
    ) +
    labs(
      x = "Proportion",
      y = y_label,
      caption = test_label 
    )
  
  return(p)
}


final_clinical_data <- merge(
  merge(phen.gr.E2.all %>% mutate(Hormone_E2 = recode_factor(Group,"NR" = "H_NR","R"  = "H_R")) %>% select(-Group),
        phen.gr.K      %>% mutate(mKI = recode_factor(Group,"NR" = "mKI_NR","R"  = "mKI_R"))%>% select(-Group),all = T),
  phen.gr.gut    %>% mutate(Gut_Type = recode_factor(CST, "GUT-P.cop"="GUT-P.cop","GUT-P.v"="GUT-P.v"))%>% select(-CST),all = T)

dim(final_clinical_data)

colors_set1 <- c("H_R" = "#2ecc71", "H_NR" = "#d35400")          
colors_set2 <- c("mKI_R" = "#1a5276", "mKI_NR" = "#bdc3c7") 
colors_set3 <- c("GUT-P.v"   = "#66c2a4", "GUT-P.cop" = "#b2e2e2")

fig_A <- plot_stacked_percentage(
  data = final_clinical_data, 
  y_var = "Gut_Type", fill_var = "Hormone_E2", 
  y_label = "Baseline Gut Type", fill_legend_title = "Hormonal Response", 
  colors = colors_set1
)
fig_A

fig_B <- plot_stacked_percentage(
  data = final_clinical_data, 
  y_var = "Gut_Type", fill_var = "mKI", 
  y_label = "Baseline Gut Type", fill_legend_title = "Clinical Response", 
  colors = colors_set2
)
fig_B

final_merged_figure <- ggarrange(
  fig_A, fig_B,
  ncol = 1, nrow = 2,align = "hv"
)

ggsave("MHT_Response_Group_Tri-Matrix.pdf", plot = final_merged_figure, 
       width = 6, height = 5)


# Change.dots.Hormone_E2 ###############
plot_trajectory_1p <- function(data, var_prefix, group_col, title_suffix, endpoint = "T24") {
  bl_col <- paste0(var_prefix, "_BL")
  ep_col <- paste0(var_prefix, "_", endpoint)
  
  plot_df <- data %>%
    select(Clinic_ID, BL = !!sym(bl_col), Endpoint = !!sym(ep_col), Group = !!sym(group_col)) %>%
    filter(!is.na(BL) & !is.na(Endpoint) & !is.na(Group))
  
  group_counts <- plot_df %>% count(Group)
  caption_text <- paste("Sample Size:", paste(group_counts$Group, "=", group_counts$n, collapse = " | "))
  
  plot_df <- plot_df %>%
    group_by(Group) %>%
    arrange(desc(BL)) %>%
    mutate(Patient_Idx = row_number()) %>%
    ungroup()
  
  is_positive_indicator <- var_prefix == "FSFI"
  
  p <- ggplot(plot_df) + 
    geom_segment(aes(x = Patient_Idx, xend = Patient_Idx, y = BL, yend = Endpoint), 
                 color = "grey70", alpha = 0.5, linewidth = 0.2,
                 arrow = arrow(length = unit(0.15, "cm"), ends = "last", type = "closed")) +
    geom_point(aes(x = Patient_Idx, y = BL), color = "#4477AA", size = 1.2, alpha = 0.8) +
    geom_point(aes(x = Patient_Idx, y = Endpoint), color = "#FF7F0E", size = 1.2, alpha = 0.8) +
    facet_grid(~ Group, scales = "free_x", space = "free") +
    labs(
      title = paste(var_prefix, "-", title_suffix),
      x = "Individual Patients (Ordered by Baseline)", 
      y = paste(var_prefix, ""),
      caption = paste(caption_text, "(Blue dot: Baseline | Orange dot:", endpoint,")")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      plot.caption = element_text(hjust = 1, face = "italic", color = "black", size = 11),  
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "white", color = "white")
    )
  return(p)
}
####

final_clinical_data <- merge(
  merge(phen.gr.E2.all %>% mutate(Hormone_E2 = recode_factor(Group,"NR" = "H_NR","R"  = "H_R")) %>% select(-Group),
        phen.gr.K      %>% mutate(mKI = recode_factor(Group,"NR" = "mKI_NR","R"  = "mKI_R"))%>% select(-Group),all = T),
  phen.gr.gut    %>% mutate(Gut_Type = recode_factor(CST, "GUT-P.cop"="GUT-P.cop","GUT-P.v"="GUT-P.v"))%>% select(-CST),all = T)

phen_merged <- phen.enroll %>% select(Clinic_ID,K_Score_BL,K_Score_T04,E2_BL,E2_T24,FSH_BL,FSH_T24)%>%
  left_join(final_clinical_data)

p_E2 <- plot_trajectory_1p(phen_merged, "E2", "Hormone_E2", "Hormonal Response")
ggsave("Change.dots.Hormone_E2.pdf", plot = p_E2, width = 12, height = 5)
write.csv(p_E2$data%>%select(-Clinic_ID), "Change.dots.Hormone_E2.csv", row.names = FALSE)


## Scatter of Delta_K_Score for Hormone_E2 #####################
my_com <- list(c("H_NR", "H_R"))

df <-  p_E2$data %>% mutate(T24 = Endpoint )%>%
  merge(Phen.Delta%>%subset(Delta_Time=="Delta_T24")%>%select(Clinic_ID,Delta_K_Score))%>%
  mutate(Delta_mKI = Delta_K_Score)

p1 <- ggplot(df, aes(x = BL, y = T24)) +
  geom_point(
    aes(
      fill = -Delta_mKI,  
      size = -Delta_mKI,  
      color = Group         
    ),
    shape = 21,             
    stroke = 1.1            
  ) +
  scale_fill_gradient(
    high = "#0072B2",          
    low = "white",          
    name = "-Delta_mKI"  
  ) +
  scale_color_manual(
    values = c("H_NR" = "#E69F00", "H_R" = "grey40"),
    name = "Group"
  ) +
  scale_size_continuous(range = c(1, 6)) +
  labs(
    title = "Scatter of Delta_mKI for E2 Group",
    x = "E2 at BL",
    y = "E2 atT24"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 12)
  )

pbl<- ggplot(df, aes(x = Group, y = BL,color=Group)) + 
  geom_boxplot() + geom_jitter(width = 0.2) +
  scale_color_manual(values = c("H_NR" = "#E69F00", "H_R" = "grey40"), name = "Group") +
  labs(x = "",y = ""  ) +
  coord_flip() + theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = my_com,
    tip.length = NA,
    label = "p.format" 
  )

p24<- ggplot(df, aes(x = Group, y = T24,color=Group)) + 
  scale_color_manual(values = c("H_NR" = "#E69F00", "H_R" = "grey40"), name = "Group") +
  labs(x = "",y = ""  ) +
  geom_boxplot() + geom_jitter(width = 0.2) + theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = my_com,
    tip.length = NA,
    label = "p.format" 
  )

ggarrange(p24, p1, NULL, pbl, align = "hv", widths = c(1.5, 3), heights = c(3, 1.2))

ggsave("Scatter of Delta_K_Score for Hormone_E2.pdf",width = 10,height = 8)

write.csv(p1$data%>%select(BL,T24,Group,Patient_Idx,Delta_mKI),"Scatter of Delta_K_Score for Hormone_E2.csv",quote = F)



# mKI.new with Hormone_E2 and mKI ###############################
library(tidyverse)

var_prefix <- "K_Score"
group_col  <- "mKI"
title_suffix <- "Clinical Response"
endpoint   <- "T04"

final_clinical_data <- merge(
  merge(
    phen.gr.E2.all %>% 
      mutate(Hormone_E2 = recode_factor(Group, "NR" = "H_NR", "R" = "H_R")) %>% 
      select(-Group),
    phen.gr.K %>% 
      mutate(mKI = recode_factor(Group, "NR" = "mKI_NR", "R" = "mKI_R")) %>% 
      select(-Group), 
    all = TRUE
  ),
  phen.gr.gut %>% 
    mutate(Gut_Type = recode_factor(CST, "GUT-P.cop" = "GUT-P.cop", "GUT-P.v" = "GUT-P.v")) %>% 
    select(-CST), 
  all = TRUE
)

phen_merged <- phen.enroll %>% 
  select(Clinic_ID, K_Score_BL, K_Score_T04, E2_BL, E2_T24, FSH_BL, FSH_T24) %>%
  left_join(final_clinical_data, by = "Clinic_ID")

bl_col <- paste0(var_prefix, "_BL")
ep_col <- paste0(var_prefix, "_", endpoint)

plot_df <- phen_merged %>%
  select(
    Clinic_ID, 
    BL         = !!sym(bl_col), 
    Endpoint   = !!sym(ep_col), 
    Group      = !!sym(group_col), 
    Hormone_E2 = Hormone_E2,
    Gut_Type   = Gut_Type
  ) %>%
  filter(!is.na(BL) & !is.na(Endpoint) & !is.na(Group)) %>%
  mutate(
    Hormone_E2_str = ifelse(is.na(Hormone_E2), "NA", as.character(Hormone_E2)),
    Gut_Type_str   = ifelse(is.na(Gut_Type), "NA", as.character(Gut_Type))
  )

group_counts <- plot_df %>% count(Group)
caption_text <- paste("Sample Size:", paste(group_counts$Group, "=", group_counts$n, collapse = " | "))

plot_df <- plot_df %>%
  group_by(Group) %>%
  arrange(desc(BL)) %>%
  mutate(Patient_Idx = row_number()) %>%
  ungroup()

e2_colors  <- c("H_R" = "#2ecc71", "H_NR" = "#d35400")
gut_colors <- c("GUT-P.v" = "#66c2a4", "GUT-P.cop" = "#b2e2e2")
all_fills  <- c(e2_colors, gut_colors, "NA" = "grey30")

p <- ggplot(plot_df) + 
  geom_segment(
    aes(
      x = Patient_Idx, xend = Patient_Idx, 
      y = BL, yend = Endpoint
    ), 
    color = "grey",
    alpha = 0.6, 
    linewidth = 0.3,
    arrow = arrow(length = unit(0.12, "cm"), ends = "last", type = "closed")
  ) +
  geom_point(aes(x = Patient_Idx, y = BL), color = "#4477AA", size = 1.2, alpha = 0.8) +
  geom_point(aes(x = Patient_Idx, y = Endpoint), color = "#FF7F0E", size = 1.2, alpha = 0.8) +
  geom_tile(
    data = plot_df,
    aes(x = Patient_Idx, y = -1, fill = Hormone_E2_str),
    height = 0.8,       
    color = "white",    
    linewidth = 0.1     
  ) +
  geom_tile(
    data = plot_df,
    aes(x = Patient_Idx, y = -2, fill = Gut_Type_str),
    height = 0.8,       
    color = "white",    
    linewidth = 0.1
  ) +
  
  facet_grid(~ Group, scales = "free_x", space = "free") +
  
  scale_color_manual(
    name = "Hormone E2 (Trajectory)", 
    values = e2_colors, 
    na.value = "grey70"
  ) +
  
  scale_fill_manual(
    name = "Annotation Tracks", 
    values = all_fills,
    breaks = c("H_R", "H_NR", "GUT-P.v", "GUT-P.cop", "NA"),
    labels = c("H_R", "H_NR", "GUT-P.v", "GUT-P.cop", "Missing (NA)")
  ) +
  labs(
    title = paste(var_prefix, "-", title_suffix),
    x = "Individual Patients (Ordered by Baseline)", 
    y = var_prefix,
    caption = paste(caption_text, "(Light dot: Baseline | Dark dot:", endpoint, ")")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(hjust = 0.5, color = "grey70"),
    plot.caption       = element_text(hjust = 1, face = "italic", color = "black", size = 11),  
    axis.text.x        = element_blank(),  
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold"),
    strip.background   = element_rect(fill = "white", color = "white"),
    legend.position    = "bottom",
    legend.box         = "vertical"
  )

print(p)
ggsave("Change.dots.mKI.new.pdf",width = 15,height = 6)
write.csv(plot_df%>%select(-Clinic_ID,-Hormone_E2_str:-Gut_Type_str),"Change.dots.mKI.new.csv",quote = F)