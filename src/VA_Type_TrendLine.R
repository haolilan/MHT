

setwd("D:/WorkProjects/Demo-MHT 2026/Results/")

stat_df <- read.csv("D:/WorkProjects/Demo-MHT 2026/Results/Responser/VA_Type/Microbe/Paired.sample.VA_Type.TimeDiff.wilcox.csv")
stat_df <- stat_df %>%subset(Site=="VA")
sig.taxa <- stat_df %>% subset(p.adj < 0.05) %>% pull("Factor") %>%unique


plot_df <- stat_df%>%subset(Factor%in% sig.taxa & Comparison%in% c("BL/T04","BL/T12","BL/T24"))


library(tidyverse)
plot_df <- stat_df%>%subset(Factor%in% sig.taxa & Comparison%in% c("BL/T04","BL/T12","BL/T24"))

bl_base  <- plot_df %>% filter(Comparison == "BL/T04") %>% select(Site, Factor, Subgroup, Value = Mean1) %>% mutate(Time_Point = "BL")
t4_base  <- plot_df %>% filter(Comparison == "BL/T04") %>% select(Site, Factor, Subgroup, Value = Mean2) %>% mutate(Time_Point = "T4")
t12_base <- plot_df %>% filter(Comparison == "BL/T12") %>% select(Site, Factor, Subgroup, Value = Mean2) %>% mutate(Time_Point = "T12")
t24_base <- plot_df %>% filter(Comparison == "BL/T24") %>% select(Site, Factor, Subgroup, Value = Mean2) %>% mutate(Time_Point = "T24")

df_points <- bind_rows(bl_base, t4_base, t12_base, t24_base) %>%
  mutate(Time_Point = factor(Time_Point, levels = c("BL", "T4", "T12", "T24")))

segment_BL_T4 <- plot_df %>%
  filter(Comparison == "BL/T04") %>%
  select(Site, Factor, Subgroup, Y_start = Mean1, Y_end = Mean2) %>%
  mutate(X_start = "BL", X_end = "T4")

segment_T4_T12 <- plot_df %>%
  filter(Comparison %in% c("BL/T04", "BL/T12")) %>%
  select(Site, Factor, Subgroup, Comparison, Mean2) %>%
  pivot_wider(names_from = Comparison, values_from = Mean2) %>%
  select(Site, Factor, Subgroup, Y_start = `BL/T04`, Y_end = `BL/T12`) %>%
  mutate(X_start = "T4", X_end = "T12")

segment_T12_T24 <- plot_df %>%
  filter(Comparison %in% c("BL/T12", "BL/T24")) %>%
  select(Site, Factor, Subgroup, Comparison, Mean2) %>%
  pivot_wider(names_from = Comparison, values_from = Mean2) %>%
  select(Site, Factor, Subgroup, Y_start = `BL/T12`, Y_end = `BL/T24`) %>%
  mutate(X_start = "T12", X_end = "T24")

df_segments <- bind_rows(segment_BL_T4, segment_T4_T12, segment_T12_T24) %>%
  drop_na(Y_start, Y_end) %>%
  mutate(
    X_start = factor(X_start, levels = c("BL", "T4", "T12", "T24")),
    X_end   = factor(X_end,   levels = c("BL", "T4", "T12", "T24"))
  )

df_text_anno <- plot_df %>%
  filter(Comparison %in% c("BL/T04", "BL/T12", "BL/T24")) %>%
  mutate(
    Clean_Label = if_else(label == "ns/ns" | is.na(label), "", as.character(label)),
    Time_Point = case_when(
      Comparison == "BL/T04" ~ "T4",
      Comparison == "BL/T12" ~ "T12",
      Comparison == "BL/T24" ~ "T24"
    ),
    Time_Point = factor(Time_Point, levels = c("BL", "T4", "T12", "T24")),
    Value = Mean2
  ) %>%
  select(Site, Factor, Subgroup, Time_Point, Value, Clean_Label) %>%
  filter(Clean_Label != "")

subgroup_palette <- c("UROG-L.c" = "#117733", "UROG-L.i" = "#74C476", "UROG-G.v" = "#FF7F0E", "UROG-Div" = "#2171B5")

p_line_trend <- ggplot() +
  geom_segment(data = df_segments, 
               aes(x = X_start, xend = X_end, y = Y_start, yend = Y_end, color = Subgroup),
               linewidth = 0.6, 
               alpha = 0.85) +
  geom_point(data = df_points, aes(x = Time_Point, y = Value, color = Subgroup), 
             size = 2.0, alpha = 0.9) +
  geom_text(data = df_text_anno, 
            aes(x = Time_Point, y = Value, label = Clean_Label, color = Subgroup),
            hjust = -0.2, 
            vjust = -0.6, 
            size = 3.2, 
            fontface = "bold", 
            show.legend = FALSE) +
  scale_color_manual(values = subgroup_palette, name = "Stratified Subgroup") +
  facet_wrap(~ Factor, scales = "free_y", nrow = 2, axes = "all") +
  theme_bw() +
  labs(
    x = "Study Longitudinal Timeline Stages",
    y = "Abundance Transformation Profile (Log Mean Value)",
    title = "Longitudinal Micro-Dynamics and Node-wise Significance Map"
  ) +
  theme(
    plot.title         = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(), 
    strip.background   = element_rect(fill = "grey95", color = "black"),
    strip.text         = element_text(size = 9, face = "bold.italic"), 
    axis.title         = element_text(face = "bold", size = 10),
    axis.text          = element_text(size = 9, color = "black"),
    legend.title       = element_text(face = "bold", size = 9),
    legend.text        = element_text(size = 8.5),
    legend.position    = "right"
  ) +
  coord_cartesian(clip = "off")

p_line_trend

unique_factors <- length(unique(df_points$Factor))
calc_height    <- max(4, ceiling(unique_factors / 2) * 2.5) 

ggsave("Microbiome_4_Points_Clean_Trend_Lines.pdf", plot = p_line_trend, width = 18, height = 6)
