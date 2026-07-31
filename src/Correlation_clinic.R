# ==========================================
# 0. Load Environment and Data Preparation
# ==========================================
library(tidyverse)
library(openxlsx)

# Load datasets (Ensure the paths are correct)
base_dir   <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

# Set working directory
target_dir <- "D:/WorkProjects/Demo-MHT 2026/Results/Association_Clinic"
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)

# ====================================================================================
# 1. Correlation Matrix Calculation Function (Batch Operation + Global FDR)
# ====================================================================================
calc_pairwise_spearman <- function(df, vars) {
  # Generate unique pairwise combinations
  grid <- expand_grid(Var1 = vars, Var2 = vars) %>%
    mutate(
      idx1 = match(Var1, vars),
      idx2 = match(Var2, vars)
    ) %>%
    filter(idx1 <= idx2) 
  
  # Perform batch correlation calculations
  results <- map2_dfr(grid$Var1, grid$Var2, ~{
    x <- as.numeric(df[[.x]])
    y <- as.numeric(df[[.y]])
    
    # Remove samples containing NA values
    valid_idx <- complete.cases(x, y)
    x <- x[valid_idx]
    y <- y[valid_idx]
    
    if (length(x) < 3) return(NULL) # Insufficient sample size to calculate
    
    ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    
    tibble(
      Test_var = .y,
      Phenotype_var = .x,
      Method = "Spearman",
      Estimate = ct$estimate,
      P_value = ct$p.value,
      N = length(x)
    )
  })
  
  # Apply global FDR correction across all generated pairs
  results <- results %>%
    mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
    relocate(P_adj, .after = P_value)
  
  return(results)
}

calc_pairwise_spearman <- function(df, vars) {
  # 🌟 优化 1：先对输入变量去重，防止重复变量造成多重计算污染
  vars <- unique(vars)
  
  # Generate unique pairwise combinations
  grid <- expand_grid(Var1 = vars, Var2 = vars) %>%
    mutate(
      idx1 = match(Var1, vars),
      idx2 = match(Var2, vars)
    ) %>%
    # 🌟 优化 2：把 <= 改为 <，严格剔除对角线元素（即剔除自己和自己的计算），且只算单侧半矩阵（避免 A-B 和 B-A 重复）
    filter(idx1 < idx2) 
  
  # Perform batch correlation calculations
  results <- map2_dfr(grid$Var1, grid$Var2, ~{
    # 转换为数值型，防止字符型向量导致的错误
    x <- as.numeric(df[[.x]])
    y <- as.numeric(df[[.y]])
    
    # Remove samples containing NA values
    valid_idx <- complete.cases(x, y)
    x <- x[valid_idx]
    y <- y[valid_idx]
    
    if (length(x) < 3) return(NULL) # Insufficient sample size to calculate
    
    # 运行 Spearman 相关分析
    ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    
    tibble(
      Test_var = .y,
      Phenotype_var = .x,
      Method = "Spearman",
      Estimate = as.numeric(ct$estimate),
      P_value = ct$p.value,
      N = length(x)
    )
  })
  
  if (is.null(results) || nrow(results) == 0) return(NULL)
  
  # Apply global FDR correction across all generated pairs
  results <- results %>%
    mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
    relocate(P_adj, .after = P_value)
  
  return(results)
}

# ==========================================
# 2. Data Filtering and Correlation Computation
# ==========================================
# Extract target variables and append the "Delta_" prefix
target_cates <- c("Therapy_Score", "Blood_Test", "Hormone")

# 🌟 在这里也加上 unique() 防御
vars.base <- phen.Cate %>% 
  filter(PhenCategory %in% target_cates, Phen_Time == "T24") %>% 
  pull(Phen) %>% 
  unique()

vars.x <- paste0("Delta_", vars.base)

# Identify qualified patients with complete data at both BL and T24 timepoints
valid_clinic_ids <- Phen.Seq %>%
  filter(Time %in% c("BL", "T24")) %>%
  count(Clinic_ID) %>%
  filter(n == 2) %>%
  pull(Clinic_ID)

# Extract the Delta_T24 data matrix for the qualified patients
out.mat <- Phen.Delta %>% 
  filter(Delta_Time == "Delta_T24", Clinic_ID %in% valid_clinic_ids)

# Execute correlation computation
final_results <- calc_pairwise_spearman(df = out.mat, vars = vars.x)

# ==========================================
# 3. Aesthetics Mapping and Attribute Assignment
# ==========================================
if (!is.null(final_results)) {
  # Create a dictionary mapping each variable to its respective category
  cate_dict <- phen.Cate %>% select(Phen, PhenCategory) %>% distinct() %>% deframe()
  
  final_results <- final_results %>%
    mutate(
      pheno.index = str_remove(Phenotype_var, "Delta_"),
      vars.index = str_remove(Test_var, "Delta_"),
      phenCate = cate_dict[pheno.index],
      varsCate = cate_dict[vars.index],
      # Generate significance star labels based on adjusted p-values
      lab.p.adj = case_when(
        P_adj < 0.001 ~ "***",
        P_adj < 0.01 ~ "**",
        P_adj < 0.05 ~ "*",
        TRUE ~ ""
      ),
      # Generate directional marks for significant correlations
      sign_mark = case_when(
        Estimate > 0 & P_adj < 0.05 ~ "+",
        Estimate < 0 & P_adj < 0.05 ~ "-",
        TRUE ~ ""
      ),
      label = ifelse(P_adj < 0.05, paste0(sign_mark, "\n", lab.p.adj), ""),
      # Scale correlation estimates using log1p transformation for enhanced color mapping
      Estimate_log = sign(Estimate) * log1p(abs(Estimate))
    )
  
  # Export comprehensive statistical results
  write.xlsx(final_results, "spearman.delta.phen2phen.xlsx", overwrite = TRUE)
}

# ==========================================
# 4. Filter Significant Features
# ==========================================
# Identify all variables that are significant (P_adj < 0.05) in at least one distinct comparison
sig_vars <- final_results %>%
  # filter(Phenotype_var != Test_var, P_adj < 0.05) %>%
  filter( phenCate == "Therapy_Score" & varsCate== "Therapy_Score")%>%
  select(Phenotype_var, Test_var) %>%
  pivot_longer(everything()) %>%
  pull(value) %>%
  unique()

# Sort significant variables by their broad category (phenCate) to create a modular heatmap structure
ordered_vars <- final_results %>%
  filter(Phenotype_var %in% sig_vars) %>%
  select(Phenotype_var, phenCate) %>%
  distinct() %>%
  arrange(phenCate, Phenotype_var) %>%
  pull(Phenotype_var)

# Mirror and complete the dataset to guarantee a solid matrix pattern
plot_data_lower_mirror <- final_results %>%
  filter(Phenotype_var %in% sig_vars, Test_var %in% sig_vars)

plot_data_lower_complete <- bind_rows(
  plot_data_lower_mirror, 
  plot_data_lower_mirror %>% rename(Phenotype_var = Test_var, Test_var = Phenotype_var)
) %>%
  distinct(Phenotype_var, Test_var, .keep_all = TRUE) %>%
  mutate(
    # Convert factors to uniform numeric coordinates based on identical orders
    Pheno_num = as.numeric(factor(Phenotype_var, levels = ordered_vars)),
    Test_num  = as.numeric(factor(Test_var, levels = ordered_vars))
  ) %>%
  # Filter strict lower triangle using numeric coordinates (Pheno_num >= Test_num)
  filter(Pheno_num >= Test_num) %>%
  mutate(
    Phenotype_var = factor(Phenotype_var, levels = ordered_vars),
    Test_var      = factor(Test_var, levels = rev(ordered_vars)) # Invert Y-axis for standard visualization
  )

# Plot clean lower-triangle heatmap without structural blank spaces
p_heat_lower <- ggplot(plot_data_lower_complete, aes(x = Phenotype_var, y = Test_var, fill = Estimate_log)) +
  geom_tile(color = "white", linewidth = 0.5) +  
  geom_text(aes(label = sign_mark), size = 3, vjust = -0.2, color = "grey20") +  
  geom_text(aes(label = lab.p.adj), size = 3, vjust = 1.2, color = "grey20") +  
  scale_fill_gradient2(
    low = "#6A0DAD",     
    mid = "#F0F0F0",       
    high = "#b35806",     
    midpoint = 0,          
    name = "Log-estimate"   
  ) +
  scale_x_discrete(position = "top") +  
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", margin = margin(b = 20))
  ) +
  labs(x = NULL, y = NULL, title = "Delta_K Score (Spearman)")

ggsave("Sfigx.Spearman.Delta_KScore.pdf", plot = p_heat_lower, width = 11, height = 9)

# ==========================================
# 5. Filter Cross-Category Significant Features
# ==========================================
# 1. Identify all variables significant in at least one cross-category comparison (phenCate != varsCate)
cross_sig_vars <- final_results %>%
  filter(phenCate != varsCate, P_adj < 0.05) %>%
  select(Phenotype_var, Test_var) %>%
  pivot_longer(everything()) %>%
  pull(value) %>%
  unique()

# 2. Sort selected variables by their broad category to structure the upper heatmap modularly
ordered_vars_upper <- final_results %>%
  filter(Phenotype_var %in% cross_sig_vars) %>%
  select(Phenotype_var, phenCate) %>%
  distinct() %>%
  arrange(phenCate, Phenotype_var) %>%
  pull(Phenotype_var)

# 3. Extract plotting data subset and create a strict upper-triangular layout (including diagonal line)
plot_data_upper <- final_results %>%
  filter(Phenotype_var %in% cross_sig_vars, Test_var %in% cross_sig_vars)

plot_data_upper_mirror <- plot_data_upper %>%
  rename(Phenotype_var = Test_var, Test_var = Phenotype_var)

plot_data_upper_complete <- bind_rows(plot_data_upper, plot_data_upper_mirror) %>%
  distinct(Phenotype_var, Test_var, .keep_all = TRUE) %>%
  mutate(
    # Use matching factor structures to compute safe coordinates
    Pheno_num = as.numeric(factor(Phenotype_var, levels = ordered_vars_upper)),
    Test_num  = as.numeric(factor(Test_var, levels = ordered_vars_upper))
  ) %>%
  # Filter strict upper triangle (Pheno_num <= Test_num) to eliminate missing steps or gaps
  filter(Pheno_num <= Test_num) %>%
  mutate(
    Phenotype_var = factor(Phenotype_var, levels = ordered_vars_upper),
    Test_var      = factor(Test_var, levels = rev(ordered_vars_upper))
  )

# 4. Plot clean upper-triangle heatmap with synchronized axis properties
p_heat_upper <- ggplot(plot_data_upper_complete, aes(x = Phenotype_var, y = Test_var, fill = Estimate_log)) +
  geom_tile(color = "white", linewidth = 0.5) +  
  geom_text(aes(label = sign_mark), size = 3, vjust = -0.2, color = "grey20") +  
  geom_text(aes(label = lab.p.adj), size = 3, vjust = 1.2, color = "grey20") +  
  scale_fill_gradient2(
    low = "#6A0DAD",     
    mid = "#F0F0F0",       
    high = "#b35806",     
    midpoint = 0,          
    name = "Log-estimate",
    na.value = "white"
  ) +
  scale_x_discrete(position = "bottom") + 
  scale_y_discrete(position = "left") +
  theme_minimal() +
  theme(
    axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y.left = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", margin = margin(b = 20))
  ) +
  labs(x = NULL, y = NULL, title = "Delta_K Score ~ Hormone ~ Metabolite (Spearman)")

ggsave("Fig1.d.Spearman.Delta_KScore_Hormone_Metabolite.intra.pdf", plot = p_heat_upper, width = 11, height = 9)
