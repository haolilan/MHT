# ==============================================================================
# Execution Master Script: Clinical and Microbiome Pipelines
# Year: 2026
# ==============================================================================

# 1. Environment & Dependencies Setup
base_dir   <- "D:/WorkProjects/Demo-MHT 2026"
utils_path <- file.path(base_dir, "/src/clinical_microbiome_longitudinal_utils.R")

# Dynamically link the architecture toolkit definitions
if (file.exists(utils_path)) {
  source(utils_path)
} else {
  stop("Missing utilities kit script file at designated base path.")
}

# Load the project dataset space
load(file.path(base_dir, "data/MHT.demo.RData"))
base_dir
# ==============================================================================
# Pipeline Run Task 1: Clinical Cohort Statistics Workflow
# ==============================================================================
clinic_dir <- file.path(base_dir, "Results/OverallEva/Clinic")
if(!dir.exists(clinic_dir)) dir.create(clinic_dir, recursive = TRUE)

# Array mappings for clinical metadata
non_k_vars  <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & !phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()
k_vars      <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()
comps_non_k <- list(c("BL", "T24"))
comb_matrix <- combn(c("BL", "T04", "T12", "T24"), 2)
comps_k     <- lapply(split(t(comb_matrix), 1:ncol(comb_matrix)), as.vector)

# Run Task Execution
run_clinic_longitudinal_pipeline(
  df = Phen.Seq, vars_k = k_vars, vars_non_k = non_k_vars,
  comps_k = comps_k, comps_non_k = comps_non_k,
  group_col = NULL, prefix = "Global", out_dir = clinic_dir, phen_cate = phen.Cate
)

# ==============================================================================
# Pipeline Run Task 2: Multi-Site Microbiome Workflow
# ==============================================================================
microbe_dir <- file.path(base_dir, "Results/OverallEva/Microbe")
if(!dir.exists(microbe_dir)) dir.create(microbe_dir, recursive = TRUE)

# Threshold mapping arrays and parameters
threshold_map  <- list("GUT" = 0.001, "VA" = 0.0001, "TO" = 0.001, "UR" = 0.0001)
my_comparisons <- list(c("BL","T04"), c("BL","T12"), c("BL","T24"))

# Run Task Execution
run_microbiome_pipeline(
  prof_list     = prof_filtered,
  meta_df       = Microbe.phen.prof,  
  threshold_map = threshold_map,
  comparisons   = my_comparisons,
  group_col     = NULL,
  prefix        = "Paired.sample.Global",
  out_dir       = microbe_dir
)