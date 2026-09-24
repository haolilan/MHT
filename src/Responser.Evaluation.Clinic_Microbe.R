# ==============================================================================
# Execution Master Script: Stratified Subgroup Analysis Pipeline
# Modules: Stratified Clinical & Stratified Microbiome Pipelines
# Year: 2026
# ==============================================================================

# 1. Environment & Utilities Initialization
base_dir   <- "D:/WorkProjects/Demo-MHT 2026"
utils_path <- file.path(base_dir, "/src/clinical_microbiome_longitudinal_utils.R")

if (file.exists(utils_path)) {
  source(utils_path)
} else {
  stop("Infrastructure toolkit script missing at target base path.")
}

load(file.path(base_dir, "data/MHT.demo.RData"))

# Define shared meta-subgroup mapping configuration matrices
subgroups_config <- list(
  "Hormone_E2" = phen.gr.E2.all %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "mKI"        = phen.gr.K      %>% select(Clinic_ID, group = Group) %>% mutate(group = factor(group, levels = c("R", "NR"))),
  "Gut_Type"   = phen.gr.gut    %>% select(Clinic_ID, group = CST)   %>% mutate(group = factor(group, levels = c("GUT-P.v", "GUT-P.cop"))),
  "To_Type"    = phen.gr.to     %>% select(Clinic_ID, group = CST)   %>% mutate(group = factor(group, levels = c("TO-P.m", "TO-N.s"))),
  "VA_Type"    = phen.gr.va     %>% select(Clinic_ID, group = CST)   %>% mutate(group = factor(group, levels = c("UROG-L.c", "UROG-L.i", "UROG-G.v", "UROG-Div")))
)

# ==============================================================================
# Task 1: Stratified Clinical Subgroup Pipeline
# ==============================================================================
# Establish required variables and parameter frameworks
non_k_vars  <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & !phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()
k_vars      <- phen.Cate$Phen[phen.Cate$VarsType == "continuous_vars" & phen.Cate$Category %in% c("K Score", "K Sub Score")] %>% unique()
comps_non_k <- list(c("BL", "T24"))
comb_matrix <- combn(c("BL", "T04", "T12", "T24"), 2)
comps_k     <- lapply(split(t(comb_matrix), 1:ncol(comb_matrix)), as.vector)


for (sg_name in names(subgroups_config)) {
  target_dir_clinic <- paste0(base_dir, "/Results/Responser/", sg_name, "/Clinic/")
  if(!dir.exists(target_dir_clinic)) dir.create(target_dir_clinic, recursive = TRUE)

  if (sg_name == "VA_Type") {
    comb_matrix_va <- combn(c("BL", "T24"), 2)
    current_comps_k <- lapply(split(t(comb_matrix_va), 1:ncol(comb_matrix_va)), as.vector)
  } else {
    current_comps_k <- comps_k
  }
  # ============================================================================
  current_gr_df <- subgroups_config[[sg_name]]
  gr_col_name   <- setdiff(colnames(current_gr_df), "Clinic_ID")[1]
  
  df_merged <- Phen.Seq %>% 
    select(-any_of(c("group", gr_col_name))) %>% 
    left_join(current_gr_df %>% select(Clinic_ID, all_of(gr_col_name)), by = "Clinic_ID") %>%
    rename(group = !!sym(gr_col_name)) %>%
    filter(!is.na(group)) 

  run_clinic_longitudinal_pipeline(
    df = df_merged, vars_k = k_vars, vars_non_k = non_k_vars,
    comps_k = current_comps_k, comps_non_k = comps_non_k, 
    group_col = "group", prefix = sg_name, out_dir = target_dir_clinic, phen_cate = phen.Cate
  )
}



# ==============================================================================
# Task 2: Stratified Microbiome Subgroup Pipeline
# ==============================================================================
threshold_map  <- list("GUT" = 0.001, "VA" = 0.0001, "TO" = 0.001, "UR" = 0.0001)
my_comparisons <- list(c("BL", "T04"), c("BL", "T12"), c("BL", "T24"))

for (sg_name in names(subgroups_config)) {
  target_dir_mg <- paste0(base_dir,"/Results/Responser/",sg_name,"/Microbe/")
  if(!dir.exists(target_dir_mg)) dir.create(target_dir_mg, recursive = TRUE)
  
  current_gr_df <- subgroups_config[[sg_name]]
  gr_col_name   <- setdiff(colnames(current_gr_df), "Clinic_ID")[1]
  
  mg_meta_merged <- Microbe.phen.prof %>%
    select(-any_of(c("group", gr_col_name))) %>% 
    left_join(current_gr_df %>% select(Clinic_ID, all_of(gr_col_name)), by = "Clinic_ID") %>%
    rename(group = !!sym(gr_col_name)) 
  
  run_microbiome_pipeline(
    prof_list = prof_filtered, meta_df = mg_meta_merged,
    threshold_map = threshold_map, comparisons = my_comparisons,
    group_col = "group", prefix = paste0("Paired.sample.", sg_name),
    out_dir = target_dir_mg, valid_4_if = FALSE
  )
}

