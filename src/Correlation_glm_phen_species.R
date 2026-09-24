# ==============================================================================
# Pipeline: Longitudinal Multi-Site Microbiome vs Clinical Phenotypes (GLM/LM Engine)
# Modules: Dynamic Sample Alignment, CLR Slicing, Covariate Adjustments, Signif-Gating
# Year: 2026
# ==============================================================================
library(tidyverse)
library(rstatix)
library(reshape2)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

clr_transform <- function(data, pseudocount = NULL) {
  if (is.null(pseudocount)) {
    min_nonzero <- min(data[data != 0], na.rm = TRUE)
    pseudocount <- min_nonzero * 0.1
  }
  data_adj <- data + pseudocount
  geo_mean <- apply(data_adj, 1, function(x) exp(mean(log(x))))
  data_clr <- log(data_adj / geo_mean)
  return(as.data.frame(data_clr))
}

get_significance_stars <- function(p_vector) {
  dplyr::case_when(
    p_vector < 0.0001 ~ "****",
    p_vector < 0.001  ~ "***",
    p_vector < 0.01   ~ "**",
    p_vector < 0.05   ~ "*",
    TRUE              ~ "ns"
  )
}

lm_taxa_calculate <- function(prof.input, groups, vars.names) {
  n_species <- ncol(prof.input)
  
  # Linear regression: V2 = scaled clinical phenotype, V3 = scaled baseline age
  model_list <- lapply(seq_len(n_species), function(j) {
    lm(prof.input[, j] ~ groups$V2 + groups$V3)
  })
  
  p_mat <- sapply(model_list, function(m) {
    if (nrow(summary(m)$coefficients) < 2) NA else summary(m)$coefficients[2, 4]
  })
  coef_mat <- sapply(model_list, function(m) {
    if (nrow(summary(m)$coefficients) < 2) NA else summary(m)$coefficients[2, 1]
  })
  
  data.frame(
    coef             = as.vector(coef_mat),
    p.value          = as.vector(p_mat),
    species          = colnames(prof.input),
    vars             = vars.names,
    stringsAsFactors = FALSE
  )
}

base_dir <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))
Phen.Seq$Age <- as.numeric(Phen.Seq$Age)

target_dir <- "D:/WorkProjects/Demo-MHT 2026/Results/Association_Clinic"
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)

# Load significant numeric factors from PERMANOVA screens
perm.files <- list.files(file.path(base_dir, "Results/Permanova/Whole"), pattern = "permanova_strata.*\\.csv", full.names = TRUE)
perm.combined <- lapply(perm.files, function(f) {
  df <- read.csv(f, header = TRUE)
  df$source_file <- sub("\\..*$", "", basename(f)) %>% gsub("permanova_strata_", "", .)
  return(df)
}) %>% bind_rows()

factor.sig <- perm.combined %>%
  dplyr::filter(var_type == "numeric", !is.na(R2))

vars.phen_list <- list(
  "VA"  = factor.sig$factor[factor.sig$source_file == "VA"],
  "UR"  = factor.sig$factor[factor.sig$source_file == "UR"],
  "TO"  = factor.sig$factor[factor.sig$source_file == "TO"],
  "GUT" = factor.sig$factor[factor.sig$source_file == "GUT"]
)

threshold_map <- list("VA" = 0.0001, "UR" = 0.0001, "TO" = 0.001, "GUT" = 0.001)

sum.res    <- list()
sum.data   <- list()
sum.data.0 <- list()

for (time.vars in c("BL", "T24")) {
  for (sites.vars in names(threshold_map)) {
    
    rel_ab.threshold <- threshold_map[[sites.vars]]
    target_phen_vars <- unique(vars.phen_list[[sites.vars]]) %>% setdiff("Age")
    if (length(target_phen_vars) == 0) next
    
    # Require complete paired sampling across BL and T24
    paired_patients <- Microbe.phen.prof %>%
      dplyr::filter(Time %in% c("BL", "T24"), Site == sites.vars) %>%
      group_by(Clinic_ID) %>%
      dplyr::filter(n() == 2) %>%
      ungroup() %>%
      pull(Clinic_ID) %>% 
      unique()
    
    if (length(paired_patients) == 0) next
    
    mapping_slice <- Microbe.phen.prof %>%
      dplyr::filter(Time == time.vars, Site == sites.vars, Clinic_ID %in% paired_patients)
    
    raw_prof_table <- prof_filtered[[sites.vars]] %>%
      as.data.frame() %>%
      rownames_to_column("SeqID") %>%
      inner_join(mapping_slice %>% select(SeqID, Clinic_ID), by = "SeqID") %>%
      dplyr::select(-SeqID) %>%
      column_to_rownames("Clinic_ID")
    
    pheno_slice <- Phen.Seq %>%
      dplyr::filter(Time == time.vars, Clinic_ID %in% rownames(raw_prof_table))
    
    baseline_age_map <- Phen.Seq %>% 
      dplyr::filter(Time == "BL") %>% 
      select(Clinic_ID, Age)
    
    for (vars.phen in target_phen_vars) {
      
      valid_clinic <- pheno_slice %>% 
        select(Clinic_ID, all_of(vars.phen)) %>% 
        drop_na() %>%
        inner_join(baseline_age_map, by = "Clinic_ID") %>%
        na.omit()
      
      common_ids <- intersect(rownames(raw_prof_table), valid_clinic$Clinic_ID)
      if (length(common_ids) < 5) next
      
      sub_prof_mat  <- raw_prof_table[common_ids, , drop = FALSE]
      sub_clinic_df <- valid_clinic %>% dplyr::filter(Clinic_ID %in% common_ids) %>% column_to_rownames("Clinic_ID")
      sub_clinic_df <- sub_clinic_df[common_ids, ]
      
      # Filter taxa with prevalence > 10% above site-specific abundance threshold
      taxa_to_keep <- sapply(sub_prof_mat, function(x) sum(x > rel_ab.threshold) > 0.1 * nrow(sub_prof_mat))
      taxa.vars    <- names(taxa_to_keep)[taxa_to_keep]
      if (length(taxa.vars) == 0) next
      
      sub_prof_filtered <- sub_prof_mat[, taxa.vars, drop = FALSE]
      data.clr          <- clr_transform(sub_prof_filtered)
      
      groups <- data.frame(
        V1 = common_ids,
        V2 = as.vector(scale(sub_clinic_df[[vars.phen]])),
        V3 = as.vector(scale(sub_clinic_df[["Age"]])),
        stringsAsFactors = FALSE
      )
      rownames(groups) <- groups$V1
      
      raw_results <- lm_taxa_calculate(data.clr, groups, vars.phen)
      
      # Multiple testing correction (BH/FDR) per phenotype feature
      results_adjusted <- raw_results %>%
        mutate(
          p.adj    = p.adjust(p.value, method = "fdr"),
          p_star   = get_significance_stars(p.value),
          q_star   = get_significance_stars(p.adj),
          raw_mark = paste(p_star, q_star, sep = "/"),
          p_mark   = if_else(p.adj > 0.05, "", raw_mark),
          sites    = sites.vars, 
          time     = time.vars
        ) %>%
        dplyr::select(-p_star, -q_star, -raw_mark)
      
      uniq_key <- paste(time.vars, sites.vars, vars.phen, sep = "_")
      sum.res[[uniq_key]] <- results_adjusted
      
      # Store CLR-transformed and raw data for downstream plotting
      data.clr$ID <- rownames(data.clr)
      data.m <- reshape2::melt(data.clr, id.vars = "ID", variable.name = "taxa", value.name = "value.clr") %>%
        mutate(value.vars = groups$V2[match(ID, groups$V1)], vars = vars.phen, sites = sites.vars, time = time.vars)
      sum.data[[uniq_key]] <- data.m
      
      sub_prof_filtered$ID <- rownames(sub_prof_filtered)
      data.m.0 <- reshape2::melt(sub_prof_filtered, id.vars = "ID", variable.name = "taxa", value.name = "value.raw") %>%
        mutate(value.vars = sub_clinic_df[[vars.phen]][match(ID, rownames(sub_clinic_df))], vars = vars.phen, sites = sites.vars, time = time.vars)
      sum.data.0[[uniq_key]] <- data.m.0
    }
  }
}

sum.res_final <- bind_rows(sum.res) %>%
  arrange(p.value) %>%
  select(species, vars, sites, time, coef, p.value, p.adj, everything())

write.csv(sum.res_final, file.path(target_dir, "glm.out.all.index.csv"), row.names = FALSE)
