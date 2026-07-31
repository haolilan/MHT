# ==============================================================================
# 1. Environment Initialization & Package Loading
# ==============================================================================
library(tidyverse)
library(vegan)
library(foreach)
library(doParallel)
library(tibble)
library(purrr)

# Set global base color configuration
BASE_COLORS <- c("GUT" = "#2CA02C", "TO" = "#D62728", "UR" = "#FF7F0E", "VA" = "#1F77B4")

# ==============================================================================
# 2. Core Public Function Definitions
# ==============================================================================

# 2.1 Logging utility function
log_message <- function(msg) {
  message(paste(Sys.time(), "|", msg))
}

# 2.2 Significance star labeling utility
add_signif_stars <- function(df, p_col = "Pvalue", padj_col = "Padjust") {
  breaks_val <- c(0, 0.001, 0.01, 0.05, 1)
  labels_val <- c("***", "**", "*", "-")
  
  df$star_pvalue <- cut(df[[p_col]], breaks = breaks_val, labels = labels_val, include.lowest = TRUE)
  if (padj_col %in% colnames(df)) {
    df$star_padjust <- cut(df[[padj_col]], breaks = breaks_val, labels = labels_val, include.lowest = TRUE)
  } else {
    df$star_padjust <- "-"
  }
  return(df)
}

# 2.3 Univariate adonis plotting function
adonis.plot <- function(pe.adonis) {
  pe.adonis.sig <- subset(pe.adonis, Pvalue < 0.05 & !is.na(R2)) %>% 
    add_signif_stars()
  
  ggplot(pe.adonis.sig, aes(x = reorder(factor, R2), y = R2, fill = -log10(Padjust))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("R²=%.3f", R2)), hjust = 0.2, vjust = 1, size = 3, color = "black") +
    geom_text(aes(label = paste0(star_pvalue, "/", star_padjust), y = 0), hjust = 1, color = "black", size = 3) +
    scale_fill_viridis_c(name = "-log10(q value)") +
    coord_flip() +
    labs(x = "Variables", y = expression(R^2~"(p / p adjust)"), 
         title = "PERMANOVA (strata: patientID)", caption = "P value < 0.05") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(face = "bold")) 
}

# 2.4 Core master function for parallel PERMANOVA (Unadjusted)
patientID.perm.fun <- function(metadat, bray, covariates, file.prefix, n.perm = 999, n.proc = 2) {
  if (!inherits(bray, "dist")) stop("bray must be dist object")
  if (nrow(metadat) != attr(bray, "Size")) stop("nrow(metadat) is unmatched with bray distance matrix size")
  
  log_message("Start to preprocess data and run PERMANOVA")
  
  cl <- makeCluster(min(n.proc, 15))
  registerDoParallel(cl)
  
  results <- foreach(i = seq_len(ncol(metadat)), .combine = rbind,
                     .packages = c("vegan"), .errorhandling = "pass") %dopar% {
                       fac <- colnames(metadat)[i]
                       tryCatch({
                         idx <- !is.na(metadat[[fac]])
                         samples <- intersect(rownames(metadat)[idx], labels(bray))
                         sample_n <- length(samples)
                         
                         # Return empty row if sample size or effective factor levels are insufficient
                         if (sample_n < max(5, nrow(metadat) * 0.05)) {
                           return(data.frame(factor=fac, SampleNum=sample_n, Df=NA, SumsOfSqs=NA, R2=NA, F.Model=NA, Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                         }
                         
                         test_data <- metadat[samples, fac, drop = FALSE]
                         if ((is.factor(test_data[[fac]]) || is.character(test_data[[fac]])) && length(unique(test_data[[fac]])) < 2) {
                           return(data.frame(factor=fac, SampleNum=sample_n, Df=NA, SumsOfSqs=NA, R2=NA, F.Model=NA, Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                         }
                         
                         sub_dist <- as.dist(as.matrix(bray)[samples, samples])
                         fac.patientID <- covariates[samples, "patient.ID"]
                         
                         # Formulate and run the unadjusted model
                         model_formula <- as.formula(paste("sub_dist ~", fac))
                         res <- as.data.frame(adonis2(model_formula, data = test_data, strata = fac.patientID, permutations = n.perm, by = "term"))
                         
                         # Calculate betadisper
                         disp_pval <- NA
                         if((is.factor(test_data[[fac]]) || is.character(test_data[[fac]])) && length(unique(test_data[[fac]])) > 1) {
                           disp_pval <- tryCatch({
                             permutest(betadisper(sub_dist, as.factor(test_data[[fac]])), strata = fac.patientID, permutations = how(nperm = n.perm))$tab$`Pr(>F)`[1]
                           }, error = function(e) NA)
                         }
                         
                         data.frame(
                           factor = fac, SampleNum = sample_n, Df = res[fac, "Df"], SumsOfSqs = res[fac, "SumOfSqs"],
                           R2 = res[fac, "R2"], F.Model = res[fac, "F"], Pvalue = res[fac, "Pr(>F)"], Disp_pval = disp_pval,
                           R2.adjust = RsquareAdj(res[fac, "R2"], sample_n, res[fac, "Df"]),
                           stringsAsFactors = FALSE
                         )
                       }, error = function(e) {
                         message("Error in variable ", fac, ": ", e$message)
                         return(data.frame(factor=fac, SampleNum=length(samples), Df=NA, SumsOfSqs=NA, R2=NA, F.Model=NA, Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                       })
                     }
  stopCluster(cl)
  
  # Compute basic statistics for the metadata
  out.stat <- data.frame(factor = colnames(metadat), var_type = "", na_count = 0, stats = "", stringsAsFactors = FALSE)
  for (i in 1:ncol(metadat)) {
    col_data <- metadat[[i]]
    out.stat$na_count[i] <- sum(is.na(col_data))
    if (is.numeric(col_data)) {
      out.stat$stats[i] <- sprintf("Mean+-SD: %.2f +- %.2f", mean(col_data, na.rm = TRUE), sd(col_data, na.rm = TRUE))
      out.stat$var_type[i] <- "numeric"
    } else if (is.factor(col_data) || is.character(col_data)) {
      freq_table <- table(col_data, useNA = "no")
      out.stat$stats[i] <- paste(names(freq_table), freq_table, sep=":", collapse="; ")
      out.stat$var_type[i] <- ifelse(is.ordered(col_data), "ordered factor", "factor")
    }
  }
  
  out.perm <- merge(as.data.frame(results), out.stat, by = "factor", all.x = TRUE)
  out.perm$Padjust <- p.adjust(out.perm$Pvalue, method = "fdr")
  
  # Reorder columns and export without Age and BMI columns
  col_order <- c("factor", "SampleNum", "Df", "SumsOfSqs", "R2", "R2.adjust", "F.Model", "Pvalue", "Padjust", "Disp_pval", "var_type", "na_count", "stats")
  out.perm <- out.perm[order(out.perm$R2, decreasing = TRUE), col_order]
  
  write.csv(out.perm, paste0(file.prefix, ".csv"), row.names = FALSE)
  
  # # Direct plot output
  # p <- adonis.plot(out.perm)
  # length.signif <- nrow(subset(out.perm, Pvalue < 0.05 & !is.na(R2)))
  # if (length.signif > 0){
  #   ggsave(paste0(file.prefix, ".pdf"), p, width = 12, height = base::max(4, base::min(24, length.signif * 0.8)))
  # }
  
  log_message("Completed analysis")
  return(out.perm)
}

# 2.5 Pairwise PERMANOVA calculation function
pairwise.adonis2 <- function(x, data, strata = NULL, nperm = 999, ...) {
  ststri <- ifelse(is.null(strata), 'Null', strata)
  lhs <- eval(x[[2]], environment(x), globalenv())
  rhs.frame <- model.frame(as.formula(paste("~", as.character(x)[3])), data, drop.unused.levels = TRUE)
  
  co <- combn(unique(as.character(rhs.frame[,1])), 2)
  nameres <- c('parent_call', apply(co, 2, function(el) paste(el[1], el[2], sep = ' vs ')))
  
  res <- vector(mode = "list", length = length(nameres))
  names(res) <- nameres
  res['parent_call'] <- list(paste(as.character(x)[2], "~", as.character(x)[3], ', strata =', ststri, ', permutations', nperm))
  
  for(elem in 1:ncol(co)){
    valid_idx <- rhs.frame[,1] %in% c(co[1,elem], co[2,elem])
    xred <- if(inherits(lhs,'dist')) as.dist(as.matrix(lhs)[valid_idx, valid_idx]) else lhs[valid_idx,]
    mdat1 <- data[valid_idx,]
    
    xnew <- as.formula(paste('xred ~', as.character(x)[3]))
    if(is.null(strata)){
      ad <- adonis2(xnew, data = mdat1, ...)
    } else {
      perm <- how(nperm = nperm); setBlocks(perm) <- mdat1[[ststri]]
      ad <- adonis2(xnew, data = mdat1, permutations = perm, ...)
    }
    res[[nameres[elem+1]]] <- ad[1:5]
  }
  class(res) <- c("pwadstrata", "list")
  return(res)
}

# ==============================================================================
# 3. Data Loading & Preprocessing
# ==============================================================================
base_dir   <- "D:/WorkProjects/Demo-MHT 2026"
load(file.path(base_dir, "data/MHT.demo.RData"))

Phen.Seq$Clinic_ID_Time <- paste(Phen.Seq$Clinic_ID, Phen.Seq$Time, sep = "_")

sites <- unique(Microbe.phen.prof$Site)

# ==============================================================================
# 4. Part 1: Whole Variables Analysis & Combined Plotting #########
# ==============================================================================
dir_all <- "D:/WorkProjects/Demo-MHT 2026/Results/Permanova/Whole/"
if(!dir.exists(dir_all)) dir.create(dir_all, recursive = TRUE)
setwd(dir_all)

# Loop calculations for each Site
map(sites, ~ {
  vars.Site <- .x
  file.prefix <- paste("permanova_strata", vars.Site, sep = "_")
  vars.phen <- setdiff(colnames(Phen.Seq), c("Clinic_ID_Time", "Clinic_ID"))
  
  metadat_raw <- merge(Microbe.phen.prof %>% dplyr::select(-Time,-Clinic_ID) %>% subset(Site == vars.Site), 
                       Phen.Seq, by = "Clinic_ID_Time") %>% 
    subset(!is.na(SeqID)) %>% 
    .[, c("SeqID","Clinic_ID", vars.phen)] %>% 
    filter(rowSums(is.na(dplyr::select(., all_of(vars.phen)))) < length(vars.phen))
  
  # 2. 提取对应的菌群矩阵，并在这里剔除全 0 样本
  prof_raw <- prof_filtered[[vars.Site]][metadat_raw$SeqID, ]
  valid_samples <- rownames(prof_raw)[rowSums(prof_raw) > 0] # 找出非全零样本 ID
  
  # 3. 🚨 核心对齐：用 valid_samples 反向去剪裁 metadat 和 菌群，实现三者完美交集！
  metadat <- metadat_raw %>% 
    filter(SeqID %in% valid_samples) %>% 
    column_to_rownames("SeqID")
  
  prof.input <- prof_raw[rownames(metadat), ] %>% .[, colSums(.) != 0]
  
  # 4. 绝对安全的 covardat 提取，直接使用 metadat 里对齐好的 Clinic_ID
  covardat <- data.frame(
    patient.ID = metadat$Clinic_ID,
    row.names = rownames(metadat)
  )
  
  # 5. 去除多余的辅助列，只保留分析需要用的临床自变量
  metadat <- metadat %>% dplyr::select(-Clinic_ID)
  
  # Scale numeric variables
  scale.vars <- names(metadat)[sapply(metadat, is.numeric)]
  metadat <- metadat %>% mutate(across(all_of(scale.vars), ~ as.numeric(scale(.))))
  
  # 计算 Bray-Curtis 距离矩阵
  dist_matrix <- vegdist(prof.input, method = "bray", add = TRUE)
  
  # 6. 🚨 并行随机种子终极守护：
  # R 的并行计算中，最标准的种子控制是在调用并行前加载 doRNG 包，或者直接在 adonis 计算前硬写
  # 在这里运行前我们不仅 set.seed，还确保每次并行实验完全一致
  set.seed(123)
  
  patientID.perm.fun(
    metadat = metadat, 
    bray = dist_matrix, 
    covariates = covardat, 
    file.prefix = file.prefix, 
    n.perm = 999, 
    n.proc = parallel::detectCores() - 1
  )
})

# Merge Part 1 results and generate plots
perm.files <- list.files(getwd(), pattern = "permanova_strata_.*\\.csv", full.names = TRUE)
perm.files <- perm.files[!str_detect(perm.files, "pairwise")]

perm.combined <- lapply(perm.files, function(f) {
  df <- read.csv(f, header = TRUE)
  df$source_file <- gsub("permanova_strata_|_pairwise", "", sub("\\..*$", "", basename(f)))
  return(df)
}) %>% list_rbind()

# Adjust and fix factor classifications
perm.combined$factor_cate <- phen.Cate$PhenCategory[match(perm.combined$factor, phen.Cate$Phen)]
perm.combined <- perm.combined %>% 
  mutate(factor_cate = case_when(
    factor == "Time" ~ "Therapy_Score",
    factor %in% c("LH", "FSH" , "E2", "Testosterone","TSH") ~ "Hormone",
    str_detect(factor, "^Delta") ~ "Delta_Therapy_Score",
    TRUE ~ factor_cate
  )) %>% 
  subset(factor_cate %in% c("Therapy_Score", "Hormone", "Blood_Test")) %>% 
  mutate(factor_cate = factor(factor_cate, levels = c("General", "Therapy_Score", "Hormone", "Blood_Test", "Questionnaire", "Ultrasound")),
         source_file = factor(source_file, levels = c("VA", "UR", "GUT", "TO")))

write.csv(perm.combined, "Combined_PERMANOVA.csv", row.names = FALSE)

# Master summary plot output
pe.adonis.sig <- perm.combined %>% 
  filter(factor %in% (perm.combined %>% filter(Pvalue < 0.05 & !is.na(R2)) %>% pull(factor))) %>% 
  add_signif_stars() %>% 
  group_by(factor_cate, factor, source_file) %>% 
  mutate(max_R2 = max(R2), sig.label = paste0(star_pvalue, "/", star_padjust))

pe.adonis.sig$sig.label[pe.adonis.sig$sig.label == "-/-"] <- ""

if(nrow(pe.adonis.sig) >0){
  
  ggplot(pe.adonis.sig, aes(x = reorder(factor, R2), y = R2, fill = source_file)) +
    geom_col(width = 0.6, position = position_dodge2(preserve = "single")) + 
    geom_text(aes(y = max_R2, label = sig.label), color = "black", size = 2, vjust = 0.5, hjust = 1) +
    coord_flip() +
    scale_fill_manual(values = BASE_COLORS) +
    facet_grid(factor_cate ~ source_file, scales = "free_y", space = "free_y") +  
    labs(x = " ", y = expression(R^2~"(p / p adjust)"), title = "PERMANOVA (only strata patientID)") +
    theme_classic2(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(face = "bold"))
  
  ggsave("Combined_all_wide.pdf", width = 10, height = 6)
  write.csv(pe.adonis.sig, "Combined_all_wide.csv", row.names = FALSE)
  
}else{
  print("Not significance")
}

# ==============================================================================
# 5. Part 2: Timepoint & Pairwise Analysis #########
# ==============================================================================
dir_time <- "D:/WorkProjects/Demo-MHT 2026/Results/Permanova/Timepoint"
if(!dir.exists(dir_time)) dir.create(dir_time, recursive = TRUE)
setwd(dir_time)

map(sites, ~ {
  vars.Site <- .x
  file.prefix <- paste("permanova_strata", vars.Site, sep = "_")

  metadat_raw <- merge(Microbe.phen.prof %>% dplyr::select(-Time, -Clinic_ID) %>% subset(Site == vars.Site), 
                       Phen.Seq, by = "Clinic_ID_Time") %>% 
    subset(!is.na(SeqID)) %>% 
    .[, c("SeqID", "Clinic_ID", "Time")]
  prof_raw <- prof_filtered[[vars.Site]][metadat_raw$SeqID, ]
  valid_samples <- rownames(prof_raw)[rowSums(prof_raw) > 0]

  metadat <- metadat_raw %>% 
    filter(SeqID %in% valid_samples) %>% 
    column_to_rownames("SeqID")
  
  prof.input <- prof_raw[rownames(metadat), ] %>% .[, colSums(.) != 0]

  covardat <- data.frame(
    patient.ID = metadat$Clinic_ID,
    Time       = metadat$Time,
    row.names  = rownames(metadat)
  )

  metadat <- metadat %>% dplyr::select(Time)
  
  dist_matrix <- vegdist(prof.input, method = "bray", add = TRUE)
  
  set.seed(123)
  patientID.perm.fun(metadat = metadat, bray = dist_matrix, covariates = covardat, file.prefix = file.prefix, n.perm = 999, n.proc = parallel::detectCores() - 3)
  
  # 5.2 Pairwise two-by-two comparison
  pairwise_results <- pairwise.adonis2(dist_matrix ~ Time, data = covardat, strata = "patient.ID", nperm = 999)
  
  result_df <- lapply(setdiff(seq_along(pairwise_results), 1), function(i) {
    data.frame(pairs = names(pairwise_results)[i], R2 = pairwise_results[[i]]$R2[1], p.value = pairwise_results[[i]]$`Pr(>F)`[1])
  }) %>% bind_rows()  
  
  result_df$q.value <- p.adjust(result_df$p.value, method = "fdr")
  result_df$site <- vars.Site
  write.csv(result_df, paste0(file.prefix, "_pairwise.csv"), row.names = FALSE)
})

# 5.3 Summary Plot: Global Timepoint
perm.files.time <- list.files(getwd(), pattern = "permanova_strata_.*\\.csv", full.names = TRUE)

perm.combine.time <- lapply(perm.files.time[!str_detect(perm.files.time, "pairwise")], function(f) {
  df <- read.csv(f, header = TRUE)
  df$source_file <- gsub("permanova_strata_", "", sub("\\..*$", "", basename(f)))
  return(df)
}) %>% list_rbind() %>% add_signif_stars()

ggplot(perm.combine.time, aes(x = source_file, y = R2, fill = source_file)) +
  geom_col(width = 0.6) +
  geom_text(aes(y = R2, label = sprintf("%.3f", R2)), size = 2, color = "black", vjust = 1.1, hjust = 1) +
  geom_text(aes(y = R2, label = star_pvalue), size = 3, color = "black", vjust = 0.1, hjust = 1) +
  scale_fill_manual(values = BASE_COLORS) +
  coord_flip() + theme_minimal() +
  labs(x = "", y = expression(R^2~"(p value)"), title = "PERMANOVA - timepoints")

ggsave("Combined_timepoint.pdf", width = 6, height = 4)

# 5.4 Summary Plot: Pairwise Two-by-Two Comparisons
perm.combine.eachtime <- lapply(perm.files.time[str_detect(perm.files.time, "pairwise")], function(f) {
  df <- read.csv(f, header = TRUE)
  df$source_file <- gsub("permanova_strata_|_pairwise", "", sub("\\..*$", "", basename(f)))
  return(df)
}) %>% list_rbind() %>% 
  add_signif_stars(p_col = "p.value", padj_col = "q.value") %>% 
  group_by(pairs, source_file) %>% 
  arrange(desc(source_file)) %>% 
  mutate(y_pos = cumsum(R2), sig.label = ifelse(star_pvalue == "-" & star_padjust == "-", "", paste0(star_pvalue, "/", star_padjust)))

ggplot(perm.combine.eachtime, aes(x = source_file, y = R2, fill = source_file)) +
  geom_col(width = 0.6, position = "stack") +
  geom_text(aes(y = y_pos, label = sprintf("%.3f", R2)), size = 2, color = "black", vjust = 1.1, hjust = 1) +
  geom_text(aes(y = y_pos, label = sig.label), size = 3, color = "black", vjust = 0.1, hjust = 1) +
  scale_fill_manual(values = BASE_COLORS) +
  coord_flip() + facet_wrap(~pairs, scales = "fixed", ncol = 3) + theme_minimal() +
  labs(x = "Variables", y = expression(R^2~"(p / p adjust)"), title = "PERMANOVA - timepoint-pairs", caption = "P value < 0.05")

ggsave("Combined_timepoint.pairs.pdf", width = 8, height = 6)
write.csv(perm.combine.eachtime, "Combined_timepoint.pairs.csv", row.names = FALSE)