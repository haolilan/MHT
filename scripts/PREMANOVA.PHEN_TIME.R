
MHT.demo <- read_rds("D:/Demo-MHT/data/MHT.demo.RDS")
## data.input 
phen.Cate <- MHT.demo$Metadata$phen.Cate
# clinic.cate.order <- MHT.demo$vars$clinic.cate.order
Phen.Seq <- MHT.demo$Metadata$Phen.Seq
prof_filtered <- MHT.demo$MG$prof_filtered
Microbe.phen.prof <- MHT.demo$MG$Microbe.phen.prof
##
Phen.Seq$Clinic_ID_Time <- paste(Phen.Seq$Clinic_ID,Phen.Seq$Time,sep = "_")

## set dir
target_dir <- paste0("D:/Demo-MHT/","Permanova")
dir.create(target_dir)
setwd(target_dir)

# function ###############################

adonis.plot <- function(pe.adonis, file.prefix) {
  pe.adonis.sig <- subset(pe.adonis, Pvalue < 0.05 & !is.na(R2)) 
  pe.adonis.sig$star_pvalue <- cut(pe.adonis.sig$Pvalue,
                                   breaks = c(0, 0.001, 0.01, 0.05, 1),
                                   labels = c("***", "**", "*", "-"),
                                   include.lowest = TRUE)
  
  pe.adonis.sig$star_padjust <- cut(pe.adonis.sig$Padjust,
                                    breaks = c(0, 0.001, 0.01, 0.05, 1),
                                    labels = c("***", "**", "*", "-"),
                                    include.lowest = TRUE)
  ##
  ggplot(pe.adonis.sig, aes(x = reorder(factor, R2), y = R2, 
                            fill = -log10(Padjust))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("R²=%.3f", R2)), 
              hjust = 0.2, vjust=1,size = 3, color = "black") +
    geom_text(aes(label = paste0(star_pvalue, "/", star_padjust),y=0), 
              hjust = 1, color = "black", size = 3) +
    scale_fill_viridis_c(name = "-log10(q value)") +
    coord_flip() +
    labs(x = "Variables", 
         y = expression(R^2~"(p / p adjust)"),  
         title = "PERMANOVA (strata: patientID)",
         caption = "P value < 0.05") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10))
    )
}
patientID.perm.fun <- function(metadat, bray, covariates, file.prefix, n.perm = 999, n.proc = 2) {
  
  require(vegan)
  require(foreach)
  require(doParallel)
  require(ggplot2)
  library(tibble)
  
  log_message <- function(msg) {
    message(paste(Sys.time(), "|", msg))
  }
  
  validate_inputs <- function(metadat, bray) {
    if (!inherits(bray, "dist")) stop("bray must be matrix")
    if (nrow(metadat) != attr(bray, "Size")) 
      stop(paste("nrow(", nrow(metadat), ") is unmatched"))
    return(TRUE)
  }
  
  preprocess_data <- function(dat) {
    dat <- dat[, colSums(!is.na(dat)) >= 3] 
    dat[] <- lapply(dat, function(col) {
      if (is.factor(col) && nlevels(col) < 8 && sum(!is.na(col)) > 0) {
        cutoff <- round(length(col) * 0.01)
        tbl <- table(col)
        rare_levels <- names(tbl)[tbl < cutoff]
        col[col %in% rare_levels] <- NA
      }
      return(col)
    })
    return(dat)
  }
  
  # main
  log_message("start to preprocess data")
  validate_inputs(metadat, bray)
  
  # 初始化并行环境
  cl <- makeCluster(min(n.proc, 15))
  registerDoParallel(cl)
  
  results <- foreach(i = seq_len(ncol(metadat)), .combine = rbind,
                     .packages = c("vegan"), 
                     .errorhandling = "pass") %dopar% {
                       tryCatch({
                         fac <- colnames(metadat)[i]
                         idx <- !is.na(metadat[[fac]])
                         samples <- intersect(rownames(metadat)[idx], labels(bray))
                         
                         if (length(samples) < max(5, nrow(metadat)*0.05)) {
                           return(data.frame(factor=fac, SampleNum=length(samples), 
                                             Df=NA, SumsOfSqs=NA, R2=NA, F.Model=NA,
                                             Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                         }
                         
                         sub_dist <- as.dist(as.matrix(bray)[samples, samples])
                         fac.dat <- metadat[samples, fac]
                         
                         fac.patientID <- covariates[samples,"patient.ID"]
                         
                         res <- adonis2(sub_dist ~ fac.dat,
                                        strata = fac.patientID,
                                        permutations = n.perm,by="term")

                         res <- as.data.frame(res)
                         
                         disp_pval <- if(length(unique(fac.dat))>1) {
                           permutest(betadisper(sub_dist, fac.dat),
                                     strata=fac.patientID,
                                     permutations = how(nperm=n.perm))$tab$`Pr(>F)`[1]
                         } else NA
                         
                         data.frame(
                           factor = fac,
                           SampleNum = length(samples),
                           Df = res["fac.dat","Df"], # res$Df[1],
                           SumsOfSqs = res["fac.dat","SumOfSqs"],# res$SumOfSqs[1],
                           R2 = res["fac.dat","R2"], #res$R2[1],
                           F.Model = res["fac.dat","F"], #res$F[1],
                           Pvalue = res["fac.dat","Pr(>F)"], # res["fac.dat","R2"],res$`Pr(>F)`[1],
                           Disp_pval = disp_pval,
                           R2.adjust = RsquareAdj(res["fac.dat","R2"], length(samples), res["fac.dat","Df"]), #RsquareAdj(res$R2[1], length(samples), res$Df[1]),
                           stringsAsFactors = FALSE
                         )
                       }, error = function(e) {
                         message("Error in variable ", fac, ": ", e$message)
                         return(data.frame(factor=fac, SampleNum=NA, Df=NA, SumsOfSqs=NA, 
                                           R2=NA, F.Model=NA, Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                       }
                       )
                     }
  
  stopCluster(cl)
  
  out.stat <- data.frame(
    factor = colnames(metadat),
    var_type = "",
    na_count = "",
    stats = "",
    stringsAsFactors = FALSE
  )
  for (i in 1:ncol(metadat)) {
    col_data <- metadat[[i]]
    na_count <- sum(is.na(col_data))

    if (is.numeric(col_data)) {
      stats_str <- sprintf("Mean+-SD: %.2f +- %.2f",
                           mean(col_data, na.rm = TRUE),
                           sd(col_data, na.rm = TRUE))
      var_type <- "numeric"
    } else if (is.factor(col_data)) {
      freq_table <- table(col_data, useNA = "no")
      stats_str <- paste(names(freq_table), freq_table, sep=":", collapse="; ")
      var_type <- ifelse(is.ordered(col_data), "ordered factor", "factor")
    } else {
      next
    }
    
    out.stat[i, c("var_type", "na_count", "stats")] <- list(
      var_type, na_count, stats_str
    )
  }
  
  # PERM
  out.perm <- as.data.frame(results)
  out.perm$Padjust <- p.adjust(out.perm$Pvalue, method = "fdr")
  out.perm <- merge(out.perm,out.stat,by="factor",all.x=T)
  out.perm <- out.perm[order(out.perm$R2, decreasing = TRUE), ]
  
  # 输出结果
  write.csv(out.perm, paste0(file.prefix, ".csv"), row.names = FALSE)
  p <- adonis.plot(out.perm, file.prefix)
  length.signif <- nrow(subset(out.perm, Pvalue < 0.05 & !is.na(R2)) )
  ggsave(paste0(file.prefix, ".pdf"), p, width = 12, 
         height = ifelse(length.signif*0.8>24,24,
                         ifelse(length.signif*0.8<4,4,length.signif*0.8)))
  
  log_message("completed analysis")
  return(list(plot = p, results = out.perm))
}

# 1. All ###################################
target_dir <- paste0("D:/Demo-MHT/","Permanova/All")
dir.create(target_dir)
setwd(target_dir)
##
sites <- unique(Microbe.phen.prof$Site)
library(purrr)
map(sites, ~ {
  vars.Site <- .x
  ##
  file.prefix <- paste("permanova_strata",vars.Site,sep = "_")
  
  ## metadata
  vars.phen <- setdiff(colnames(Phen.Seq),c("Clinic_ID_Time","Clinic_ID"))
  metadat <-
    merge(Microbe.phen.prof %>% dplyr::select(-Time)%>% subset(Site==vars.Site),  
          Phen.Seq, 
          by="Clinic_ID_Time")%>% subset(!is.na(SeqID))%>%
    .[,c("SeqID",vars.phen)]  
  
  # 
  metadat <- metadat %>%
    filter(rowSums(is.na(dplyr::select(., all_of(vars.phen)))) < length(vars.phen))
  rownames(metadat) <- NULL
  metadat <- metadat %>% column_to_rownames("SeqID")
  
  ##  COV, CLINIC_ID，Clinic_ID,Time
  covardat <- data.frame(
    patient.ID = Microbe.phen.prof$Clinic_ID[match(rownames(metadat),Microbe.phen.prof$SeqID)] 
  ) 
  rownames(covardat) <- rownames(metadat)
  
  #  
  identical(rownames(covardat),rownames(metadat))
  str(covardat)
  str(metadat[,1:10])
  
  # 对部分变量进行 scale ###
  scale.vars <- names(metadat)[sapply(metadat, is.numeric)]
  metadat <- metadat %>%
    mutate(across(all_of(scale.vars), ~ scale(.)))
  
  ## prof.input
  prof.input <-
    prof_filtered[[vars.Site]]%>% .[rownames(metadat),] %>% .[rowSums(.) != 0,colSums(.) != 0 ]
  identical(rownames(metadat),rownames(prof.input))
  ##
  perm <- patientID.perm.fun(
    metadat = metadat,
    bray = vegdist(prof.input, method = "bray",add = TRUE) , 
    covariates = covardat,
    file.prefix = file.prefix,
    n.perm = 9, ## 999
    n.proc = parallel::detectCores() - 1
  )
}
)

# 2. Timepoint ###################################
target_dir <- paste0("D:/Demo-MHT/","Permanova/Timepoint")
dir.create(target_dir)
setwd(target_dir)

patientID.perm.fun <- function(metadat, bray, covariates, file.prefix, n.perm = 999, n.proc = 2) {
  
  require(vegan)
  require(foreach)
  require(doParallel)
  require(ggplot2)
  
  log_message <- function(msg) {
    message(paste(Sys.time(), "|", msg))
  }
  
  validate_inputs <- function(metadat, bray) {
    if (!inherits(bray, "dist")) stop("bray must be matrix")
    if (nrow(metadat) != attr(bray, "Size")) 
      stop(paste("nrow(", nrow(metadat), ") is unmatched"))
    return(TRUE)
  }
  
  preprocess_data <- function(dat) {
    dat <- dat[, colSums(!is.na(dat)) >= 3]  
    dat[] <- lapply(dat, function(col) {
      if (is.factor(col) && nlevels(col) < 8 && sum(!is.na(col)) > 0) {
        cutoff <- round(length(col) * 0.01)
        tbl <- table(col)
        rare_levels <- names(tbl)[tbl < cutoff]
        col[col %in% rare_levels] <- NA
      }
      return(col)
    })
    return(dat)
  }
  
  
  # 主流程
  log_message("start to ")
  validate_inputs(metadat, bray)
  
  # 初始化并行环境
  cl <- makeCluster(min(n.proc, 15))
  registerDoParallel(cl)
  
  results <- foreach(i = seq_len(ncol(metadat)), .combine = rbind,
                     .packages = c("vegan"), 
                     .errorhandling = "pass") %dopar% {
                       tryCatch({
                         fac <- colnames(metadat)[i]
                         idx <- !is.na(metadat[[fac]])
                         samples <- intersect(rownames(metadat)[idx], labels(bray))
                         
                         
                         if (length(samples) < max(5, nrow(metadat)*0.05)) {
                           return(data.frame(factor=fac, SampleNum=length(samples), 
                                             Df=NA, SumsOfSqs=NA, R2=NA, F.Model=NA,
                                             Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                         }
                         
                         sub_dist <- as.dist(as.matrix(bray)[samples, samples])
                         fac.dat <- metadat[samples, fac]
                         
                         fac.patientID <- covariates[samples,"patient.ID"]
                         
                         res <- adonis2(sub_dist ~ fac.dat,
                                        strata = fac.patientID,
                                        permutations = n.perm,by="term")

                         res <- as.data.frame(res)
                         
                         disp_pval <- if(length(unique(fac.dat))>1) {
                           permutest(betadisper(sub_dist, fac.dat),
                                     strata=fac.patientID,
                                     permutations = how(nperm=n.perm))$tab$`Pr(>F)`[1]
                         } else NA
                         
                         # 
                         data.frame(
                           factor = fac,
                           SampleNum = length(samples),
                           Df = res["fac.dat","Df"], # res$Df[1],
                           SumsOfSqs = res["fac.dat","SumOfSqs"],# res$SumOfSqs[1],
                           R2 = res["fac.dat","R2"], #res$R2[1],
                           F.Model = res["fac.dat","F"], #res$F[1],
                           Pvalue = res["fac.dat","Pr(>F)"], # res["fac.dat","R2"],res$`Pr(>F)`[1],
                           Disp_pval = disp_pval,
                           R2.adjust = RsquareAdj(res["fac.dat","R2"], length(samples), res["fac.dat","Df"]), #RsquareAdj(res$R2[1], length(samples), res$Df[1]),
                           stringsAsFactors = FALSE
                         )
                       }, error = function(e) {
                         message("Error in variable ", fac, ": ", e$message)
                         return(data.frame(factor=fac, SampleNum=NA, Df=NA, SumsOfSqs=NA, 
                                           R2=NA, F.Model=NA, Pvalue=NA, Disp_pval=NA, R2.adjust=NA))
                       }
                       )
                     }
  
  stopCluster(cl)
  results
  
  out.stat <- data.frame(
    factor = colnames(metadat),
    var_type = "",
    na_count = "",
    stats = "",
    stringsAsFactors = FALSE
  )
  for (i in 1:ncol(metadat)) {
    col_data <- metadat[[i]]
    na_count <- sum(is.na(col_data))
    
    # 统计量计算
    if (is.numeric(col_data)) {
      stats_str <- sprintf("Mean+-SD: %.2f +- %.2f",
                           mean(col_data, na.rm = TRUE),
                           sd(col_data, na.rm = TRUE))
      var_type <- "numeric"
    } else if (is.factor(col_data)) {
      freq_table <- table(col_data, useNA = "no")
      stats_str <- paste(names(freq_table), freq_table, sep=":", collapse="; ")
      var_type <- ifelse(is.ordered(col_data), "ordered factor", "factor")
    } else {
      next
    }
    
    # 填充统计列
    out.stat[i, c("var_type", "na_count", "stats")] <- list(
      var_type, na_count, stats_str
    )
  }
  
  # PERM结果后处理
  out.perm <- as.data.frame(results)
  out.perm$Padjust <- p.adjust(out.perm$Pvalue, method = "fdr")
  out.perm <- merge(out.perm,out.stat,by="factor",all.x=T)
  out.perm <- out.perm[order(out.perm$R2, decreasing = TRUE), ]
  
  # 输出结果
  write.csv(out.perm, paste0(file.prefix, ".csv"), row.names = FALSE)
  p <- adonis.plot(out.perm, file.prefix)
  length.signif <- nrow(subset(out.perm, Pvalue < 0.05 & !is.na(R2)) )
  ggsave(paste0(file.prefix, ".pdf"), p, width = 12, 
         height = ifelse(length.signif*0.8>24,24,
                         ifelse(length.signif*0.8<4,4,length.signif*0.8)))
  
  log_message("Completed analysis")
  return(list(plot = p, results = out.perm))
}

sites <- unique(Microbe.phen.prof$Site)
vars.phen <- "Time"

map(sites, ~ {
  vars.Site <- .x
  ##
  file.prefix <- paste("permanova_strata",vars.Site,sep = "_")
  
  ## metadat
  metadat <-
    merge(Microbe.phen.prof%>% select(-Time)%>% subset(Site==vars.Site),  
          Phen.Seq, 
          by="Clinic_ID_Time")%>% subset(!is.na(SeqID))%>%
    .[,c("SeqID",vars.phen)] 
  rownames(metadat) <- NULL
  metadat <- metadat %>% column_to_rownames("SeqID")
  
  ## 
  covardat <- data.frame(
    patient.ID = Microbe.phen.prof$Clinic_ID[match(rownames(metadat),Microbe.phen.prof$SeqID)],
    Time = Microbe.phen.prof$Time[match(rownames(metadat),Microbe.phen.prof$SeqID)]
  ) 
  rownames(covardat) <- rownames(metadat)

  identical(rownames(covardat),rownames(metadat))
  str(covardat)
  str(metadat)
  
  ## prof.input
  prof.input <-
    prof_filtered[[vars.Site]]%>% .[rownames(metadat),] %>% .[rowSums(.) != 0,colSums(.) != 0 ]
  identical(rownames(metadat),rownames(prof.input))
  
  ##
  perm <- patientID.perm.fun(
    metadat = metadat,
    bray = vegdist(prof.input, method = "bray",add = TRUE) , 
    covariates = covardat,
    file.prefix = file.prefix,
    n.perm = 9, #999,
    n.proc = parallel::detectCores() - 3
  )
  

  pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
    
    #describe parent call function
    ststri <- ifelse(is.null(strata),'Null',strata)
    fostri <- as.character(x)
    #list to store results
    
    #copy model formula
    x1 <- x
    # extract left hand side of formula
    lhs <- eval(x1[[2]], environment(x1), globalenv())
    environment(x1) <- environment()
    # extract factors on right hand side of formula
    rhs <- x1[[3]]
    # create model.frame matrix
    x1[[2]] <- NULL
    rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
    
    # create unique pairwise combination of factors
    co <- combn(unique(as.character(rhs.frame[,1])),2)
    
    # create names vector
    nameres <- c('parent_call')
    for (elem in 1:ncol(co)){
      nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
    }
    #create results list
    res <- vector(mode="list", length=length(nameres))
    names(res) <- nameres
    
    #add parent call to res
    res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
    
    
    #start iteration trough pairwise combination of factors
    for(elem in 1:ncol(co)){
      
      #reduce model elements
      if(inherits(eval(lhs),'dist')){
        xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                             rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
      }else{
        xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
      }
      
      mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
      
      # redefine formula
      if(length(rhs) == 1){
        xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))
      }else{
        xnew <- as.formula(paste('xred' ,
                                 paste(rhs[-1],collapse= as.character(rhs[1])),
                                 sep='~'))}
      
      #pass new formula to adonis
      if(is.null(strata)){
        ad <- adonis2(xnew,data=mdat1, ... )
      }else{
        perm <- how(nperm = nperm)
        setBlocks(perm) <- with(mdat1, mdat1[,ststri])
        ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
      
      res[nameres[elem+1]] <- list(ad[1:5])
    }
    #names(res) <- names
    class(res) <- c("pwadstrata", "list")
    return(res)
  }

  # 
  pairwise_results <- pairwise.adonis2(
    vegdist(prof.input, method = "bray",add = TRUE)~Time, 
    data = covardat, 
    strata = "patient.ID",
    nperm = 999
  )
  
  # 
  result_df <- lapply(setdiff(seq_along(pairwise_results),1), function(i) {
    data.frame(
      pairs = names(pairwise_results)[i], 
      R2 = pairwise_results[[i]]$R2[1],
      p.value = pairwise_results[[i]]$`Pr(>F)`[1]
    )
  }) %>% bind_rows()  
  result_df$q.value <- p.adjust(result_df$p.value,method = "fdr")
  result_df$site <- vars.Site
  write.csv(result_df, paste0(file.prefix, "_pairwise.csv"), row.names = FALSE)
}
)



# 3. Combine ################
## 3.0 all  ################
target_dir <- paste0("D:/Demo-MHT/","Permanova/All")
setwd(target_dir)
perm.files <- list.files(getwd(), 
                         pattern = "permanova_strata.*\\.csv", 
                         full.names = TRUE)
perm.combined <- lapply(perm.files, function(f) {
  df <- read.csv(f, header = TRUE) 
  filename <- sub("\\..*$", "", basename(f))  
  df$source_file <- filename%>%gsub("permanova_strata_","",.)    
  return(df)
})%>% list_rbind() 

## match category
perm.combined$factor_cate <- phen.Cate$PhenCategory[match(perm.combined$factor,phen.Cate$Phen)]
perm.combined<- perm.combined%>% mutate(factor_cate=case_when(
  factor%in%"Time" ~ "Therapy_Score",
  factor%in%c("LH","FSH" ,"E2","Testosterone") ~ "Hormone",
  str_detect(factor, "^Delta") ~ "Delta_Therapy_Score", 
  TRUE  ~ factor_cate)
)
perm.combined$factor_cate%>%unique()

## p ajust redo
perm.combined <- perm.combined%>%subset(factor_cate%in%c("Therapy_Score","Hormone","Blood_Test"))
perm.combined$Padjust <- p.adjust(perm.combined$Pvalue, method = "fdr")

### order
perm.combined$factor_cate <- factor(perm.combined$factor_cate,levels = c("General","Therapy_Score","Hormone","Blood_Test","Questionnaire","Ultrasound"))
perm.combined$source_file <- factor(perm.combined$source_file,levels = c("VA" ,"UR","GUT", "TO"  ))

##### figure output #####
factor.sig <- perm.combined%>%subset(Pvalue<0.05 & !is.na(R2))%>%pull(factor)
pe.adonis.sig <- subset(perm.combined, factor %in% factor.sig) 

# 创建显著性星号标签（根据P值范围）
base_colors <- c("#2CA02C", "#D62728", "#FF7F0E","#1F77B4" ) ## gut, to, ur , va
names(base_colors) <- c("GUT","TO","UR","VA")
pe.adonis.sig$star_pvalue <- cut(pe.adonis.sig$Pvalue,
                                 breaks = c(0, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", "-"),
                                 include.lowest = TRUE)
pe.adonis.sig$star_padjust <- cut(pe.adonis.sig$Padjust,
                                  breaks = c(0, 0.001, 0.01, 0.05, 1),
                                  labels = c("***", "**", "*", "-"),
                                  include.lowest = TRUE)

###### plots ####################
pe.adonis.sig <- pe.adonis.sig %>%
  group_by(factor_cate,factor, source_file) %>%
  mutate(
    max_R2 = max(R2),
    meadian_R2 = max_R2 - R2/2,
    sig.label = paste0(star_pvalue, "/", star_padjust)
  )
pe.adonis.sig$sig.label[pe.adonis.sig$sig.label=="-/-"] <- ""
##
ggplot(pe.adonis.sig, aes(x = reorder(factor, R2), y = R2,fill = source_file)) +
  geom_col(width = 0.6,position = position_dodge2(preserve = "single") ) + 
  geom_text(aes(y = max_R2, label = sig.label),
            color = "black", size = 2, vjust = 0.5, hjust=1) +
  coord_flip() +
  scale_fill_manual(values = base_colors)+
  facet_grid(factor_cate~source_file,  scales = "free_y", space = "free_y")+  
  labs(x = " ",
       y = expression(R^2~"(p / p adjust)"),  
       title = "PERMANOVA (only strata patientID)",
       caption = " ") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10))
  )
ggsave("Combined_all_wide.pdf",width=10,height=6)
write.csv(pe.adonis.sig, "Combined_all_wide.csv", row.names = FALSE)




## 3.1 timepoints #####
target_dir <- paste0("D:/Demo-MHT/","Permanova/Timepoint")
setwd(target_dir)
perm.files <- list.files(getwd(), 
                         pattern = "permanova_strata.*\\.csv", 
                         full.names = TRUE)
## whole
perm.files1 <- perm.files[grep("pairwise",perm.files,invert = T)]
perm.combine.time <- lapply(perm.files1, function(f) {
  df <- read.csv(f, header = TRUE) 
  filename <- sub("\\..*$", "", basename(f))  
  df$source_file <- filename%>%gsub("permanova_strata_","",.)               
  return(df)
})%>% list_rbind() 

perm.combine.time$star_pvalue <- cut(perm.combine.time$Pvalue,
                                     breaks = c(0, 0.001, 0.01, 0.05, 1),
                                     labels = c("***", "**", "*", "-"),
                                     include.lowest = TRUE)
perm.combine.time$star_padjust <- cut(perm.combine.time$Padjust,
                                      breaks = c(0, 0.001, 0.01, 0.05, 1),
                                      labels = c("***", "**", "*", "-"),
                                      include.lowest = TRUE)
base_colors <- c("#2CA02C", "#D62728", "#FF7F0E","#1F77B4" ) ## gut, to, ur , va
perm.combine.time <- perm.combine.time %>%
  arrange(desc(source_file)) %>%     
  mutate(
    y_pos = R2,         
    label_y = y_pos - R2/2 ,
    sig.label = paste0(star_pvalue, "/", star_padjust)
  )
perm.combine.time$sig.label[perm.combine.time$sig.label=="-/-"] <- ""

ggplot(perm.combine.time,
       aes(x = source_file, y = R2, fill = source_file)) +
  geom_col(width = 0.6, position = "stack") +
  geom_text(
    aes(y = y_pos, label = sprintf("%.3f", R2)), 
    size = 2, color = "black", vjust = 1.1,hjust=1,
  ) +
  geom_text(
    aes(y = y_pos, label = star_pvalue), 
    size = 3, color = "black", vjust = 0.1, hjust=1,
  ) +
  scale_fill_manual(values = base_colors)+
  coord_flip() +
  theme_minimal() +
  labs(x = "", 
       y = expression(R^2~"(p value)"), 
       title = "PERMANOVA - timepoints")
ggsave("Combined_timepoint.pdf",width=6,height=4)


## pairs 
perm.files2 <- perm.files[grep("pairwise",perm.files,invert = F)]
perm.combine.eachtime <- lapply(perm.files2, function(f) {
  df <- read.csv(f, header = TRUE) 
  filename <- sub("\\..*$", "", basename(f))  
  df$source_file <- filename%>%gsub("permanova_strata_|_pairwise","",.)                 
  return(df)
})%>% list_rbind() 
perm.combine.eachtime$star_pvalue <- cut(perm.combine.eachtime$p.value,
                                         breaks = c(0, 0.001, 0.01, 0.05, 1),
                                         labels = c("***", "**", "*", "-"),
                                         include.lowest = TRUE)
perm.combine.eachtime$star_padjust <- cut(perm.combine.eachtime$q.value,
                                          breaks = c(0, 0.001, 0.01, 0.05, 1),
                                          labels = c("***", "**", "*", "-"),
                                          include.lowest = TRUE)
# 定义原始颜色
base_colors <- c("#2CA02C", "#D62728", "#FF7F0E","#1F77B4" ) ## gut, to, ur , va
names(base_colors) <- c("GUT","TO","UR","VA")
perm.combine.eachtime <- perm.combine.eachtime %>%
  group_by(pairs,source_file)%>%
  arrange(desc(source_file)) %>%     
  mutate(
    y_pos = cumsum(R2),              
    label_y = y_pos - R2/2,
    sig.label = paste0(star_pvalue, "/", star_padjust)
  )

ggplot(perm.combine.eachtime,
       aes(x = source_file, y = R2, fill = source_file)) +
  geom_col(width = 0.6, position = "stack") +
  geom_text(
    aes(y = y_pos, label = sprintf("%.3f", R2)), 
    size = 2, color = "black", vjust = 1.1,hjust=1,
  ) +
  geom_text(
    aes(y = y_pos, label = sig.label), 
    size = 3, color = "black", vjust = 0.1, hjust=1,
  ) +
  scale_fill_manual(values = base_colors)+
  coord_flip() +
  facet_wrap(~pairs,nrow = 2)+
  theme_minimal() +
  labs(x = "Variables", 
       y = expression(R^2~"(p / p adjust)"),  
       title = "PERMANOVA - timepoint-pairs",
       caption = "P value < 0.05")
ggsave("Combined_timepoint.pairs.pdf",width=12,height=8)






