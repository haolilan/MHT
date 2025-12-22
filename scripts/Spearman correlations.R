library(tidyverse)
library(readr)
library(purrr)
## install.packages("openxlsx") # if not installed
MHT.demo <- read_rds("D:/Demo-MHT/data/MHT.demo.RDS")
## data.input 
phen.Cate <- MHT.demo$Metadata$phen.Cate
Phen.Seq <- MHT.demo$Metadata$Phen.Seq
Phen.Delta <- MHT.demo$Metadata$Phen.Delta

## set dir
target_dir <- paste0("D:/Demo-MHT/","Clinic_association")
dir.create(target_dir)
setwd(target_dir)

############################################################################################################
analyze_spearman <- function(Test_var, df, 
                             continuous_vars) {
  
 
  cor_continuous <- map_dfr(continuous_vars, ~{
 
    x_vec <- as.numeric(df[[Test_var]])
    y_vec <- as.numeric(df[[.x]])
  
    complete_cases <- complete.cases(x_vec, y_vec)
    x_vec <- x_vec[complete_cases]
    y_vec <- y_vec[complete_cases]
    
    # 计算相关系数
    ct <- cor.test(x_vec, y_vec, method = "spearman")
    
    tibble(
      Test_var = Test_var,
      Phenotype_var = .x,
      Phenotype_vartype = "Continuous",
      Method = "Spearman",
      Statistic = ct$statistic,
      Estimate = ct$estimate,
      P_value = ct$p.value,
      N = length(x_vec)   
    )
  })
  
  # 合并结果并校正
  cor_continuous%>%
    mutate(P_adj = p.adjust(P_value, "fdr"))%>%
    relocate(P_adj, .after = P_value)
}

## 1.4 delta_K ~ delta E2 ~ delta Metabolite############################
vars.x <- phen.Cate%>%subset(PhenCategory%in%c("Therapy_Score","Blood_Test","Hormone"))%>%
  subset(Phen_Time=="T24")%>%pull(Phen)
vars.x <- paste0("Delta_",vars.x)
vars.x
##
clinic_ID_2 <- Phen.Seq%>%
  subset(Time%in%c("BL","T24") )%>% group_by(Clinic_ID)%>%summarise(pairs = n())%>%filter(pairs==2)%>%pull("Clinic_ID")
##
out.mat <- Phen.Delta%>% subset(Delta_Time == "Delta_T24")  
out.mat <- out.mat[out.mat$Clinic_ID%in%clinic_ID_2,]

###
final_results <- map_dfr(
  vars.x, 
  analyze_spearman,
  df = out.mat,
  continuous_vars = vars.x
)

### select significance ###########
final_results$pheno.index <- gsub("Delta_","",final_results$Phenotype_var)
final_results$vars.index <- gsub("Delta_","",final_results$Test_var)
##
final_results$phenCate <- phen.Cate$PhenCategory[match(final_results$pheno.index,phen.Cate$Phen)]
final_results$varsCate <- phen.Cate$PhenCategory[match(final_results$vars.index,phen.Cate$Phen)]
final_results <- final_results %>%
  mutate(lab.p.adj = case_when(                
    P_adj < 0.001 ~ "***",
    P_adj < 0.01 ~ "**",
    P_adj < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  mutate(sign = case_when(                
    sign(Estimate) > 0 & P_adj < 0.05 ~ "+",
    sign(Estimate) < 0 & P_adj < 0.05 ~ "-",
    TRUE ~ ""
  ))%>%
  mutate(label = paste0(sign,"\n",lab.p.adj))

final_results <- final_results %>%
  mutate(Estimate_log = sign(Estimate) * log1p(abs(Estimate)))
final_results$p.adj.all = p.adjust(final_results$P_value,n=length(final_results$P_value), method = "fdr")
final_results$label
write.xlsx(final_results, "spearman.delta.phen2phen.xlsx", quote = FALSE)

## 3. all signatures #####################################
final_results.padj <- final_results %>%reshape2::dcast(.,Phenotype_var~Test_var,value.var = "P_adj")%>%column_to_rownames("Phenotype_var")

## sig.vars 
sig.vars <- c(rownames(final_results.padj)[final_results.padj%>%apply(.,1,function(x)any(x<0.05, na.rm = TRUE))],
              colnames(final_results.padj)[final_results.padj%>%apply(.,2,function(x)any(x<0.05, na.rm = TRUE))])
final_results.sig <- final_results%>%subset(Phenotype_var%in%sig.vars&Test_var%in%sig.vars)

## order 
sig.vars.order <- c(final_results.sig%>%arrange(phenCate)%>%pull(Phenotype_var),final_results.sig%>%arrange(phenCate)%>%pull(Test_var))%>%
  unique()%>%rev()
final_results.sig$Phenotype_var <- factor(final_results.sig$Phenotype_var,levels = sig.vars.order)
final_results.sig$Test_var <- factor(final_results.sig$Test_var,levels = sig.vars.order)
final_results.sig%>%ggplot(aes(Phenotype_var,Test_var,fill=Estimate))+geom_tile()+facet_grid()
n.line <- table(final_results$phenCate[match(sig.vars.order,final_results$Phenotype_var)])%>%rev()%>%.[-length(.)]%>%cumsum()%>%as.vector()+0.5

## plot
plot_data <- final_results.sig
n <- length(unique(plot_data$Phenotype_var))

plot_data_modified <- plot_data %>%
  group_by(Phenotype_var, Test_var) %>% 
  mutate(
    pheno_index = as.numeric(Phenotype_var),
    vars_index = as.numeric(Test_var)
  ) %>%
  ungroup() %>%
  mutate(
    Estimate_log = ifelse(pheno_index <= vars_index, Estimate_log, NA),  
    sign = ifelse(pheno_index <= vars_index, sign, NA),         
    lab.p.adj = ifelse(pheno_index <= vars_index, lab.p.adj, NA)
  ) %>%
  select(-pheno_index, -vars_index) %>%
  subset(!is.na(Estimate_log))

plot_data_modified %>%
  ggplot(aes(x = Phenotype_var, y = Test_var, fill = Estimate_log)) +
  geom_tile( color = "white", linewidth = 0.5) + 
  geom_text(aes(label = sign), size = 3, vjust = 0.2,color = "grey30") +  
  geom_text(aes(label = lab.p.adj), size = 2, vjust = 1.5, color = "grey30") +  
  scale_fill_gradient2(
    low = "#6A0DAD",     
    mid = "#F0F0F0",       
    high = "#b35806",     
    midpoint = 0,          
    name = "log(Estimate)"   
  ) +
  scale_x_discrete(position = "top") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 0,hjust = 0),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = -90),   
    panel.spacing.y = unit(0.5, "lines")    
  ) +
  labs(x = " ", y = " ", title = "Delta_K Score ~ Hormone ~ Metabolite")


ggsave("Spearman.Delta_KScore_Hormone_Metabolite.all.pdf",width = 11,height = 9)

