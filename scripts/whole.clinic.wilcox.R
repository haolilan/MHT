
library(reshape2)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(data.table)
library(stringr)
library(forcats)
library(rstatix)

#############################################
MHT.demo <- read_rds("D:/Demo-MHT/data/MHT.demo.RDS")

## data.input 
phen.Cate <- MHT.demo$Metadata$phen.Cate
phen.input <- MHT.demo$Metadata$Phen.Seq
clinic.cate.order <- MHT.demo$vars$clinic.cate.order

## set dir
target_dir <- paste0("D:/Demo-MHT/","Clinic_whole_feature")
dir.create(target_dir)
setwd(target_dir)

# 1.0 clinic index compared: as T24 vs BL ############################################################################################################
## 1.1 clinic only, BL&T24  ############
phen.vars <- phen.Cate$Phen[phen.Cate$VarsType=="continuous_vars" &
                              phen.Cate$Category%in%setdiff(clinic.cate.order,c("K Score","K Sub Score"))]%>%
  unique()
####
prefix <- "Phen.nonK" 
my_comparisons <- list(c("BL","T24"))
top_factors <- "BL" 

####
data <- phen.input[,c("Time","Clinic_ID",phen.vars)]%>%
  rename(Response=Time)
xtabs(~Response,data)

sum_stats_results <- c()
sum_plot_data_results <- c()
for (phen.var in phen.vars){    
  stats_results <- c()
  plot_data_results <- c()
  
  for (i in 1:length(my_comparisons)) {
    df <- data[,c("Clinic_ID","Response",phen.var)]%>%na.omit()
    pairs.samp <- df %>%subset(Response%in%my_comparisons[[i]])%>% group_by(Clinic_ID)%>%summarise(pairs = n())%>%filter(pairs==2)%>%pull(Clinic_ID)
    plot_data <- df %>%
      subset(Clinic_ID %in% pairs.samp & Response%in%my_comparisons[[i]])%>%
      arrange(Clinic_ID,Response)%>%
      select(Response,Clinic_ID,phen.var) %>%
      reshape2::melt(id.vars = c("Response","Clinic_ID"),variable.name = "Factor",value.name = "Value")%>%
      mutate(Response = fct_relevel(Response, top_factors)) 
    head(plot_data )
    
    # 4. 批量计算组间差异p值（校正多重检验）
    stat_data <- plot_data %>%
      group_by(Factor) %>%
      wilcox_test(
        Value ~ Response,
        paired = TRUE,
        detailed = TRUE) %>%  
      adjust_pvalue(method = "fdr")
    stat_plot <- plot_data %>%
      group_by(Factor, Response) %>%
      summarise(Mean_Value = mean(Value, na.rm = TRUE))%>%
      reshape2::dcast(.,Factor~Response,value.var="Mean_Value")%>%
      mutate(Mean1=.[,2],Mean2=.[,3],
             F2C=log2(Mean2/Mean1),
             Mean.diff=Mean2-Mean1)%>%
      mutate(Mean.diff.log=sign(Mean.diff) * log1p(abs(Mean.diff)))%>%
      select(Factor,Mean1,Mean2,F2C,Mean.diff.log)
    stat_data <- merge(stat_data,stat_plot,by="Factor",all.x = T)
    
    stats_results <- rbind(stats_results,stat_data%>%mutate(comparisons=paste(my_comparisons[[i]],collapse = "/")))
    plot_data_results <- rbind(plot_data_results,plot_data%>%mutate(comparisons=paste(my_comparisons[[i]],collapse = "/")))
  }
  
  sum_stats_results <- rbind(sum_stats_results,stats_results)
  sum_plot_data_results <- rbind(sum_plot_data_results,plot_data_results)
}
sum_stats_results <- sum_stats_results %>%
  mutate(p.adj = p.adjust(p,method="fdr")) %>%  
  mutate(lab.p.adj = case_when(               
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%  
  mutate(lab.p = case_when(              
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%
  mutate(label=paste(lab.p,lab.p.adj,sep = "/"))

write.csv(sum_stats_results,paste0(prefix,".phen.wilcox.csv"),quote = F, row.names = F)
write.csv(sum_plot_data_results,paste0(prefix,".phen.wilcox.box.csv"),quote = F, row.names = F)

## 1.2 k only T04 - T24 ############
phen.vars <-  phen.Cate$Phen[phen.Cate$VarsType=="continuous_vars" &
                               phen.Cate$Category%in%c("K Score","K Sub Score")]%>%
  unique()
phen.vars
####
prefix <- "Phen.ONLYK"  # Paired.4sample

comb_matrix <- combn(c("BL","T04","T12","T24"), 2)
my_comparisons <- lapply(split(t(comb_matrix), 1:ncol(comb_matrix)), as.vector)
top_factors <- "BL" # group1

####
data <- phen.input[,c("Time","Clinic_ID",phen.vars)]%>%
  rename(Response=Time)
xtabs(~Response,data)

sum_stats_results <- c()
sum_plot_data_results <- c()
for (phen.var in phen.vars){    
  stats_results <- c()
  plot_data_results <- c()
  
  for (i in 1:length(my_comparisons)) {
    df <- data[,c("Clinic_ID","Response",phen.var)]%>%na.omit()
    pairs.samp <- df %>%subset(Response%in%my_comparisons[[i]])%>% group_by(Clinic_ID)%>%summarise(pairs = n())%>%filter(pairs==2)%>%pull(Clinic_ID)
    plot_data <- df %>%
      subset(Clinic_ID %in% pairs.samp & Response%in%my_comparisons[[i]])%>%
      arrange(Clinic_ID,Response)%>%
      select(Response,Clinic_ID,phen.var) %>%
      reshape2::melt(id.vars = c("Response","Clinic_ID"),variable.name = "Factor",value.name = "Value")%>%
      mutate(Response = fct_relevel(Response, top_factors)) 
    head(plot_data )
    
    # 4. 批量计算组间差异p值（校正多重检验）
    stat_data <- plot_data %>%
      group_by(Factor) %>%
      wilcox_test(
        Value ~ Response,
        paired = TRUE,
        detailed = TRUE) %>%  
      adjust_pvalue(method = "fdr")
    stat_plot <- plot_data %>%
      group_by(Factor, Response) %>%
      summarise(Mean_Value = mean(Value, na.rm = TRUE))%>%
      reshape2::dcast(.,Factor~Response,value.var="Mean_Value")%>%
      mutate(Mean1=.[,2],Mean2=.[,3],
             F2C=log2(Mean2/Mean1),
             Mean.diff=Mean2-Mean1)%>%
      mutate(Mean.diff.log=sign(Mean.diff) * log1p(abs(Mean.diff)))%>%
      select(Factor,Mean1,Mean2,F2C,Mean.diff.log)
    stat_data <- merge(stat_data,stat_plot,by="Factor",all.x = T)
    
    stats_results <- rbind(stats_results,stat_data%>%mutate(comparisons=paste(my_comparisons[[i]],collapse = "/")))
    plot_data_results <- rbind(plot_data_results,plot_data%>%mutate(comparisons=paste(my_comparisons[[i]],collapse = "/")))
  }
  
  
  sum_stats_results <- rbind(sum_stats_results,stats_results)
  sum_plot_data_results <- rbind(sum_plot_data_results,plot_data_results)
}
sum_stats_results <- sum_stats_results %>%
  mutate(p.adj = p.adjust(p,method="fdr")) %>%  
  mutate(lab.p.adj = case_when(               
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%  
  mutate(lab.p = case_when(              
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%
  mutate(label=paste(lab.p,lab.p.adj,sep = "/"))

write.csv(sum_stats_results,paste0(prefix,".phen.wilcox.csv"),quote = F, row.names = F)
write.csv(sum_plot_data_results,paste0(prefix,".phen.wilcox.box.csv"),quote = F, row.names = F)

## 1.3 combined ######################################
prefix <- "whole" 
sum_stats_results1 <- read.csv(paste0(target_dir,"/","Phen.nonK.phen.wilcox.csv"))
sum_stats_results2 <- read.csv(paste0(target_dir,"/","Phen.ONLYK.phen.wilcox.csv"))
sum_stats_results <- rbind(sum_stats_results1%>%mutate(Group="nonK"),
                           sum_stats_results2%>%mutate(Group="onlyK"))
write.csv(sum_stats_results,paste0("merged.",prefix,".phen.wilcox.csv"),quote = F, row.names = F)
#### heatmap
draw <- sum_stats_results %>%
  mutate(p.adj = p.adjust(p,method="fdr")) %>%   
  mutate(lab.p.adj = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%
  mutate(lab.p = case_when(   
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))%>%
  mutate(label=paste(lab.p,lab.p.adj,sep = "/"))

draw$label <- ifelse(draw$Mean.diff.log>0&draw$label!="ns/ns",paste0(draw$label,"\n","+"),
                     ifelse(draw$Mean.diff.log<0&draw$label!="ns/ns",paste0(draw$label,"\n","-"),as.vector(draw$label)))
draw$label <- gsub("ns/ns","",draw$label) 
draw$Cate <- phen.Cate$Category[match(draw$Factor,phen.Cate$Phen)]
draw$Cate <- factor(draw$Cate,levels = rev(c("K Score","K Sub Score","Sex Hormone","Glucose", "Lipid","Inflammation_Thyroxione", "Liver_Renal","Coagulation","CBC","Ultrasound")))
#
ggplot(draw,aes(comparisons,Factor,fill=Mean.diff.log))+geom_tile(color="lightgrey")+
  geom_text(aes(label = label),angle=90, size = 2, color = "black") +  
  facet_grid(
    Cate ~ ., 
    scales = "free_y", 
    space = "free_y" 
  ) +
  scale_fill_gradient2(  
    low =  "#0072B2", 
    mid = "white",  
    high = "#FF7F0E",  
    midpoint = 0   
  ) +
  theme_minimal() +   
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_text(angle = -180,hjust = 0,size=12),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 90,size=8), 
    panel.spacing.y = unit(0.5, "lines")     
  )+labs(x="",y="")
##
ggsave(paste0("Merged.wilcox.heamap.pdf"),width =7,height = 16)
