library(forcats)
library(rstatix)
library(ggbeeswarm)

MHT.demo <- read_rds("D:/Demo-MHT/data/MHT.demo.RDS")
prof_filtered <- MHT.demo$MG$prof_filtered
Microbe.phen.prof <- MHT.demo$MG$Microbe.phen.prof

## set dir
target_dir <- paste0("D:/Demo-MHT/","Microbe_whole_feature")
dir.create(target_dir)
setwd(target_dir)

#  4 sites to compare BL ############################################################################################################
threshold_map <- list(
  "VA" = 0.0001,    # 阴道
  "UR" = 0.0001,    # 尿液
  "GUT" = 0.001,  # 肠道
  "TO" = 0.001    # 口腔
)
prefix <- "Paired.4sample"   
sum_stats_results <- c()
sum_plot_data_results <- c()

for (sites.vars in names(threshold_map)) {
  rel_ab.threshold <- threshold_map[[sites.vars]]
  xtabs(~Site+Time,Microbe.phen.prof)
  
  ####
  data <- prof_filtered[[sites.vars]] 
  ##
  taxa_to_keep <- sapply(data, function(x) {
    sum(x > rel_ab.threshold) > 0.1 * nrow(data)  
  })
  #  
  taxa.vars <- names(taxa_to_keep)[taxa_to_keep]
  data <- data[,taxa.vars]
  pseu_value <- min(data[data!=0])*0.1
  data <- log10(data+pseu_value)
  ##
  df <- data%>% 
    mutate(Response=Microbe.phen.prof$Time[match(rownames(data),Microbe.phen.prof$SeqID)],
           Clinic_ID=Microbe.phen.prof$Clinic_ID[match(rownames(data),Microbe.phen.prof$SeqID)])
  top_factors <- levels(df$Response)
  
  if (prefix == "Paired.4sample") {
    clinic_ids_4 <- df %>% group_by(Clinic_ID,Response)%>% group_by(Clinic_ID)%>%summarise(pairs = n())%>%filter(pairs==4)%>%pull(Clinic_ID)
    df <- df %>% subset(Clinic_ID %in% clinic_ids_4)
  }
  
  xtabs(~Response,df)
  ##
  my_comparisons <- list(c("BL","T04"),c("BL","T12"),c("BL","T24"))
  stats_results <- c()
  plot_data_results <- c()
  
  for (i in 1:length(my_comparisons)) {
    pairs.samp <- df %>%subset(Response%in%my_comparisons[[i]])%>% group_by(Clinic_ID)%>%summarise(pairs = n())%>%filter(pairs==2)%>%pull(Clinic_ID)
    plot_data <- df %>%
      subset(Clinic_ID %in% pairs.samp & Response%in%my_comparisons[[i]])%>%
      arrange(Clinic_ID,Response)%>%
      select(Response,Clinic_ID,taxa.vars) %>%
      reshape2::melt(id.vars = c("Response","Clinic_ID"),variable.name = "Factor",value.name = "Value")%>%
      mutate(Response = fct_relevel(Response, top_factors)) 
    head(plot_data )
    
    
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
             Mean.diff=Mean2-Mean1)%>%
      mutate(Mean.diff.log=sign(Mean.diff) * log1p(abs(Mean.diff)))%>%
      dplyr::select(Factor,Mean1,Mean2,Mean.diff.log)
    
    stat_data <- merge(stat_data,stat_plot,by="Factor",all.x = T)
    
    stats_results <- rbind(stats_results,stat_data%>%mutate(comparisons=paste(my_comparisons[[i]],collapse = "/")))
    plot_data_results <- rbind(plot_data_results,plot_data%>%mutate(comparisons=paste(my_comparisons[[i]],collapse = "/")))
  }
  stats_results <- stats_results %>%
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
  ##
  stats_results$Sites <- sites.vars
  
  sum_stats_results <- rbind(sum_stats_results,stats_results)
  sum_plot_data_results <- rbind(sum_plot_data_results,plot_data_results%>%mutate(Sites=sites.vars))
}

write.csv(sum_stats_results,paste0(prefix,".TimeDiff.wilcox.csv"),quote = F, row.names = F)
# write.csv(sum_plot_data_results,paste0(prefix,".TimeDiff.wilcox.box.csv"),quote = F, row.names = F)

## heatmap ###################
draw <- sum_stats_results%>%subset(p.adj<0.05)%>%select(Factor,Sites)%>%
  merge(sum_stats_results)%>%.[!duplicated(.),]

draw$groups <- paste(draw$group1,draw$group2,sep = "/")
draw$label <- gsub("ns/ns","",draw$label) 
draw%>%select(Sites,n1:n2,groups)%>%.[!duplicated(.),]
write.csv(draw,paste0(prefix,".TimeDiff.wilcox.padj.0.05.csv"),quote = F, row.names = F)

##
prefix <- "Paired.4sample" 
draw <-  read.csv(paste0(prefix,".TimeDiff.wilcox.padj.0.05.csv"))
ggplot(draw,aes(groups,Factor,fill=Mean.diff.log))+geom_tile()+
  geom_text(aes(label = label), size = 2.5, color = "black") +  
  facet_grid(
    Sites ~ ., 
    scales = "free_y", 
    space = "free_y"  
  ) +
  scale_fill_gradient2(  
    low = "#5e3c99",    
    mid = "white",    
    high = "#b35806",  
    midpoint = 0       
  ) +
  theme_minimal() +     
  theme(
    axis.text.y = element_text(face = "bold.italic"),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = -90),  
    panel.spacing.y = unit(0.5, "lines")    
  )
##
ggsave(paste0(prefix,".TimeDiff.wilcox.heatmap.pdf"),width = 6,height = 9)


# whole shannon ##########################################################################################################
phen.vars <- "Shannon"
phen.input <- MHT.demo$MG$diversity
####
prefix <- "Shannon" 
comb_matrix <- combn(c("BL","T04","T12","T24"), 2)
my_comparisons <- lapply(split(t(comb_matrix), 1:ncol(comb_matrix)), as.vector)
top_factors <- "BL" 

####
data <- phen.input[,c("Time","Clinic_ID","Site",phen.vars)]%>%
  rename(Response=Time)
xtabs(~Response,data)

sum_stats_results <- c()
sum_plot_data_results <- c()
for (sites.vars in names(threshold_map)) {   
  stats_results <- c()
  plot_data_results <- c()
  
  for (i in 1:length(my_comparisons)) {
    df <- data[,c("Clinic_ID","Response","Site",phen.var)]%>%subset(Site==sites.vars)%>%na.omit()
    pairs.samp <- df %>%subset(Response%in%my_comparisons[[i]])%>% group_by(Clinic_ID)%>%summarise(pairs = n())%>%filter(pairs==2)%>%pull(Clinic_ID)
    plot_data <- df %>%
      subset(Clinic_ID %in% pairs.samp & Response%in%my_comparisons[[i]])%>%
      arrange(Clinic_ID,Response)%>%
      select(Response,Clinic_ID,phen.var) %>%
      reshape2::melt(id.vars = c("Response","Clinic_ID"),variable.name = "Factor",value.name = "Value")%>%
      mutate(Response = fct_relevel(Response, top_factors)) 
    head(plot_data )
    
    
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
  stats_results$Sites <- sites.vars
  
  sum_stats_results <- rbind(sum_stats_results,stats_results)
  sum_plot_data_results <- rbind(sum_plot_data_results,plot_data_results)
}
sum_stats_results <- sum_stats_results %>%
  mutate(p.adj = p.adjust(p,method="fdr")) %>%   # P adj
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

write.csv(sum_stats_results,paste0(prefix,".TimeDiff.wilcox.csv"),quote = F, row.names = F)
