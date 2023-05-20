rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(survival)
library(dplyr)
library(survminer)

opt1 <- paste0("Results/03-survival/")#输出位置
ifelse(dir.exists(opt1),FALSE,dir.create(opt1,recursive = T))

  proj <- "UCEC"
  opt <- paste0("./Results/02-ConsensusClusterPlus/",proj,"_Clouster.txt")
  claa <- read.table(file = opt,row.names = 1)
  claa <- cbind(rownames(claa),claa)
  colnames(claa) <- c("cc","group")

  #加载生存数据
  cla <- read.table("./01-data/TCGA-UCEC.survival.tsv",row.names = 1)
  cla <- cbind(rownames(cla),cla)
  colnames(cla) <- cla[1,]
  cla <- cla[-1,]
  colnames(cla)[3] <- "PATIENT"
  cla <- cla[!duplicated(cla$PATIENT),]
  GC <- merge(claa,cla,by.x="cc",by.y="PATIENT")# 数据融合
  GC$OS<-as.numeric(GC$OS)
  GC$OS.time<-as.numeric(GC$OS.time)
  
  GC$group <- sub('1','Cluster1',GC$group)
  GC$group <- sub('2','Cluster2',GC$group)
  GC$group <- sub('3','Cluster3',GC$group)

  fit <- survfit(Surv(OS.time,OS)~group,data=GC)
  png(file = paste0(opt1,proj,"_survival.png"),width = 500,height = 600,)
  p <- ggsurvplot(fit,
                  risk.table = TRUE,
                  pval = TRUE,)
  p$plot <- p$plot + labs(title = proj)
  print(p,newpage = FALSE)
  dev.off()
  





