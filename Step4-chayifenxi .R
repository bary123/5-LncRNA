#Step3-1  ---提取Counts肿瘤样本------###
###----------------------------------###
{
  rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(TCGAbiolinks)
library(dplyr)
library(WGCNA)

opt <- paste0("Results/04-DEseq2/Step1/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

   ##加载表达数据
    load("./01-data/TCGA-UCEC_RNASeq_Counts.RData")#加载后的变量名称为：expDataCounts
    head(expDataCounts)[,1:3]

    ###删除正常样本的的列
      ##正常组织样本
        SamN <- TCGAquery_SampleTypes(barcode = colnames(expDataCounts),
                              typesample = c("NT","NB","NBC","NEBV","NBM"))
    ##肿瘤组织样本,并去除重复病人的数据
    SamT <- setdiff(colnames(expDataCounts),SamN)
    tur_expDataCounts <- expDataCounts[,SamT]
    dim(tur_expDataCounts)

    source("lcnRNAfenxi/00-fun/del_dup_sample.R")
    tur_expDataCounts <- del_dup_sample(tur_expDataCounts)
    dim(tur_expDataCounts)
    ##过滤不表达的基因在少数样本中表达的基因
    sat <- lapply(rownames(tur_expDataCounts) ,function(x){
       as.numeric(tur_expDataCounts[x,] != 0)
    })

    length(sat)
    dim(tur_expDataCounts)
    tur_expDataCounts<- tur_expDataCounts[lapply(sat,sum) > ncol(tur_expDataCounts)/2,] 
    dim(tur_expDataCounts)
    ##保存肿瘤样本
    save(tur_expDataCounts,file = paste0(opt,"UCEC_Step3_output.Rdata"))
}




#Step3-2  -----DESeq2方法差异分析-两两分组进行比较---###
###------------- Cluster1和Cluster2------------------###
{rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(dplyr)

opt <- paste0("Results/04-DEseq2/Cluster1-2/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

  proj <- "UCEC" #输入癌症类型
  message(proj)
  ##加载癌症样本分型
  subtype <- read.table(file=paste0("./Results/02-ConsensusClusterPlus/",
                                  proj,"_Clouster.txt"),row.names = 1)
  table(subtype)
  subtype  <- as.matrix(subtype)       
  subtype1 <- subtype  %>% as.data.frame()
  row.names(subtype1) <- rownames(subtype)
  colnames(subtype1)[1] <- "group"
  dim(subtype1)


  #subtype1 <-factor(subtype[,1],levels=c(1,2,3),labels = c("x","y","z")) %>% as.data.frame()
  
  
  ##分离出三个亚型的表达数据
  load("./Results/04-DEseq2/Step1/UCEC_Step3_output.Rdata")
  dim(tur_expDataCounts)
  gsExp0 <- tur_expDataCounts %>% t()  %>% as.data.frame()
  gsExp1 <- as.data.frame(lapply(gsExp0,as.integer))
  rownames(gsExp1) <-rownames(gsExp0)  
  gsExp1 <- gsExp1  %>% t()  %>% as.data.frame()
  dim(subtype1)
   
  
  #两两分组进行比较
  
  ###Cluster1和Cluster2
  subtype1group1 <- filter(subtype1,subtype1$group <3 )  
  table(subtype1group1)
  gsExp1group1 <-  gsExp1[,rownames(subtype1group1)]%>% as.data.frame()
  #subtype1group1 <-factor(subtype1group1[,1],levels=c(1,2),labels = c("x","y"))

  
  # 包的安装和加载
  library("DESeq2")
  #构建DESeqDataSet对象
  dds <- DESeqDataSetFromMatrix(countData =gsExp1group1, colData =subtype1group1,design = ~ group)
  dds <- DESeq(dds)# 函数分析差异
  
  sizeFactors(dds)# 计算标准化因子
  
  res <- results(dds)#提取差异表达结果
  class(res)
  res <- as.data.frame(res)

  ##----结果处理与保存
    res <- cbind(rownames(res), res)
  # 重命名列名
  colnames(res) <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                   "pval", "padj")
  dim(res)
  #保存文件到本地
  write.table(res,file=paste0(opt,proj,"-1_2-DESeq2.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
            na = "")

  ##--------差异基因筛选
  # 获取差异基因
  resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 1),]
  resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
  resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
  dim(resSig)
  write.table(resSig,file=paste0(opt,proj,"-1_2-DESeq2-1.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
  
  rm(resSig)
  #绘制火山图
  resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 0),]
  resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
  resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
  dim(resSig)
  write.table(resSig,file=paste0(opt,proj,"-1_2-DESeq2-2.gene.txt"),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
  
}
# ------火山图 -----###
###----------------------------------###
{
rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(dplyr)
source("lcnRNAfenxi/00-fun/del_dup_sample.R")
#创建输出地址
opt <- paste0("Results/04-DEseq2/Cluster1-2/")
#ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#加载差异分析结果
cla<- read.table("./Results/04-DEseq2/Cluster1-2/UCEC-1_2-DESeq2-2.gene.txt",row.names = 1)
cla <- cbind(rownames(cla),cla)
colnames(cla) <- cla[1,]
cla <- cla[-1,]

cla$pval <- as.numeric(as.character(cla$pval))
cla$log2FoldChange <- as.numeric(as.character(cla$log2FoldChange))

resSig <- cla
resSig$up_down <- "No"
resSig[which(abs(resSig$log2FoldChange) > 2), "up_down"] <- "Down"
resSig[which(resSig$log2FoldChange >2), "up_down"] <- "Up"
table(resSig$up_down)
dim(resSig)

#which(rownames(resSig)=='H19')
#head(resSig[2994,])
data = resSig
head(data)
data =data [!duplicated(data$gene_id),]
data <- data[,c(1,3,7,8)]
colnames(data) <- c("Gene","Log2FC","padj","Group")
data$padj <- as.numeric(as.character(data$padj))
data$Log2FC <- as.numeric(as.character(data$Log2FC))
head(data)


#找出纵轴的虚线位置
AAA <- filter(data,data$padj<0.05 & abs(data$Log2FC) > 2)
AAA$padj <- (-log10(AAA$padj))
Ylin <- min(AAA$padj)

colr <- brewer.pal(9,"Set1")[c(3,9,1)]    #设置颜色
png(file = paste0(opt,"UCEC_volcano plot.png"),width = 500,height = 600,)
 p <-ggplot(data = data,
            aes(x = Log2FC, 
                y = -log10(padj), 
                color = Group)) +  
  geom_point(alpha=0.8, size = 2) +  
  theme_bw(base_size = 20) +       #调整坐标轴大小
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +  
  scale_color_manual(values=colr)+
  scale_fill_manual(values =colr)+ 
  labs(x="Log2FC",y="-log10(padj)")+
  geom_vline(xintercept=c(-2,2), linetype="dotted")+   #画横坐标的虚线
  geom_hline(yintercept=Ylin,linetype="dotted")+            #画纵坐标的虚线
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),  #调整图注大小
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.direction = "horizontal",
        legend.position = c(0.5,0.96),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_blank())
    print(p,newpage = FALSE)
dev.off()

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggrepel)) install.packages("ggrepel")
if(!require(dplyr))install.packages("dplyr")
library(ggplot2)
library(ggrepel)
library(dplyr)

data$label=ifelse( -log10(data$padj) > 10 & abs(data$Log2FC) >=2,data$Gene,"")
table(data$label)

AB <- data[which( -log10(data$padj) > 10 & abs(data$Log2FC) >=2),]
table(AB$Group)
write.table(AB,file=paste0(opt,"AB.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

png(file = paste0(opt,"UCEC_volcano plot1.png"),width = 500,height = 600,)
ff <- p+geom_text_repel(data = data, aes(x = Log2FC, 
                                    y =-log10(padj), 
                                   label =label),
                  max.overlaps =100,
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
print(ff,newpage = FALSE)
dev.off()


png(file = paste0(opt,"UCEC_volcano plot2.png"),width = 500,height = 600,)
f1 <- p+geom_text_repel(data = data, aes(x = Log2FC, 
                                         y =-log10(padj), 
                                         label =label),
                        max.overlaps =100,
                        size = 3,box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE)+
                        ylim(0,80)
print(f1,newpage = FALSE)
dev.off()
}

###------------- Cluster1和Cluster3------------------###
{
  rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(dplyr)

opt <- paste0("Results/04-DEseq2/Cluster1-3/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

proj <- "UCEC" #输入癌症类型
message(proj)
##加载癌症样本分型
subtype <- read.table(file=paste0("./Results/02-ConsensusClusterPlus/",
                                  proj,"_Clouster.txt"),row.names = 1)
table(subtype)
subtype  <- as.matrix(subtype)       
subtype1 <- subtype  %>% as.data.frame()
row.names(subtype1) <- rownames(subtype)
colnames(subtype1)[1] <- "group"
dim(subtype1)


#subtype1 <-factor(subtype[,1],levels=c(1,2,3),labels = c("x","y","z")) %>% as.data.frame()


##分离出三个亚型的表达数据
load("./Results/04-DEseq2/Step1/UCEC_Step3_output.Rdata")
dim(tur_expDataCounts)
gsExp0 <- tur_expDataCounts %>% t()  %>% as.data.frame()
gsExp1 <- as.data.frame(lapply(gsExp0,as.integer))
rownames(gsExp1) <-rownames(gsExp0)  
gsExp1 <- gsExp1  %>% t()  %>% as.data.frame()
dim(subtype1)


#两两分组进行比较

###Cluster1和Cluster3
subtype1group2 <- filter(subtype1,subtype1$group != 2 )
table(subtype1group2)
gsExp1group2 <-  gsExp1[,rownames(subtype1group2)]%>% as.data.frame()
#subtype1group2 <-factor(subtype1group2[,1],levels=c(1,3),labels = c("x","z"))

# 包的安装和加载
library("DESeq2")
#构建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(countData =gsExp1group2, colData =subtype1group2,design = ~ group)
dds <- DESeq(dds)# 函数分析差异

sizeFactors(dds)# 计算标准化因子

res <- results(dds)#提取差异表达结果
class(res)
res <- as.data.frame(res)

##----结果处理与保存
res <- cbind(rownames(res), res)
# 重命名列名
colnames(res) <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                   "pval", "padj")
dim(res)
#保存文件到本地
write.table(res,file=paste0(opt,proj,"-1_3-DESeq2.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
            na = "")

##--------差异基因筛选
# 获取差异基因
resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 1),]
resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
dim(resSig)
write.table(resSig,file=paste0(opt,proj,"-1_3-DESeq2-1.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

rm(resSig)
#绘制火山图
resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 0),]
resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
dim(resSig)
#resSig <- resSig[!duplicated(resSig$padj),]

resSig <- na.omit(resSig)
head(resSig[1034,])
write.table(resSig,file=paste0(opt,proj,"-1_3-DESeq2-2.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

}
# ------火山图 -----###
###----------------------------------###

{
rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(dplyr)
source("lcnRNAfenxi/00-fun/del_dup_sample.R")
#创建输出地址
opt <- paste0("Results/04-DEseq2/Cluster1-3/")
#ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#加载差异分析结果
cla<- read.table("./Results/04-DEseq2/Cluster1-3/UCEC-1_3-DESeq2-2.gene.txt",row.names = 1)
cla <- cbind(rownames(cla),cla)
colnames(cla) <- cla[1,]
cla <- cla[-1,]

cla$pval <- as.numeric(as.character(cla$pval))
cla$log2FoldChange <- as.numeric(as.character(cla$log2FoldChange))

resSig <- cla
resSig$up_down <- "No"
resSig[which(abs(resSig$log2FoldChange) > 2), "up_down"] <- "Down"
resSig[which(resSig$log2FoldChange >2), "up_down"] <- "Up"
table(resSig$up_down)
dim(resSig)

which(rownames(resSig)=='H19')
head(resSig[2994,])
data = resSig
head(data)
data =data [!duplicated(data$gene_id),]
data <- data[,c(1,3,7,8)]
colnames(data) <- c("Gene","Log2FC","padj","Group")
data$padj <- as.numeric(as.character(data$padj))
data$Log2FC <- as.numeric(as.character(data$Log2FC))
head(data)


#找出纵轴的虚线位置
AAA <- filter(data,data$padj<0.05 & abs(data$Log2FC) > 2)
AAA$padj <- (-log10(AAA$padj))
Ylin <- min(AAA$padj)

colr <- brewer.pal(9,"Set1")[c(3,9,1)]    #设置颜色
png(file = paste0(opt,"UCEC_volcano plot.png"),width = 500,height = 600,)
p <-ggplot(data = data,
           aes(x = Log2FC, 
               y = -log10(padj), 
               color = Group)) +  
  geom_point(alpha=0.8, size = 2) +  
  theme_bw(base_size = 20) +       #调整坐标轴大小
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +  
  scale_color_manual(values=colr)+
  scale_fill_manual(values =colr)+ 
  labs(x="Log2FC",y="-log10(padj)")+
  geom_vline(xintercept=c(-2,2), linetype="dotted")+   #画横坐标的虚线
  geom_hline(yintercept=Ylin,linetype="dotted")+            #画纵坐标的虚线
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),  #调整图注大小
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.direction = "horizontal",
        legend.position = c(0.5,0.96),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_blank())
print(p,newpage = FALSE)
dev.off()

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggrepel)) install.packages("ggrepel")
if(!require(dplyr))install.packages("dplyr")
library(ggplot2)
library(ggrepel)
library(dplyr)



data$label=ifelse( -log10(data$padj) > 10 & abs(data$Log2FC) >=2,data$Gene,"")
table(data$label)

AC <- data[which( -log10(data$padj) > 10 & abs(data$Log2FC) >=2),]
table(AC$Group)
write.table(AC,file=paste0(opt,"AC.gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

png(file = paste0(opt,"UCEC_volcano plot1.png"),width = 500,height = 600,)
ff <- p+geom_text_repel(data = data, aes(x = Log2FC, 
                                         y =-log10(padj), 
                                         label =label),
                        max.overlaps =100,
                        size = 3,box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE)
print(ff,newpage = FALSE)
dev.off()


png(file = paste0(opt,"UCEC_volcano plot2.png"),width = 500,height = 600,)
f1 <- p+geom_text_repel(data = data, aes(x = Log2FC, 
                                         y =-log10(padj), 
                                         label =label),
                        max.overlaps =100,
                        size = 3,box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE)+
  ylim(0,80)
print(f1,newpage = FALSE)
dev.off()
}

###------------- Cluster2和Cluster3------------------###
{
  rm(list = ls())
  options(stringsAsFactors = F)
  #setwd("C:\\Users\\86147\\Desktop\\SingleGene")
  library(dplyr)
  
  opt <- paste0("Results/04-DEseq2/Cluster2-3/")#创建输出地址
  ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))
  
  proj <- "UCEC" #输入癌症类型
  message(proj)
  ##加载癌症样本分型
  subtype <- read.table(file=paste0("./Results/02-ConsensusClusterPlus/",
                                    proj,"_Clouster.txt"),row.names = 1)
  table(subtype)
  subtype  <- as.matrix(subtype)       
  subtype1 <- subtype  %>% as.data.frame()
  row.names(subtype1) <- rownames(subtype)
  colnames(subtype1)[1] <- "group"
  dim(subtype1)
  
  
  #subtype1 <-factor(subtype[,1],levels=c(1,2,3),labels = c("x","y","z")) %>% as.data.frame()
  
  
  ##分离出三个亚型的表达数据
  load("./Results/04-DEseq2/Step1/UCEC_Step3_output.Rdata")
  dim(tur_expDataCounts)
  gsExp0 <- tur_expDataCounts %>% t()  %>% as.data.frame()
  gsExp1 <- as.data.frame(lapply(gsExp0,as.integer))
  rownames(gsExp1) <-rownames(gsExp0)  
  gsExp1 <- gsExp1  %>% t()  %>% as.data.frame()
  dim(subtype1)
  
  
  #两两分组进行比较
  
  ###Cluster2和Cluster3
  subtype1group3 <- filter(subtype1,subtype1$group >1 )
  table(subtype1group3)
  gsExp1group3 <-  gsExp1[,rownames(subtype1group3)]%>% as.data.frame()
  #subtype1group3 <-factor(subtype1group3[,1],levels=c(2,3),labels = c("y","z"))
  
  
  # 包的安装和加载
  library("DESeq2")
  #构建DESeqDataSet对象
  dds <- DESeqDataSetFromMatrix(countData =gsExp1group3, colData =subtype1group3,design = ~ group)
  dds <- DESeq(dds)# 函数分析差异
  
  sizeFactors(dds)# 计算标准化因子
  
  res <- results(dds)#提取差异表达结果
  class(res)
  res <- as.data.frame(res)
  
  ##----结果处理与保存
  res <- cbind(rownames(res), res)
  # 重命名列名
  colnames(res) <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                     "pval", "padj")
  dim(res)
  #保存文件到本地
  write.table(res,file=paste0(opt,proj,"-2_3-DESeq2.gene.txt"),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
              na = "")
  
  ##--------差异基因筛选
  # 获取差异基因
  resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 1),]
  resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
  resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
  dim(resSig)
  write.table(resSig,file=paste0(opt,proj,"-2_3-DESeq2-1.gene.txt"),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
  
  rm(resSig)
  #绘制火山图
  resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 0),]
  resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
  resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
  dim(resSig)
  write.table(resSig,file=paste0(opt,proj,"-2_3-DESeq2-2.gene.txt"),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
  
}
  # ------火山图 -----###
  ###----------------------------------###
  
{
rm(list = ls())
  options(stringsAsFactors = F)
  library(ggplot2)
  library(RColorBrewer)
  library(ggbeeswarm)
  library(tidyverse)
  library(ggpubr)
  library(ggsci)
  library(dplyr)
  source("lcnRNAfenxi/00-fun/del_dup_sample.R")
  #创建输出地址
  opt <- paste0("Results/04-DEseq2/Cluster2-3/")
  #ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))
  
  #加载差异分析结果
  cla<- read.table("./Results/04-DEseq2/Cluster2-3/UCEC-2_3-DESeq2-2.gene.txt",row.names = 1)
  cla <- cbind(rownames(cla),cla)
  colnames(cla) <- cla[1,]
  cla <- cla[-1,]
  
  cla$pval <- as.numeric(as.character(cla$pval))
  cla$log2FoldChange <- as.numeric(as.character(cla$log2FoldChange))
  
  resSig <- cla
  resSig$up_down <- "No"
  resSig[which(abs(resSig$log2FoldChange) > 2), "up_down"] <- "Down"
  resSig[which(resSig$log2FoldChange >2), "up_down"] <- "Up"
  table(resSig$up_down)
  dim(resSig)
  
  which(rownames(resSig)=='H19')
  head(resSig[2994,])
  data = resSig
  head(data)
  data =data [!duplicated(data$gene_id),]
  data <- data[,c(1,3,7,8)]
  colnames(data) <- c("Gene","Log2FC","padj","Group")
  data$padj <- as.numeric(as.character(data$padj))
  data$Log2FC <- as.numeric(as.character(data$Log2FC))
  head(data)
  
  
  #找出纵轴的虚线位置
  AAA <- filter(data,data$padj<0.05 & abs(data$Log2FC) > 2)
  AAA$padj <- (-log10(AAA$padj))
  Ylin <- min(AAA$padj)
  
  colr <- brewer.pal(9,"Set1")[c(3,9,1)]    #设置颜色
  png(file = paste0(opt,"UCEC_volcano plot.png"),width = 500,height = 600,)
  p <-ggplot(data = data,
             aes(x = Log2FC, 
                 y = -log10(padj), 
                 color = Group)) +  
    geom_point(alpha=0.8, size = 2) +  
    theme_bw(base_size = 20) +       #调整坐标轴大小
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +  
    scale_color_manual(values=colr)+
    scale_fill_manual(values =colr)+ 
    labs(x="Log2FC",y="-log10(padj)")+
    geom_vline(xintercept=c(-2,2), linetype="dotted")+   #画横坐标的虚线
    geom_hline(yintercept=Ylin,linetype="dotted")+            #画纵坐标的虚线
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold"),    #调整图注大小
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.direction = "horizontal",
          legend.position = c(0.5,0.96),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = "black"),
          plot.background = element_blank())
  print(p,newpage = FALSE)
  dev.off()
  
  if(!require(ggplot2)) install.packages("ggplot2")
  if(!require(ggrepel)) install.packages("ggrepel")
  if(!require(dplyr))install.packages("dplyr")
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  data$label=ifelse( -log10(data$padj) > 10 & abs(data$Log2FC) >=2,data$Gene,"")
  table(data$label)
  BC <- data[which( -log10(data$padj) > 10 & abs(data$Log2FC) >=2),]
  table(BC$Group)
  write.table(BC,file=paste0(opt,"BC.gene.txt"),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
  
  
  png(file = paste0(opt,"UCEC_volcano plot1.png"),width = 500,height = 600,)
  ff <- p+geom_text_repel(data = data, aes(x = Log2FC, 
                                           y =-log10(padj), 
                                           label =label),
                          max.overlaps =100,
                          size = 3,box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.8, "lines"), 
                          segment.color = "black", 
                          show.legend = FALSE)
  print(ff,newpage = FALSE)
  dev.off()
  
  
  png(file = paste0(opt,"UCEC_volcano plot2.png"),width = 500,height = 600,)
  f1 <- p+geom_text_repel(data = data, aes(x = Log2FC, 
                                           y =-log10(padj), 
                                           label =label),
                          max.overlaps =100,
                          size = 3,box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.8, "lines"), 
                          segment.color = "black", 
                          show.legend = FALSE)+
                          ylim(0,80)
  print(f1,newpage = FALSE)
  dev.off()
}


###画3 Venn 图####
{
  rm(list = ls())

AB<- read.table("./Results/04-DEseq2/Cluster1-2/AB.gene.txt",row.names = 1)
AB <- cbind(rownames(AB),AB)
colnames(AB) <- AB[1,]
AB <- AB[-1,]
A <- AB$Gene

AC<- read.table("./Results/04-DEseq2/Cluster1-3/AC.gene.txt",row.names = 1)
AC <- cbind(rownames(AC),AC)
colnames(AC) <- AC[1,]
AC <- AC[-1,]
B <- AC$Gene

BC<- read.table("./Results/04-DEseq2/Cluster2-3/BC.gene.txt",row.names = 1)
BC <- cbind(rownames(BC),BC)
colnames(BC) <- BC[1,]
BC <- BC[-1,]
C <- BC$Gene

#Final_genes <- intersect(lasso_geneids,randomForestSRC_geneids)#取交集
library(VennDiagram)

Length_A<-length(A)
Length_B<-length(B)
Length_C<-length(C)
Length_AB<-length(intersect(A,B))
Length_BC<-length(intersect(B,C))
Length_AC<-length(intersect(A,C))
Length_ABC<-length(intersect(intersect(A,B),C))

T<-venn.diagram(list(Cluster1_2=A,## Cluster1和Cluster2筛选的结果
                     Cluster1_3=B,## Cluster1和Cluster3筛选的结果
                     Cluster2_3=C),## Cluster2和Cluster3筛选的结果
                
                filename = paste0('./Results/04-DEseq2/AB_AC_BC_venn.jpg'), ## 输出图片的名字
                  # 数字 number
                  cex = 2,  # 数字大小
                  fontface = "bold",  # 加粗
                  fontfamily = "sans",  # 字体
                  # 标签 category
                  cat.cex = 1.5,  # 字体大小
                  #cat.fontface = "bold",  # 加粗
                  cat.default.pos = "outer",  # 位置, outer 内 text 外
                  cat.pos = c(-27, 27, 180),  # 位置，用圆的度数
                  #cat.dist = c(0.055, 0.055, 0.085),  # 位置，离圆的距离
                  cat.fontfamily = "sans",  # 字体
                 # 圈
                lwd=1,           # 圈线条粗细 1 2 3 4 5
                lty=2,           # 线条类型, 1 实线, 2 虚线, blank 无线条
                col=c('red','green','blue'),     # 线条色
                fill=c('red','green','blue'),    # 填充色
                cat.col=c('red','green','blue'), # 字体色
                reverse=TRUE
              
                
        )
grid.draw(T)

}


