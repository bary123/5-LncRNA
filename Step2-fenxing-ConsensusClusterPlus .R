###ConsensusClusterPlus进行分型
rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(TCGAbiolinks)
library(dplyr)
library(WGCNA)

proj <- "UCEC"
opts <- paste0("Results/02-ConsensusClusterPlus/")#输出位置
ifelse(dir.exists(opts),FALSE,dir.create(opts,recursive = T))

load("./Results/01-getSig-lnRNA/UCEC_Step1_output.Rdata")#加载AAM相关LncRNA

#准备工作
library(ConsensusClusterPlus)#执行样本分型
results <- ConsensusClusterPlus(as.matrix(siglncExp),
                                maxK=9,
                                reps= 100,
                                pItem= 0.8,
                                pFeature=1,
                                title= opts,
                                clusterAlg="hc",
                                distance="pearson",
                                innerLinkage = "ward.D2",
                                seed=123456,
                                plot="pdf")
#提取每一个cluster
icl <- calcICL(results, title =opts ,plot = "png")

#输出结果
clusterNum <- 3 #分几类，根据判断标准判断
cluster <- results[[clusterNum]][["consensusClass"]]
write.table(cluster,file=paste(opts,"/",proj,"_Clouster.txt",sep = ""),
            sep="\t",quote=F,col.names=F)



