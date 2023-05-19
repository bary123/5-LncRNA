
rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(TCGAbiolinks)
library(dplyr)
library(WGCNA)
source("lcnRNAfenxi/00-fun/del_dup_sample.R")

#创建输出地址
opt <- paste0("Results/01-getSig-lnRNA/")
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

###读入基因集
gss <- read.csv("lcnRNAfenxi/geneset/AA.csv",header = T)[,3] %>% unique()


##加载基因注释信息，加载后的变量名为：hsaGeneInfo
load("./02-base_files/hsaGeneInfo.Rdata")

###TCGA数据库中UCEC癌症数据
  load("./01-data/02-TCGA-RNASeq-TPM/TCGA-UCEC_RNASeq_TPM.RData") #加载表达数据expDataTPM
  dim(expDataTPM)
  head(expDataTPM)[,1:3]
  SamN <- TCGAquery_SampleTypes(barcode = colnames(expDataTPM),
                                typesample = c("NT","NB","NBC","NEBV","NBM"))##正常组织样本
  ####肿瘤组织样本,并去除重复病人的数据
  SamT <- setdiff(colnames(expDataTPM),SamN)
  tur_expDataTPM = expDataTPM[,SamT]
  source("lcnRNAfenxi/00-fun/del_dup_sample.R")
  tur_expDataTPM <- del_dup_sample(tur_expDataTPM)
  dim(tur_expDataTPM)
  ##过滤不表达的基因在少数样本中表达的基因
  sat <- lapply(rownames(tur_expDataTPM) ,function(x){
  as.numeric(tur_expDataTPM[x,] != 0)
  })
  
  length(sat)
  dim(tur_expDataTPM)
  tur_expDataTPM <- tur_expDataTPM[lapply(sat,sum) > ncol(tur_expDataTPM)/2,] 
  dim(tur_expDataTPM)

  ##读入的基因集gs与表达数据中的基因取交集
  gs <- intersect(gss,rownames(tur_expDataTPM))

  ##基因集的表达数据
  gsExp <- tur_expDataTPM[gs,] %>% t() %>% as.data.frame()
  head(gsExp)[,1:5]

##找出表达矩阵中的lnRNA
table(hsaGeneInfo$gene_type)
lncs <- hsaGeneInfo$symbol[hsaGeneInfo$gene_type == "lncRNA"]
lnc <- intersect(lncs,rownames(tur_expDataTPM))

##lncRNA的表达数据
lncExp <- tur_expDataTPM[lnc,] %>% t() %>% as.data.frame()
head(lncExp)[,1:5]
dim(lncExp)

###相关性分析
corcef <- cor(gsExp,lncExp,method = "pearson")
corp <- corPvalueStudent(corcef,nrow(lncExp))

corpval <-corp
##筛选出显著性lncRNA
cutPval <- 0.05
cutCef <- 0.4 ##绝对值

siglnc <- lapply(gs, function(x){
  colnames(corcef)[abs(corcef[x,]) > cutCef  & corpval[x,] < cutCef ] %>% na.omit()
}) %>% unlist() %>% unique()

siglncExp <- lncExp[,siglnc] %>% t() %>% as.data.frame()
dim(siglncExp)
save(gsExp,siglncExp,file =paste0(opt,"UCEC_Step1_output.Rdata"))
AAMlncRNA <- row.names(siglncExp)
write.table(AAMlncRNA,file=paste0(opt,"AAM-lncRNA_gene.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

