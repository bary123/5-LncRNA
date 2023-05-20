###--------------------------------------------------###
###------------- Cluster1和Cluster2------------------###
rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(clusterProfiler) #用来做富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)#这个包里存有人的注释文件
library(BiocManager)

opt <- paste0("Results/05-KEGG/Cluster1-2/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#读入差异基因
  res<- read.table("./Results/04-DEseq2/Cluster1-2/UCEC-1_2-DESeq2-1.gene.txt",row.names = 1)
  res <- cbind(rownames(res), res)
  colnames(res) <-res[1,] 
  res <- res[-1,] 
  names(res)

gene<- as.character(res$gene_id)  #获得基因 symbol ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =gene,
                       keytype = "SYMBOL",
                       column = "ENSEMBL") ##需要将symbolID转换成ENSEMBL
res<- mutate(res,gene_id=DEG.entrez_id)
ccc <- bitr(geneID=res$gene_id,fromType = "ENSEMBL",
                        toType = c("ENSEMBL","SYMBOL","ENTREZID"),
                        OrgDb = org.Hs.eg.db)
res <- merge(ccc,res,by.x="ENSEMBL",by.y="gene_id")# 数据融合
res=res[is.na(res[,"ENTREZID"])==F,]#去除基因id为NA的基因
ENTREZ_ID <- res$ENTREZID

kegg <- enrichKEGG(gene = res$ENTREZID,  #指定基因名称，例如 ENTREZ ID
                   keyType = 'kegg',  #指定基因名称类型，kegg 代表了使用 ENTREZ ID
                   organism = 'hsa',  #hsa代表人类，其它物种更改这行即可，可在 KEGG 官网查询物种名缩写
                   pvalueCutoff = 1,  #p 值 p 调整值之类的，指定 1 为阈值，也就是将所有的条目都输出了
                   qvalueCutoff = 1, 
                   pAdjustMethod = 'BH'  #p 值调整方法
)

pdf(file=paste0(opt,"/1_2-KEGG柱状图.pdf"),width = 10,height = 7)
barplot(kegg,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_2-KEGG点图.pdf"),width = 10,height = 7)
dotplot(kegg,showCategory = 12)
dev.off()

write.table(kegg,file=paste0(opt,"KEGGId.txt"),sep="\t",quote=F,row.names = F)#保存富集结果

###------------- Cluster1和Cluster3------------------###

rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(clusterProfiler) #用来做富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)#这个包里存有人的注释文件
library(BiocManager)

opt <- paste0("Results/05-KEGG/Cluster1-3/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#读入差异基因
res<- read.table("./Results/04-DEseq2/Cluster1-3/UCEC-1_3-DESeq2-1.gene.txt",row.names = 1)
res <- cbind(rownames(res), res)
colnames(res) <-res[1,] 
res <- res[-1,] 
names(res)

gene<- as.character(res$gene_id)  #获得基因 symbol ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =gene,
                       keytype = "SYMBOL",
                       column = "ENSEMBL") ##需要将symbolID转换成ENSEMBL
res<- mutate(res,gene_id=DEG.entrez_id)
ccc <- bitr(geneID=res$gene_id,fromType = "ENSEMBL",
            toType = c("ENSEMBL","SYMBOL","ENTREZID"),
            OrgDb = org.Hs.eg.db)
res <- merge(ccc,res,by.x="ENSEMBL",by.y="gene_id")# 数据融合
res=res[is.na(res[,"ENTREZID"])==F,]#去除基因id为NA的基因
ENTREZ_ID <- res$ENTREZID

kegg <- enrichKEGG(gene = res$ENTREZID,  #指定基因名称，例如 ENTREZ ID
                   keyType = 'kegg',  #指定基因名称类型，kegg 代表了使用 ENTREZ ID
                   organism = 'hsa',  #hsa代表人类，其它物种更改这行即可，可在 KEGG 官网查询物种名缩写
                   pvalueCutoff = 1,  #p 值 p 调整值之类的，指定 1 为阈值，也就是将所有的条目都输出了
                   qvalueCutoff = 1, 
                   pAdjustMethod = 'BH'  #p 值调整方法
)

pdf(file=paste0(opt,"/1_3-KEGG柱状图.pdf"),width = 10,height = 7)
barplot(kegg,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_3-KEGG点图.pdf"),width = 10,height = 7)
dotplot(kegg,showCategory = 12)
dev.off()

write.table(kegg,file=paste0(opt,"KEGGId.txt"),sep="\t",quote=F,row.names = F)#保存富集结果

###------------- Cluster2和Cluster3------------------###
rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(clusterProfiler) #用来做富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)#这个包里存有人的注释文件
library(BiocManager)

opt <- paste0("Results/05-KEGG/Cluster2-3/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#读入差异基因
res<- read.table("./Results/04-DEseq2/Cluster2-3/UCEC-2_3-DESeq2-1.gene.txt",row.names = 1)
res <- cbind(rownames(res), res)
colnames(res) <-res[1,] 
res <- res[-1,] 
names(res)

gene<- as.character(res$gene_id)  #获得基因 symbol ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =gene,
                       keytype = "SYMBOL",
                       column = "ENSEMBL") ##需要将symbolID转换成ENSEMBL
res<- mutate(res,gene_id=DEG.entrez_id)
ccc <- bitr(geneID=res$gene_id,fromType = "ENSEMBL",
            toType = c("ENSEMBL","SYMBOL","ENTREZID"),
            OrgDb = org.Hs.eg.db)
res <- merge(ccc,res,by.x="ENSEMBL",by.y="gene_id")# 数据融合
res=res[is.na(res[,"ENTREZID"])==F,]#去除基因id为NA的基因
ENTREZ_ID <- res$ENTREZID

kegg <- enrichKEGG(gene = res$ENTREZID,  #指定基因名称，例如 ENTREZ ID
                   keyType = 'kegg',  #指定基因名称类型，kegg 代表了使用 ENTREZ ID
                   organism = 'hsa',  #hsa代表人类，其它物种更改这行即可，可在 KEGG 官网查询物种名缩写
                   pvalueCutoff = 1,  #p 值 p 调整值之类的，指定 1 为阈值，也就是将所有的条目都输出了
                   qvalueCutoff = 1, 
                   pAdjustMethod = 'BH'  #p 值调整方法
)

pdf(file=paste0(opt,"2_3-KEGG柱状图.pdf"),width = 10,height = 7)
barplot(kegg,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_3-KEGG点图.pdf"),width = 10,height = 7)
dotplot(kegg,showCategory = 12)
dev.off()

write.table(kegg,file=paste0(opt,"KEGGId.txt"),sep="\t",quote=F,row.names = F)#保存富集结果



