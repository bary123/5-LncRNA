###------------- Cluster1和Cluster2------------------###
{
  rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(clusterProfiler) #用来做富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)#这个包里存有人的注释文件
library(dplyr)

opt <- paste0("Results/05-GO/Cluster1-2/")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#读入差异基因
lsss<- read.table("./Results/04-DEseq2/Cluster1-2/UCEC-1_2-DESeq2-1.gene.txt",row.names = 1)
res <- cbind(rownames(lsss), lsss)
res <- na.omit(res) ## 去掉NA
colnames(res) <-res[1,] 
res <- res[-1,] 
names(res)

gene<- as.character(res$gene_id)  #获得基因 symbol ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = gene,
                       keytype = "SYMBOL",
                       column = "ENSEMBL") ##需要将symbolID转换成ENSEMBL
res<- mutate(res,gene_id=DEG.entrez_id)
ccc <- bitr(geneID=res$gene_id,fromType = "ENSEMBL",
                        toType = c("ENSEMBL","SYMBOL","ENTREZID"),
                        OrgDb = org.Hs.eg.db)
res <- merge(ccc,res,by.x="ENSEMBL",by.y="gene_id")# 数据融合

write.table(res,file=paste0(opt,"/id.xlsx"),sep="\t",quote=F,row.names=F)
ENTREZ_ID <- res$ENTREZID
go_BP <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="BP",pvalueCutoff = 0.01)
go_MF <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="MF",pvalueCutoff = 0.01)
go_CC <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="CC",pvalueCutoff = 0.01)

go_BPR <- go_BP@result
write.table(go_BPR,file=paste0(opt,"/1_2-BP.xlsx"),sep="\t",quote=F,row.names=F)
go_MFR <- go_MF@result
write.table(go_MFR,file=paste0(opt,"/1_2-MF.xlsx"),sep="\t",quote=F,row.names=F)
go_CCR <- go_CC@result
write.table(go_CCR,file=paste0(opt,"/1_2-CC.xlsx"),sep="\t",quote=F,row.names=F)

pdf(file=paste0(opt,"/1_2-BP图.pdf"),width = 10,height = 7)
barplot(go_BP,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_2-MF图.pdf"),width = 10,height = 7)
barplot(go_MF,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_2-CC图.pdf"),width = 10,height = 7)
barplot(go_CC,drop = TRUE,showCategory = 12)
dev.off()
}


###------------- Cluster1和Cluster3------------------###
{
  rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(clusterProfiler) #用来做富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)#这个包里存有人的注释文件
library(dplyr)

opt <- paste0("Results/05-GO/Cluster1-3")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#读入差异基因
lsss<- read.table("./Results/04-DEseq2/Cluster1-3/UCEC-1_3-DESeq2-1.gene.txt",row.names = 1)
res <- cbind(rownames(lsss), lsss)
res <- na.omit(res) ## 去掉NA
colnames(res) <-res[1,] 
res <- res[-1,] 
names(res)

gene<- as.character(res$gene_id)  #获得基因 symbol ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = gene,
                       keytype = "SYMBOL",
                       column = "ENSEMBL") ##需要将symbolID转换成ENSEMBL
res<- mutate(res,gene_id=DEG.entrez_id)
ccc <- bitr(geneID=res$gene_id,fromType = "ENSEMBL",
            toType = c("ENSEMBL","SYMBOL","ENTREZID"),
            OrgDb = org.Hs.eg.db)
res <- merge(ccc,res,by.x="ENSEMBL",by.y="gene_id")# 数据融合

write.table(res,file=paste0(opt,"/id.xlsx"),sep="\t",quote=F,row.names=F)
ENTREZ_ID <- res$ENTREZID
go_BP <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="BP",pvalueCutoff = 0.01)
go_MF <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="MF",pvalueCutoff = 0.01)
go_CC <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="CC",pvalueCutoff = 0.01)

go_BPR <- go_BP@result
write.table(go_BPR,file=paste0(opt,"/1_3-BP.xlsx"),sep="\t",quote=F,row.names=F)
go_MFR <- go_MF@result
write.table(go_MFR,file=paste0(opt,"/1_3-MF.xlsx"),sep="\t",quote=F,row.names=F)
go_CCR <- go_CC@result
write.table(go_CCR,file=paste0(opt,"/1_3-CC.xlsx"),sep="\t",quote=F,row.names=F)

pdf(file=paste0(opt,"/1_3-BP图.pdf"),width = 10,height = 7)
barplot(go_BP,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_3-MF图.pdf"),width = 10,height = 7)
barplot(go_MF,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/1_3-CC图.pdf"),width = 10,height = 7)
barplot(go_CC,drop = TRUE,showCategory = 12)
dev.off()


}

###------------- Cluster2和Cluster3------------------###
{
rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(clusterProfiler) #用来做富集分析
library(topGO)#画GO图用的
library(Rgraphviz)
library(pathview)#看KEGG pathway的
library(org.Hs.eg.db)#这个包里存有人的注释文件
library(dplyr)

opt <- paste0("Results/05-GO/Cluster2-3")#创建输出地址
ifelse(dir.exists(opt),FALSE,dir.create(opt,recursive = T))

#读入差异基因
lsss<- read.table("./Results/04-DEseq2/Cluster2-3/UCEC-2_3-DESeq2-1.gene.txt",row.names = 1)
res <- cbind(rownames(lsss), lsss)
res <- na.omit(res) ## 去掉NA
colnames(res) <-res[1,] 
res <- res[-1,] 
names(res)

gene<- as.character(res$gene_id)  #获得基因 symbol ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = gene,
                       keytype = "SYMBOL",
                       column = "ENSEMBL") ##需要将symbolID转换成ENSEMBL
res<- mutate(res,gene_id=DEG.entrez_id)
ccc <- bitr(geneID=res$gene_id,fromType = "ENSEMBL",
            toType = c("ENSEMBL","SYMBOL","ENTREZID"),
            OrgDb = org.Hs.eg.db)
res <- merge(ccc,res,by.x="ENSEMBL",by.y="gene_id")# 数据融合

write.table(res,file=paste0(opt,"/id.xlsx"),sep="\t",quote=F,row.names=F)
ENTREZ_ID <- res$ENTREZID
go_BP <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="BP",pvalueCutoff = 0.01)
go_MF <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="MF",pvalueCutoff = 0.01)
go_CC <- enrichGO(ENTREZ_ID,"org.Hs.eg.db",ont="CC",pvalueCutoff = 0.01)

go_BPR <- go_BP@result
write.table(go_BPR,file=paste0(opt,"/2_3-BP.xlsx"),sep="\t",quote=F,row.names=F)
go_MFR <- go_MF@result
write.table(go_MFR,file=paste0(opt,"/2_3-MF.xlsx"),sep="\t",quote=F,row.names=F)
go_CCR <- go_CC@result
write.table(go_CCR,file=paste0(opt,"/2_3-CC.xlsx"),sep="\t",quote=F,row.names=F)

pdf(file=paste0(opt,"/2_3-BP图.pdf"),width = 10,height = 7)
barplot(go_BP,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/2_3-MF图.pdf"),width = 10,height = 7)
barplot(go_MF,drop = TRUE,showCategory = 12)
dev.off()
pdf(file=paste0(opt,"/2_3-CC图.pdf"),width = 10,height = 7)
barplot(go_CC,drop = TRUE,showCategory = 12)
dev.off()

}




