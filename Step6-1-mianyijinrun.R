#Step3------免疫浸润
#----------------
rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
##===========
library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(tibble)
library(ggpubr)
library(pheatmap)
library(ggrepel)
#source("./daima/00_function.R")

###加载所有癌症样本数据
timer = read.table("./01-data/04-infiltrationTIMER2/infiltration_estimation_for_tcga.csv",
                   sep = ",",
                   header = T,
                   check.names = F)

#读取输入文件，并对输入文件处理

proj = "UCEC"
cancer = "UCEC"
all_immue_eatimate = data.frame()
allcancer_ssGSEA_Score = data.frame()
gene_immCell_cor = data.frame()#ssGSEA
gene_immCell_cor2 = data.frame()#TIMER2中的几种方法

###====================ssGSEA，TIMER
#文件输出路径
outdir = paste0("Results/08-mianyifenxi/Step3/")
ifelse(dir.exists(outdir),print("文件夹已经存在"),dir.create(outdir))


color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")

  
  ssGSEA_Score = read.table(paste0("01-data/TCGA-UCEC-norm_turmor_ssGSEA_Score.txt"),
                            header = T,sep = "\t",check.names = F)
  rownames(ssGSEA_Score) = ssGSEA_Score$id
  colnames(ssGSEA_Score) <- substring(colnames(ssGSEA_Score),1,12)
  
  
  load("./Results/08-mianyifenxi/Step1/TCGA-UCEC_gene_output.Rdata")
  tur_expDataTPM <- siglncExp %>% t() %>% as.data.frame()
  rm(gsExp,siglncExp)
  tur_expDataTPM <- cbind(ID=rownames(tur_expDataTPM),tur_expDataTPM)
  #加载癌症样本分型
  subtype <- read.table(file="./Results/08-mianyifenxi/Step2/UCEC_Clouster.txt",row.names = 1)
  table(subtype)
  subtype  <- as.matrix(subtype)       
  subtype1 <-factor(subtype[,1],levels=c(1,2,3),labels = c("Cluster1","Cluster2","Cluster3")) %>% as.data.frame()
  row.names(subtype1) <- rownames(subtype)
  colnames(subtype1)[1] <- "group"
  subtype1 <- cbind(ID=rownames(subtype1),subtype1)
  gene_meta  <- merge(tur_expDataTPM,subtype1,by="ID")# 数据融合
  rm(subtype,subtype1)
  rownames(gene_meta) = gene_meta[,1]
  gene_meta <- gene_meta[,-1]
  
  coid = intersect(rownames(gene_meta), colnames(ssGSEA_Score))
  gene_meta = gene_meta[coid,]
  rm(coid)
  gene_meta = arrange(gene_meta,group)
  ssGSEA_Score = ssGSEA_Score[,rownames(gene_meta)]
  #outdirp = paste0(outdir)   #输出位置
  #ifelse(dir.exists(outdirp),
         #print("文件夹已经存在"),
         #dir.create(outdirp))
  
  
  which(colnames(gene_meta)=='group') #查看分组信息
  samples_meta = data.frame(group = gene_meta[,61])# 获取分组信息
  rownames(samples_meta) = rownames(gene_meta)
  ann_colors = list(group = c( Cluster1= "#E41A1C",Cluster2 = "#4DAF4A",Cluster3 ="#3300CC"))
  
  
  f1 = paste0(outdir,"/",proj,"-ssGSEA-heatmap.pdf")
  pdf(f1,height=3.5,width=9)
  pheatmap(ssGSEA_Score,
           annotation_col = samples_meta,
           color =colorRampPalette(color.key)(50),
           cluster_cols =F,
           fontsize=8,
           fontsize_row=8,
           scale="row",
           show_colnames=F,
           annotation_colors = ann_colors,
           fontsize_col=3)
  dev.off()
  norm_ssGSEA_Score2 =  ssGSEA_Score %>% t() %>% as.data.frame()
  norm_ssGSEA_Score2 = cbind(samples_meta,norm_ssGSEA_Score2)
  
  ldata <- melt(norm_ssGSEA_Score2,id.vars = "group")
  head(ldata)
  
  write.csv(ldata,file = paste0(outdir,proj,"-eatimate.csv"))
  
  
  
  p = ggplot(ldata, aes(x=variable, y=value,color = group)) +
    #geom_jitter(size = 0.5,aes(x=variable, y=value,color = Group))+
    geom_boxplot(aes(color = group),alpha =1,
                 lwd = 0.5, outlier.size = 1,
                 outlier.colour = "white")+ #color = c("red", "blue"),
    theme_bw()+
    stat_compare_means(label = "p.signif") +
    #labs(title = 'ImmuneCellAI') +
    theme(legend.position = "top",
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.title =element_blank(),
          legend.text = element_text(size = 12, face = "bold",colour = "#1A1A1A")
    ) +
    ylab("ssGSEA score") +
    xlab("") +
    #ylim(c(0.45,1.6))+
    scale_color_manual(values = c("#E41A1C","#3300CC","#EE7700")) #+ ,
  
  
  pdf(paste0(outdir,"/",proj,"-ssGSEA-boxplot.pdf"),height=5.5,width=9)
  print(p)
  dev.off()
  rm(ldata)
  
  
  ssGSEA_Score = ssGSEA_Score %>% t() %>% as.data.frame()
  ssGSEA_Score = rownames_to_column(ssGSEA_Score,var = "id")
  ssGSEA_Score = add_column(ssGSEA_Score,cancer = proj,.before = 1)
  gene_meta = rownames_to_column(gene_meta,var = "id")
  ssGSEA_Score = merge(gene_meta,ssGSEA_Score,by = "id")
  allcancer_ssGSEA_Score = rbind(allcancer_ssGSEA_Score, ssGSEA_Score)
  
  
  which(colnames(gene_meta)=='group') #查看分组信息
  tpm = gene_meta[,c(-1,-62)]  #去掉基因行名，分组信息，只留基因列
  one_cell_df = data.frame() 
  which(colnames(ssGSEA_Score)=='Activated B cell')  #查看去掉Activated B cell前面的列：cell in colnames(ssGSEA_Score)[-c(1:732)]
  for(t in colnames(tpm)){
      tt = tpm[,t]
    for(cell in colnames(ssGSEA_Score)[-c(1:63)]){
        immScore = ssGSEA_Score[,cell]
        cor.spearman = cor.test(tt,immScore,method ="spearman")
        df = data.frame(gene = t,cell = cell,
                    spearman= cor.spearman[["estimate"]],
                    p.value = cor.spearman[["p.value"]])
    one_cell_df = rbind(one_cell_df,df)
    }
  }
  
  
  gene_immCell_cor = rbind(gene_immCell_cor,one_cell_df)
  
  ###ssGSEA方法基因分组
  gene_ssGSEA<- gene_immCell_cor
  ssGSEA <- gene_ssGSEA[!duplicated(gene_ssGSEA$cell),]  #去重
  ssGSEA$method <- "ssGSEA"
  ssGSEA <- ssGSEA[,c(2,5)]
  
  
  cell_subd<- as.data.frame(p$data)
  #unique(cell_subd)
  head(cell_subd[1:5,1:3])
  colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
  head(cell_subd[1:5,1:3])
  compaired <- list(c("Cluster1", "Cluster2"),
                    c("Cluster1", "Cluster3"),
                    c("Cluster2", "Cluster3"))
 
  #Method$cell<- sub('Macrophage/Monocyte','Macrophage-Monocyte',Method$cell)
  
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/ssGSEA/")
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))

    for(tt in unique(ssGSEA$cell)){
      #  x<- Activated B cell
      message(tt)
      cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
      fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
        # geom_boxplot(aes(color = Sample_Origin),alpha =1,
        #              lwd = 0.5, outlier.size = 1,
        #              outlier.colour = "white")+ #color = c("red", "blue"),
        geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
        geom_boxplot(width = 0.1,fill = "white")+
        theme_bw()+
        scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
        labs(title = tt,x ="ssGSEA",y = "Immunoinfiltration scores")+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test = wilcox.test)+
        theme(
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.position="none"
        )
      #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
      pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
      print(fig)
      dev.off()
    }
    
  
    
##-----------------------------
  class(samples_meta$group)
  #samples_meta$group = factor(samples_meta$group,levels = c("Cluster1","Cluster2","Cluster3"))
  id =rownames(samples_meta)
  samples_meta$cell_type = id
  samples_meta = samples_meta[!duplicated(samples_meta$cell_type),]
  
  timer$cell_type <- substring(timer$cell_type,1,12)
  conp = intersect(samples_meta$cell_type,timer$cell_type)
  onetimer = timer[timer$cell_type %in% conp,]
  onetimer = merge(samples_meta,onetimer,by = "cell_type")
  onetimer = onetimer[!duplicated(onetimer$cell_type),]
  rownames(onetimer) = onetimer$cell_type
  immue_eatimate = add_column(onetimer,cancer,.before = 1)
  colnames(immue_eatimate)[2] = "id"
  
  which(colnames(immue_eatimate)=='group')#查看并去掉分组信息
  immue_eatimate = immue_eatimate[,-3]
  immue_eatimate = merge(gene_meta,immue_eatimate,by = "id")
  all_immue_eatimate = rbind(all_immue_eatimate,immue_eatimate)
  
  
  iedf = data.frame()
  immue_eatimate = immue_eatimate %>% t() %>% na.omit() %>% t() %>% as.data.frame()
  which(colnames(immue_eatimate)=='B cell_TIMER')  #查看并去掉B cell_TIMER前面的列
  
  for(t in colnames(tpm)){
    tt = tpm[,t]
    for(cell_m in colnames(immue_eatimate)[-c(1:63)]){
        imm = immue_eatimate[,cell_m]
        cor.spearman = cor.test(tt,
                                as.numeric(imm),
                                method ="spearman")
        na.omit(c(NA,12,431))
          df = data.frame(gene=t ,
                    cell = unlist(strsplit(cell_m,"_" ))[1],
                    method = unlist(strsplit(cell_m,"_" ))[2],
                    spearman= cor.spearman[["estimate"]],
                    p.value = cor.spearman[["p.value"]])
        iedf = rbind(iedf,df)
    }
  }
  
  gene_immCell_cor2 =rbind(gene_immCell_cor2,iedf)
  onetimer = onetimer[,-1]
  colnames(onetimer)[1] = "gene"
  
  ###TIMER方法基因分组
  gene_TIMER<- gene_immCell_cor2
  TIMER <- gene_TIMER[,c(2,3)]
  
  ##根据基因高低表达分组，用于绘制从TIMER2下载的数据
  plotbox_TIMER = function(data,meth,proj){
    coln = colnames(data)[grep(paste0("_",meth,"$"),colnames(data))]
    onetimer_meth = data[,c("gene",coln)]
    colnames(onetimer_meth) = gsub(paste0("_",meth),"",colnames(onetimer_meth))
    onetimer_meth <- melt(onetimer_meth,id.vars = "gene")
    dim(onetimer_meth)
    onetimer_meth = onetimer_meth[onetimer_meth$variable !="uncharacterized cell",]
    onetimer_meth = onetimer_meth[onetimer_meth$variable !="Cancer associated fibroblast",]
    fig = ggplot(onetimer_meth, aes(x=variable, y=value,color = gene)) +
      #geom_jitter(size = 0.5,aes(x=variable, y=value,color = Group))+
      geom_boxplot(aes(color = gene),alpha =1,
                   lwd = 0.6, outlier.size = 1,
                   outlier.colour = "white")+ #color = c("red", "blue"),
      theme_bw()+
      stat_compare_means(label = "p.signif") +
      ylab("Cell immune infltration") +
      xlab(paste0(proj," dataset from TCGA ","(",meth,")")) +
      scale_color_manual(values = c("#E41A1C","#4DAF4A","#3300CC"))+ #+ ,"#EAE838"
      theme(legend.position = "top",
            plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
            axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
            axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
            legend.text = element_text(size = 12, face = "bold",colour = "#1A1A1A")
      )
    return(fig)
  }
  
  
  pdf(paste0(outdir,"/",proj,"-TIMER2-CIBERSORT-ABS-箱形图.pdf"),height=5.5,width=9)
  p <- plotbox_TIMER(data = onetimer,meth = "CIBERSORT-ABS",proj = proj)+ylim(0,0.2) 
  print(p)
  dev.off()
 
  {  
  cell_subd<- as.data.frame(p$data)
  #unique(cell_subd)
  head(cell_subd[1:5,1:3])
  colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
  head(cell_subd[1:5,1:3])
  
  #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
  #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
  TIMER1=TIMER[TIMER$method == "CIBERSORT-ABS",]
    AA <- TIMER1[1,2]
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
  for(tt in unique(TIMER1$cell)){
    #  x<- Activated B cell
    
    message(tt)
    cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
    fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
      # geom_boxplot(aes(color = Sample_Origin),alpha =1,
      #              lwd = 0.5, outlier.size = 1,
      #              outlier.colour = "white")+ #color = c("red", "blue"),
      geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
      geom_boxplot(width = 0.1,fill = "white")+
      theme_bw()+
      scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
      labs(title = tt,x ="CIBERSORT-ABS",y = "Immunoinfiltration scores")+
      geom_signif(comparisons = compaired,
                  step_increase = 0.1,
                  map_signif_level = T,
                  test = wilcox.test)+
      theme(
        plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
        axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
        axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
        legend.position="none"
      )
    #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
    pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
    print(fig)
    dev.off()
    
}
  }  
    
    
  pdf(paste0(outdir,"/",proj,"-TIMER2-CIBERSORT-箱形图.pdf"),height=5.5,width=9)
  p <- plotbox_TIMER(data = onetimer,meth = "CIBERSORT",proj = proj)+ylim(0,0.35)
  print(p)
  dev.off()
  
  
 {
   cell_subd<- as.data.frame(p$data)
  #unique(cell_subd)
  head(cell_subd[1:5,1:3])
  colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
  head(cell_subd[1:5,1:3])
  
  #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
  #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
  TIMER1=TIMER[TIMER$method == "CIBERSORT",]
  AA <- TIMER1[1,2]
  outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
  ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
  for(tt in unique(TIMER1$cell)){
    #  x<- Activated B cell
    
    message(tt)
    cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
    fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
      # geom_boxplot(aes(color = Sample_Origin),alpha =1,
      #              lwd = 0.5, outlier.size = 1,
      #              outlier.colour = "white")+ #color = c("red", "blue"),
      geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
      geom_boxplot(width = 0.1,fill = "white")+
      theme_bw()+
      scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
      labs(title = tt,x = AA,y = "Immunoinfiltration scores")+
      geom_signif(comparisons = compaired,
                  step_increase = 0.1,
                  map_signif_level = T,
                  test = wilcox.test)+
      theme(
        plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
        axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
        axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
        legend.position="none"
      )
    #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
    pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
    print(fig)
    dev.off()
    
  }
  
  }
  
  
  
  pdf(paste0(outdir,"/",proj,"-TIMER2-EPIC箱形图.pdf"),height=5.5,width=6)
  p <- plotbox_TIMER(data = onetimer,meth = "EPIC",proj = proj)+ylim(0,0.1)
  print(p)
  dev.off()
  
  {
    cell_subd<- as.data.frame(p$data)
    #unique(cell_subd)
    head(cell_subd[1:5,1:3])
    colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
    head(cell_subd[1:5,1:3])
    
    #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
    #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    TIMER1=TIMER[TIMER$method == "EPIC",]
    AA <- TIMER1[1,2]
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    for(tt in unique(TIMER1$cell)){
      #  x<- Activated B cell
      
      message(tt)
      cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
      fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
        # geom_boxplot(aes(color = Sample_Origin),alpha =1,
        #              lwd = 0.5, outlier.size = 1,
        #              outlier.colour = "white")+ #color = c("red", "blue"),
        geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
        geom_boxplot(width = 0.1,fill = "white")+
        theme_bw()+
        scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
        labs(title = tt,x = AA,y = "Immunoinfiltration scores")+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test = wilcox.test)+
        theme(
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.position="none"
        )
      #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
      pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
      print(fig)
      dev.off()
      
    }
    
  }
  
  
  pdf(paste0(outdir,"/",proj,"-TIMER2-MCPCOUNTER箱形图.pdf"),height=5.5,width=7)
  p <- plotbox_TIMER(data = onetimer,meth = "MCPCOUNTER",proj = proj)+ylim(0,20)
  print(p)
  dev.off()
  
  {
    cell_subd<- as.data.frame(p$data)
    #unique(cell_subd)
    head(cell_subd[1:5,1:3])
    colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
    head(cell_subd[1:5,1:3])
    
    #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
    #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    TIMER1=TIMER[TIMER$method == "MCPCOUNTE",]
    AA <- TIMER1[1,2]
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    for(tt in unique(TIMER1$cell)){
      #  x<- Activated B cell
      
      message(tt)
      cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
      fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
        # geom_boxplot(aes(color = Sample_Origin),alpha =1,
        #              lwd = 0.5, outlier.size = 1,
        #              outlier.colour = "white")+ #color = c("red", "blue"),
        geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
        geom_boxplot(width = 0.1,fill = "white")+
        theme_bw()+
        scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
        labs(title = tt,x = AA,y = "Immunoinfiltration scores")+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test = wilcox.test)+
        theme(
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.position="none"
        )
      #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
      pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
      print(fig)
      dev.off()
      
    }
    
  } 
  
  pdf(paste0(outdir,"/",proj,"-TIMER2-QUANTISEQ箱形图.pdf"),height=5.5,width=7)
  p <- plotbox_TIMER(data = onetimer,meth = "QUANTISEQ",proj = proj)+ylim(0,0.075)
  print(p)
  dev.off()
  
  
  {
    cell_subd<- as.data.frame(p$data)
    #unique(cell_subd)
    head(cell_subd[1:5,1:3])
    colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
    head(cell_subd[1:5,1:3])
    
    #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
    #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    TIMER1=TIMER[TIMER$method == "QUANTISEQ",]
    AA <- TIMER1[1,2]
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    for(tt in unique(TIMER1$cell)){
      #  x<- Activated B cell
      
      message(tt)
      cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
      fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
        # geom_boxplot(aes(color = Sample_Origin),alpha =1,
        #              lwd = 0.5, outlier.size = 1,
        #              outlier.colour = "white")+ #color = c("red", "blue"),
        geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
        geom_boxplot(width = 0.1,fill = "white")+
        theme_bw()+
        scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
        labs(title = tt,x = AA,y = "Immunoinfiltration scores")+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test = wilcox.test)+
        theme(
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.position="none"
        )
      #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
      pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
      print(fig)
      dev.off()
      
    }
    
  } 
  
  
  pdf(paste0(outdir,"/",proj,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
  p <- plotbox_TIMER(data = onetimer,meth = "TIMER",proj = proj)+ylim(0,0.75)
  print(p)
  dev.off()
  
  {
    cell_subd<- as.data.frame(p$data)
    #unique(cell_subd)
    head(cell_subd[1:5,1:3])
    colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
    head(cell_subd[1:5,1:3])
    
    #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
    #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    TIMER1=TIMER[TIMER$method == "TIMER",]
    AA <- TIMER1[1,2]
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    for(tt in unique(TIMER1$cell)){
      #  x<- Activated B cell
      
      message(tt)
      cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
      fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
        # geom_boxplot(aes(color = Sample_Origin),alpha =1,
        #              lwd = 0.5, outlier.size = 1,
        #              outlier.colour = "white")+ #color = c("red", "blue"),
        geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
        geom_boxplot(width = 0.1,fill = "white")+
        theme_bw()+
        scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
        labs(title = tt,x = AA,y = "Immunoinfiltration scores")+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test = wilcox.test)+
        theme(
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.position="none"
        )
      #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
      pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
      print(fig)
      dev.off()
      
    }
    
  } 
  
  pdf(paste0(outdir,"/",proj,"-TIMER2-XCELL-箱形图.pdf"),height=5.5,width=15)
  p <- plotbox_TIMER(data = onetimer,meth = "XCELL",proj = proj)+ylim(0,0.3)
  print(p)
  dev.off()
  
  {
    cell_subd<- as.data.frame(p$data)
    #unique(cell_subd)
    head(cell_subd[1:5,1:3])
    colnames(cell_subd) <- c("Sample_Origin","Cell_subtype","Value")
    head(cell_subd[1:5,1:3])
    
    #outdirP2= paste0("Results/08-mianyifenxi/Step4-2/CIBERSORT-ABS")
    #ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    TIMER1=TIMER[TIMER$method == "XCELL",]
    AA <- TIMER1[1,2]
    outdirP2= paste0("Results/08-mianyifenxi/Step4-2/",AA)
    ifelse(dir.exists(outdirP2),print("文件夹已经存在"),dir.create(outdirP2))
    for(tt in unique(TIMER1$cell)){
      #  x<- Activated B cell
      
      message(tt)
      cell_subd1=cell_subd[cell_subd$Cell_subtype == tt,]
      fig = ggplot(cell_subd1, aes(x= Sample_Origin, y= Value,color = Sample_Origin)) +
        # geom_boxplot(aes(color = Sample_Origin),alpha =1,
        #              lwd = 0.5, outlier.size = 1,
        #              outlier.colour = "white")+ #color = c("red", "blue"),
        geom_violin(aes(fill = Sample_Origin),trim = FALSE)+
        geom_boxplot(width = 0.1,fill = "white")+
        theme_bw()+
        scale_color_manual(values = c("#E41A1C","#3300CC","#4DAF4A"))+ #+ ,"#EAE838"
        labs(title = tt,x = AA,y = "Immunoinfiltration scores")+
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test = wilcox.test)+
        theme(
          plot.title = element_text(color = 'black', size = 12, hjust = 0.5),
          axis.text.x = element_text(angle = 45,face = "bold",colour = "#1A1A1A",hjust=1,vjust=1),
          axis.text.y = element_text(face = "bold",colour = "#1A1A1A"),
          axis.title = element_text(size = 12,face = "bold", colour = "#1A1A1A"),
          legend.position="none"
        )
      #pdf(paste0(outdirP,"/",tt,"-TIMER2-TIMER箱形图.pdf"),height=5.5,width=6)
      pdf(paste0(outdirP2,"/",tt,"-图.pdf"),height=5.5,width=6)
      print(fig)
      dev.off()
      
    }
    
  } 
  
write.csv(all_immue_eatimate,
          file = paste0(outdir,"/all_immue_eatimate.csv"))
write.csv(allcancer_ssGSEA_Score,
          file = paste0(outdir,"/allcancer_ssGSEA_Score.csv"))
write.csv(gene_immCell_cor,
          file = paste0(outdir,"/gene_immCell_ssGSEA_cor.csv"))
write.csv(gene_immCell_cor2,
          file = paste0(outdir,"/gene_immCell_TIMER_cor.csv"))


