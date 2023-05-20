rm(list = ls())
options(stringsAsFactors = F)
#setwd("C:\\Users\\86147\\Desktop\\SingleGene")
library(dplyr)

    #选择癌症类型加载并加载临床表型数据
    proj <- "UCEC"
    opts <- paste0("./Results/07-phenotype/")#输出位置
    ifelse(dir.exists(opts),FALSE,dir.create(opts,recursive = T))
    
    cl <- read.table("./01-data/TCGA-UCEC.GDC_phenotype.tsv",header = TRUE,sep = "\t",quote = "\"",)
    tmp=as.data.frame(colnames(cl))
    which(colnames(cl)=='submitter_id.samples')     #样本名
    which(colnames(cl)=='submitter_id')
    which(colnames(cl)=='vital_status.demographic')  #是否存活
    table(cl[,'vital_status.demographic'])
    colnames(cl)[82] <- "vital_status"
    which(colnames(cl)=='clinical_stage')    #分期
    table(cl[,'clinical_stage'])
    colnames(cl)[15] <- "stage"
    
    #生存时间
    which(colnames(cl)=='day_of_dcc_upload')
    table(cl[,'day_of_dcc_upload'])
    which(colnames(cl)=='day_of_form_completion')
    table(cl[,'day_of_form_completion'])

    which(colnames(cl)=='age_at_initial_pathologic_diagnosis')#年龄
    table(cl[,'age_at_initial_pathologic_diagnosis'])
    colnames(cl)[6] <- "age"
    which(colnames(cl)=='lost_follow_up')
    which(colnames(cl)=='days_to_death.demographic')
    table(cl[,'days_to_death.demographic'])
    which(colnames(cl)=='race.demographic')
    which(colnames(cl)=='gender.demographic')   #性别
    table(cl[,'gender.demographic'])
    colnames(cl)[80] <- "gender"
    
    colnames(cl[,c(1,13,82,15,31,6,80,17,18)])
    meta=as.data.frame(cl[,c(1,13,82,15,31,6,80,17,18)])
    colnames(meta)
    meta <- meta[!duplicated(meta$submitter_id),]  #去重
    meta[1:4,1:4]
    #去掉lost follow up的样本
    phe <- meta
    table(phe[,5])
    phe <- phe[!phe$lost_follow_up == 'YES',]
    
    #划分age_group
    phe$age=as.numeric(phe$age)
    phe$age_group=ifelse(phe$age>median(phe$age,na.rm=T),'older','younger')
    table(phe$age)
    table(phe$age_group)
    
    #stage
    table(phe$stage)
    phe <- phe[phe$stage != 'not reported',]
    phe$stage <- gsub('[1-2]','',phe$stage)
    phe$stage <- gsub('[A-C]','',phe$stage)
    table(phe$stage)

    clinical_info <- phe
    row.names(clinical_info) <- c(1:length(clinical_info[,1]))
    library(tableone)
    ###对需要观测的临床特质值进行重新编码
    colnames(clinical_info)
    class(clinical_info$age)
    clinical_info$age<-as.numeric(clinical_info$age)
    clinical_info$AGE<-factor(ifelse(clinical_info$age>60,'>60','<=60'),ordered = T)
    clinical_info$Gender<-factor(toupper(clinical_info$gender),ordered = T)
    table(clinical_info$stage)
    clinical_info$stage<-factor(clinical_info$stage,ordered = T) 
    
    save(phe,clinical_info,file=paste0(opts,proj,"_Clinical phenotypic11.Rdata"))
    
    