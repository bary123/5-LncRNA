#### Step I 设置环境变量 ####
rm(list = ls()) ## 清空环境变量
#setwd('F:\\Other_analysis\\JLX_Nomogram\\Lesson3\\Lesson3_Results')
if (!requireNamespace("BiocManager", quietly = TRUE)) ## 安装R包
  install.packages("BiocManager") ## 安装R包Method I
# BiocManager::install('ResourceSelection') ## 安装R包Method II

library(tidyverse) ## 加载R包
library(rms)
library(Hmisc)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(foreign)
library(survMisc)
library(plyr)
library(timeROC)
library(survminer)
library(ResourceSelection)
options(stringsAsFactors = FALSE) #禁止chr转成factor

dir_name <-paste0("Results/07-cox/")#输出地址
ifelse(dir.exists(dir_name),FALSE,dir.create(dir_name,recursive = T))
## load data ##
load("./Results/01-getSig-lnRNA/UCEC_Step1_output.RData")
load("./Results/07-07-phenotype/UCEC_Clinical phenotypic11.Rdata")
df_omit_clear <-clinical_info
GEO_rt <- siglncExp
rm(gsExp,clinical_info,phe,siglncExp)
tmp <- intersect(df_omit_clear$submitter_id,colnames(GEO_rt))#取交集，一共是258个交集的
df <- subset(GEO_rt,select = tmp)#取tmp这些列
rm(GEO_rt,tmp)
df <- data.frame(t(df))#转换成行名是样本名
save(df,df_omit_clear,file=paste0(dir_name,"StepI-phenotypic.Rdata" ))

#####----加载生存时间
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)   
load("./Results/07-cox/StepI-phenotypic.Rdata")  
cla <- read.table(file="./01-data/TCGA-UCEC.survival.tsv",row.names = 1 )
cla <- cbind(rownames(cla),cla)
colnames(cla) <- cla[1,]
cla <- cla[-1,]
colnames(cla)[3] <- "PATIENT"
class(cla$OS.time)
cla <- cla[!duplicated(cla$PATIENT),]
tmp <- merge(df_omit_clear,cla,by.x="submitter_id",by.y="PATIENT")#
class(df_omit_clear$time)
tmp$time<-tmp$OS.time
tmp$time <- as.numeric(tmp$time)
tmp =tmp[,-c(13:15)] 
df_omit_clear=tmp
save(df,df_omit_clear,file="./Results/07-cox/StepI-phenotypic1.Rdata" )

# head(df)
# head(df_omit_clear)
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr) 
load("./Results/07-cox/StepI-phenotypic1.Rdata")
df$Samples <- rownames(df)
df_omit_clear$Samples <- df_omit_clear$submitter_id
df <- merge(df_omit_clear, df ,by = 'Samples')#通过Samples来合并
## 为了加快运算速度，我们取前100个样本，样本是随机取的，就是为了运算速度 ！！！
# df <- df[sample(nrow(df),100), ]
# head(df)

rm(df_omit_clear)

colnames(df)[4]<-"fustatus"
table(df$fustatus)
df$fustatus <- sub('Alive','0',df$fustatus)
df$fustatus <- sub('Dead','1',df$fustatus)
table(df$fustatus)
df =df[,-2]
set.seed(2023000)## 给一个随机种子，随机种子的功能：似每一次的随机结果相同，结果可以重复，因为下一步我们要随机7:3得到训练集和验证集，为了每次随机结果一致，我们需要给你一个随机种子.
tmp <- sort(sample(nrow(df),nrow(df)*.7),)##将数据集随机按照8:2分开,把不是训练集中数据，作为验证集。
df_train <-df[tmp,]
df_validation <-df[-tmp,]
rm(df,tmp)

save.image(file = './Results/07-cox/StepI-Chr2_tidyverse_rawdata.RData')


#### Step II 迭代特征选取 ####

## Chr I 随机森林（生存分析) ##

rm(list=ls())
options(stringsAsFactors = FALSE) #禁止chr转成factor
if(!"survivalsvm"%in%installed.packages()) 
  install.packages("survivalsvm")

if(!"randomForestSRC"%in%installed.packages()) 
  install.packages("randomForestSRC")
library(survival)
library(survivalsvm) ## 加载R包
library(randomForestSRC)
load("./Results/07-cox/StepI-Chr2_tidyverse_rawdata.RData") ## 加载数据
df_train$fustatus <- as.numeric(as.character(df_train$fustatus)) ## 调整数据类型,因子型转变为数值型得先变成字符型
class(df_train$fustatus)

df_train <- df_train[,-c(1:2,4:12)] ## 只留下生存状态和生存时间和基因表达量信息，进行特征选取

colnames(df_train)[ncol(df_train)]#看一下最后一列是不是基因名

out.rsf.1 <- rfsrc(Surv(time, fustatus) ~ . ,
                   data = df_train,
                   ntree = 100, ## tree的个数
                   nsplit = 1) ## 随机拆分数（非负整数），可提高运算速度，默认值为0

print.rfsrc(out.rsf.1) ## 输出结果信息

pred <- predict(out.rsf.1,
                df_train,
                OOB=TRUE,
                type="response")

gene.vs<-var.select(object = out.rsf.1,
                    method="vh") ## two options : 'vh':速度慢，准确度高；'md':速度快，准确度相对较低。

randomForestSRC_geneids <- gene.vs$topvars ## 提取分析结果中方差较大的变量

save(gene.vs,randomForestSRC_geneids,file = './Results/07-cox/StepII-rfsrc_genes.RData')


## Chr II Lasso回归筛选变量 ##

rm(list = ls()) ## 清空环境变量

if(!"glmnet"%in%installed.packages()) 
  install.packages("glmnet")

library("glmnet") ## 加载R包
library(tidyverse)
library(survival)
load("./Results/07-cox/StepI-Chr2_tidyverse_rawdata.RData") ## load 数据

df_train$fustatus <- as.numeric(as.character(df_train$fustatus)) ## 设置数据格式

expr_df <- df_train[,-c(1:13)] %>% as.matrix() ## 根据R包的要求，将数据需要筛选的部分提取转换为矩阵,lasso不需要生存时间和生存状态，所以多删除两列
class(expr_df)
# Lasso筛选变量动态过程图
cvfit = cv.glmnet(expr_df,
                  Surv(df_train$time,df_train$fustatus), 
                  nfold=10,#10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提
                  family = "cox"
) 
pdf(paste0('./Results/07-cox/Plot1_lasso_.pdf'),onefile = F)
plot(cvfit) ## 画图
dev.off()


fit <- glmnet(expr_df, Surv(df_train$time,df_train$fustatus), 
              family = "cox") 
pdf(paste0('./Results/07-cox/Plot2_lasso.pdf'),onefile = F)
plot(fit, label = TRUE)
dev.off()
#变量筛选
coef.min = coef(cvfit,s ="lambda.min")  ## lambda.min & lambda.1se 取一个
active.min = which(coef.min != 0 ) ## 找出那些回归系数没有被惩罚为0的

(lasso_geneids <- colnames(expr_df)[active.min]) ## 提取基因名称

save(lasso_geneids,file = './Results/07-cox/StepII-lasso_geneids.RData')



## Chr 3 Venn plot ##
rm(list = ls())
load("./Results/07-cox/StepII-rfsrc_genes.RData")
load("./Results/07-cox/StepII-lasso_geneids.RData")

Final_genes <- intersect(lasso_geneids,randomForestSRC_geneids)#取交集

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    LASSO = lasso_geneids, ## lasso筛选的结果
    RandomForestSRC = randomForestSRC_geneids ## 随机森林筛选的结果
  ),
  filename = paste0('./Results/07-cox/Plot1_lasso_RFSRC_sur_venn.jpg'), ## 输出图片的名字
  lwd = 3, ## 线的粗细
  fill = c("cornflowerblue", "darkorchid1"), ## 输出图像的颜色
  alpha = 0.6, ## 透明度
  label.col = "black", ## label的颜色
  cex = 1.5, ## 字符大小
  fontfamily = "serif", ## 字符格式
  fontface = "bold", ## 是否加粗
  cat.col = c("cornflowerblue", "darkorchid1"), ## 图例的颜色
  cat.cex = 2, ## 图例字体大小
  cat.fontfamily = "serif", ## 字体格式
  cat.fontface = "bold", ## 是否加粗
  margin = 0.05,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 20)
)

save(Final_genes,file = './Results/07-cox/StepII-Final_genes.RData')

#### Step III 单因素&多因素Cox分析 ####

rm(list = ls()) ## 清空环境变量
library(survival)
library(dplyr)
library(plyr)
load("./Results/07-cox/StepI-Chr2_tidyverse_rawdata.RData") ## 加载数据（基因表达）
load("./Results/07-cox/StepII-Final_genes.RData") ## 加载数据（临床病理学参数）
load("./Results/07-cox/StepII-lasso_geneids.RData") ## 加载数据（临床病理学参数）
load("./Results/07-cox/StepII-rfsrc_genes.RData") ## 加载数据（临床病理学参数）
df_train$fustatus <- as.numeric(as.character(df_train$fustatus))
df_train <- df_train[,-2]
df_validation <- df_validation[,-2]

BaSurv <- Surv(time = df_train$time, ## 生存时间
               event = df_train$fustatus) ## 生存状态

UniCox <- function(x){ ## 构建一个R function 便于后期调用
  FML <- as.formula(paste0('BaSurv~',x)) ## 构建生存分析公式
  GCox <- coxph(FML, data = df_train) ## Cox分析
  GSum <- summary(GCox) ## 输出结果
  HR <- round(GSum$coefficients[,2],2) ## 输出HR值
  PValue <- round(GSum$coefficients[,5],3) ## 输出P值
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") ## 输出HR的执行区间
  Unicox <- data.frame("characteristics" = x, ## 返回结果，并构建数据框
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       "P Value" = PValue)
  return(Unicox)
}

UniCox(colnames(df_train)[235]) ## 试验一下 function是否存在错误？

VarNames <- lasso_geneids
#VarNames <- randomForestSRC_geneids
#VarNames <- Final_genes ## 输出需要分析的变量名字
UniVar <- lapply(VarNames,UniCox) ## 批量做Cox分析
UniVar <- ldply(UniVar,data.frame) ## 将结果整理为一个数据框
(GetFactors <- UniVar$characteristics[which(UniVar$P.Value < 0.2)] %>% as.character()) ## 筛选其中P值<0.2的变量纳入多因素cox分析。

## 多因素分析 ##
fml <-  as.formula(paste0('BaSurv~',paste0(GetFactors,collapse = '+'))) ## 将上述单因素Cox分析中P值<0.2的变量纳入多因素cox分析
MultiCox <- coxph(fml, data = df_train) ## 多因素Cox回归
MultiSum <- summary(MultiCox) ## 提取分析结果

MultiName <- as.character(GetFactors)
MHR <- round(MultiSum$coefficients[,2],2)
MPV <- round(MultiSum$coefficients[,5],3)
MCIL <- round(MultiSum$conf.int[,3],2)
MCIU <- round(MultiSum$conf.int[,4],2)
MCI <- paste0(MCIL,'-',MCIU)
MulCox <- data.frame('characteristics' = MultiName, ## 构建多因素cox分析结果数据框
                     'Hazard Ratio' = MHR,
                     'CI95' = MCI,
                     'P Value' = MPV)

Final <- merge.data.frame(UniVar,MulCox, by = "characteristics", all = T, sodf_train = T) ## 将单因素Cox分析的结果和多因素Cox分析的结果合并成一个数据框
#View(Final)
    #write.csv(Final,'./Results/07-cox/CoxReg_new.csv')

(Final_GetFactors <- Final$characteristics[which(Final$P.Value.y < 0.05)] %>% as.character()) ## 确定最后有生存意义的因子
    #save(GetFactors,file = './Results/07-cox/Chr3_Univariate_Cox.RData')

fml <-  as.formula(paste0('BaSurv~',paste0(Final_GetFactors,collapse = '+'))) ## 通过多因素Cox回归构建一个评分
MultiCox <- coxph(fml, data = df_train)
MultiSum <- summary(MultiCox)

index.min <- MultiSum$coefficients[,1]
(index.min <- as.numeric(index.min))

signature <- as.matrix(subset(df_train,select = Final_GetFactors)) %*% as.matrix(exp(index.min))
df_train$signature <- signature

getwd()

save(df_train,df_validation,file = './Results/07-cox/StepIII-Signature_rawdata.RData')

## KM分析进行验证 ##
df_train$signature_by2 <- ifelse(df_train$signature > median(df_train$signature),'High-Score','Low-Score') ## 这里的命名要注意！！！

fit <- survfit(Surv(time, fustatus)~signature_by2,data = df_train)
library(survminer)
p <- ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c('High-Score','Low-Score'), 
                legend.title='Risk Score', palette = c("dodgerblue2","orchid2"),
                risk.table.height = 0.3)

pdf(paste0('./Results/07-cox/Plot2_Signature_KMplot.pdf'),onefile = F)
print(p)
dev.off()





#### Step VI Time-ROC  ####
## Method I ##
rm(list = ls())
library(timeROC)
load("./Results/07-cox//StepIII-Signature_rawdata.RData")
bc <- df_train
bc$signature_by2 <- ifelse(bc$signature >= median(bc$signature),'High-Risk','Low-Risk')
cox.timp2 <- coxph(Surv(bc$time, bc$fustatus == 1) ~ bc$signature_by2)
lpFit <- cox.timp2$linear.predictors
roc.fit <-timeROC(T = bc$time, # 结局时间
                  delta = bc$fustatus, # 生存结局
                  marker = lpFit, # 预测变量
                  cause = 1, # 阳性结局赋值，比如死亡，复发的赋值
                  weighting = "marginal", ## 'marginal' = Kaplan-Meier estimator of the censoring distribution 
                  ## or 'cox' = censoring by the Cox model
                  times = c(365*c(1,3,5)), ## 时间点，选取1年、3年和5年生存率
                  ROC = T,
                  iid = TRUE) 
library(tidyverse)

CI95 <- confint(roc.fit)[1] %>% unlist()

ROC_1 <- paste0(' : ',roc.fit$AUC[1] %>% as.numeric() %>% round(.,2),'(',CI95[1],'-',CI95[4],')')
ROC_3 <- paste0(' : ',roc.fit$AUC[2] %>% as.numeric() %>% round(.,2),'(',CI95[2],'-',CI95[5],')')
ROC_5 <- paste0(' : ',roc.fit$AUC[3] %>% as.numeric() %>% round(.,2),'(',CI95[3],'-',CI95[6],')')

filename <- paste0(paste0('./Results/07-cox/Plot4_OS_AUC_t_group_TypeI.pdf'))

pdf(filename,onefile = F)

plot(roc.fit,time=365,col = "blue",add =FALSE,title = F)#time 是时间点，col是线条颜色，add指是否添加在上一张图中
plot(roc.fit,time=365*3,col = "red",add = TRUE)
plot(roc.fit,time=365*5,col = "green",add = TRUE)

title(main = 'Time-dependent ROC curve')

legend("bottomright",
       legend = c(paste0('AUC at 1 year',ROC_1), paste0('AUC at 3 years',ROC_3),paste0('AUC at 5 years',ROC_5)),
       lty = c("solid","solid","solid"),
       col = c("blue","red","green"),
       bty = 'n')

dev.off()

## Method II ##
filename2 <- paste0(paste0('./Results/07-cox/Plot5_OS_AUC_t_group_TypeII.pdf'))
pdf(filename2,onefile = F)
plotAUCcurve(roc.fit, col = "red")

legend("topright",
       legend=c("Time-dependent ROC curve"),
       lty=c("solid"),
       col=c("red"),
       bty = 'n')

dev.off()






