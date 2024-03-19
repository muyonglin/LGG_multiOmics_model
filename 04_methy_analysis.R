.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.2.0/lib/R/library","/thinker/storage/udata/sunhr/R/R4.2.0/lib/R/library/", .libPaths()))
.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.1.3/lib/R/library","/thinker/aid2/udata/all/muyl/others/software/condaEnvs/methyArray3/lib/R/library/", .libPaths()))
"/thinker/aid2/udata/all/muyl/others/software/condaEnvs/methyArray3/lib/R/library"

# .libPaths(c('/thinker/storage/biostaff/Biosoft/gao/R_lib161'
#             ,'/thinker/storage/udata/bing/RPackages'
#             ,'/thinker/storage/udata/bing/anaconda3/lib/R/library/'
#             ,'/thinker/storage/udata/muyl/software/R_lib161', .libPaths()))
.libPaths()
# .libPaths(c('/thinker/storage/biostaff/Biosoft/gao/R_lib161'
#             ,'/thinker/storage/udata/bing/RPackages'
#             ,'/thinker/storage/udata/bing/anaconda3/lib/R/library/'
#             ,'/thinker/storage/udata/muyl/software/R_lib161', .libPaths()))
getwd()
rm(list = ls())
#

# 安装ChAMP -----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChAMP")

library(ChAMP)


# 安装失败，直接引用穆博的环境 ----------------------------------------------------------
# conda activate /thinker/aid2/udata/all/muyl/others/software/condaEnvs/methyArray3

# testDir=system.file("extdata",package="ChAMPdata")
# myLoad <- champ.load(testDir,arraytype="450K")

# sample_info_use=read.csv("/thinker/aid2/udata/all/gao/others/tcga/lgg/file_use_info_table.xls",sep = "\t")
# sample_info_use
# rownames(sample_info_use)=sample_info_use$sample.submitter_id
# write.table(sample_info_use$sample.submitter_id, file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/samples.txt",quote = F, row.names = F, col.names = F)

library(data.table)
library(Hmisc)
#表型数据
pd.all <- read.delim("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA.LGG.sampleMap%2FLGG_clinicalMatrix", header = T, stringsAsFactors = F)
pd.all
# pd <- pd.all[,c("sampleID","sample_type","bcr_sample_barcode","bcr_patient_barcode","X_GENOMIC_ID_TCGA_LGG_hMethyl450","X_GENOMIC_ID_TCGA_LGG_hMethyl450_MethylMix")]
pd <- pd.all[,c("sampleID","sample_type","bcr_sample_barcode","bcr_patient_barcode")]
pd$sample_type <- ifelse(pd$sample_type=="Primary Tumor","Primary","Recurrent")
# pd$patient <- substr(pd$sampleID,1,12)
# pd <- pd[pd$patient %in% sample_info_use$cases.submitter_id,]
recurrent_sample <- pd[pd$sample_type=="Recurrent",]$bcr_patient_barcode
pd01 <- pd[pd$bcr_patient_barcode %in% recurrent_sample,]
rownames(pd01)=pd01$sampleID
pd=na.omit(pd01)
# 保存配对样本（此时有48对）
# write.table(pd$sampleID, file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/samples01.txt",quote = F, row.names = F, col.names = F)
# write.table(sample_info_use$sample.submitter_id, file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/samples.txt",quote = F, row.names = F, col.names = F)
write.table(pd01$sampleID, file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_sampleMap_samples01.txt",quote = F, row.names = F, col.names = F)

# 在/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao这个目录下：
## 筛选配对的LGG配对Recurrent_Primary样本
cat TCGA.LGG.sampleMapFHumanMethylation450 | cut -f 1-6 |head -5
# 1.获得LGG样本所在的列
# $ cols=($(sed '1!d;s/\t/\n/g' TCGA-LGG.methylation450.tsv | grep -nf samples.txt | sed 's/:.*$//'))
$ cols=($(sed '1!d;s/\t/\n/g' TCGA.LGG.sampleMapFHumanMethylation450 | grep -nf TCGA-LGG_Recurrent_Primary_sampleMap_samples01.txt | sed 's/:.*$//'))
# sed '1!d;s/\t/\n/g' TCGA.LGG.sampleMapFHumanMethylation450
# # 2.筛选出LGG样本Recurrent_Primary的甲基化信号矩阵
$ cut -f 1$(printf ",%s" "${cols[@]}") TCGA-LGG.methylation450.tsv > TCGA-LGG_Recurrent_Primary_methy450k.txt
$ cut -f 1$(printf ",%s" "${cols[@]}") TCGA.LGG.sampleMapFHumanMethylation450 > TCGA-LGG_Recurrent_Primary_sampleMap_methy450k.txt
## 读取具有配对Recurrent_Primary样本配对的LGG甲基化信号矩阵
getwd()
# a <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_methy450k.txt", data.table = F )
# a02 <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG.methylation450.tsv",data.table = F )
# a02[1:4,1:4]
# nrow(a02)
#
# 如果要读取所有的样本的话 ------------------------------------------------------------
# a01 <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA.LGG.sampleMapFHumanMethylation450", data.table = F )



# 如果要配对样本的话 ---------------------------------------------------------------

a01 <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_sampleMap_methy450k.txt", data.table = F )
# a
# a[1:4,1:4]
a01[1:4,1:4]
colnames(a01)
colnames(a02)
substr(colnames(a02), 1,15)
substr(colnames(a02), 1,15) %in% colnames(a01)
colnames()

sum(a01$sample %in% a02$`Composite Element REF`)
# head(a01$sample)

nrow(a01)
rownames(a01)=a01[,1]
a=a01[,-1]
a=na.omit(a)
# b=a[complete.cases(a),]
# beta=as.matrix(b)
library(impute)
# beta信号值矩阵里面不能有NA值
# beta=impute.knn(beta)
# betaData=beta$data
# betaData=betaData+0.00001
# sum(is.na(betaData))
# load("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/01.RData")

# #表型数据
# pd <- pd[colnames(betaData),]
# table(pd$patient) #此时有些normal样本没有配对的tumor
# NT.s <- pd$patient[duplicated(pd$patient)] #得到32对样本，和文章一致了
# pd <- pd[pd$patient %in% NT.s,]
# betaData <- betaData[,pd$sampleID]
# dim(betaData)
# 载入为ChAMP对象
library(ChAMP)
betaData <- as.matrix(a)
# 匹配临床数据和methy数据
pd <- pd[match(colnames(a),row.names(pd)),]
nrow(pd)
dim(betaData)
identical(colnames(betaData),row.names(pd))
colnames(betaData)[1:10]
row.names(pd)[1:10]
myLoad=champ.filter(beta = betaData ,pd = pd) #这一步已经自动完成了过滤
# 然后出现如下：
# > myLoad=champ.filter(beta = betaData ,pd = pd)
# [===========================]
# [<<<< ChAMP.FILTER START >>>>>]
# -----------------------------
#
#   In New version ChAMP, champ.filter() function has been set to do filtering on the result of champ.import(). You can use champ.import() + champ.filter() to do Data Loading, or set "method" parameter in champ.load() as "ChAMP" to get the same effect.
#
# This function is provided for user need to do filtering on some beta (or M) matrix, which contained most filtering system in champ.load except beadcount. User need to input beta matrix, pd file themselves. If you want to do filterintg on detP matrix and Bead Count, you also need to input a detected P matrix and Bead Count information.
#
# Note that if you want to filter more data matrix, say beta, M, intensity... please make sure they have exactly the same rownames and colnames.
#
#
# [ Section 1:  Check Input Start ]
# You have inputed beta for Analysis.
#
# pd file provided, checking if it's in accord with Data Matrix...
#     !!! Your pd file does not have Sample_Name column, we can not check your Sample_Name, please make sure the pd file is correct.
#
#   Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...
#     !!! Parameter detP is not found, filterDetP is reset FALSE now.
#
#   Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...
#     !!! Parameter beadcount is not found, filterBeads is reset FALSE now.
#
#   parameter autoimpute is TRUE. Checking if the conditions are fulfilled...
#     !!! ProbeCutoff is 0, which means you have no needs to do imputation. autoimpute has been reset FALSE.
#
#   Checking Finished :filterMultiHit,filterSNPs,filterNoCG,filterXY would be done on beta.
# [ Section 1: Check Input Done ]
#
#
# [ Section 2: Filtering Start >>
#
#     Filtering NoCG Start
#   Only Keep CpGs, removing 978 probes from the analysis.
#
#   Filtering SNPs Start
#   Using general 450K SNP list for filtering.
#   Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
#     Removing 22439 probes from the analysis.
#
#   Filtering MultiHit Start
#     Filtering probes that align to multiple locations as identified in Nordlund et al
#     Removing 4 probes from the analysis.
#
#   Filtering XY Start
#     Filtering probes located on X,Y chromosome, removing 7899 probes from the analysis.
#
#   Updating PD file
#     filterDetP parameter is FALSE, so no Sample Would be removed.
#
#   Fixing Outliers Start
#     Replacing all value smaller/equal to 0 with smallest positive value.
#     Replacing all value greater/equal to 1 with largest value below 1..
# [ Section 2: Filtering Done ]
#
#  All filterings are Done, now you have 338881 probes and 530 samples.
#
# [<<<<< ChAMP.FILTER END >>>>>>]
# [===========================]
# [You may want to process champ.QC() next.]
#


dim(myLoad$beta)
save(myLoad,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/step1-output.Rdata')


# 在R终端中：
rm(list = ls())
library(ChAMP)
load('/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/step1-output.Rdata')
# 数据归一化
getwd()
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,plotBMIQ=FALSE,resultsDir="./CHAMP_Normalization_paired/")
# 19：16分，16个核的时候，内存230G
# > myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,plotBMIQ=FALSE,resultsDir="./CHAMP_Normalization/")
# 结果成功了
# [===========================]
# [>>>>> ChAMP.NORM START <<<<<<]
# -----------------------------
#   champ.norm Results will be saved in ./CHAMP_Normalization/
#   [ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]
#
# << Normalizing data with BMIQ Method >>
#   Note that,BMIQ function may fail for bad quality samples (Samples did not even show beta distribution).
# 16 cores will be used to do parallel BMIQ computing.
# [>>>>> ChAMP.NORM END <<<<<<]
# [===========================]
# [You may want to process champ.SVD() next.]

# myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K") # 显示有NA值
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na)
table(num.na)
# myNorm <- myNorm[,which(num.na < 250000)]
save(myNorm,file = './Rdata/step2-champ_myNorm.Rdata')
# 所有样本的myNorm在：/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/
save(myNorm,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/step2-champ_myNorm.Rdata')
# 这里是配对样本的myNorm
save(myNorm,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/step2-champ_myNorm.Rdata')

# 主成分分析
library("FactoMineR")
library("factoextra")

dat <- t(myNorm)
group_list=pd$sample_type
pd <- myLoad$pd[colnames(myNorm),] #去掉异常样本
group_list=pd$sample_type
table(group_list)
pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/01_methy_PCA.pdf")
pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/01_methy_PCA.pdf")
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups")
dev.off()

# 热图
group_list=pd$sample_type
# 第一个参数是指要参与计算的矩阵；
# 第二个参数是指按行计算还是按列计算，1——表示按行计算，2——按列计算；
# 第三个参数是指具体的运算参数。
# sd函数计算数据列或者向量标准差（standard deviation）
#  tail() 函数用于获取向量、矩阵、表、 DataFrame 或函数的最后部分
cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)

rownames(ac)=colnames(myNorm)
pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/02_heatmap_top1000_sd.pdf')


# 相关关系矩阵热图
# 组内的样本的相似性应该是要高于组间的！

pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F,
                   filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/02_heatmap_top1000_corresponding.pdf')

# 去除异常样本
# 图中看出TCGA-CV-5971-01，TCGA-CV-6953-11，TCGA-CV-6955-11这3个样本有点异常
# myNorm <- myNorm[,!(colnames(myNorm) %in% c("TCGA-CV-5971-01","TCGA-CV-6953-11","TCGA-CV-6955-11"))]
pd <- pd[colnames(myNorm),]
# save(pd,myNorm,file = "./Rdata/filtered.Rdata")
save(pd,myNorm,file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/filtered.Rdata")


# 直接差异分析 ------------------------------------------------------------------
# 差异分析
group_list <- pd$sample_type
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
myDMP_01 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "none",adjPVal = 0.001)
myDMP_02 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "none",adjPVal = 1)
# myDMP_01 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "BH",adjPVal = 0.1)
# myDMP_01 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "BY",adjPVal = 0.05)
# myDMP_01 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "holm",adjPVal = 0.05)

# > myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
# [===========================]
# [<<<<< ChAMP.DMP START >>>>>]
# -----------------------------
#   !!! Important !!! New Modification has been made on champ.DMP():
#
#   (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.
#
# (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted "pheno" parameter is "numeric" type.
#
# --------------------------------
#
#   [ Section 1:  Check Input Pheno Start ]
#
# You pheno is character type.
# Your pheno information contains following groups. >>
#   <Primary>:516 samples.
# <Recurrent>:14 samples.
# [The power of statistics analysis on groups contain very few samples may not strong.]
# pheno contains only 2 phenotypes
# compare.group parameter is NULL, two pheno types will be added into Compare List.
# Primary_to_Recurrent compare group : Primary, Recurrent
#
# [ Section 1:  Check Input Pheno Done ]
#
#
# [ Section 2:  Find Differential Methylated CpGs Start ]
#
# -----------------------------
#   Start to Compare : Primary, Recurrent
# Contrast Matrix
# Contrasts
# Levels       pRecurrent-pPrimary
# pPrimary                    -1
# pRecurrent                   1
# You have found 11894 significant MVPs with a BH adjusted P-value below 0.05.
# Calculate DMP for Primary and Recurrent done.
#
# [ Section 2:  Find Numeric Vector Related CpGs Done ]
#
#
# [ Section 3:  Match Annotation Start ]
#
#
# [ Section 3:  Match Annotation Done ]
#
# [<<<<<< ChAMP.DMP END >>>>>>]
# [===========================]
# [You may want to process DMP.GUI() or champ.GSEA() next.]


head(myDMP[[1]])
save(myDMP,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/step3-output-myDMP.Rdata')
save(myDMP_01,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/step3-output-myDMP_001.Rdata')
save(myDMP_02,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/step3-output-myDMP_002.Rdata')

hmc <- myDMP[[1]][myDMP[[1]]$deltaBeta>0,]

hmc_01 <- myDMP_01[[1]][myDMP_01[[1]]$deltaBeta>0,]
hmc_02 <- myDMP_02[[1]][myDMP_02[[1]]$deltaBeta>0,]


write.table(hmc,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/hydroxymethylation_CpGs.xls",sep = "\t")
write.table(hmc_01,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/hydroxymethylation_CpGs.xls",sep = "\t")
write.table(hmc_02,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/hydroxymethylation_CpGs02.xls",sep = "\t")
#根据差异分析结果画火山图，因为作者的差异分析用的是paired t test，自己选了差异甲基化的cutoff。

# df_DMP <- myDMP$Primary_to_Recurrent[,1:5]
df_DMP <- myDMP$Recurrent_to_Primary
df_DMP_01 <- myDMP_01$Recurrent_to_Primary
df_DMP_02 <- myDMP_02$Recurrent_to_Primary

logFC_cutoff <- 0.45
df_DMP$change <- ifelse(df_DMP$adj.P.Val < 10^-15 & abs(df_DMP$logFC) > logFC_cutoff,
                        ifelse(df_DMP$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
table(df_DMP$change)

logFC_cutoff <- 0.05
df_DMP$change <- ifelse(df_DMP$adj.P.Val < 10^-7 & abs(df_DMP$logFC) > logFC_cutoff,
                          ifelse(df_DMP$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
write.table(df_DMP,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/df_DMP.xls",quote = F, row.names = F, col.names = F)

logFC_cutoff <- 0.05
df_DMP_01$change <- ifelse(df_DMP_01$adj.P.Val < 10^-6 & abs(df_DMP_01$logFC) > logFC_cutoff,
                        ifelse(df_DMP_01$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')

logFC_cutoff <- 0.0001

logFC_cutoff <- 0.001

# df_DMP_02$change <- ifelse(df_DMP_02$adj.P.Val < 10^-5 & abs(df_DMP_02$logFC) > logFC_cutoff,
#                            ifelse(df_DMP_02$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')

df_DMP_02$change <- ifelse(df_DMP_02$adj.P.Val < 10^-3 & abs(df_DMP_02$logFC) > logFC_cutoff,
                           ifelse(df_DMP_02$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
table(df_DMP_02$change)
write.table(df_DMP_02,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_002.xls",quote = F, row.names = T, col.names = T,sep = "\t")

# this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
#                     '\nThe number of up gene is ',nrow(df_DMP[df_DMP$change =='UP',]) ,
#                     '\nThe number of down gene is ',nrow(df_DMP[df_DMP$change =='DOWN',]))

this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(df_DMP_02[df_DMP_02$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(df_DMP_02[df_DMP_02$change =='DOWN',]))


library(ggplot2)
pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/volcanic_DMP.pdf")
g <- ggplot(data=df_DMP_02,
            aes(x=logFC, y=-log10(adj.P.Val),
                color=change)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g)
dev.off()

g <- ggplot(data=df_DMP_02,
            aes(x=logFC, y=-log10(adj.P.Val),
                color=change)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g)
# 上面的火山图暂时还画不了,由于甲基化的beta值不同于表达量，实际上用logFC也不是很对。曾老师课上推荐用deltabeta值替代logFC，就是甲基化信号值的差值。
# 参考：https://www.jianshu.com/p/8a2ad6d92790?u_atoken=f986bec4-9070-4ca7-a780-e02fbdef6ad5&u_asession=01NK0hARm3L8se2pg864I-q5NFWGXhIf4UKySmmEaKyyk_v60KgWEiglnrTUYbHwr_X0KNBwm7Lovlpxjd_P_q4JsKWYrT3W_NKPr8w6oU7K_UBmg6fyrsuo00ft-5CneIhUF3o-sVtq6Wun3JL3SJe2BkFo3NEHBv0PZUm6pbxQU&u_asig=05kTpl17A0roe5q06Aol0_P0GL_ZzE-SvDI99KQCrehsOn-olEqMNE1XdtrmQSA9VqF3ddxk-PleNEcQhraz2ywBdypBTmNUGOvXilDbYRBxLmn_lIGkDOn3c-i_OhmXuqde6AtnBSEteGtOzYL1A1ZlpHLNeVRbNjw1qKTeQ_qXr9JS7q8ZD7Xtz2Ly-b0kmuyAKRFSVJkkdwVUnyHAIJzbPtesF_23FoQf5wlVCYKsbvKFUWwEuGKX2ovy3mtYWgom7nzSzR1LP16f45fIKp-e3h9VXwMyh6PgyDIVSG1W_AK00aI3WFh7mc8vZyx2bF6n2obtm-JMAyZ9HwsFQtaDx0TBEOkQThS0-RPetNR3JGFMbCpYbPTweDIz8OwNuzmWspDxyAEEo4kbsryBKb9Q&u_aref=VewX2hHzy1Ei9m1a4kptgqUPMYg%3D
QC.GUI(beta=myNorm,arraytype="450K") # 显示有NA值
df_DMP_01 <- myDMP$Primary_to_Recurrent
df_DMP_02=df_DMP_01[df_DMP_01$gene!="",]
deltabeta_t <- 0.0001

P.Value_t <- 10^-10

df_DMP_02$change <- ifelse(df_DMP_02$adj.P.Val < P.Value_t & abs(df_DMP_02$deltaBeta) > deltabeta_t,
                        ifelse(df_DMP_02$deltaBeta > deltabeta_t ,'UP','DOWN'),'NOT')
table(df_DMP_02$change)
#>
#>   DOWN    NOT     UP
#>    345 108379    814
save(df_DMP_02,file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/step3.df_DMP_02.Rdata")

library(dplyr)
library(ggplot2)

this_tile <- paste0('Cutoff for deltabeta is ',round(deltabeta_t,3),
                    '\nThe number of up gene is ',nrow(df_DMP_02[df_DMP_02$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(df_DMP_02[df_DMP_02$change =='DOWN',]))

dat = rownames_to_column(df_DMP_02)
for_label <- dat %>% head(3)
pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/volcanic_DMP.pdf")
p <- ggplot(data=df_DMP_02,
            aes(x=deltaBeta,
                y= -log10(adj.P.Val)))+
  geom_point(alpha=0.4,size=3.5,aes(color=change))+ylab("-log10(Pvalue")+
  scale_color_manual(values = c("green","grey","red"))+
  geom_vline(xintercept = c(-deltabeta_t,deltabeta_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
dev.off()

library(ggplot2)

pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/volcanic_DMP002.pdf")
g <- ggplot(data=df_DMP_02,
            aes(x=deltaBeta, y=-log10(adj.P.Val),
                color=change)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("deltaBeta") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g)
dev.off()


# 热图
df_DMP_02 <- df_DMP
cg <- rownames(df_DMP_02[df_DMP_02$change != "NOT",])
plot_matrix <- myNorm[cg,]
annotation_col <- data.frame(Sample=pd$sample_type)
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors = list(Sample = c(Primary="#4DAF4A",Recurrent="#E41A1C"))
library(pheatmap)

pheatmap::pheatmap(plot_matrix,show_colnames = T,
                   annotation_col = annotation_col,
                   border_color = NA,
                   color = colorRampPalette(colors = c("white","navy"))(50),
                   annotation_colors = ann_colors,show_rownames = F,
                   filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/02_881_heatmap.pdf')

pheatmap::pheatmap(plot_matrix,show_colnames = T,
                   annotation_col = annotation_col,
                   border_color = NA,
                   color = colorRampPalette(colors = c("white","navy"))(50),
                   annotation_colors = ann_colors,show_rownames = F,
                   filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/02_881_heatmap.pdf')


# DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=group_list,resultsDir="/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/CHAMP_QCimages")

# 直接差异基因画KEGG
# df_DMP_02

# QC.GUI(beta=myNorm,arraytype="450K") # 显示有NA值
# QC =champ.QC(beta=myLoad$beta,pheno=myLoad$pd$sample_type)
# > QC =champ.QC(beta=myLoad$beta,pheno=myLoad$pd$sample_type)
# [===========================]
# [<<<<< ChAMP.QC START >>>>>>]
# -----------------------------
#   champ.QC Results will be saved in ./CHAMP_QCimages/
#   [QC plots will be proceed with 338881 probes and 530 samples.]
#
# << Prepare Data Over. >>
#   << plot mdsPlot Done. >>
#
#   << Plot densityPlot Done. >>
#
#   < Dendrogram Plot Feature Selection Method >: No Selection, directly use all CpGs to calculate distance matrix.
# << Plot dendrogram Done. >>
#
#   [<<<<<< ChAMP.QC END >>>>>>>]
# [===========================]
# [You may want to process champ.norm() next.]

# df_DMP <- myDMP$Primary_to_Recurrent[,1:5]
# logFC_cutoff <- 0.45
# df_DMP$change <- ifelse(df_DMP$adj.P.Val < 10^-15 & abs(df_DMP$logFC) > logFC_cutoff,
#                         ifelse(df_DMP$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
# table(df_DMP$change)
#
# this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
#                     '\nThe number of up gene is ',nrow(df_DMP[df_DMP$change =='UP',]) ,
#                     '\nThe number of down gene is ',nrow(df_DMP[df_DMP$change =='DOWN',]))
#
#

# 区域：
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$sample_type,method="Bumphunter")
write.table(myDMR$BumphunterDMR,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/myDMR.xls",sep = "\t")
DMR.GUI(DMR=myDMR)
# [===========================]
# [<<<<< ChAMP.DMR START >>>>>]
# -----------------------------
#   !!! important !!! We just upgrate champ.DMR() function, since now champ.DMP() could works on multiple phenotypes, but ProbeLasso can only works on one DMP result, so if your pheno parameter contains more than 2 phenotypes, and you want to use ProbeLasso function, you MUST specify compare.group=c("A","B"). Bumphunter and DMRcate should not be influenced.
#
# [ Section 1:  Check Input Pheno Start ]
#
# You pheno is character type.
# Your pheno information contains following groups. >>
#   <Primary Tumor>:516 samples.
# <Recurrent Tumor>:14 samples.
#
# [ Section 1:  Check Input Pheno Done ]
#
#
# [ Section 2:  Run DMR Algorithm Start ]
#
# Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
# << Find DMR with Bumphunter Method >>
#   3 cores will be used to do parallel Bumphunter computing.
# According to your data set, champ.DMR() detected 6373 clusters contains MORE THAN 7 probes within300 maxGap. These clusters will be used to find DMR.
#
# [bumphunterEngine] Parallelizing using 3 workers/cores (backend: doParallelMC, version: 1.0.17).
# [bumphunterEngine] Computing coefficients.
# [bumphunterEngine] Smoothing coefficients.
# Loading required package: rngtools
# [bumphunterEngine] Performing 250 bootstraps.
# [bumphunterEngine] Computing marginal bootstrap p-values.
# [bumphunterEngine] Smoothing bootstrap coefficients.
# [bumphunterEngine] cutoff: 0.895
# [bumphunterEngine] Finding regions.
# [bumphunterEngine] Found 319 bumps.
# [bumphunterEngine] Computing regions for each bootstrap.
# [bumphunterEngine] Estimating p-values and FWER.
# << Calculate DMR success. >>
#   Bumphunter detected 4 DMRs with P value <= 0.05.
#
# [ Section 2:  Run DMR Algorithm Done ]

# [<<<<<< ChAMP.DMR END >>>>>>]
# [===========================]
# [You may want to process DMR.GUI() or champ.GSEA() next.]
ulimit::memory_limit(60000)
# KEGG

library(org.Hs.eg.db)
library(clusterProfiler)

library(stringr)
library(BSgenome)
#library(RIdeogram)
library(circlize)
gene.v <- unique(as.character(df_DMP$gene))
write.table(gene.v,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/df_DMP_gene.list",quote = F, row.names = F, col.names = F)
gene.v <- unique(as.character(df_DMP_02[df_DMP_02$change!="NOT",]$gene))
write.table(gene.v,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_gene.list",quote = F, row.names = F, col.names = F)
gene.v <- read.csv("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_gene.list",header = F)
gene_L <- gene.v$V
# gene_L <- rownames(all.N.VS.T.deg.res.pcut.0.05)
# gene_L <- str_remove(gene_L, '\\.\\d+')
# oL <- bitr(gene_L, fromType='ENSEMBL', toType='UNIPROT', OrgDb = org.Hs.eg.db)
oL <- bitr(gene_L, fromType='SYMBOL', toType='UNIPROT', OrgDb = org.Hs.eg.db)
# KEGG富集的数据库
ken <- enrichKEGG(oL$UNIPROT,organism = 'hsa',keyType = 'uniprot', pAdjustMethod='none', pvalueCutoff=0.05)
# GO三个数据库
goCC <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='CC', pAdjustMethod='none', pvalueCutoff=0.05)
goBP <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='BP', pAdjustMethod='none', pvalueCutoff=0.05)
goMF <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='MF', pAdjustMethod='none', pvalueCutoff=0.05)


pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/samples.barplot.DMP.all.gene.functional.analysis.pdf')
p <- dotplot(ken,title='KEGG enrich of mutation gene')
print(p)
p <- dotplot(goCC,title='GO CC enrich of mutation gene')
print(p)
p <- dotplot(goBP,title='GO BP enrich of mutation gene')
print(p)
p <- dotplot(goMF,title='GO MF enrich of mutation gene')
print(p)

dev.off()
