# 0. set environment ----------
ph <- '/thinker/aid2/udata/all/gao/others/tcga/lgg/'
setwd(ph)


.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.2.0/lib/R/library","/thinker/storage/udata/sunhr/R/R4.2.0/lib/R/library/", .libPaths()))
.libPaths()
library(data.table)
library(stringr)

library(vidger)
library(clusterProfiler)
library(TCGAbiolinks)

library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(stringr)

library(pheatmap)
library(org.Hs.eg.db)
library(BiocParallel)
##new
library(RCircos)
library(enrichplot)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(limma)
library(WGCNA)
library(ggpubr)

rm(list = ls())

# 1. prepare mutation survival p value --------
load("/thinker/aid2/udata/all/gao/others/tcga/lgg.RData")
# load("/thinker/aid2/udata/all/gao/others/tcga/lgg/RNA_star_count_14_sample.RData")
sampleTable$sampleName
lgg_maf


maf_df <- lgg_maf@data
maf_df$bcr_patient_barcode <- substr(maf_df$Tumor_Sample_Barcode,1,12)
maf_df$bcr_patient_barcode
# sampleTable$bcr_patient_barcode <- substr(sampleTable$sampleName,1,12)


# paired_sample_14_maf_df <- maf_df[maf_df$bcr_patient_barcode %in% sampleTable$bcr_patient_barcode,]
# unique(paired_sample_14$Tumor_Sample_Barcode)
# lgg_maf@variants.per.sample$Tumor_Sample_Barcode

gene_V <- unique(maf_df$Hugo_Symbol)

row_n_surv_cut <- length(unique(maf_df$Tumor_Sample_Barcode)) * 0.04
patient_num <- length(unique(maf_df$Tumor_Sample_Barcode))
mut_gene_coad_V <- c()
mut_freq_coad <- c()
high_mut_gene_coad_V <- c()
for (gene in gene_V){
  tmp <- maf_df[maf_df$Hugo_Symbol == gene, c('Tumor_Sample_Barcode')]
  mut_sample <- unique(tmp$Tumor_Sample_Barcode)
  mut_gene_coad_V <- c(mut_gene_coad_V,gene)
  mut_freq_coad <- c(mut_freq_coad, length(mut_sample)/patient_num)
  # print(length(tmp))
  # tmp_gene_mut <- c(tmp_gene_mut, nrow(tmp))
  if (length(mut_sample) >= row_n_surv_cut){
    # print(gene)
    # surv_gene_df[, gene] <- 0
    # surv_gene_df[surv_gene_df$sample_name %in% tmp$Tumor_Sample_Barcode, gene] <- 1
    high_mut_gene_coad_V <- c(high_mut_gene_coad_V,gene)
  }
}
high_mut_gene_coad_V

# pd.all <- read.delim("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA.LGG.sampleMap%2FLGG_clinicalMatrix", header = T, stringsAsFactors = F)
# pd.all

# maf_df_clinic_df <- pd[pd.all$X_GENOMIC_ID_TCGA_LGG_mutation_ucsc_maf_gene %in% unique(maf_df$Tumor_Sample_Barcode),]
# paired_sample_clinic_df <- pd.all[pd.all$bcr_patient_barcode %in% sampleTable$bcr_patient_barcode ,]
# paired_sample_clinic_df

# write.table(paired_sample_clinic_df,"/thinker/aid2/udata/all/gao/others/tcga/lgg/model/paired_sample_clinic_df.xls",row.names = F,sep = "\t",quote = F)
# surv_df_paired <- lgg_clinical_index[lgg_clinical_index$bcr_patient_barcode %in% sampleTable$bcr_patient_barcode,]
# surv_df <- lgg_clinical_index
surv_info <- lgg_clinical_index
surv_info
# write.table(surv_df,"/thinker/aid2/udata/all/gao/others/tcga/lgg/model/paired_sample_clinic_df_os_time.xls",row.names = F,sep = "\t",quote = F)

# Overall Survival (OS) The event call is derived from “vital status” parameter.
# The time_to_event is in days, equals to days_to_death if patient deceased;
# in the case of a patient is still living, the time variable is the maximum(days_to_last_known_alive, days_to_last_followup).
#This pair of clinical parameters are called _EVENT and _TIME_TO_EVENT on the cancer browser.

surv_df <- lgg_clinical_index[, c('bcr_patient_barcode', 'days_to_death', 'days_to_last_follow_up', 'vital_status')]
nrow(surv_df)
unique(surv_df$vital_status)
surv_df <- surv_df[surv_df$vital_status %in% c("Alive","Dead"),]
nrow(surv_df)
# 因为dataframe各列为factor，需转换成character'
# surv_df <- surv_info
surv_df$bcr_patient_barcode <- as.character(surv_df$bcr_patient_barcode)
surv_df$days_to_death <- as.character(surv_df$days_to_death)
surv_df$days_to_last_follow_up <- as.character(surv_df$days_to_last_follow_up)
surv_df$vital_status <- as.character(surv_df$vital_status)

# 因为死亡和存活患者的生存时间在表中两个不同的列，需要将两列时间合并
bcr_barcode <- c()
surv.status <- c()
surv.days <- c()
for (i in 1:nrow(surv_df)){
  bcr_barcode <- c(bcr_barcode, surv_df$bcr_patient_barcode[i])
  if (surv_df$vital_status[i] == 'Dead'){
    surv.days <- c(surv.days, surv_df$days_to_death[i])
    surv.status <- c(surv.status, 1)
  }else{
    if (surv_df$vital_status[i] == 'Alive'){
      surv.days <- c(surv.days, surv_df$days_to_last_follow_up[i])
      surv.status <- c(surv.status, 0)
    }else{
      surv.days <- c(surv.days, NA)
      surv.status <- c(surv.status, NA)
    }
  }
}

surv_df_2 <- data.frame(bcr_patient_barcode=bcr_barcode, time=surv.days, status=surv.status, stringsAsFactors = F)
surv_df_2$time <- as.numeric(surv_df_2$time)
# 有的项目会出现一个患者，两个生存信息的情况，如果有这种情况可以使用下一步的代码，提取较长的follow up time
# 可以使用下述代码判断一下是否出现患者重复的情况
# 如果有重复出现，则会有FALSE
surv_df_2$bcr_patient_barcode %in% unique(surv_df_2$bcr_patient_barcode)
# 或者
length(surv_df_2$bcr_patient_barcode) == length(unique(surv_df_2$bcr_patient_barcode))
surv_df_2$time[is.na(surv_df_2$time)]
surv_df_2 <- surv_df_2[complete.cases(surv_df_2),]
surv_df_2$time[is.na(surv_df_2$time)]
nrow(surv_df_2)

# select survival longest follow up time
barcode_V <- c()
surv.days <- c()
surv.status <- c()
for (i in unique(surv_df_2$bcr_patient_barcode)){
  tmpdf <- surv_df_2[surv_df_2$bcr_patient_barcode == i, ]
  tmpdf <- unique(tmpdf)
  if (nrow(tmpdf) == 1){
    barcode_V <- c(barcode_V, i)
    surv.days <- c(surv.days, tmpdf$time)
    surv.status <- c(surv.status, tmpdf$status)
  }else{
    max_d <- max(tmpdf$time)
    tmpdf_2 <- tmpdf[tmpdf$time == max_d, ]
    barcode_V <- c(barcode_V, i)
    surv.days <- c(surv.days, unique(tmpdf_2$time))
    surv.status <- c(surv.status, unique(tmpdf_2$status))
  }
}

surv_df_3 <- data.frame(patient=barcode_V, time=surv.days, status=surv.status, stringsAsFactors = F)
rownames(surv_df_3) <- surv_df_3$bcr_patient_barcode
surv_df_3 <- surv_df_3[complete.cases(surv_df_3),]
nrow(surv_df_3)



# surv_df_paired$OS_time_days <- surv_df_paired$days_to_death
# surv_df_paired$OS_time_days[is.na(surv_df_paired$OS_time_days)] <- 0
# surv_df_paired$days_to_last_follow_up_01 <- surv_df_paired$days_to_last_follow_up
# surv_df_paired$days_to_last_follow_up_01[is.na(surv_df_paired$days_to_last_follow_up_01)] <- 0
# surv_df_paired$days_to_last_follow_up_01
# surv_df_paired$OS_time_days[surv_df_paired$OS_time_days==0] <- surv_df_paired$days_to_last_follow_up_01
# surv_df_paired$OS_time_days <- ifelse(surv_df_paired$vital_status=="Alive",surv_df_paired$days_to_last_follow_up_01,surv_df_paired$OS_time_days)
# surv_df_paired[c("OS_time_days","vital_status","days_to_death","days_to_last_follow_up")]
# surv_info_df <- as.data.frame(surv_info)
# surv_df$OS_time_days <- surv_df$days_to_death
# surv_df$OS_time_days[is.na(surv_df$OS_time_days)] <- 0
# surv_df$days_to_last_follow_up_01 <- surv_df$days_to_last_follow_up
# surv_df$days_to_last_follow_up_01[is.na(surv_df$days_to_last_follow_up_01)|] <- 0
# surv_df$days_to_last_follow_up_01
# surv_df$OS_time_days[surv_df$OS_time_days==0] <- surv_df$days_to_last_follow_up_01
# surv_df$OS_time_days <- ifelse(surv_df$vital_status=="Alive",surv_df$days_to_last_follow_up_01,surv_df$OS_time_days)
# surv_df[c("OS_time_days","vital_status","days_to_death","days_to_last_follow_up")][255,]
#
#
# surv_df$os <- surv_df$OS_time_days/30
# surv_df$os
# surv_df$status <-  ifelse(surv_df$vital_status=="Alive",0,1)
#
# for (gene in high_mut_gene_coad_V){
#   tmp <- paired_sample_14_maf_df[paired_sample_14_maf_df$Hugo_Symbol == gene, 'bcr_patient_barcode']
#   surv_df_paired[, gene] <- 0
#   surv_df_paired[surv_df_paired$bcr_patient_barcode %in% tmp$bcr_patient_barcode, gene] <- 1
# }
# surv_df_paired$survival <- surv_df_paired$status
# surv_df_paired$OS <- surv_df_paired$os
# # survival
# os_df <-  surv_df_paired[!is.na( surv_df_paired$os),]

maf.df <- maf_df
head(surv_df_3)
length(unique(maf.df$bcr_patient_barcode))
nrow(surv_df_3)
unique(maf.df$bcr_patient_barcode) %in% surv_df_3$patient
unique(maf.df$bcr_patient_barcode)[!(unique(maf.df$bcr_patient_barcode) %in% surv_df_3$patient)]
tmp <- unique(maf.df$bcr_patient_barcode)[!(unique(maf.df$bcr_patient_barcode) %in% surv_df_3$patient)]
tmp %in% surv_info$bcr_patient_barcode
surv_info[surv_info$bcr_patient_barcode %in% tmp, ]

unique(surv_df_3$patient)[!unique(surv_df_3$patient) %in% unique(maf.df$bcr_patient_barcode)]

# surv_df_4 <- surv
surv_df_4 <- surv_df_3[surv_df_3$patient %in% unique(maf.df$bcr_patient_barcode) ,]

# surv_df_3 <- surv_df_3[surv_df_3$Tumor_Sample_Barcode %in% unique(maf.df$Tumor_Sample_Barcode),]
# surv_df_3 <- surv_df_3[unique(maf.df$bcr_patient_barcode)%in% surv_df_3$patient ,]
nrow(surv_df_4)


surv_info_df <- as.data.frame(surv_df_4)
colnames(surv_info_df) <- c('Tumor_Sample_Barcode', 'OS', 'survival')
# surv_info_df <- surv_info_df[, c('PFS', 'OS', 'reoccurence', 'survival', 'Tumor_Sample_Barcode')]
for (gene in high_mut_gene_coad_V){
  tmp <- maf_df[maf_df$Hugo_Symbol == gene, 'bcr_patient_barcode']
  surv_info_df[, gene] <- 0
  surv_info_df[surv_info_df$Tumor_Sample_Barcode %in% tmp$bcr_patient_barcode, gene] <- 1
}

mut_gene_surv_df <- surv_info_df
write.table(mut_gene_surv_df,"./ana_20221121/mut_gene_surv_df.xls",row.names = F,sep = "\t",quote = F)
# 2. prepare RNA deg survival p value ---------
getwd()
setwd("/thinker/aid2/udata/all/gao/others/tcga/lgg")
# 用所有的样本的star_count
htseq_l <- list.files('./star_count_bak', pattern = '*star_count', full.names = T)
htseq_l
sample.info.df <- data.frame(htseq=htseq_l, stringsAsFactors = F)

sample.info.df[, 'condition'] <- basename(sample.info.df$htseq)
sample.info.df$condition <- str_remove(sample.info.df$condition, '_star_count')
sample.info.df$condition <- str_remove(sample.info.df$condition, '.+\\d_')
sample.info.df$condition
# sample.info.df$condition <- str_remove(sample.info.df$condition, 'P\\d+')
# sample.info.df$condition <- str_remove(sample.info.df$condition, '\\.htseqCount')
# sample.info.df$condition <- str_replace(sample.info.df$condition, 'L', 'N')
sample.info.df[, 'sample_name'] <- basename(sample.info.df$htseq)

sample.info.df$sample_name <- str_remove(sample.info.df$sample_name, '_.+')
sample.info.df
sampleTable <- data.frame(sampleName = sample.info.df$sample_name
                          , fileName = sample.info.df$htseq
                          , condition = sample.info.df$condition)
table(sampleTable$condition)
ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable
                                       , directory = './'
                                       , design = ~condition)
register(MulticoreParam(20))
keep <- rowSums(counts(ddsHTseq)) >= min(1 * ncol(ddsHTseq), 10)
dds <- ddsHTseq[keep, ]
dds <- DESeq(dds, parallel = T)
dds

dds_all_sample <- dds

# # 1.3.1 DEGs between Primary_Tumor and Recurrent_Tumor
res <- results(dds_all_sample, contrast=c('condition', 'Primary_Tumor', 'Recurrent_Tumor'), parallel=T)
resOrdered <- res[order(res$padj),]
resCut <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)##新增
#resCut <- subset(resOrdered, padj < 0.001)
all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05 <- resCut##新增
#all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.001 <- resCut
N.vs.T.table.for.volcano <- resOrdered
all_sample_deg <- all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05
#
# # 1.3.2 deg expression heatmap
# vsd <- vst(dds)
# countMat <- assay(vsd)
#
# degMat <- countMat[rownames(all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05),]##新增
# #degMat <- countMat[rownames(all.deg.control.VS.T1.pcut.0.001),]
#
# tmp_table <- sampleTable
# rownames(tmp_table) <- sampleTable$sampleName
# tmp_table <- tmp_table[colnames(degMat),]
# annot_df <- data.frame(condition=tmp_table$condition)
# rownames(annot_df) <- colnames(degMat)
# ##show top 100 gene ordered by p.adj values.
# degMat_top100<-degMat[1:min(nrow(degMat),100),]
# # rownames(annot_df) <- rownames(inF2)
# # pheatmap(degMat, show_rownames = F, border_color = NA, show_colnames = F, annotation_col = annot_df)
# # pheatmap(degMat_top100, show_rownames = F, border_color = NA, show_colnames = F, annotation_col = annot_df
# #          , main = 'all DEGs between Recurrent_Tumor and Primary_Tumor'
# #          , filename = './plots/rnaseq/1.degs/deg.all.unique.heatmap.normalized.count.mat.padj.cut.0.05.pdf')
#
# ## 1.3.3 deg function enrichment
# gene_V <- c(rownames(all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05))
# #gene_V <- c(rownames(all.deg.control.VS.T1.pcut.0.001))
# gene_V <- unique(gene_V)
# gene_V <- as.character(gene_V)
# gene_V <- str_remove(gene_V, '\\.\\d+')
# oL <- bitr(gene_V, fromType=toupper('ENSEMBL'), toType='ENTREZID', OrgDb = org.Hs.eg.db)
# #oL <- bitr(gene_V, fromType=toupper('ENSEMBL'), toType='UNIPROT', OrgDb = org.Hs.eg.db)
# inV <- oL$ENTREZID
# #inV <- oL$UNIPROT
# # 1.3.4 deg output
# # prepare ensg data frame
# # gtf <- '/mnt/phoenix/bio-web/pipeline/rnaseq/reference/reference/gencode.v34.annotation.gtf'
# gtf <- '/thinker/storage/udata/muyl/genomes/human/GRCH38/gencode.v34.annotation.gtf'
# inT <- fread(gtf, skip = 5)
#
# inT <- as.data.frame(inT)
# tmp <- inT[inT$V3 == 'gene',]
#
# tmp[,'gene_id'] <- str_extract(tmp$V9, 'ENSG\\d+\\.\\d+')
# # tmp[,'transcript_name'] <- str_extract(tmp$V9, 'ENST\\d+\\.\\d+')
# tmp[, 'gene_name'] <- str_extract(tmp$V9, '(gene_name ").*?("; level)')
# tmp$gene_name <- str_remove(tmp$gene_name, 'gene_name "')
# tmp$gene_name <- str_remove(tmp$gene_name, '"; level')
# # tmp[, 'exon_length'] <- tmp$V5 - tmp$V4 + 1
#
# tmp[,"chr"]<-tmp$V1
# tmp[,'start']<-tmp$V4
# tmp[,'end']<-tmp$V5
#
# # tmp[,'gene_type'] <- str_extract(tmp$V9, 'transcript_type ".*?"; transcript_status')
# tmp[,'gene_type'] <- str_extract(tmp$V9, 'gene_type ".*?"; gene_name')
# tmp$gene_type <- str_remove(tmp$gene_type, 'gene_type "')
# tmp$gene_type <- str_remove(tmp$gene_type, '"; gene_name')
#
# tmp[,'gene_id_trim'] <- str_remove(tmp$gene_id, '\\.\\d+')
#
# #gene_ensg_df <- tmp[,c('gene_name', 'gene_id', 'gene_type', 'gene_id_trim')]
# gene_ensg_df <- tmp[,c('chr','gene_name','start','end', 'gene_id', 'gene_type', 'gene_id_trim')]
# gene_ensg_df <- unique(gene_ensg_df)
#
#
# tmp <- volcano
# tmp[, 'threshold'] <- threshold
# tmp1 <- tmp[tmp$threshold == TRUE, ]
# print(str_c('N vs T1 deg number is: ', as.character(nrow(tmp1))))
# print(str_c('N vs T1 up regulated deg number: ', as.character(nrow(tmp1[tmp1$log2FoldChange > 0, ]))))
# print(str_c('N vs T1 down regulated deg number: ', as.character(nrow(tmp1[tmp1$log2FoldChange < 0, ]))))
#
# tmp_gene_ensg_df <- gene_ensg_df[gene_ensg_df$gene_id %in% rownames(tmp1),]
# # rownames(tmp_gene_ensg_df) <- tmp_gene_ensg_df$gene_id
#
# tmp1[, 'gene_id'] <- rownames(tmp1)
# tmp2 <- left_join(tmp_gene_ensg_df,tmp1)
# # 18个的RNA
# htseq_l <- list.files('./star_count', pattern = '*star_count', full.names = T)
# htseq_l
# sample.info.df <- data.frame(htseq=htseq_l, stringsAsFactors = F)
#
# sample.info.df[, 'condition'] <- basename(sample.info.df$htseq)
# sample.info.df$condition <- str_remove(sample.info.df$condition, '_star_count')
# sample.info.df$condition <- str_remove(sample.info.df$condition, '.+\\d_')
# sample.info.df$condition
# # sample.info.df$condition <- str_remove(sample.info.df$condition, 'P\\d+')
# # sample.info.df$condition <- str_remove(sample.info.df$condition, '\\.htseqCount')
# # sample.info.df$condition <- str_replace(sample.info.df$condition, 'L', 'N')
# sample.info.df[, 'sample_name'] <- basename(sample.info.df$htseq)
#
# sample.info.df$sample_name <- str_remove(sample.info.df$sample_name, '_.+')
# sample.info.df
# sampleTable <- data.frame(sampleName = sample.info.df$sample_name
#                           , fileName = sample.info.df$htseq
#                           , condition = sample.info.df$condition)
# table(sampleTable$condition)
# ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable
#                                        , directory = './'
#                                        , design = ~condition)
# register(MulticoreParam(20))
# keep <- rowSums(counts(ddsHTseq)) >= min(1 * ncol(ddsHTseq), 10)
# dds <- ddsHTseq[keep, ]
# dds <- DESeq(dds, parallel = T)
# dds
#
# dds_14_pair <- dds

# res <- results(dds_14_pair, contrast=c('condition', 'Primary_Tumor', 'Recurrent_Tumor'), parallel=T)
# resOrdered <- res[order(res$padj),]
# resCut <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)##新增
# #resCut <- subset(resOrdered, padj < 0.001)
# all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05 <- resCut##新增
# #all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.001 <- resCut
# N.vs.T.table.for.volcano <- resOrdered

# 2.1 prepare foldchange data frame --------
gtf <- '/thinker/storage/udata/muyl/genomes/human/GRCH38/gencode.v34.annotation.gtf'
inT <- fread(gtf, skip = 5)

inT <- as.data.frame(inT)
tmp <- inT[inT$V3 == 'gene',]

tmp[,'gene_id'] <- str_extract(tmp$V9, 'ENSG\\d+\\.\\d+')
# tmp[,'transcript_name'] <- str_extract(tmp$V9, 'ENST\\d+\\.\\d+')
tmp[, 'gene_name'] <- str_extract(tmp$V9, '(gene_name ").*?("; level)')
tmp$gene_name <- str_remove(tmp$gene_name, 'gene_name "')
tmp$gene_name <- str_remove(tmp$gene_name, '"; level')
# tmp[, 'exon_length'] <- tmp$V5 - tmp$V4 + 1

# tmp[,'gene_type'] <- str_extract(tmp$V9, 'transcript_type ".*?"; transcript_status')
tmp[,'gene_type'] <- str_extract(tmp$V9, 'gene_type ".*?"; gene_name')
tmp$gene_type <- str_remove(tmp$gene_type, 'gene_type "')
tmp$gene_type <- str_remove(tmp$gene_type, '"; gene_name')

tmp[,'gene_id_trim'] <- str_remove(tmp$gene_id, '\\.\\d+')

gene_ensg_df <- tmp[,c('gene_name', 'gene_id', 'gene_type', 'gene_id_trim')]
gene_ensg_df <- unique(gene_ensg_df)


# deg gene name
tmp_deg <- all_sample_deg
tmp_deg <- as.data.frame(tmp_deg)

deg_ensg <- rownames(tmp_deg)
deg_ensg <- str_remove(deg_ensg, '\\.\\d+')

deg_ensg %in% gene_ensg_df$gene_id_trim
sum(deg_ensg %in% gene_ensg_df$gene_id_trim)
nrow(tmp_deg)

deg_symbol <- gene_ensg_df[gene_ensg_df$gene_id_trim %in% deg_ensg, 'gene_name']
deg_symbol <- unique(deg_symbol)
length(deg_symbol)


# count matrix
# in 16.f20603  # 9.1.2.2 model deg
vsd <- vst(dds_all_sample)
countMat <- assay(vsd)

degM <- countMat[rownames(tmp_deg),]
degM[1:5,1:5]
deg_mat <- degM

head(sampleTable)
sampleTable$condition
tumor_coln <- sampleTable[sampleTable$condition == 'Primary_Tumor', 'sampleName']
control_coln <- sampleTable[sampleTable$condition == "Recurrent_Tumor", 'sampleName']


# control_coln <- coln[grepl('_control', coln)]
control_mat <- deg_mat[,control_coln]
tmp <- apply(control_mat, 1, mean)
head(tmp)

deg_mean_control <- data.frame(gene_id=names(tmp), mean_vst=tmp)

# tumor_coln <- coln[grepl('_T1', coln)]
tumor_mat <- deg_mat[, tumor_coln]

# check gene order
sum(rownames(deg_mean_control) == rownames(tumor_mat))

fc_df <- data.frame(gene_id=rownames(tumor_mat))

for(tc in tumor_coln){
  fc_df[, tc] <- tumor_mat[, tc] - deg_mean_control$mean_vst
}
# fpkmM <- as.data.frame(fpkmM)
# coln <- colnames(fpkmM)
# coln <- str_remove(coln, 'cufflinksOut.aws.')
# coln <- str_remove(coln, 'cufflinksOut.d2.')
# coln <- str_remove(coln, 'cufflinksOut.67.')
# colnames(fpkmM) <- coln
# rownames(fpkmM) <- fpkmM$V1
# fpkmM1 <- fpkmM[,c(2:ncol(fpkmM))]
# deg_fpkm_mat <- fpkmM1[deg_symbol,]
# deg_fpkm_mat <- log2(deg_fpkm_mat + 1 )
fc_df[1:5,1:5]

deg_fc_df <- t(fc_df)
deg_fc_df <- as.data.frame(deg_fc_df)
colnames(deg_fc_df) <- as.character(fc_df$gene_id)
deg_fc_df <- deg_fc_df[2:nrow(deg_fc_df),]
deg_fc_df[1:5,1:5]

deg_fc_df[,'patient'] <- rownames(deg_fc_df)
# deg_fc_df[,'patient'] <- str_remove(deg_fc_df$patient, '_T1')
deg_fc_df[,'patient'] <- substr(deg_fc_df$patient,1,12)
deg_fc_df[1:5,c(1:5, ncol(deg_fc_df))]


rna_surv_info <- surv_df_3

# surv_info <- fread('./ana_20211105/update_table_20211117.csv')
tmp_surv_df <- as.data.frame(rna_surv_info)
head(tmp_surv_df)
colnames(tmp_surv_df) <- c('patient', 'OS', 'survival')
# tmp_surv_df[, 'patient'] <- str_remove(tmp_surv_df$Tumor_Sample_Barcode, '_T1')
# tmp_surv_df <- tmp_surv_df[,  c('PFS', 'OS', 'reoccurence', 'survival', 'patient')]

tmp_surv_df$patient %in% deg_fc_df$patient
sum(tmp_surv_df$patient %in% deg_fc_df$patient)
length(tmp_surv_df$patient %in% deg_fc_df$patient)
deg_fc_df$patient %in% tmp_surv_df$patient
sum(deg_fc_df$patient %in% tmp_surv_df$patient)

tmp_patient <- intersect(deg_fc_df$patient, tmp_surv_df$patient)

tmp_deg_fc_df <- deg_fc_df[deg_fc_df$patient %in% tmp_patient,]
tmp_surv_df <- tmp_surv_df[tmp_surv_df$patient %in% tmp_patient,]

tmp_deg_fc_surv_df <- left_join(tmp_surv_df, deg_fc_df, by='patient')
tmp_deg_fc_surv_df[1:5,1:10]
write.table(tmp_dmg_delta_beta_surv_df,"./ana_20221121/tmp_deg_fc_surv_df.xls",row.names = F,sep = "\t",quote = F)


# 3. 所有样本chAMP DMG survival p value -------
# conda activate /thinker/aid2/udata/all/muyl/others/software/condaEnvs/methyArray3
.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.1.3/lib/R/library","/thinker/aid2/udata/all/muyl/others/software/condaEnvs/methyArray3/lib/R/library/", .libPaths()))
library(data.table)
library(Hmisc)
library(ChAMP)
# write.table(all_sample_name,file="/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_all_samples01.txt",quote = F, row.names = F, col.names = F)
# 1.获得LGG样本所在的列
#表型数据
pd.all <- read.delim("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA.LGG.sampleMap%2FLGG_clinicalMatrix", header = T, stringsAsFactors = F)
pd.all
pd <- pd.all[,c("sampleID","sample_type","bcr_sample_barcode","bcr_patient_barcode")]
pd$sample_type <- ifelse(pd$sample_type=="Primary Tumor","Primary","Recurrent")
head(pd)
pd
# recurrent_sample <- pd[pd$sample_type=="Recurrent",]$bcr_patient_barcode
# pd01 <- pd[pd$bcr_patient_barcode %in% recurrent_sample,]
# rownames(pd01)=pd01$sampleID
rownames(pd)=pd$sampleID
pd=na.omit(pd)
# pd
nrow(pd)
# 保存配对样本
# write.table(pd01$sampleID, file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_sampleMap_samples01.txt",quote = F, row.names = F, col.names = F)
#$ cols=($(sed '1!d;s/\t/\n/g' TCGA.LGG.sampleMapFHumanMethylation450 | grep -nf TCGA-LGG_Recurrent_Primary_sampleMap_samples01.txt | sed 's/:.*$//'))

# sed '1!d;s/\t/\n/g' TCGA.LGG.sampleMapFHumanMethylation450
# $ cut -f 1$(printf ",%s" "${cols[@]}") TCGA.LGG.sampleMapFHumanMethylation450 > TCGA-LGG_Recurrent_Primary_sampleMap_methy450k.txt

# a01 <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_sampleMap_methy450k.txt", data.table = F )
a01 <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA.LGG.sampleMapFHumanMethylation450", data.table = F )
a01[1:4,1:4]
colnames(a01)
# colnames(a02)
# substr(colnames(a02), 1,15)
# substr(colnames(a02), 1,15) %in% colnames(a01)
# colnames()

# sum(a01$sample %in% a02$`Composite Element REF`)
# head(a01$sample)
nrow(a01)
ncol(a01)
rownames(a01)=a01[,1]
a=a01[,-1]
a=na.omit(a)
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
dim(myLoad$beta)
# save(myLoad,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/step1-output.Rdata')
# library(ChAMP)
# load('/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/step1-output.Rdata')
# 数据归一化
# getwd()
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,plotBMIQ=FALSE,resultsDir="./ana_20221116/CHAMP_Normalization_all/")
all.sample.myNorm <- myNorm
# saveRDS(myNorm, './ana_20221116/all.sample.myNorm.champ.rds')

group_list <- pd$sample_type
myDMP01 <- champ.DMP(beta = myNorm,pheno=group_list, adjust.method = "none",adjPVal = 0.05)
myDMP02 <- champ.DMP(beta = myNorm,pheno=group_list)
# 8000 dmp in myDMP02
# saveRDS(myDMP02, './ana_20221116/all.sample.myDMP02.champ.rds')

head(myDMP02$Primary_to_Recurrent)
all.sample.dmp.gene <- unique(as.character(myDMP02$Primary_to_Recurrent$gene))
all.sample.dmp.gene[is.na(all.sample.dmp.gene)]
all.sample.dmp.gene[all.sample.dmp.gene == '']
all.sample.dmp.gene <- all.sample.dmp.gene[all.sample.dmp.gene != '']
write.csv(all.sample.dmp.gene,"/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/all.sample.dmp.gene.list",row.names = F )

tmp.dmp <- rownames(myDMP02$Primary_to_Recurrent)
# calculate delta beta as foldchange
# https://rdrr.io/bioc/ChAMP/src/R/champ.DMP.R
all.sample.dmp.mynorm <- all.sample.myNorm[tmp.dmp,]
ncol(all.sample.dmp.mynorm)

tmp_dmp_df <- data.frame(probe=rownames(myDMP02$Primary_to_Recurrent), gene=myDMP02$Primary_to_Recurrent$gene)


# gene norm
gene_prob_list <- list()
for (gene in all.sample.dmp.gene){
  tmp <- tmp_dmp_df[tmp_dmp_df$gene == gene,]
  probe <- tmp$probe
  probe <- probe[probe!='']
  gene_prob_list[[gene]] <- probe
}


gene_norm_df <- data.frame(patient=colnames(myNorm))
for(i in names(gene_prob_list)){
  if (length(gene_prob_list[[i]]) > 1){
    tmp <- myNorm[gene_prob_list[[i]],]
    gene_norm_df[,i] <- colMeans(tmp)
  }else{
    gene_norm_df[,i] <- myNorm[gene_prob_list[[i]],]
  }
}

gene_norm_df[1:5,1:5]
gene_norm_df[(nrow(gene_norm_df) - 5): nrow(gene_norm_df),1:5]
gene_norm_df_1 <- gene_norm_df
rownames(gene_norm_df_1) <- gene_norm_df_1$patient
gene_norm_df_1 <- gene_norm_df_1[,2:ncol(gene_norm_df_1)]
gene_norm_df_1 <- as.data.frame(t(gene_norm_df_1))

pd
primary_sample <- pd[pd$sample_type=="Primary",]$sampleID
recurrent_sample <- pd[pd$sample_type=="Recurrent",]$sampleID
# primary.sample.dmp.mynorm <- all.sample.myNorm[,primary_sample]
# recurrent.sample.dmp.mynorm <- all.sample.myNorm[,recurrent_sample]
primary.sample.gene_norm_df_1 <- gene_norm_df_1[,primary_sample]
recurrent.sample.gene_norm_df_1 <- gene_norm_df_1[,recurrent_sample]

mean_beta_recurrent_sample <- rowMeans(recurrent.sample.gene_norm_df_1)


# check gene order
sum(rownames(deg_mean_control) == rownames(tumor_mat))
sum(rownames(primary.sample.gene_norm_df_1) == rownames(recurrent.sample.gene_norm_df_1))

# delta_beta_df <- data.frame(gene_id=rownames(tumor_mat))
delta_beta_df <- data.frame(gene_id=rownames(primary.sample.gene_norm_df_1))
head(delta_beta_df)
# tumor_coln
# for(tc in tumor_coln){
#   delta_beta_df[, tc] <- primary.sample.gene_norm_df_1[, tc] - mean_beta_recurrent_sample
# }
primary_sample
for(tc in primary_sample){
  delta_beta_df[, tc] <- primary.sample.gene_norm_df_1[, tc] - mean_beta_recurrent_sample
}
# fpkmM <- as.data.frame(fpkmM)
# coln <- colnames(fpkmM)
# coln <- str_remove(coln, 'cufflinksOut.aws.')
# coln <- str_remove(coln, 'cufflinksOut.d2.')
# coln <- str_remove(coln, 'cufflinksOut.67.')
# colnames(fpkmM) <- coln
# rownames(fpkmM) <- fpkmM$V1
# fpkmM1 <- fpkmM[,c(2:ncol(fpkmM))]
# deg_fpkm_mat <- fpkmM1[deg_symbol,]
# deg_fpkm_mat <- log2(deg_fpkm_mat + 1 )
delta_beta_df[1:5,1:5]

dmg_delta_beta_df <- t(delta_beta_df)
dmg_delta_beta_df <- as.data.frame(dmg_delta_beta_df)
colnames(dmg_delta_beta_df) <- as.character(delta_beta_df$gene_id)
dmg_delta_beta_df <- dmg_delta_beta_df[2:nrow(dmg_delta_beta_df),]
dmg_delta_beta_df[1:5,1:5]

dmg_delta_beta_df[,'patient01'] <-rownames(dmg_delta_beta_df)
# dmg_delta_beta_df[,'patient'] <- str_remove(dmg_delta_beta_df$patient, '_T1')
dmg_delta_beta_df[,'patient'] <- substr(dmg_delta_beta_df$patient01,1,12)
dmg_delta_beta_df[1:5,c(1:5, ncol(dmg_delta_beta_df))]
# saveRDS(dmg_delta_beta_df, './ana_20221116/all.sample.dmg_delta_beta_df.champ.rds')
#
dmg_delta_beta_df <- readRDS("/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/all.sample.dmg_delta_beta_df.champ.rds")
dmg_delta_beta_df[is.na(dmg_delta_beta_df)]
dmg_delta_beta_df[1:5,1:5]
rna_surv_info <- surv_df_3
rna_surv_info <- surv_df_3[!surv_df_3$time<0,]

# surv_info <- fread('./ana_20211105/update_table_20211117.csv')
tmp_surv_df <- as.data.frame(rna_surv_info)
head(tmp_surv_df)
colnames(tmp_surv_df) <- c('patient', 'OS', 'survival')
# tmp_surv_df[, 'patient'] <- str_remove(tmp_surv_df$Tumor_Sample_Barcode, '_T1')
# tmp_surv_df <- tmp_surv_df[,  c('PFS', 'OS', 'reoccurence', 'survival', 'patient')]

tmp_surv_df$patient %in% dmg_delta_beta_df$patient
sum(tmp_surv_df$patient %in% dmg_delta_beta_df$patient)
length(tmp_surv_df$patient %in% dmg_delta_beta_df$patient)
dmg_delta_beta_df$patient %in% tmp_surv_df$patient
sum(dmg_delta_beta_df$patient %in% tmp_surv_df$patient)

tmp_patient <- intersect(dmg_delta_beta_df$patient, tmp_surv_df$patient)
length(tmp_patient)
tmp_dmg_delta_beta_df <- dmg_delta_beta_df[dmg_delta_beta_df$patient %in% tmp_patient,]
tmp_surv_df <- tmp_surv_df[tmp_surv_df$patient %in% tmp_patient,]
tmp_surv_df[tmp_surv_df$OS<0,]
tmp_dmg_delta_beta_surv_df <- left_join(tmp_surv_df, dmg_delta_beta_df, by='patient')
tmp_dmg_delta_beta_surv_df[1:5,1:10]
write.table(tmp_dmg_delta_beta_surv_df,"./ana_20221121/tmp_dmg_delta_beta_surv_df.xls",row.names = F,sep = "\t",quote = F)

# 4 select survival significant genes --------
# 4.1mutation-----------------------
tmp_table <- mut_gene_surv_df
tmp_table$OS <- as.numeric(tmp_table$OS)
tmp_table$survival <- as.numeric(tmp_table$survival)
tmp_table_1 <- tmp_table[!is.na(tmp_table$OS) & !is.na(tmp_table$survival),]

tmp_table[is.na(tmp_table$OS),]
tmp_table[is.na(tmp_table$survival),]

# gene_v <- c()
# pvalue_v <- c()
# for (i in high_mut_gene_coad_V){
#   tmp <- coxph(Surv(tmp_table_1$OS, tmp_table_1$survival) ~ tmp_table_1[,i], data = tmp_table_1)
#   tmp1 <- summary(tmp)
#   gene_v <- c(gene_v, i)
#   pvalue_v <- c(pvalue_v, tmp1$waldtest[3])
# }
pvalL <- list()
for (i in high_mut_gene_coad_V){
  fit <- survfit(Surv(tmp_table_1$OS, tmp_table_1$survival) ~ tmp_table_1[,i], data = tmp_table_1)
  pva <- surv_pvalue(fit, data = tmp_table_1)
  pva[1,1] <- i
  pvalL[[i]] <- pva
  # sfn <- paste('./plots/7.cnv.selected.segment.survival/survival.plots/', i , '.survival.curve.pdf', sep='')
  # pdf(sfn)
  # p <- ggsurvplot(fit, pval = T)
  # print(p)
  # dev.off()
}

p.value.df <- data.table::rbindlist(pvalL)
p.value.df_bak <- p.value.df
p.value.df.order.mut <- p.value.df[order(p.value.df$pval),]

# 4.2deg--------------------------
tmp_table <- tmp_deg_fc_surv_df
tmp_table$OS <- as.numeric(tmp_table$OS)
tmp_table$survival <- as.numeric(tmp_table$survival)
tmp_table_1 <- tmp_table[!is.na(tmp_table$OS) & !is.na(tmp_table$survival),]
ncol(tmp_table_1)
tmp_table_1[1:5,1:5]
# tmp_table[is.na(tmp_table$OS),]
# tmp_table[is.na(tmp_table$survival),]

# coln[5898] is NA
coln <- colnames(tmp_table_1)

gene_v <- c()
pvalue_v <- c()
for (a in c(4:length(coln))){
  # for (a in c(5899:length(coln))){
  # print(c(a, coln[a]))
  # print()
  i <- coln[a]
  tmp <- coxph(Surv(tmp_table_1$OS, tmp_table_1$survival) ~ as.numeric(as.character(tmp_table_1[,i])), data = tmp_table_1)
  tmp1 <- summary(tmp)
  gene_v <- c(gene_v, i)
  pvalue_v <- c(pvalue_v, tmp1$waldtest[3])
}

p_df_deg <- data.frame(gene=gene_v, pval=pvalue_v)
p_df_deg_bak <- p_df_deg
p_df_deg_order <- p_df_deg[order(p_df_deg$pval),]

# 4.3dmp----------------------

# tmp_table <- tmp_deg_fc_surv_df
tmp_table <- tmp_dmg_delta_beta_surv_df
tmp_table$OS <- as.numeric(tmp_table$OS)
tmp_table$survival <- as.numeric(tmp_table$survival)
tmp_table_1 <- tmp_table[!is.na(tmp_table$OS) & !is.na(tmp_table$survival),]
ncol(tmp_table_1)
tmp_table_1[1:5,1:5]
# tmp_table[is.na(tmp_table$OS),]
# tmp_table[is.na(tmp_table$survival),]

# coln[5898] is NA
coln <- colnames(tmp_table_1)

gene_v <- c()
pvalue_v <- c()
for (a in c(4:length(coln)-1)){
  # for (a in c(5899:length(coln))){
  # print(c(a, coln[a]))
  # print()
  i <- coln[a]
  tmp <- coxph(Surv(tmp_table_1$OS, tmp_table_1$survival) ~ as.numeric(as.character(tmp_table_1[,i])), data = tmp_table_1)
  tmp1 <- summary(tmp)
  gene_v <- c(gene_v, i)
  pvalue_v <- c(pvalue_v, tmp1$waldtest[3])
}

p_df_dmp <- data.frame(gene=gene_v, pval=pvalue_v)
p_df_dmp_order <- p_df_dmp[order(p_df_dmp$pval),]


# deg_gene <- p_df_deg_order[complete.cases(p_df_deg_order),]

# 5.0 mutaition Lasso-----------------------------
library(glmnet)
library(survival)
library(survivalROC)

set.seed(111)
data=mut_gene_surv_df
data$OS.time <- data$OS
data$OS <- data$survival
data <- data[data$OS.time > 0,]
x=as.matrix(data[,c(4:(ncol(data) -1))])
y=data.matrix(Surv(data$OS.time,data$OS))
fit_mut_5.0 <- glmnet(x, y, family = "cox",alpha=1, maxit = 1000)
pdf('./ana_20221116/test.lasso.fit.parameter.plot.mut.model.pdf')
plot(fit_mut_5.0, xvar = "lambda", label = TRUE)
plot(fit_mut_5.0)
dev.off()

set.seed(111)
cvfit_mut_5.0 <- cv.glmnet(x, y, family="cox",alpha=1, maxit = 1000)
pdf('./ana_20221116/test.lasso.cv.fit.parameter.plot.mut.model.pdf')
plot(cvfit_mut_5.0)
abline(v=log(c(cvfit_mut_5.0$lambda.min,cvfit_mut_5.0$lambda.1se)),lty="dashed")
dev.off()
cvfit$lambda.min
cvfit$lambda.1se

# fit_mut_5.0 <- fit
# cvfit_mut_5.0 <- cvfit
# saveRDS(c('fit_mut_5.0', 'cvfit_mut_5.0'), './ana_20221116/fit.cvfit.for.mut.gene.lasso.selection.5.0.rds')
Coeff <- coef(fit_mut_5.0, s = cvfit_mut_5.0$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])

write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/lasso.feature.selection.mutation.res.xls",
            quote = F, sep = '\t')

mut_mat <- as.matrix(data[,colnames(data)%in%rownames(coef_df)])
# mut_mat1 <- data[,colnames(data)%in%rownames(coef_df)]
# mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat %*% as.numeric(coef_df$coef))
score_df$status <- data$OS
score_df$time <-data$OS.time

risk_surv_df_mut <- score_df
risk_surv_df_mut$group[risk_surv_df_mut$score >= median(risk_surv_df_mut$score)] = "high"
risk_surv_df_mut$group[risk_surv_df_mut$score < median(risk_surv_df_mut$score)] = "low"

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.mut.DNA_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_mut$time, risk_surv_df_mut$status) ~ risk_surv_df_mut$group, data = risk_surv_df_mut)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff mut DNA linear model')
                , data = risk_surv_df_mut)
print(p)
dev.off()

mut_gene_select <- rownames(coef_df)
mut_gene_select
mut_gene_select_surv_df <-  data[,c("Tumor_Sample_Barcode", "OS", "survival",mut_gene_select,"OS.time")]

# 5.0.1 train test------------------
# mut_mat <- as.matrix(as.numeric(data[,colnames(data)%in%rownames(coef_df)]))
data_train <- data[sample(1:nrow(data),round(nrow(data)*0.7)),c("Tumor_Sample_Barcode", "OS", "survival",mut_gene_select,"OS.time")]
data_test <- data[!(data$Tumor_Sample_Barcode%in% data_train$Tumor_Sample_Barcode),c("Tumor_Sample_Barcode", "OS", "survival",mut_gene_select,"OS.time")]
tmp_x_mut <- data_train[,c(4:(ncol(data_train) -1))]
# tmp_x_deg <- apply(tmp_x_deg, 2, as.character)
tmp_x_mut <- apply(tmp_x_mut, 2, as.numeric)
tmp_x_mut[1:5,1:5]

x=as.matrix(tmp_x_mut)
dim(x)
# x=as.matrix(data[,c(4:(ncol(data) -1))])
y=data.matrix(Surv(data_train$OS.time,data_train$OS))
dim(y)
x[is.na(x)]

set.seed(111)
cvfit_mut_5.1_train <- cv.glmnet(x, y, family="cox",alpha=0, maxit = 1000,nfolds = 10)
pdf('./ana_20221116/mut_train.lasso.cv.fit.parameter.plot.mut.DNA.model.pdf')
plot(cvfit_mut_5.1_train)
abline(v=log(c(cvfit_mut_5.1_train$lambda.min,cvfit_mut_5.1_train$lambda.1se)),lty="dashed")
dev.off()
cvfit_mut_5.1_train$lambda.min
cvfit_mut_5.1_train$lambda.1se
# fit_deg_5.1 <- fit
# cvfit_deg_5.1 <- cvfit
Coeff <- coef(cvfit_mut_5.1_train, s = cvfit_mut_5.1_train$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])

write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/mut_train.cffit.1.10.crossvalid.model.parmas.xls",
            quote = F, sep = '\t')

# intersect(data_train$patient,data_test$patient)
mut_mat <- as.matrix(data_train[,colnames(data_train)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_train$OS
score_df$time <-data_train$OS.time
risk_surv_df_train <- score_df
risk_surv_df_train$group[risk_surv_df_train$score >= median(risk_surv_df_train$score)] = "high"
risk_surv_df_train$group[risk_surv_df_train$score < median(risk_surv_df_train$score)] = "low"

write.table(risk_surv_df_train,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/mut_train.cvfit.1.10.crossvalid.risk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.mut.DNA._train_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_train$time, risk_surv_df_train$status) ~ risk_surv_df_train$group, data = risk_surv_df_train)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff mut DNA train linear model')
                , data = risk_surv_df_train)
print(p)
dev.off()

library(pROC)
roc_mut_train<-roc(risk_surv_df_train$status,risk_surv_df_train$score)
pdf("./ana_20221116/ROC.curve.train.dataset.MUT_DNA.model.pdf")
plot(roc_mut_train,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 5.0.1.2 DNA mut data_test ---------------------------------------------------------
mut_mat <- as.matrix(data_test[,colnames(data_test)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
# mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat %*% as.numeric(coef_df$coef))
score_df$status <- data_test$OS
score_df$time <-data_test$OS.time

risk_surv_df_test <- score_df
risk_surv_df_test$group[risk_surv_df_test$score >= median(risk_surv_df_test$score)] = "high"
risk_surv_df_test$group[risk_surv_df_test$score < median(risk_surv_df_test$score)] = "low"

write.table(risk_surv_df_test,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/mut_test.cvfit.1.10.crossvalid.risk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.mut.DNA._test_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_test$time, risk_surv_df_test$status) ~ risk_surv_df_test$group, data = risk_surv_df_test)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff mut DNA test linear model')
                , data = risk_surv_df_test)
print(p)
dev.off()

library(pROC)
roc_mut_test<-roc(risk_surv_df_test$status,risk_surv_df_test$score)
pdf("./ana_20221116/ROC.curve.test.dataset.MUT_DNA.model.pdf")
plot(roc_mut_test,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 5.1 RNA deg  Lasso-----------------------------
library(glmnet)
library(survival)
library(survivalROC)


data=tmp_deg_fc_surv_df
data$OS.time <- data$OS
data$OS <- data$survival
data <- data[data$OS.time > 0,]
data[1:5,1:5]
data[1:5,(ncol(data) - 5):ncol(data)]
data$OS
data$OS.time
tmp_x_deg <- data[,c(4:(ncol(data) -1))]
# tmp_x_deg <- apply(tmp_x_deg, 2, as.character)
tmp_x_deg <- apply(tmp_x_deg, 2, as.numeric)

tmp_x_deg[1:5,1:5]

x=as.matrix(tmp_x_deg)
# x=as.matrix(data[,c(4:(ncol(data) -1))])
y=data.matrix(Surv(data$OS.time,data$OS))
x[is.na(x)]
set.seed(111)
fit_deg_5.1 <- glmnet(x, y, family = "cox",alpha=1, maxit = 10000)
pdf('./ana_20221116/test.lasso.fit.parameter.plot.deg.model.pdf')
plot(fit_deg_5.1, xvar = "lambda", label = TRUE)
plot(fit_deg_5.1)
dev.off()


set.seed(111)
cvfit_deg_5.1 <- cv.glmnet(x, y, family="cox",alpha=1)
pdf('./ana_20221116/test.lasso.cv.fit.parameter.plot.deg.model.pdf')
plot(cvfit_deg_5.1)
abline(v=log(c(cvfit_deg_5.1$lambda.min,cvfit_deg_5.1$lambda.1se)),lty="dashed")
dev.off()
cvfit$lambda.min
cvfit$lambda.1se

# fit_deg_5.1 <- fit
# cvfit_deg_5.1 <- cvfit


Coeff <- coef(fit_deg_5.1, s = cvfit_deg_5.1$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])
coef_df

head(gene_ensg_df)
rownames(coef_df) %in% gene_ensg_df$gene_id
coef_df[,'gene_id_trim'] <- str_remove(rownames(coef_df), '\\.\\d+')

rownames(coef_df) %in% gene_ensg_df$gene_id
coef_df$gene_id %in% gene_ensg_df$gene_id_trim
rownames(coef_df)[!(rownames(coef_df) %in% gene_ensg_df$gene_id)]
'ENSG00000171540' %in% gene_ensg_df$gene_id_trim
gene_ensg_df[gene_ensg_df$gene_id_trim == 'ENSG00000171540',]

tmp_df <- gene_ensg_df[gene_ensg_df$gene_id_trim %in% coef_df$gene_id, ]
rownames(tmp_df) <- tmp_df$gene_id_trim
# coef_df[, 'gene_name'] <- tmp_df[coef_df$gene_id,'gene_name']
coef_df <- left_join(coef_df, tmp_df)

write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/lasso.feature.selection.RNA.deg.res.xls",
            quote = F, sep = '\t')

mut_mat <- as.matrix(data[,colnames(data)%in%rownames(coef_df)])
mut_mat1 <- data[,colnames(data)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data$OS
score_df$time <-data$OS.time

risk_surv_df_deg <- score_df
risk_surv_df_deg$group[risk_surv_df_deg$score >= median(risk_surv_df_deg$score)] = "high"
risk_surv_df_deg$group[risk_surv_df_deg$score < median(risk_surv_df_deg$score)] = "low"

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.deg.RNA.high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_deg$time, risk_surv_df_deg$status) ~ risk_surv_df_deg$group, data = risk_surv_df_deg)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff deg RNA linear model')
                , data = risk_surv_df_deg)
print(p)
dev.off()


deg_select <- rownames(coef_df)
deg_select

deg_select_surv_df <-  data[,c("patient", "OS", "survival",deg_select,"OS.time")]
deg_select_surv_df

# 5.1.1 train test------------------
ph <- '/thinker/aid2/udata/all/gao/others/tcga/lgg/'
setwd(ph)

# mut_mat <- as.matrix(as.numeric(data[,colnames(data)%in%rownames(coef_df)]))
data_train <- data[sample(1:nrow(data),round(nrow(data)*0.7)),c("patient", "OS", "survival",deg_select,"OS.time")]
data_test <- data[!(data$patient%in% data_train$patient),c("patient", "OS", "survival",deg_select,"OS.time")]
tmp_x_deg <- data_train[,c(4:(ncol(data_train) -1))]
# tmp_x_deg <- apply(tmp_x_deg, 2, as.character)
tmp_x_deg <- apply(tmp_x_deg, 2, as.numeric)

tmp_x_deg[1:5,1:5]

x=as.matrix(tmp_x_deg)
dim(x)
# x=as.matrix(data[,c(4:(ncol(data) -1))])
y=data.matrix(Surv(data_train$OS.time,data_train$OS))
x[is.na(x)]
set.seed(111)
cvfit_deg_5.1_train <- cv.glmnet(x, y, family="cox",alpha=0, nfolds = 10)
pdf('./ana_20221116/deg_train.lasso.cv.fit.parameter.plot.deg.model.pdf')
plot(cvfit_deg_5.1_train)
abline(v=log(c(cvfit_deg_5.1_train$lambda.min,cvfit_deg_5.1_train$lambda.1se)),lty="dashed")
dev.off()
cvfit_deg_5.1_train$lambda.min
cvfit_deg_5.1_train$lambda.1se
# fit_deg_5.1 <- fit
# cvfit_deg_5.1 <- cvfit
Coeff <- coef(cvfit_deg_5.1_train, s = cvfit_deg_5.1_train$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])
coef_df
write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/RNA.deg_train.cffit.1.10.crossvalid.model.parmas.xls",
            quote = F, sep = '\t')


# intersect(data_train$patient,data_test$patient)
mut_mat <- as.matrix(data_train[,colnames(data_train)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_train$OS
score_df$time <-data_train$OS.time

risk_surv_df_train <- score_df
risk_surv_df_train$group[risk_surv_df_train$score >= median(risk_surv_df_train$score)] = "high"
risk_surv_df_train$group[risk_surv_df_train$score < median(risk_surv_df_train$score)] = "low"

write.table(risk_surv_df_train,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/RNA.deg_train.cvffit.1.10.crossvalid.ridk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.deg.RNA._train_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_train$time, risk_surv_df_train$status) ~ risk_surv_df_train$group, data = risk_surv_df_train)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff deg RNA train linear model')
                , data = risk_surv_df_train)
print(p)
dev.off()

library(pROC)
roc_deg<-roc(risk_surv_df_train$status,risk_surv_df_train$score)
pdf("./ana_20221116/ROC.curve.train.dataset.DEG_rna.model.pdf")
plot(roc_deg,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 5.1.1.2 RNA deg data_test ---------------------------------------------------------

mut_mat <- as.matrix(data_test[,colnames(data_test)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_test$OS
score_df$time <-data_test$OS.time

risk_surv_df <- score_df
risk_surv_df$group[risk_surv_df$score >= median(risk_surv_df$score)] = "high"
risk_surv_df$group[risk_surv_df$score < median(risk_surv_df$score)] = "low"

write.table(risk_surv_df_test,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/RNA.deg_test.cvfit.1.10.crossvalid.risk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.deg.RNA._test_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df$time, risk_surv_df$status) ~ risk_surv_df$group, data = risk_surv_df)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff deg RNA test linear model')
                , data = risk_surv_df)
print(p)
dev.off()



library(pROC)
roc_deg<-roc(risk_surv_df$status,risk_surv_df$score)
pdf("./ana_20221116/ROC.curve.test.dataset.DEG_rna.model.pdf")
plot(roc_deg,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 5.3 dmp lasso selection --------
library(glmnet)
library(survival)
library(survivalROC)

set.seed(111)
data=tmp_dmg_delta_beta_surv_df
data$OS.time <- data$OS
data$OS <- data$survival
data <- data[data$OS.time > 0,]
data[1:5,1:5]
data[1:5,(ncol(data) - 5):ncol(data)]

tmp_x_dmp <- data[,c(4:(ncol(data) -2))]
tmp_x_dmp <- apply(tmp_x_dmp, 2, as.character)
tmp_x_dmp <- apply(tmp_x_dmp, 2, as.numeric)

tmp_x_dmp[1:5,1:5]

x=as.matrix(tmp_x_dmp)
y=data.matrix(Surv(data$OS.time,data$OS))
fit_dmp_5.3 <- glmnet(x, y, family = "cox",alpha=1)
pdf('./ana_20221116/test.lasso.fit.parameter.plot.dmp.model.pdf')
plot(fit_dmp_5.3, xvar = "lambda", label = TRUE)
plot(fit_dmp_5.3)
dev.off()

set.seed(111)
cvfit_dmp_5.3 <- cv.glmnet(x, y, family="cox",alpha=1)
pdf('./ana_20221116/test.lasso.cv.fit.parameter.plot.dmp.model.pdf')
plot(cvfit_dmp_5.3)
abline(v=log(c(cvfit_dmp_5.3$lambda.min,cvfit_dmp_5.3$lambda.1se)),lty="dashed")
dev.off()
cvfit$lambda.min
cvfit$lambda.1se

# fit_deg_5.1 <- fit
# cvfit_deg_5.1 <- cvfit


Coeff <- coef(fit_dmp_5.3, s = cvfit_dmp_5.3$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])
coef_df
write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/lasso.feature.selection.RNA.dmg.res.xls",
            quote = F, sep = '\t')

mut_mat <- as.matrix(data[,colnames(data)%in%rownames(coef_df)])
mut_mat1 <- data[,colnames(data)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat1, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data$OS
score_df$time <-data$OS.time

risk_surv_df_dmp <- score_df
risk_surv_df_dmp$group[risk_surv_df_dmp$score >= median(risk_surv_df_dmp$score)] = "high"
risk_surv_df_dmp$group[risk_surv_df_dmp$score < median(risk_surv_df_dmp$score)] = "low"

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.gene.DMP_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_dmp$time, risk_surv_df_dmp$status) ~ risk_surv_df_dmp$group, data = risk_surv_df_dmp)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff delta DMP linear model')
                , data = risk_surv_df_dmp)
print(p)
dev.off()

dmp_gene_select <- rownames(coef_df)
dmp_gene_select
dmp_gene_select_surv_df <-  data[,c("patient", "OS", "survival",dmp_gene_select,"OS.time")]
dmp_gene_select_surv_df
# 5.3.1 DMP train test------------------
# mut_mat <- as.matrix(as.numeric(data[,colnames(data)%in%rownames(coef_df)]))
data_train <- data[sample(1:nrow(data),round(nrow(data)*0.7)),c("patient", "OS", "survival",dmp_gene_select,"OS.time")]
data_test <- data[!(data$patient%in% data_train$patient),c("patient", "OS", "survival",dmp_gene_select,"OS.time")]
tmp_x_mut <- data_train[,c(4:(ncol(data_train) -1))]
# tmp_x_deg <- apply(tmp_x_deg, 2, as.character)
tmp_x_mut <- apply(tmp_x_mut, 2, as.numeric)
tmp_x_mut[1:5,1:5]

x=as.matrix(tmp_x_mut)
dim(x)
# x=as.matrix(data[,c(4:(ncol(data) -1))])
y=data.matrix(Surv(data_train$OS.time,data_train$OS))
dim(y)
x[is.na(x)]

set.seed(111)
cvfit_dmp_5.1_train <- cv.glmnet(x, y, family="cox",alpha=0, maxit = 1000,nfolds = 10)
pdf('./ana_20221116/dmp_train.lasso.cv.fit.parameter.plot.gene.DMP.model.pdf')
plot(cvfit_dmp_5.1_train)
abline(v=log(c(cvfit_dmp_5.1_train$lambda.min,cvfit_dmp_5.1_train$lambda.1se)),lty="dashed")
dev.off()
cvfit_dmp_5.1_train$lambda.min
cvfit_dmp_5.1_train$lambda.1se
# fit_deg_5.1 <- fit
# cvfit_deg_5.1 <- cvfit
Coeff <- coef(cvfit_dmp_5.1_train, s = cvfit_dmp_5.1_train$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])
coef_df
write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/methy_train.cffit.1.10.crossvalid.model.parmas.xls",
            quote = F, sep = '\t')
# intersect(data_train$patient,data_test$patient)
mut_mat <- as.matrix(data_train[,colnames(data_train)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_train$OS
score_df$time <-data_train$OS.time
risk_surv_df_train <- score_df
risk_surv_df_train$group[risk_surv_df_train$score >= median(risk_surv_df_train$score)] = "high"
risk_surv_df_train$group[risk_surv_df_train$score < median(risk_surv_df_train$score)] = "low"

write.table(risk_surv_df_train,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/methy_train.cvffit.1.10.crossvalid.ridk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.gene.DMP._train_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_train$time, risk_surv_df_train$status) ~ risk_surv_df_train$group, data = risk_surv_df_train)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff mut DNA train linear model')
                , data = risk_surv_df_train)
print(p)
dev.off()

library(pROC)
roc_dmp_train<-roc(risk_surv_df_train$status,risk_surv_df_train$score)
pdf("./ana_20221116/ROC.curve.train.dataset.gene.DMP.model.pdf")
plot(roc_dmp_train,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 5.3.1.2 DMP data_test ---------------------------------------------------------
mut_mat <- as.matrix(data_test[,colnames(data_test)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_test$OS
score_df$time <-data_test$OS.time

risk_surv_df_test <- score_df
risk_surv_df_test$group[risk_surv_df_test$score >= median(risk_surv_df_test$score)] = "high"
risk_surv_df_test$group[risk_surv_df_test$score < median(risk_surv_df_test$score)] = "low"

write.table(risk_surv_df_test,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/methy_test.cvfit.1.10.crossvalid.risk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.gene.DMP._test_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_test$time, risk_surv_df_test$status) ~ risk_surv_df_test$group, data = risk_surv_df_test)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff gene DMP test linear model')
                , data = risk_surv_df_test)
print(p)
dev.off()
getwd()

library(pROC)
roc_dmp_test<-roc(risk_surv_df_test$status,risk_surv_df_test$score)
pdf("./ana_20221116/ROC.curve.test.dataset.gene.DMP.model.pdf")
plot(roc_dmp_test,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 6.0 three merge-----------------------------

temp_x_mut <- mut_gene_select_surv_df
temp_x_deg <- deg_select_surv_df
temp_x_dmp <- dmp_gene_select_surv_df
write.table(mut_gene_select_surv_df,"./ana_20221121/mut_gene_select_surv_df.xls",row.names = F,sep = "\t",col.names = F,quote = F)
write.table(mut_gene_select_surv_df,"./ana_20221121/deg_select_surv_df.xls",row.names = F,sep = "\t",col.names = F,quote = F)
write.table(mut_gene_select_surv_df,"./ana_20221121/dmp_gene_select_surv_df.xls",row.names = F,sep = "\t",col.names = F,quote = F)

sample_mut <- temp_x_mut$Tumor_Sample_Barcode
sample_deg <- temp_x_deg$patient
sample_dmp <- temp_x_dmp$patient
mut_deg_sample <- intersect(sample_mut ,sample_deg)
all_three_sample <- intersect(mut_deg_sample,sample_dmp)
all_three_sample

need_temp_x_mut <- temp_x_mut[temp_x_mut$Tumor_Sample_Barcode%in%all_three_sample,]
colnames(need_temp_x_mut)
need_temp_x_deg <- temp_x_deg[temp_x_deg$patient%in%all_three_sample,]
need_temp_x_dmp <- temp_x_dmp[temp_x_dmp$patient%in%all_three_sample,]


tmp_all_mat1 <- left_join(need_temp_x_mut, need_temp_x_deg, by=c("Tumor_Sample_Barcode"="patient"))
tmp_all_mat1$OS.x==tmp_all_mat1$OS.y
tmp_all_mat1

tmp_all_mat2 <- left_join(tmp_all_mat1, need_temp_x_dmp, by=c("Tumor_Sample_Barcode"='patient'))
tmp_all_mat2

set.seed(111)
deg_select
mut_gene_select
dmp_gene_select

data=tmp_all_mat2[,c("Tumor_Sample_Barcode","OS","survival",mut_gene_select,deg_select,dmp_gene_select,"OS.time")]
# data$OS.time <- data$OS
# data$OS <- data$survival
data <- data[data$OS.time > 0,]
data[1:5,1:5]
data[1:5,(ncol(data) - 5):ncol(data)]

# tmp_x_three <- data[,c(4:(ncol(data) -1))]
# tmp_x_three <- apply(tmp_x_three, 2, as.character)
# tmp_x_three <- apply(tmp_x_three, 2, as.numeric)

# tmp_x_three[1:5,1:5]
# 6.3.1 three train test------------------
# mut_mat <- as.matrix(as.numeric(data[,colnames(data)%in%rownames(coef_df)]))
data_train <- data[sample(1:nrow(data),round(nrow(data)*0.7)),]
data_test <- data[!(data$Tumor_Sample_Barcode%in% data_train$Tumor_Sample_Barcode),]
tmp_x_mut <- data_train[,c(4:(ncol(data_train) -1))]
tmp_x_mut <- apply(tmp_x_mut, 2, as.character)
tmp_x_mut <- apply(tmp_x_mut, 2, as.numeric)
tmp_x_mut[1:5,1:5]

x=as.matrix(tmp_x_mut)
dim(x)
# x=as.matrix(data[,c(4:(ncol(data) -1))])
y=data.matrix(Surv(data_train$OS.time,data_train$OS))
dim(y)
x[is.na(x)]

set.seed(111)
cvfit_three_5.1_train <- cv.glmnet(x, y, family="cox",alpha=0, maxit = 2000,nfolds = 10)
pdf('./ana_20221116/three_train.lasso.cv.fit.parameter.plot.gene.three.model.pdf')
plot(cvfit_three_5.1_train)
abline(v=log(c(cvfit_three_5.1_train$lambda.min,cvfit_three_5.1_train$lambda.1se)),lty="dashed")
dev.off()
cvfit_three_5.1_train$lambda.min
cvfit_three_5.1_train$lambda.1se
# fit_deg_5.1 <- fit
# cvfit_deg_5.1 <- cvfit
Coeff <- coef(cvfit_three_5.1_train, s = cvfit_three_5.1_train$lambda.min)
Active.Index <- which(Coeff != 0)
Active.Coefficients <- Coeff[Active.Index]
Active.Index
Active.Coefficients

coef_df <- as.matrix(Coeff)
coef_df <- data.frame(coef=coef_df[coef_df[,1]!=0,])
coef_df
write.table(coef_df,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/three_train.cffit.1.10.crossvalid.model.parmas.xls",
            quote = F, sep = '\t')

# intersect(data_train$patient,data_test$patient)
mut_mat <- as.matrix(data_train[,colnames(data_train)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_train$OS
score_df$time <-data_train$OS.time
risk_surv_df_train <- score_df
risk_surv_df_train$group[risk_surv_df_train$score >= median(risk_surv_df_train$score)] = "high"
risk_surv_df_train$group[risk_surv_df_train$score < median(risk_surv_df_train$score)] = "low"
dim(risk_surv_df_train)

write.table(risk_surv_df_train,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/three_train.cvffit.1.10.crossvalid.ridk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.gene.three._train_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_train$time, risk_surv_df_train$status) ~ risk_surv_df_train$group, data = risk_surv_df_train)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff train linear model')
                , data = risk_surv_df_train)
print(p)
dev.off()

library(pROC)
roc_three_train<-roc(risk_surv_df_train$status,risk_surv_df_train$score)
pdf("./ana_20221116/ROC.curve.train.dataset.gene.three.model.pdf")
plot(roc_three_train,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()

# 6.3.1.2 three data_test ---------------------------------------------------------
mut_mat <- as.matrix(data_test[,colnames(data_test)%in%rownames(coef_df)])
#mut_mat1 <- data[,colnames(data_train)%in%rownames(coef_df)]
mut_mat2 <- apply(mut_mat, 2, as.numeric)
score_df <- data.frame(score = mut_mat2 %*% as.numeric(coef_df$coef))
score_df$status <- data_test$OS
score_df$time <-data_test$OS.time

risk_surv_df_test <- score_df
risk_surv_df_test$group[risk_surv_df_test$score >= median(risk_surv_df_test$score)] = "high"
risk_surv_df_test$group[risk_surv_df_test$score < median(risk_surv_df_test$score)] = "low"
dim(risk_surv_df_test)
write.table(risk_surv_df_test,
            file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/three_test.cvfit.1.10.crossvalid.risk.df.xls",
            quote = F, sep = '\t')

pdf('./ana_20221116/elastic.net.model.group.by.pcutoff.gene.three._test_high.fc.featrue.median.survival.curve.pdf',onefile = F)
fit <- survfit(Surv(risk_surv_df_test$time, risk_surv_df_test$status) ~ risk_surv_df_test$group, data = risk_surv_df_test)
p <- ggsurvplot(fit, pval = T, risk.table = T
                , title=paste('survival curve of  pcutoff three test linear model')
                , data = risk_surv_df_test)
print(p)
dev.off()

library(pROC)
roc_three_test<-roc(risk_surv_df_test$status,risk_surv_df_test$score)
pdf("./ana_20221116/ROC.curve.test.dataset.three.model.pdf")
plot(roc_three_test,print.auc=TRUE,plot=TRUE,print.thres=TRUE,col="red")
dev.off()
# # 7.5.3.1 select deg and demir using lasso --------
# x_deg <- tmp_deg_fc_surv_df
# tmp_x_deg <- x_deg[,c('OS', 'survival', deg_gene_v)]
# tmp_x_deg <- apply(tmp_x_deg, 2, as.character)
# tmp_x_deg <- apply(tmp_x_deg, 2, as.numeric)
# tmp_x_deg <- tmp_x_deg[complete.cases(tmp_x_deg),]
#
# x_deg_mat <- tmp_x_deg[,deg_gene_v]
#
# y_deg <- tmp_x_deg[, c('OS', 'survival')]
#
# # x <- apply(x_deg_mat, 2, as.character)
# # x <- apply(x, 2, as.numeric)
# #
# # y <- apply(y_deg, 2, as.character)
# # y <- apply(y_deg, 2, as.numeric)
# # x <- as.matrix(x_deg_mat)
# cv.fit <- cv.glmnet(x_deg_mat, Surv(y_deg[,'OS'], y_deg[,'survival'])
#                     , family='cox', maxit=1000, alpha=1, nfolds = 4)
# # fit <- glmnet(x, Surv(y$time, y$status), family='cox', maxit=1000, alpha=1)
#
# plot(cv.fit)
# coefficience_cv.fit <- coef(cv.fit, s = cv.fit$lambda.min)
# active.index <- which(coefficience_cv.fit != 0)
# active.coef <- coefficience_cv.fit[active.index,]
# x_mir <- tmp_mir_fc_surv_df


# 3.1 配对样本chAMP DMG survival p value -------
# .libPaths(c("/thinker/aid2/udata/all/gao/R/R4.1.3/lib/R/library","/thinker/aid2/udata/all/muyl/others/software/condaEnvs/methyArray3/lib/R/library/", .libPaths()))
# library(data.table)
# library(Hmisc)
# #表型数据
# pd.all <- read.delim("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA.LGG.sampleMap%2FLGG_clinicalMatrix", header = T, stringsAsFactors = F)
# pd.all
# pd <- pd.all[,c("sampleID","sample_type","bcr_sample_barcode","bcr_patient_barcode")]
# pd$sample_type <- ifelse(pd$sample_type=="Primary Tumor","Primary","Recurrent")
# head(pd)
# pd
# recurrent_sample <- pd[pd$sample_type=="Recurrent",]$bcr_patient_barcode
# pd01 <- pd[pd$bcr_patient_barcode %in% recurrent_sample,]
# rownames(pd01)=pd01$sampleID
# pd=na.omit(pd01)
# pd
# nrow(pd)
# # 保存配对样本
# # write.table(pd01$sampleID, file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_sampleMap_samples01.txt",quote = F, row.names = F, col.names = F)
# #$ cols=($(sed '1!d;s/\t/\n/g' TCGA.LGG.sampleMapFHumanMethylation450 | grep -nf TCGA-LGG_Recurrent_Primary_sampleMap_samples01.txt | sed 's/:.*$//'))
#
# # sed '1!d;s/\t/\n/g' TCGA.LGG.sampleMapFHumanMethylation450
# # $ cut -f 1$(printf ",%s" "${cols[@]}") TCGA.LGG.sampleMapFHumanMethylation450 > TCGA-LGG_Recurrent_Primary_sampleMap_methy450k.txt
#
# a01 <- fread("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/TCGA-LGG_Recurrent_Primary_sampleMap_methy450k.txt", data.table = F )
# a01[1:4,1:4]
# colnames(a01)
# # colnames(a02)
# # substr(colnames(a02), 1,15)
# # substr(colnames(a02), 1,15) %in% colnames(a01)
# # colnames()
#
# # sum(a01$sample %in% a02$`Composite Element REF`)
# # head(a01$sample)
# nrow(a01)
# rownames(a01)=a01[,1]
# a=a01[,-1]
# a=na.omit(a)
# library(ChAMP)
# betaData <- as.matrix(a)
# # 匹配临床数据和methy数据
# pd <- pd[match(colnames(a),row.names(pd)),]
# nrow(pd)
# dim(betaData)
# identical(colnames(betaData),row.names(pd))
# colnames(betaData)[1:10]
# row.names(pd)[1:10]
# myLoad=champ.filter(beta = betaData ,pd = pd) #这一步已经自动完成了过滤
# dim(myLoad$beta)
# # save(myLoad,file = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/step1-output.Rdata')
# # library(ChAMP)
# # load('/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/step1-output.Rdata')
# # 数据归一化
# # getwd()
# myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,plotBMIQ=FALSE,resultsDir="./ana_20221116/CHAMP_Normalization_paired/")
# # saveRDS(myNorm, './ana_20221116/paired.sample.myNorm.champ.rds')
# # 主成分分析
# # library("FactoMineR")
# # library("factoextra")
# # dat <- t(myNorm)
# # group_list=pd$sample_type
# # pd <- myLoad$pd[colnames(myNorm),] #去掉异常样本
# # group_list=pd$sample_type
# # table(group_list)
# # pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/01_methy_PCA.pdf")
# # pdf("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/01_methy_PCA.pdf")
# # dat.pca <- PCA(dat, graph = FALSE)
# # fviz_pca_ind(dat.pca,
# #              geom.ind = "point",
# #              col.ind = group_list,
# #              addEllipses = TRUE,
# #              legend.title = "Groups")
# # dev.off()
# # # 热图
# # group_list=pd$sample_type
# # 第一个参数是指要参与计算的矩阵；
# # 第二个参数是指按行计算还是按列计算，1——表示按行计算，2——按列计算；
# # 第三个参数是指具体的运算参数。
# # sd函数计算数据列或者向量标准差（standard deviation）
# #  tail() 函数用于获取向量、矩阵、表、 DataFrame 或函数的最后部分
# # cg=names(tail(sort(apply(myNorm,1,sd)),1000))
# # library(pheatmap)
# # ac=data.frame(group=group_list)
# #
# # rownames(ac)=colnames(myNorm)
# # pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
# #          annotation_col=ac,filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/02_heatmap_top1000_sd.pdf')
# # # 相关关系矩阵热图
# # # 组内的样本的相似性应该是要高于组间的！
# # pheatmap::pheatmap(cor(myNorm[cg,]),
# #                    annotation_col = ac,
# #                    show_rownames = F,
# #                    show_colnames = F,
# #                    filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/02_heatmap_top1000_corresponding.pdf')
# # pd <- pd[colnames(myNorm),]
# # # save(pd,myNorm,file = "./Rdata/filtered.Rdata")
# # save(pd,myNorm,file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/filtered.Rdata")
# group_list <- pd$sample_type
# myDMP.pair.01 <- champ.DMP(beta = myNorm,pheno=group_list)
# # only 51 dmp between paired samples
# nrow(myDMP.pair.01$Recurrent_to_Primary)
#
#
# # myDMP_01 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "none",adjPVal = 0.001)
# myDMP_02 <- champ.DMP(beta = myNorm,pheno=group_list,adjust.method = "none",adjPVal = 1)
# hmc_02 <- myDMP_02[[1]][myDMP_02[[1]]$deltaBeta>0,]
# write.table(hmc_02,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/hydroxymethylation_CpGs02.xls",sep = "\t")
# logFC_cutoff <- 0.001
#
# # df_DMP_02$change <- ifelse(df_DMP_02$adj.P.Val < 10^-5 & abs(df_DMP_02$logFC) > logFC_cutoff,
# #                            ifelse(df_DMP_02$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
#
# df_DMP_02$change <- ifelse(df_DMP_02$adj.P.Val < 10^-3 & abs(df_DMP_02$logFC) > logFC_cutoff,
#                            ifelse(df_DMP_02$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
# table(df_DMP_02$change)
# write.table(df_DMP_02,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_002.xls",quote = F, row.names = T, col.names = T,sep = "\t")
#
#
# df_DMP_02=df_DMP_01[df_DMP_01$gene!="",]
# deltabeta_t <- 0.0001
#
# P.Value_t <- 10^-10
#
# df_DMP_02$change <- ifelse(df_DMP_02$adj.P.Val < P.Value_t & abs(df_DMP_02$deltaBeta) > deltabeta_t,
#                            ifelse(df_DMP_02$deltaBeta > deltabeta_t ,'UP','DOWN'),'NOT')
# table(df_DMP_02$change)
# #>
# #>   DOWN    NOT     UP
# #>    345 108379    814
# save(df_DMP_02,file = "/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out/step3.df_DMP_02.Rdata")
# dat = rownames_to_column(df_DMP_02)
# for_label <- dat %>% head(3)
#
#
# # 区域：
# myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$sample_type,method="Bumphunter")
# write.table(myDMR$BumphunterDMR,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/myDMR.xls",sep = "\t")
# gene.v <- unique(as.character(df_DMP_02[df_DMP_02$change!="NOT",]$gene))
# write.table(gene.v,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_gene.list",quote = F, row.names = F, col.names = F)
# gene.v <- read.csv("/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_gene.list",header = F)
#
#


library(survival)
library(survminer)

pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/model/os.high.mut.gene.all.available.patient.survival.curve.pdf')
for (gene in high_mut_gene_coad_V){
  group <- os_df[,gene]
  fit <- survfit(Surv(os_df$OS, os_df$survival) ~ group, data = os_df)
  pva <- surv_pvalue(fit, data = os_df)
  pva$variable <-  gene
  print(pva)
  p <- ggsurvplot(fit, pval = T, risk.table = T
                  , title=paste('survival curve of ', gene, ' mutattion status in LGG_paired')
                  , data = os_df)
  print(p)
}
dev.off()

mut_gene_surv_df <- surv_df_paired
# degMat
# all.deg.N.VS.T.pcut
# 7.5.2.2 prepare fold change matrix for deg--------
# deg gene name
gene_ensg_df
tmp_deg <- all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05
tmp_deg <- as.data.frame(tmp_deg)

# vst reaults is log2 transformed
vsd <- vst(dds)
# vsd <- vst(dds)
countMat <- assay(vsd)
degM <- countMat[rownames(tmp_deg),]
coln <- colnames(degM)
# 这里是要干啥？
coln <- str_remove(coln, 'P')



all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05

deg_ensg <- rownames(tmp_deg)
deg_ensg <- str_remove(deg_ensg, '\\.\\d+')

deg_ensg %in% gene_ensg_df$gene_id_trim
deg_symbol <- gene_ensg_df[gene_ensg_df$gene_id_trim %in% deg_ensg, 'gene_name']
deg_symbol <- unique(deg_symbol)






# 8.0 kegg----------------------------------------------

library(org.Hs.eg.db)
library(clusterProfiler)

library(stringr)
library(BSgenome)
#library(RIdeogram)
library(circlize)
R.utils::setOption( "clusterProfiler.download.method",'wget')
# gene.v <- unique(as.character(df_DMP$gene))
# write.table(gene.v,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired/df_DMP_gene.list",quote = F, row.names = F, col.names = F)
gene.v <- unique(as.character(df_DMP_02[df_DMP_02$change!="NOT",]$gene))
# write.table(gene.v,"/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg_meth_gao/methy_out_paired_001/df_DMP_gene.list",quote = F, row.names = F, col.names = F)
write.table(gene.v,"/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116//df_DMP_gene.list",quote = F, row.names = F, col.names = F)
# 先在终端产生这个gene.v

gene.v <- read.csv("/mnt/phoenix/gao/projects/test_kegg/df_DMP_gene.list",header = F)
gene.v <- read.csv("/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116//df_DMP_gene.list",header = F)
gene_L <- gene.v$V1
oL <- bitr(gene_L, fromType='SYMBOL', toType='UNIPROT', OrgDb = org.Hs.eg.db)

# gene_L <- rownames(all.N.VS.T.deg.res.pcut.0.05)
# gene_L <- str_remove(gene_L, '\\.\\d+')
# oL <- bitr(gene_L, fromType='ENSEMBL', toType='UNIPROT', OrgDb = org.Hs.eg.db)
oL <- bitr(gene_L, fromType='SYMBOL', toType='UNIPROT', OrgDb = org.Hs.eg.db)
# KEGG富集的数据库
ken <- enrichKEGG(oL$UNIPROT,organism = 'hsa',keyType = 'uniprot', pAdjustMethod='none', pvalueCutoff=0.05)
ken
# GO三个数据库
goCC <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='CC', pAdjustMethod='none', pvalueCutoff=0.05)
goBP <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='BP', pAdjustMethod='none', pvalueCutoff=0.05)
goMF <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='MF', pAdjustMethod='none', pvalueCutoff=0.05)


pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/ana_20221116/samples.barplot.DMP.all.gene.functional.analysis.pdf')
pdf('/mnt/phoenix/gao/projects/test_kegg/samples.barplot.DMP.all.gene.functional.analysis.pdf')
p <- dotplot(ken,title='KEGG enrich of mutation gene')
print(p)
p <- dotplot(goCC,title='GO CC enrich of mutation gene')
print(p)
p <- dotplot(goBP,title='GO BP enrich of mutation gene')
print(p)
p <- dotplot(goMF,title='GO MF enrich of mutation gene')
print(p)

dev.off()
save.image("/thinker/aid2/udata/all/gao/others/tcga/lgg/05.model_20221121.RData")

# 9.0 kegg_of_all_sample----------------------------------------------
# 在终端：

ph <- '/thinker/aid2/udata/all/gao/others/tcga/lgg/'
setwd(ph)
library(org.Hs.eg.db)
library(clusterProfiler)

library(stringr)
library(BSgenome)
#library(RIdeogram)
library(circlize)
R.utils::setOption( "clusterProfiler.download.method",'wget')

# 先在终端产生这个gene.v,根据3.所有样本chAMP产生这个基因。
write.csv(all.sample.dmp.gene,"/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/all.sample.dmp.gene.list",row.names = F )
gene.v <- read.csv("/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/all.sample.dmp.gene.list")
gene_L <- gene.v$x
oL <- bitr(gene_L, fromType='SYMBOL', toType='UNIPROT', OrgDb = org.Hs.eg.db)

# gene_L <- rownames(all.N.VS.T.deg.res.pcut.0.05)
# gene_L <- str_remove(gene_L, '\\.\\d+')
# oL <- bitr(gene_L, fromType='ENSEMBL', toType='UNIPROT', OrgDb = org.Hs.eg.db)
oL <- bitr(gene_L, fromType='SYMBOL', toType='UNIPROT', OrgDb = org.Hs.eg.db)
# KEGG富集的数据库
ken <- enrichKEGG(oL$UNIPROT,organism = 'hsa',keyType = 'uniprot', pAdjustMethod='none', pvalueCutoff=0.05)
ken
# GO三个数据库
goCC <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='CC', pAdjustMethod='none', pvalueCutoff=0.05)
goBP <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='BP', pAdjustMethod='none', pvalueCutoff=0.05)
goMF <- enrichGO(gene_L,org.Hs.eg.db,keyType='SYMBOL', ont='MF', pAdjustMethod='none', pvalueCutoff=0.05)



pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/all_samples.barplot.DMP.all.gene.functional_KEGG_GO.analysis.pdf')
p <- dotplot(ken,title='KEGG enrich of differentially methylated gene')
print(p)
p <- dotplot(goCC,title='GO CC enrich of differentially methylated gene')
print(p)
p <- dotplot(goBP,title='GO BP enrich of differentially methylated gene')
print(p)
p <- dotplot(goMF,title='GO MF enrich of differentially methylated gene')
print(p)

dev.off()
write.csv(ken@result, '/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/KEGG.enrichment.all.deg.pvalue.table.csv')
write.csv(goCC@result, '/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/GO.CC.enrichment.all.deg.pvalue.table.csv')
write.csv(goBP@result, '/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/GO.BP.enrichment.all.deg.pvalue.table.csv')
write.csv(goMF@result, '/thinker/aid2/udata/all/gao/others/tcga/lgg/all_sample_kegg/GO.MF.enrichment.all.deg.pvalue.table.csv')
save.image("/thinker/aid2/udata/all/gao/others/tcga/lgg/06.all_sample_kegg_go20230310.RData")
