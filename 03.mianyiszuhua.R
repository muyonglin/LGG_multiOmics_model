# 0. setup environment  ---------------------------------------------------
.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.2.0/lib/R/library","/thinker/storage/udata/sunhr/R/R4.2.0/lib/R/library/", .libPaths()))
.libPaths()
library(data.table)
library(stringr)

library(vidger)
library(clusterProfiler)

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
# 将reads的count值转换为FPKM值和TPM值
# 该代码可以实现count转为FPKM
##首先需要在gtf文件中得到外显子的'chr','start','end', 'gene_id', 'gene_type', 'gene_id_trim'这些列
#注意gtf文件位置，下面的是Phoenix服务器的
# gtf <- '/mnt/phoenix/bio-web/pipeline/rnaseq/reference/reference/gencode.v34.annotation.gtf'
rm(list = ls())
gtf <- '/thinker/storage/udata/muyl/genomes/human/GRCH38/gencode.v34.annotation.gtf'
inT <- data.table::fread(gtf, skip = 5)
inT <- as.data.frame(inT)
inT
tmp <- inT[inT$V3 == 'exon',]
tmp[,"chr"]<-tmp$V1
tmp[,'start']<-tmp$V4
tmp[,'end']<-tmp$V5
tmp[,'gene_id'] <- stringr::str_extract(tmp$V9, 'ENSG\\d+\\.\\d+')
tmp[,'gene_type'] <- stringr::str_extract(tmp$V9, 'gene_type ".*?"; gene_name')
tmp$gene_type <- stringr::str_remove(tmp$gene_type, 'gene_type "')
tmp$gene_type <- stringr::str_remove(tmp$gene_type, '"; gene_name')
tmp[,'gene_id_trim'] <- str_remove(tmp$gene_id, '\\.\\d+')
exon_ensg_df <- tmp[,c('chr','start','end', 'gene_id', 'gene_type', 'gene_id_trim')]
tmp
tmp[,'gene_name'] <- stringr::str_extract(tmp$V9, 'gene_name ".*?"; transcript_type')

tmp$gene_name <- stringr::str_remove(tmp$gene_name, 'gene_name "')
tmp$gene_name <- stringr::str_remove(tmp$gene_name, '"; transcript_type')
exon_ensg_df02 <- tmp[,c('chr','start','end', 'gene_id', 'gene_type', 'gene_id_trim','gene_name')]
exon_ensg_df02
# 计算外显子长度
g_l = lapply(split(exon_ensg_df,exon_ensg_df$gene_id),function(x){
  # x=split(t1,t1$geneid)[[1]]
  tmp=apply(x,1,function(y){
    y[2]:y[3]
  })
  length(unique(unlist(tmp)))
  # sum(x[,4])
})
g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))
g_l

# 1. differential expression analysis --------
getwd()
htseq_l <- list.files('./star_count', pattern = '*star_count', full.names = T)
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
ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable
                                       , directory = './'
                                       , design = ~condition)
register(MulticoreParam(20))
keep <- rowSums(counts(ddsHTseq)) >= min(1 * ncol(ddsHTseq), 10)
dds <- ddsHTseq[keep, ]
dds <- DESeq(dds, parallel = T)
dds@colData
dds@assays
##获取count数文件，并得到各个样本的reads总数
dds_count<-assay(dds)
dds_count
ng=intersect(rownames(dds_count),g_l$gene_id)

exprSet=dds_count[ng,]
lengths=g_l[match(ng,g_l$gene_id),2]
#head(lengths)
#head(rownames(exprSet))
#exprSet[1:4,1:4]
total_count<- colSums(exprSet)

head(total_count)

#head(lengths)
#total_count[1]
#lengths[2]
#53*10^9/(1476*3978357)

#计算FPKM
#计算公式FPKM = reads(基因reads数)*10^9/exon(外显子长度)*all reads(每个样本的reads总数)
fpkm <- t(do.call( rbind,
                   lapply(1:length(total_count),
                          function(i){
                            10^9*exprSet[,i]/lengths/total_count[i]
                          }) ))
colnames(fpkm)<-colnames(dds_count)
fpkm
# total_fpkm<- colSums(fpkm)
# total_fpkm
getwd()
write.table(fpkm,"/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/star_count_to_fpkm.csv",sep = "\t",quote = F)

# 使用FPKM数据作为输入时的前处理过程
# 将各样本cufflinks结果文件夹中的FPKM文件，合并为一个FPKM矩阵，用于后续处理
# fpkm_01path <- '/thinker/storage/biostaff/bioinformatics.code.tutorial.attchment/muyl/cibersort/cufflinks.fpkm.del.dup.del.zero.mat'
# fpkm.01mat <- fread(fpkm_01path)
# head(fpkm.01mat)
# fpkm.01mat <- as.data.frame(fpkm.01mat)
# fpkm.01mat
# matrix_fpkm <- as.matrix(total_fpkm)
# matrix_fpkm
fpkm.mat <- as.data.frame(fpkm)
fpkm.mat

fpkm.mat <- fpkm.mat[complete.cases(fpkm.mat),]
coln <- colnames(fpkm.mat)
coln
# fpkm.mat <- fpkm.mat[,c(1:(ncol(fpkm.mat) - 1 ))]
# colnames(fpkm.mat) <- c('gene_name', coln[2:length(coln)])
#
tmp <- fpkm.mat
tmp[,'rowsum'] <- apply(tmp[,c(2:ncol(fpkm.mat))], 1, sum)
tmp
# any
# any value bigger than 5
# big5 <- c()
# for (i in 1:nrow(tmp)){
#   tmp_vec <- tmp[i, c(2:ncol(fpkm.mat))]
#   big5 <- c(big5, any(tmp_vec > 5))
# }
# tmp[,'bigger5'] <- big5

# tmp_mat <- tmp[tmp$bigger5,]
# 此处需要将一些在各样本中FPKM值均太小的基因去除，还要将一些小数点后位数太多的FPKM值改成整数，不然在cibersort运行过程中会报错
tmp_mat <- tmp[tmp$rowsum > 10,]
tmp_mat[tmp_mat < 0.000001] <- 0
tmp_mat <- tmp_mat[complete.cases(tmp_mat),]
# ph <- '/thinker/storage/udata/muyl/projects/18.k19607.rnaseq.kangliang/cufflinks.fpkm.bigger5.mat'
# write.table(tmp_mat[,c(1:ncol(fpkm.mat))], ph, row.names = FALSE, quote = F)
write.table(tmp_mat[,c(1:(ncol(tmp_mat) - 1))], '/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/fpkm.del.dup.del.zero.for.cibersort.mat'
            , sep = '\t', row.names = F, quote = F)
write.table(tmp_mat[,c(1:(ncol(tmp_mat) - 1))], '/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/fpkm.del.dup.del.zero.for.cibersort.mat.txt'
            , sep = '\t', row.names = F, quote = F)
write.table(tmp_mat[,c(1:(ncol(tmp_mat) - 1))], '/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/01fpkm.del.dup.del.zero.for.cibersort.mat'
            , sep = '\t', quote = F)
# fpkm.mat

need_exon_ensg_df02 <- exon_ensg_df02[,c("gene_id","gene_name")]
need_exon_ensg_df03 <- distinct(need_exon_ensg_df02)
need_exon_ensg_df03
head(tmp_mat)
tmp_mat02 <- tmp_mat
tmp_mat02
tmp_mat02$gene_id <- row.names(tmp_mat)
tmp_mat02
tmp_mat03 <- left_join(tmp_mat02,need_exon_ensg_df03,by="gene_id")
tmp_mat03
# row.names(tmp_mat03) <- tmp_mat03$gene_name
rownames(tmp_mat03) <- tmp_mat03$gene_name
tmp_mat03$`Gene symbol` <- tmp_mat03$gene_name
tmp_mat03
rownames(tmp_mat03) <- tmp_mat03[,c("gene_name")]
tmp_mat03[,c(1:2)]
write.table(tmp_mat03[,c(1:(ncol(tmp_mat) - 1))], '/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/01fpkm.del.dup.del.zero.for.cibersort.mat03.txt'
            , sep = '\t', quote = F,row.names = F)
library(dplyr)
# Gene symbol
new_tmp_mat04 <- tmp_mat03 %>% dplyr::select(`Gene symbol`, everything())
new_tmp_mat04
tmp_mat05 <- new_tmp_mat04[!duplicated(new_tmp_mat04[1]),]

write.table(new_tmp_mat04[,c(1:(ncol(tmp_mat) - 4))], '/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/01fpkm.del.dup.del.zero.for.cibersort04.mat.txt'
            , sep = '\t', quote = F,row.names = F)
write.table(tmp_mat05[,c(1:(ncol(tmp_mat05) - 4))], '/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/05fpkm.del.dup.del.zero.for.cibersort.mat.txt'
            , sep = '\t', quote = F,row.names = F)
save.image('/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/01_fpkm04.RData')



# # 自己通过R源码产生 ---------------------------------------------------------------
# source("/thinker/aid2/udata/all/gao/others/tcga/lgg/CIBERSORT/CIBERSORT.R")
#
# # Define LM22 file
# LM22.file <- "/thinker/aid2/udata/all/gao/others/tcga/lgg/CIBERSORT/LM22.txt"
# temp_lm22 <- read.table(LM22.file,header=T,sep="\t",row.names=1,check.names=F)
# temp_lm22
# exp.file <- "/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/01fpkm.del.dup.del.zero.for.cibersort04.mat.txt"
# exp.file <-"/thinker/aid2/udata/all/gao/others/tcga/lgg/fpkm/05fpkm.del.dup.del.zero.for.cibersort.mat.txt"
#
# TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000)

# 这个代码不行，还是通过网页来生成results -------------------------------------------------
# 本分析基于cibersort的结果进行，cibersort文件位于 /thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersortx/
inF <- read.table('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersortx/CIBERSORTx_Job13_Results.csv'
                  , sep=',', header = T, stringsAsFactors = F)
rownames(inF) <- inF[,1]

inF1 <- inF[,2:23]

rown <- rownames(inF1)
# rown <- str_remove(rown, 'cufflinksOut.aws.')
# rown <- str_remove(rown, 'cufflinksOut.67.')
# rown <- str_remove(rown, 'cufflinksOut.d2.')
# rown <- str_replace(rown, '\\.', '_')

tmp_table <- sampleTable
sampleTable
rownames(tmp_table) <- sampleTable$sampleName

# tmp[,'p_n'] <- str_c(tmp$pateint_id_new, '_', tmp$location, '_', tmp$condition, '_', tmp$lym_meta)
inF2 <- inF1
inF2$condition <- tmp_table[rown, 'condition']
inF2
rownames(inF2) <- rown
# write.csv(inF2, './lncrna_ana/plots/3.cibersort/sample.immune.cell.fraction.matrix.table.with.condition.csv')

frac <- c()
cell <- c()
patient <- c()
condition <- c()
loc_V <- c()
for (i in 1:22){
  for (j in 1:nrow(inF2)){
    frac <- c(frac, inF2[j, i])
    cell <- c(cell, colnames(inF2)[i])
    condition <- c(condition, as.character(inF2$condition)[j])
    patient <- c(patient, as.character(rown)[j])
  }
}

# longdf <- data.frame(fraction=frac, cellType=cell, ps_name=patient, condition=condition, location=loc_V)
longdf <- data.frame(fraction=frac, cellType=cell, condition=condition, ps_name=patient)

pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/fpkm.stacked.barplot.immune.cell.fraction.per.cell.pdf', width = 10)
# ggplot() + geom_bar(aes(y = fraction, x = ps_name, fill = cellType), data = longdf,
#                     stat="identity", width=1) + xlab('') +
#   theme(legend.position = 'bottom', axis.text.x = element_text(angle = 90)) + ggtitle('all.samples.immune.cell.fraction')

ggplot() + geom_bar(aes(y = fraction, x = ps_name, fill = cellType), data = longdf,
                    stat="identity", width=1) + xlab('') +
  theme(legend.position = 'bottom', axis.text.x =  element_blank()) + ggtitle('all.samples.immune.cell.fraction')


dev.off()

# c. cibersort immune cell fraction heatmap


inF3 <- inF2[,c(1:22)]
inMat <- t(inF3)
pheatmap(inMat, cluster_cols = T)

# add annotation
annot_df <- data.frame(condition=inF2$condition)
rownames(annot_df) <- rownames(inF2)
# annot_df <- annot_df[, c('condition')]
# pheatmap(inMat, cluster_cols = T, annotation_col = annot_df, cellwidth = 8)
pheatmap(inMat, cluster_cols = T, annotation_col = annot_df, cellwidth = 12
         , filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.sample.immune.cell.fraction.heatmap.annotated.NTM.wide.pdf')
dev.off()
pheatmap(inMat, cluster_cols = T, annotation_col = annot_df, show_colnames = F
         , filename = '/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.sample.immune.cell.fraction.heatmap.annotated.NTM.pdf')
dev.off()

# 1.6 update information for manuscript RNAseq 20210901 --------
# 1.6.1 cibersort violin plot --------
ciber_long_df <- longdf
ggplot(ciber_long_df, aes(x=cellType, y=fraction, fill=condition)) +
  geom_violin()

pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.violin.plot.multiple.cell.type.two.condition.pdf', width = 30)
ggplot(ciber_long_df, aes(x=cellType, y=fraction, fill=condition)) +
  geom_violin(position=position_dodge(1))
dev.off()


tmp <- ciber_long_df[ciber_long_df$cellType == 'B.cells.naive',]
# ggplot(tmp, aes(x=cellType, y=fraction, fill=condition)) +
#   geom_violin() + geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))

ggplot(tmp, aes(x=condition, y=fraction, fill=condition)) +
  geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(width=0.03)

ctp_v <- as.character(unique(ciber_long_df$cellType))
pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.violin.plot.two.conditon.one.cell.type.per.plot.pdf')
for (i in ctp_v){
  tmp <- ciber_long_df[ciber_long_df$cellType == i, ]
  p <- ggplot(tmp, aes(x=condition, y=fraction, fill=condition)) +
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    geom_boxplot(width=0.03, fill='white') +
    labs(title=paste0("Distribution of ", i,  " cell fraction in normal and tumor samples")) +
    stat_compare_means()
  print(p)
}
dev.off()


# 1.6.2 correlation matrix plot --------
# install.packages('corrplot')
library(corrplot)
inMat
cell_fraction_mat <- t(inMat)
cell_fraction_mat <- as.matrix(cell_fraction_mat)
mydata.cor = cor(cell_fraction_mat)

pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/correlation.matrix.bubble.plot.cibersort.cell.types.pdf')
p <- corrplot(mydata.cor, type = 'upper'
              , title = 'correlation matrix among reletive fraction of immune cell types'
              , mar = c(1,0,3,0))
print(p)
dev.off()

# tcga.lgg.clinic.information.table.csv
# 1.6.4 cibersort consensus cluster and survival analysis --------
# BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
# 2 consensus cluster
inF_tmp <- inF2[,c(1:22)]
inF_tmp01 <- inF_tmp
rown02 <- str_remove(rownames(inF_tmp01),pattern = "-0.+")
rown02
inF_tmp01$bcr_patient_barcode <- str_remove(rownames(inF_tmp01),pattern = "-0.+")
inF_tmp01
rownames(inF_tmp01) <- inF_tmp01$bcr_patient_barcode

inF_tmp02 <- inF_tmp01[!duplicated(inF_tmp01["bcr_patient_barcode"]),]
# inF_tmp02 <- inF_tmp01
rownames(inF_tmp02) <- inF_tmp02$bcr_patient_barcode
inF_tmp03 <-inF_tmp02[,c(1:ncol(inF_tmp02)-1)]

#
# [526] "TCGA-WY-A858-01A-11R-A36H-07" "TCGA-WY-A859-01A-12R-A36H-07" "TCGA-WY-A85A-01A-21R-A36H-07"
# [529] "TCGA-WY-A85B-01A-11R-A36H-07" "TCGA-WY-A85C-01A-11R-A36H-07" "TCGA-WY-A85D-01A-11R-A36H-07"
# [532] "TCGA-WY-A85E-01A-11R-A36H-07"

# bcr_patient_barcode
# 1:        TCGA-FG-7643
# 2:        TCGA-DU-8158
# 3:        TCGA-DU-5854


# inF_tmp <- as.matrix(inF_tmp)
inF_tmp <- as.matrix(inF_tmp03)
# inF_tmp02 <- as.matrix(inF_tmp03)
# 总生存期（O：
surv_info <- fread('/thinker/aid2/udata/all/gao/others/tcga/lgg/tcga.lgg.clinic.information.table.csv')
surv_info
surv_info01 <- surv_info[surv_info$vital_status !="Not Reported",]

surv_info01$os_day <- ifelse(is.na(surv_info01$days_to_death),surv_info01$days_to_last_follow_up ,surv_info01$days_to_death)
surv_info01$OS <- surv_info01$os_day/30
surv_info01$OS
surv_info01$survival <- ifelse(surv_info01$vital_status=="Alive",1,0)
colnames(surv_info)

tmp_surv_info <- as.data.frame(surv_info01)
tmp_surv_info <- tmp_surv_info[, c('bcr_patient_barcode',  'OS', 'survival')]
# tmp_surv_info[, 'sample_code'] <- str_remove(tmp_surv_info$Tumor_Sample_Barcode, '_.+')

# 3.5 survival difference between clusters
results <- ConsensusClusterPlus(t(inF_tmp), maxK = 10, reps = 1000, pItem = 0.8
                                , pFeature = 1, title = '/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx'
                                , clusterAlg = 'km', distance = 'euclidean'
                                , seed=1111111, plot = 'png')


tmpdf1 <- results[[3]]$consensusClass
tmp_group_df <- as.data.frame(tmpdf1)
tmp_group_df[, 'bcr_patient_barcode'] <- rownames(tmp_group_df)
tmp_group_df
# tmp_group_df[, 'sample_code'] <- str_remove(rownames(tmp_group_df), 'P')
# tmp_group_df[, 'sample_code'] <- str_remove(tmp_group_df$sample_code, '_.+')
# tmp_group_df[, 'sample_code'] <- str_remove(tmp_group_df$sample_code, '[N,T]')
# tmp_group_df[, 'Tumor_Sample_Barcode'] <- rownames(tmp_group_df)
# tmp_group_df[, 'condition'] <- 'N'
# tmp_group_df[grepl('T', tmp_group_df$Tumor_Sample_Barcode), 'condition'] <- 'T'
tmp_surv_group_df <- tmp_group_df

tmp_surv_info_df <- left_join(tmp_surv_info, tmp_surv_group_df, by='bcr_patient_barcode')
tmp_surv_info_df <- tmp_surv_info_df[,c('bcr_patient_barcode', 'OS',  'survival', 'tmpdf1')]
colnames(tmp_surv_info_df) <- c('sample_code',  'OS', 'survival', 'cluster_label')
library(survival)
library(survminer)
# install.packages("survminer")
fit <- survfit(Surv(OS, survival) ~ cluster_label, data = tmp_surv_info_df)
pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.consensus.cluster.OS.diff.cluster.knum.3.pdf')
ggsurvplot(fit, pval = T, risk.table = T)
dev.off()

# fit <- survfit(Surv(PFS, reoccurence) ~ cluster_label, data = tmp_surv_info_df)
# pdf('./ana_20210901/cibersort.consensus.cluster.PFS.diff.cluster.knum.3.pdf')
# ggsurvplot(fit, pval = T, risk.table = T)
# dev.off()
# 1.6.6 relationship between immune cell fraction and survival --------
icf_tmp <- inF_tmp
icf_tmp <- as.data.frame(icf_tmp)
icf_tmp
icf_tmp[, 'bcr_patient_barcode'] <- rownames(icf_tmp)
tmp_icf_surv_df <- left_join(tmp_surv_info, icf_tmp, by='bcr_patient_barcode')
tmp_icf_surv_df
tmp_icf_surv_df01 <- tmp_icf_surv_df[complete.cases(tmp_icf_surv_df),]
tmp_icf_surv_df <- tmp_icf_surv_df01
cell_v <- colnames(inF_tmp)
cell_v
# gene_v <- deg_ol_df$ENSEMBL
pvaldf <- list()
pval_df_coln <- c('immune cell type', 'logrank test score', 'df', 'pvalue')
log_score_v <- c()
df_value_v <- c()
pvalue_value_v <- c()

pdf('/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.immune.cell.types.median.group.survival.curve.pdf',width = 10)
for (i in cell_v){
  tmp_pval_v <- c()
  cox_pvalue <- summary(coxph(Surv(tmp_icf_surv_df$OS, tmp_icf_surv_df$survival) ~ tmp_icf_surv_df[, i]
                              , data = tmp_icf_surv_df))
  # pvalL[[i]] <- c(i, cox_pvalue$sctest[1], cox_pvalue$sctest[2], cox_pvalue$sctest[3])
  log_score_v <- c(log_score_v, as.vector(cox_pvalue$sctest)[1])
  df_value_v <- c(df_value_v, as.vector(cox_pvalue$sctest)[2])
  pvalue_value_v <- c(pvalue_value_v, as.vector(cox_pvalue$sctest)[3])

  tmp <- tmp_icf_surv_df[,c(i, 'OS', 'survival')]
  tmp[, 'group'] <- 1
  tmp[tmp[,i] > median(tmp[,i]), 'group'] <- 2

  fit <- survfit(Surv(OS, survival) ~ group, data = tmp)
  pva <- surv_pvalue(fit, data = tmp)
  pva[1,1] <- i
  pvaldf[[i]] <- pva
  # sfn <- paste('./plots/7.cnv.selected.segment.survival/survival.plots.cnv.DFS.updated/', i , '.DFS.curve.cnv.updated.pdf', sep='')
  # pdf(sfn)
  p <- ggsurvplot(fit, pval = T, risk.table = T
                  , title=paste('survival curve of ', i, ' cell grouped by fraction median')
                  , data = tmp)
  print(p)
  # dev.off()
}
dev.off()
p.value.df <- rbindlist(pvaldf)
p.value.df <- as.data.frame(p.value.df)
# p.value.df[, 'cell'] <- p.value.df$variable
# median.group.surv.curve.p.value.df <- left_join(p.value.df, deg_ol_df)
cox.p.value.df <- data.frame(ENSEMBL=cell_v, score=log_score_v, method='log-rand', pval=pvalue_value_v)
# cox.p.value.df.with.gene.name <- left_join(cox.p.value.df, deg_ol_df)
write.table(p.value.df
            , '/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.immune.cell.types.median.group.logrank.test.p.value.csv'
            , quote = F, row.names = F, sep = ',')
write.table(cox.p.value.df
            , '/thinker/aid2/udata/all/gao/others/tcga/lgg/cibersortx/cibersort.immune.cell.types.coxph.logrank.test.p.value.csv'
            , quote = F, row.names = F, sep = ',')

save.image("/thinker/aid2/udata/all/gao/others/tcga/lgg/20221018mianyizuhua.RData")
