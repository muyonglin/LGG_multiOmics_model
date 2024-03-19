rm(list = ls());
# 0. setup environment  ---------------------------------------------------
.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.2.0/lib/R/library","/thinker/storage/udata/sunhr/R/R4.2.0/lib/R/library/", .libPaths()))

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

library(vidger)
library(clusterProfiler)
library(stringr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
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

# mRNA   fpkm
query.star.count <- GDCquery(project = "TCGA-LGG",
                             experimental.strategy = "RNA-Seq",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "STAR - Counts"
                             # , legacy = TRUE
)
# GDCdownload(query.star.count)
# # hg19
# query.star.count19 <- GDCquery(project = "TCGA-LGG",
#                                experimental.strategy = "RNA-Seq",
#                                data.category = "Gene expression",
#                                # workflow.type = "STAR - Counts",
#                                legacy = TRUE
# )
# unique(query.star.count19$results[[1]]$sample_type =="Recurrent Tumor")
# lihc_Rnaseq <- GDCprepare(query, summarizedExperiment = FALSE)
# lihc_Rnaseq <- GDCprepare(query)
# prepare STAR results
load("/thinker/aid2/udata/all/gao/others/tcga/lgg/RNA_star_count.RData")
file_table <- query.star.count$results[[1]]
file_table
file_table[file_table$cases.submitter_id == "TCGA-DU-6404",]

duplicated_name <- file_table %>% group_by(cases.submitter_id) %>% summarise(freq=n()) %>% filter(freq >1)
Recurrent_name <- file_table[file_table$sample_type=="Recurrent Tumor",]$cases.submitter_id
second_need_table <- file_table[file_table$cases.submitter_id %in% Recurrent_name,]
# A：Vial, 在一系列患者组织中的顺序，绝大多数样本该位置编码都是A; 很少数的是B，表示福尔马林固定石蜡包埋组织，已被证明用于测序分析的效果不佳，所以不建议使用-01B的样本数据：
second_need02_table <- second_need_table[order(second_need_table$cases.submitter_id,decreasing = T),]
second_need02_table$standard <- str_remove(second_need02_table$sample.submitter_id,".+\\d")
second_need02_table$standard02 <- str_remove(second_need02_table$sample.submitter_id,".+-")
second_need02_table[second_need02_table$cases.submitter_id=="TCGA-FG-5965",]
second_need03_table <- second_need02_table[(second_need02_table$standard =="A")|(second_need02_table$standard02 =="01B"),]
write.table(second_need03_table,"/thinker/aid2/udata/all/gao/others/tcga/lgg/file_use_info_table.xls",row.names = F, quote = F, sep="\t")
second_need03_table$tumor_type <- str_replace(second_need03_table$sample_type," ", "_")
second_need03_table$star_file_name <- paste0(second_need03_table$cases,"_",second_need03_table$tumor_type,"_star_count")
second_need03_table$star_file_name
write.csv(second_need03_table[c("star_file_name"),],"/thinker/aid2/udata/all/gao/others/tcga/lgg/file_use_info_table.csv")

# annovar.tmb.df[order(annovar.tmb.df$TMB, decreasing = T),]
nrow(second_need_table)
duplicated_second_need_table_name <- second_need_table %>% group_by(cases.submitter_id) %>% summarise(freq=n()) %>% filter(freq >1)

unique(file_table$sample_type)
write.table(file_table,"/thinker/aid2/udata/all/gao/others/tcga/lgg/file_table.xls",row.names = F, quote = F, sep="\t")

# 只
# file_table[file_table$file_id == 'ff12abd3-0f45-4063-afa7-fa5cad973159',]
fn <- paste0('./GDCdata/TCGA-LGG/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/'
             , file_table[i, 'file_id'], '/', file_table[i, 'file_name'])
tmp_rc <- fread(fn)
tmp_rc <- tmp_rc[5:nrow(tmp_rc),]
tmp_rc
lihc_star_count <- data.frame(gene_name=tmp_rc$gene_name, stringsAsFactors = FALSE)
head(lihc_star_count)

for (i in 1:nrow(file_table)){
  fn <- paste0('./GDCdata/TCGA-LGG/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/'
               , file_table[i, 'file_id'], '/', file_table[i, 'file_name'])
  tmp_rc <- fread(fn)
  tmp_rc <- tmp_rc[5:nrow(tmp_rc),]
  lihc_star_count[, file_table[i, 'cases']] <- tmp_rc$unstranded
}

for (i in 1:nrow(file_table)){
  fn <- paste0('./GDCdata/TCGA-LGG/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/'
               , file_table[i, 'file_id'], '/', file_table[i, 'file_name'])
  tmp_rc <- fread(fn)
  tmp_rc <- tmp_rc[5:nrow(tmp_rc),]
  to_tmp_rc <- tmp_rc[,c("gene_id","unstranded")]
  type = str_replace(file_table[i, 'sample_type']," ", "_")
  file_case_name=file_table[i, 'cases']
  write.table(to_tmp_rc,file.path("/thinker/aid2/udata/all/gao/others/tcga/lgg/star_count/",stringr::str_interp("${file_case_name}_${type}_star_count")),
            row.names = F,col.names = F,quote = F, sep = "\t")
  # lihc_star_count[, file_table[i, 'cases']] <- tmp_rc$unstranded
}

for (i in 1:nrow(second_need03_table)){
  fn <- paste0('./GDCdata/TCGA-LGG/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/'
               , second_need03_table[i, 'file_id'], '/', second_need03_table[i, 'file_name'])
  tmp_rc <- fread(fn)
  tmp_rc <- tmp_rc[5:nrow(tmp_rc),]
  to_tmp_rc <- tmp_rc[,c("gene_id","unstranded")]
  type = str_replace(second_need03_table[i, 'sample_type']," ", "_")
  file_case_name=second_need03_table[i, 'cases']
  write.table(to_tmp_rc,file.path("/thinker/aid2/udata/all/gao/others/tcga/lgg/star_count/",stringr::str_interp("${file_case_name}_${type}_star_count")),
              row.names = F,col.names = F,quote = F, sep = "\t")
  # lihc_star_count[, file_table[i, 'cases']] <- tmp_rc$unstranded
}

# lihc_rna <- lihc_star_count
# lihc_rna
# lihc_rna
# # rownames(lihc_rna) <- lihc_rna$gene_name
# # 9.1.2.2 model deg  --------
# counts_df_lihc <- lihc_rna
# tmp_name_lihc_rna <- data.frame(index_num=rownames(counts_df_lihc)
#                                 , gene=counts_df_lihc$gene_name)
# # rownames(counts_df_lihc) <- counts_df_lihc$gene_name
# counts_df_lihc <- counts_df_lihc[,2:ncol(counts_df_lihc)]
# # counts_df_lihc <- t(counts_df_lihc)
# # colnames(counts_df_lihc) <- lihc_rna$gene_name
# counts_df_lihc <- as.matrix(counts_df_lihc)
#
# # vsd_lihc_rna <- varianceStabilizingTransformation(counts_df_lihc)
# vsd_lihc_rna <- vst(counts_df_lihc)
#
# # deg_v_tmp <- param_df[grepl('ENSG', param_df$gene), 'gene_name']
#
# # deg_v_tmp <- deg_gene_v
# # deg_v_tmp %in% tmp_name_lihc_rna$gene
# #
# # deg_index_tmp <- tmp_name_lihc_rna[tmp_name_lihc_rna$gene %in% deg_v_tmp, 'index_num']
# #
# # vsd_lihc_rna_deg <- vsd_lihc_rna[deg_index_tmp,]
#
#
# # calculate fold change
# nrow(query.star.count$results[[1]][query.star.count$results[[1]]$sample_type =="Primary Tumor",])
# nrow(query.star.count$results[[1]][query.star.count$results[[1]]$sample_type =="Recurrent Tumor",])
#
#
# deg_file_info <- query.star.count$results[[1]]
# unique(deg_file_info$sample_type)
#
# normal_col <- deg_file_info[deg_file_info$sample_type == "Primary Tumor", 'cases']
# tumor_col <- deg_file_info[deg_file_info$sample_type != "Recurrent Tumor", 'cases']
#
#
# vsd_normal_df <- vsd_lihc_rna_deg[, normal_col]
# rownames(vsd_normal_df) <- deg_v_tmp
# tmp <- apply(vsd_normal_df, 1, mean)
# head(tmp)
#
# deg_lihc_mean_Recurrent_Tumor <- data.frame(gene_id=names(tmp), mean_vst=tmp)

# 开始RNA的分析：
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
table(sampleTable$condition)
ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable
                                       , directory = './'
                                       , design = ~condition)
register(MulticoreParam(20))
keep <- rowSums(counts(ddsHTseq)) >= min(1 * ncol(ddsHTseq), 10)
dds <- ddsHTseq[keep, ]
dds <- DESeq(dds, parallel = T)
dds

# 1.3.0.0 pca analysis for all samples --------
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"))
# pcaData <- plotPCA(vst_sample)
tmpdf <- pcaData$data
# tmpdf[, 'barcode'] <- str_remove(basename(as.character(tmpdf$name)), '_combined_R1.f.htseqCount')
# tmpdf[,'Tumor_Sample_Barcode'] <- rename_df[tmpdf$barcode, 'Tumor_Sample_Barcode']
# p <- ggplot(tmpdf, aes(PC1, PC2, colour = condition)) + geom_point(size=3) + geom_text_repel(label=tmpdf$name)
p <- ggplot(tmpdf, aes(PC1, PC2, colour = condition)) + geom_point(size=3)
dir.create('plots/rnaseq/1.degs/', recursive = TRUE)
pdf('./plots/rnaseq/1.degs/pca.plot.all.samples.all.starOut.pdf')
print(p)
dev.off()


# 1.3.0.1 boxplot for FPM distribution unber different conditions --------
# manually modified
.getDeseqBox <- function(data, d.factor) {
  if(is.null(d.factor)) {
    stop(
      'This appears to be a DESeq object.
            Please state d.factor variable.'
    )
  }
  dat1 <- as.data.frame(colData(data))
  dat2 <- fpm(data)
  nam <- as.vector(unique(dat1[[d.factor]]))
  ls.nam <- list()
  ls.mean <- list()
  for (i in nam) {
    ls.nam[[i]] <- row.names(dat1[which(dat1[d.factor] == i), ])
    for (j in seq_along(ls.nam)) {
      ls.mean[[j]] <- rowMeans(dat2[, ls.nam[[j]]])
    }
  }
  names(ls.mean) <- sapply(nam, paste)
  dat3 <- as.data.frame(ls.mean)
  dat3 <- tidyr::gather(as.data.frame(dat3))
  dat3$key <- as.factor(dat3$key)
  return(dat3)
}

dat <- .getDeseqBox(dds, d.factor='condition')

pdf('./plots/rnaseq/1.degs/jittered.boxplot.all.gene.fpm.distribution.pdf')
ggplot(
  dat, aes(x = key, y = log10(value + 1), fill = key)
) + geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16,
                                                   color = "grey14",
                                                   alpha = 0.7,
                                                   size = 0.1,
                                                   position = position_jitter(0.375)) +
  xlab("Condition") +
  ylab(paste("log","10", ' (FPM)')) +
  guides(fill = guide_legend(title = "Condition")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("FPM distribution") + theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))
dev.off()


# vidger original
pdf('./plots/rnaseq/1.degs/boxplot.all.gene.fpm.distribution.pdf')
vsBoxPlot(
  data = dds, d.factor = 'condition', type = 'deseq',
  title = TRUE, legend = TRUE, grid = TRUE, aes = 'box'
)
dev.off()

# 1.3.1 DEGs between Primary_Tumor and Recurrent_Tumor --------
res <- results(dds, contrast=c('condition', 'Primary_Tumor', 'Recurrent_Tumor'), parallel=T)
resOrdered <- res[order(res$padj),]
resCut <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)##新增
#resCut <- subset(resOrdered, padj < 0.001)
all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05 <- resCut##新增
#all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.001 <- resCut
N.vs.T.table.for.volcano <- resOrdered


# 1.3.1.1 volcano plot for all degs --------
#p.adj=0.05
pdf('./plots/rnaseq/1.degs/vidger.volcano.plot.deseq2.pdf')
vsVolcano(
  x = 'Primary_Tumor', y = 'Recurrent_Tumor',
  data = dds, d.factor = 'condition', type = 'deseq',
  padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE,
  legend = TRUE, grid = TRUE, data.return = FALSE)
dev.off()

# manual volcano plot
n_t_volcano <- as.data.frame(N.vs.T.table.for.volcano)
n_t_volcano$log2FoldChange <- as.numeric(as.character(n_t_volcano$log2FoldChange))
n_t_volcano$pvalue <- as.numeric(as.character(n_t_volcano$pvalue))
n_t_volcano$padj <- as.numeric(as.character(n_t_volcano$padj))

pdf('./plots/rnaseq/1.degs/ggplot2.Primary_Tumor.vs.Recurrent_Tumor.deg.volcano.plot.pdf')
volcano<-subset(n_t_volcano,select = c(pvalue,padj,log2FoldChange))
#volcano<-subset(n_t_volcano,select = c(pvalue,log2FoldChange))
volcano <- volcano[complete.cases(volcano),]
threshold<-as.factor((volcano$log2FoldChange>1|volcano$log2FoldChange<(-1))&volcano$padj<0.05)
#threshold<-as.factor((volcano$log2FoldChange>1|volcano$log2FoldChange<(-1))&volcano$pvalue<0.01)
r03=ggplot(volcano,aes(log2FoldChange,-log2(pvalue),colour=threshold))+geom_point(size=0.5)
r04=r03+labs(title="Volcanoplot of Primary_Tumor VS Recurrent_Tumor DEGs")+theme(plot.title = element_text(hjust = 0.5))+xlim(-6,22)
r05=r04+geom_vline(xintercept=c(-1,1),col='blue',size=0.5)+geom_hline(yintercept=-log2(0.01),col="blue")
r05
dev.off()

tmp <- volcano
tmp[, 'threshold'] <- threshold
tmp1 <- tmp[tmp$threshold == TRUE, ]
print(str_c('N vs T deg number is: ', as.character(nrow(tmp1))))
print(str_c('N vs T up regulated deg number: ', as.character(nrow(tmp1[tmp1$log2FoldChange > 0, ]))))
print(str_c('N vs T down regulated deg number: ', as.character(nrow(tmp1[tmp1$log2FoldChange < 0, ]))))
write.csv(tmp1[tmp1$log2FoldChange > 0, ], './plots/rnaseq/1.degs/Primary_Tumor VS Recurrent_Tumor.deg.volcano.plot.up.regulated.csv')
write.csv(tmp1[tmp1$log2FoldChange < 0, ], './plots/rnaseq/1.degs/Primary_Tumor VS Recurrent_Tumor.deg.volcano.plot.down.regulated.csv')

# 1.3.2 deg expression heatmap ---------
vsd <- vst(dds)
countMat <- assay(vsd)

degMat <- countMat[rownames(all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05),]##新增
#degMat <- countMat[rownames(all.deg.control.VS.T1.pcut.0.001),]

tmp_table <- sampleTable
rownames(tmp_table) <- sampleTable$sampleName
tmp_table <- tmp_table[colnames(degMat),]
annot_df <- data.frame(condition=tmp_table$condition)
rownames(annot_df) <- colnames(degMat)
##show top 100 gene ordered by p.adj values.
degMat_top100<-degMat[1:min(nrow(degMat),100),]
# rownames(annot_df) <- rownames(inF2)
# pheatmap(degMat, show_rownames = F, border_color = NA, show_colnames = F, annotation_col = annot_df)
pheatmap(degMat_top100, show_rownames = F, border_color = NA, show_colnames = F, annotation_col = annot_df
         , main = 'all DEGs between Recurrent_Tumor and Primary_Tumor'
         , filename = './plots/rnaseq/1.degs/deg.all.unique.heatmap.normalized.count.mat.padj.cut.0.05.pdf')

## 1.3.3 deg function enrichment --------
gene_V <- c(rownames(all.deg.Recurrent_Tumor.VS.Primary_Tumor.pcut.0.05))
#gene_V <- c(rownames(all.deg.control.VS.T1.pcut.0.001))
gene_V <- unique(gene_V)
gene_V <- as.character(gene_V)
gene_V <- str_remove(gene_V, '\\.\\d+')
oL <- bitr(gene_V, fromType=toupper('ENSEMBL'), toType='ENTREZID', OrgDb = org.Hs.eg.db)
#oL <- bitr(gene_V, fromType=toupper('ENSEMBL'), toType='UNIPROT', OrgDb = org.Hs.eg.db)
inV <- oL$ENTREZID
#inV <- oL$UNIPROT

R.utils::setOption( "clusterProfiler.download.method",'wget')
ken <- enrichKEGG(inV,organism = 'hsa', keyType = 'kegg', pAdjustMethod='none', pvalueCutoff=1,qvalueCutoff = 1)
goCC <- enrichGO(inV,org.Hs.eg.db, ont='CC', pAdjustMethod='none', pvalueCutoff=1,qvalueCutoff = 1)
goBP <- enrichGO(inV,org.Hs.eg.db, ont='BP', pAdjustMethod='none', pvalueCutoff=1,qvalueCutoff = 1)
goMF <- enrichGO(inV,org.Hs.eg.db, ont='MF', pAdjustMethod='none', pvalueCutoff=1,qvalueCutoff = 1)

ken <- setReadable(ken, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
goCC <- setReadable(goCC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
goBP <- setReadable(goBP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
goMF <- setReadable(goMF, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

##Barplot grouped by GOterms.
pdf('./plots/rnaseq/1.degs/deg.Recurrent_Tumor_and_Primary_Tumor.GO.enrichment.barplot.pdf',height = 12)
go_all<-rbind(goBP[1:10,] ,goCC[1:10,], goMF[1:10,])
go_all$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
go_all$Description <- factor(go_all$Description,levels = rev(go_all$Description))
go_all$group<-factor(go_all$group,levels=c("BP","CC","MF"))
ggplot(data = go_all, aes(x = Description, y = -log10(pvalue), fill = group))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "-log10(p-value)",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11))# 图例标签大小
dev.off()

##Cnetplot of KEGG pathway
pdf('./plots/rnaseq/1.degs/deg.Recurrent_Tumor_and_Primary_Tumor.KEGG.enrichment.cnetplot.pdf',width = 10)
cnetplot(ken, circular = TRUE, colorEdge = TRUE,)
dev.off()

##Dotplot of KEGG and GO terms.
pdf('./plots/rnaseq/1.degs/deg.Recurrent_Tumor_and_Primary_Tumor.GO.KEGG.enrichment.dotplot.pdf'
    , width = 10)
p1 <- dotplot(ken,title='KEGG enrich of degs')
p2 <- dotplot(goCC,title='GO CC enrich of degs')
p3 <- dotplot(goBP,title='GO BP enrich of degs')
p4 <- dotplot(goMF,title='GO MF enrich of degs')
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# 1.3.4 deg output --------
# prepare ensg data frame
# gtf <- '/mnt/phoenix/bio-web/pipeline/rnaseq/reference/reference/gencode.v34.annotation.gtf'
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

tmp[,"chr"]<-tmp$V1
tmp[,'start']<-tmp$V4
tmp[,'end']<-tmp$V5

# tmp[,'gene_type'] <- str_extract(tmp$V9, 'transcript_type ".*?"; transcript_status')
tmp[,'gene_type'] <- str_extract(tmp$V9, 'gene_type ".*?"; gene_name')
tmp$gene_type <- str_remove(tmp$gene_type, 'gene_type "')
tmp$gene_type <- str_remove(tmp$gene_type, '"; gene_name')

tmp[,'gene_id_trim'] <- str_remove(tmp$gene_id, '\\.\\d+')

#gene_ensg_df <- tmp[,c('gene_name', 'gene_id', 'gene_type', 'gene_id_trim')]
gene_ensg_df <- tmp[,c('chr','gene_name','start','end', 'gene_id', 'gene_type', 'gene_id_trim')]
gene_ensg_df <- unique(gene_ensg_df)


tmp <- volcano
tmp[, 'threshold'] <- threshold
tmp1 <- tmp[tmp$threshold == TRUE, ]
print(str_c('N vs T1 deg number is: ', as.character(nrow(tmp1))))
print(str_c('N vs T1 up regulated deg number: ', as.character(nrow(tmp1[tmp1$log2FoldChange > 0, ]))))
print(str_c('N vs T1 down regulated deg number: ', as.character(nrow(tmp1[tmp1$log2FoldChange < 0, ]))))

tmp_gene_ensg_df <- gene_ensg_df[gene_ensg_df$gene_id %in% rownames(tmp1),]
# rownames(tmp_gene_ensg_df) <- tmp_gene_ensg_df$gene_id

tmp1[, 'gene_id'] <- rownames(tmp1)
tmp2 <- left_join(tmp_gene_ensg_df,tmp1)
#tmp2 <- left_join(tmp1, tmp_gene_ensg_df)

##circle plot of DEGs
# Human chromosome data(hg38)
data(UCSC.HG38.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(cyto.info = UCSC.HG38.Human.CytoBandIdeogram,
                            chr.exclude = NULL, tracks.inside = 10, tracks.outside = 0)

pdf(file="./plots/rnaseq/1.degs/deg.Recurrent_Tumor_and_Primary_Tumor.circles.pdf", height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
if(length(which(tmp2$chr=="chrM"))>0){
  circosdata<-tmp2[-which(tmp2$chr=="chrM"),]
}else{circosdata <- tmp2}
circosdata_pro<-subset(circosdata, gene_type=="protein_coding") %>%
  .[order(.$padj),] %>% .[1:(nrow(.)*0.3),c("chr","start","end","gene_name")]
RCircos.Gene.Connector.Plot(circosdata_pro, track.num = 1, side = "in")

RCircos.Gene.Name.Plot(circosdata_pro, name.col = 4,track.num = 2, side = "in");

circosdata_pro<-subset(circosdata, gene_type=="protein_coding") %>%
  .[order(.$padj),c("chr","start","end","gene_name","padj")]
circosdata_pro$padj<-(-log10(circosdata_pro$padj))
RCircos.Heatmap.Plot(circosdata_pro, data.col =5 , track.num = 5 , side = "in");

circosdata_pro<-subset(circosdata, gene_type=="protein_coding") %>%
  .[order(.$padj),c("chr","start","end","gene_name","log2FoldChange")]
RCircos.Scatter.Plot(circosdata_pro, data.col = 5,track.num = 6 , side = "in", by.fold = 1);
dev.off()

write.csv(tmp2[tmp2$log2FoldChange > 0, ], './plots/rnaseq/1.degs/Recurrent_Tumor_and_Primary_Tumor.deg.volcano.plot.up.regulated.with.gene.name.csv', row.names = F)
write.csv(tmp2[tmp2$log2FoldChange < 0, ], './plots/rnaseq/1.degs/Recurrent_Tumor_and_Primary_Tumor.deg.volcano.plot.down.regulated.with.gene.name.csv', row.names = F)

#1.3.5 all gene tpm matrix output
tmp <- inT[inT$V3 == 'exon',]


tmp[,"chr"]<-tmp$V1
tmp[,'start']<-tmp$V4
tmp[,'end']<-tmp$V5
tmp[,'gene_id'] <- str_extract(tmp$V9, 'ENSG\\d+\\.\\d+')
tmp[,'gene_type'] <- str_extract(tmp$V9, 'gene_type ".*?"; gene_name')
tmp$gene_type <- str_remove(tmp$gene_type, 'gene_type "')
tmp$gene_type <- str_remove(tmp$gene_type, '"; gene_name')

tmp[,'gene_id_trim'] <- str_remove(tmp$gene_id, '\\.\\d+')

exon_ensg_df <- tmp[,c('chr','start','end', 'gene_id', 'gene_type', 'gene_id_trim')]

# calculate tpm value
g_l = lapply(split(exon_ensg_df,exon_ensg_df$gene_id),function(x){
  # x=split(t1,t1$geneid)[[1]]
  tmp=apply(x,1,function(y){
    y[2]:y[3]
  })
  length(unique(unlist(tmp)))
  # sum(x[,4])
})
g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))

dds_count<-assay(dds)
ng=intersect(rownames(dds_count),g_l$gene_id)

exprSet=dds_count[ng,]
lengths=g_l[match(ng,g_l$gene_id),2]
#head(lengths)
#head(rownames(exprSet))
#exprSet[1:4,1:4]
total_count<- colSums(exprSet)
#head(total_count)
#head(lengths)
#total_count[1]
#lengths[2]
#53*10^9/(1476*3978357)
fpkm <- t(do.call( rbind,
                   lapply(1:length(total_count),
                          function(i){
                            10^9*exprSet[,i]/lengths/total_count[i]
                          }) ))
colnames(fpkm)<-colnames(dds_count)
total_fpkm<- colSums(fpkm)
tpm<- t(do.call( rbind,
                 lapply(1:length(total_fpkm),
                        function(i){
                          10^6*fpkm[,i]/total_fpkm[i]
                        }) ))
colnames(tpm)<-colnames(dds_count)
#total_fpkm[1]
#fpkm[2,1]*10^6/total_fpkm[1]
fpkm_log<-log2(fpkm+1)
tpm_log<-log2(tpm+1)

.getTPM <- function(data, d.factor) {
  if(is.null(d.factor)) {
    stop(
      'This appears to be a DESeq object.
            Please state d.factor variable.'
    )
  }
  dat1 <- as.data.frame(colData(dds))
  dat2 <- data
  nam <- as.vector(unique(dat1[[d.factor]]))
  ls.nam <- list()
  ls.mean <- list()
  for (i in nam) {
    ls.nam[[i]] <- row.names(dat1[which(dat1[d.factor] == i), ])
    for (j in seq_along(ls.nam)) {
      ls.mean[[j]] <- rowMeans(dat2[, ls.nam[[j]]])
    }
  }
  names(ls.mean) <- sapply(nam, paste)
  dat3 <- as.data.frame(ls.mean)
  dat3 <- tidyr::gather(as.data.frame(dat3))
  dat3$key <- as.factor(dat3$key)
  return(dat3)
}
dat_tpm <- .getTPM(tpm_log, d.factor='condition')

#1.3.5.1 density plot for TPM distribution under different condition
pdf('./plots/rnaseq/1.degs/density.all.gene.tpm.distribution.pdf')
ggplot(dat_tpm,aes(value,fill=as.factor(key)))+
  geom_density(alpha = 0.3)+theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank())+
  labs(fill = "conditions")##表达值密度曲线
dev.off()

#1.3.5.2 boxplot for TPM distribution under different condition
pdf('./plots/rnaseq/1.degs/boxplot.all.gene.tpm.distribution.pdf')
ggplot(dat_tpm, aes(x = key, y = value, fill = key)
) + geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16,
                                                   color = "grey14",
                                                   alpha = 0.7,
                                                   size = 0.1,
                                                   position = position_jitter(0.375)) +
  xlab("Condition") +
  ylab(paste("log","2", ' (TPM)')) +
  guides(fill = guide_legend(title = "Condition")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("TPM distribution") + theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()
write.csv(tpm_log, './plots/rnaseq/1.degs/all.gene.TPM.matrix.csv', row.names = T)

#1.3.5.3 Sample correlation matrix diagram
pdf('./plots/rnaseq/1.degs/corrplot.all.sample04.pdf')
tmp_cor<-cor(tpm_log)
sample_cor<-round(tmp_cor,digits = 2)
annotation_row=data.frame(sample=str_remove(rownames(sample_cor),".+_"))##设置注释分组
rownames(annotation_row)<-rownames(sample_cor)##注释起名
annotation_col=data.frame(sample=str_remove(colnames(sample_cor),".+_"))
rownames(annotation_col)<-colnames(sample_cor)
pheatmap(sample_cor,border_color=F,annotation_row=annotation_row,annotation_col = annotation_col)
pheatmap(sample_cor,border_color=F,annotation_row = NA, annotation_col = NA,annotation_names_row = F,
         annotation_names_col = T,cluster_rows = T,
         cluster_cols = T,show_rownames = F, show_colnames = F)
dev.off()
getwd()
pdf('./plots/rnaseq/1.degs/corrplot.all.sample05.pdf',width = 25,height = 25)
pheatmap(sample_cor,border_color=F,annotation_row=annotation_row,annotation_col = annotation_col)
dev.off()
#1.3.6.1 GSEA analysis
##Hallmark
hallmark_df<-  msigdbr(species = "Homo sapiens",
                       category = "H") %>% dplyr::select(.,gs_name, gene_symbol) %>% as.data.frame
##KEGG
kegg_df<-  msigdbr(species = "Homo sapiens",
                   category = "C2",subcategory = "CP:KEGG") %>%
  dplyr::select(.,gs_name, gene_symbol) %>% as.data.frame

##GO
go_df<-  msigdbr(species = "Homo sapiens",
                 category = "C5") %>%
  dplyr::select(.,gs_name, gene_symbol,gs_subcat) %>% .[.$gs_subcat=="GO:BP",] %>%
  as.data.frame
##REACTOME
reactome_df<-  msigdbr(species = "Homo sapiens",
                       category = "C2",subcategory = "CP:REACTOME") %>%
  dplyr::select(.,gs_name, gene_symbol) %>% as.data.frame

n_t_gsea<-n_t_volcano[order(n_t_volcano$log2FoldChange,decreasing = T),]
gene_G<-n_t_gsea$log2FoldChange
names(gene_G)<-str_remove(rownames(n_t_gsea),'\\.\\d+')
gene_G1 <- bitr(names(gene_G), fromType=toupper('ENSEMBL'), toType='SYMBOL', OrgDb = org.Hs.eg.db)
gene_GSEA<-data.frame(log2fc = gene_G)
gene_GSEA$ensmbl<-rownames(gene_GSEA)
gene_GSEA1<-merge(gene_G1,gene_GSEA,by.x="ENSEMBL",by.y="ensmbl") %>% .[order(.$log2fc,decreasing = T),]
inG <- gene_GSEA1$log2fc
names(inG) <- gene_GSEA1$SYMBOL##ordered.by.log2FC.gene.list

gsea_hallmark<-GSEA(inG,TERM2GENE = hallmark_df, pvalueCutoff = 1,pAdjustMethod = "none")
gsea_kegg<-GSEA(inG,TERM2GENE = kegg_df, pvalueCutoff = 1,pAdjustMethod = "none")
gsea_gobp<-GSEA(inG,TERM2GENE = go_df, pvalueCutoff = 1,pAdjustMethod = "none")
gsea_reactome<-GSEA(inG,TERM2GENE = reactome_df, pvalueCutoff = 1,pAdjustMethod = "none")
pdf('./plots/rnaseq/1.degs/GSEA.all.gene.GO.KEGG.enrichment.pdf')
p1 <- gseaplot2(gsea_hallmark,geneSetID=c(1:5),pvalue_table=F,title = "GSEA_HALLMARK")
p2 <- gseaplot2(gsea_kegg,geneSetID=c(1:5),pvalue_table=F,title = "GSEA_KEGG")
p3 <- gseaplot2(gsea_gobp,geneSetID=c(1:5),pvalue_table=F,title = "GSEA_GOBP")
p4 <- gseaplot2(gsea_reactome,geneSetID=c(1:5),pvalue_table=F,title = "GSEA_REACTOME")
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

#1.3.6.2 GSVA analysis
##Hallmark
hallmark_df_all<-  msigdbr(species = "Homo sapiens",
                           category = "H") %>% dplyr::select(.,gs_name, gs_exact_source, gene_symbol)
hallmark_list <- split(hallmark_df_all$gene_symbol, hallmark_df_all$gs_name)

##KEGG
kegg_df_all<-  msigdbr(species = "Homo sapiens",
                       category = "C2",subcategory = "CP:KEGG") %>%
  dplyr::select(.,gs_name, gs_exact_source, gene_symbol)
kegg_list <- split(kegg_df_all$gene_symbol, kegg_df_all$gs_name)

##GO
go_df_all<-  msigdbr(species = "Homo sapiens",
                     category = "C5") %>%
  dplyr::select(.,gs_name, gs_exact_source, gene_symbol,gs_subcat) %>% .[.$gs_subcat=="GO:BP",]
go_list <- split(go_df_all$gene_symbol, go_df_all$gs_name)

save.image("./RNA_lgg_tcga_star20220813.RData")
##REACTOME
reactome_df_all<-  msigdbr(species = "Homo sapiens",
                           category = "C2",subcategory = "CP:REACTOME") %>%
  dplyr::select(.,gs_name, gs_exact_source, gene_symbol)
reactome_list <- split(reactome_df_all$gene_symbol, reactome_df_all$gs_name)

tpm_log<-as.data.frame(tpm_log)
tpm_log$name<-str_remove(rownames(tpm_log),'\\.\\d+')
gene_df <- bitr(tpm_log$name, fromType=toupper('ENSEMBL'), toType='SYMBOL', OrgDb = org.Hs.eg.db)
symbol_tpm_log<-merge(tpm_log,gene_df,by.x="name",by.y="ENSEMBL") %>%
  .[!.$SYMBOL%in%names(which(table(.$SYMBOL)>=2)),]
rownames(symbol_tpm_log)<-symbol_tpm_log$SYMBOL
symbol_tpm_log<-symbol_tpm_log[,-c(which(colnames(symbol_tpm_log)=="name" |
                                           colnames(symbol_tpm_log)=="SYMBOL"))]
symbol_tpm_log<-as.matrix(symbol_tpm_log)

gsva_hallmark<-gsva(symbol_tpm_log, gset.idx.list=hallmark_list,
                    method = "gsva",kcdf="Gaussian", verbose=TRUE)
gsva_kegg<-gsva(symbol_tpm_log, gset.idx.list=kegg_list,
                method = "gsva",kcdf="Gaussian", verbose=TRUE)
gsva_go<-gsva(symbol_tpm_log, gset.idx.list=go_list,
              method = "gsva",kcdf="Gaussian", verbose=TRUE)
gsva_reactome<-gsva(symbol_tpm_log, gset.idx.list=reactome_list,
                    method = "gsva",kcdf="Gaussian", verbose=TRUE)

##N.vs T1.limma.gsva.pathway.differential.analysis
.gsva_diff<-function(gsva_data,control,tumor){
  tmp_l<-gsva_data[,c(grep(paste("*_",control,sep=""),colnames(gsva_hallmark)),
                      grep(paste("*_",tumor,sep=""),colnames(gsva_hallmark)))]
  group_list<-str_remove(colnames(tmp_l),".+_")
  design <- model.matrix(~0+factor(group_list))
  colnames(design) <- levels(factor(group_list))
  rownames(design) <- colnames(tmp_l)
  contrast.matrix <- makeContrasts(contrasts=paste0(control,'-',tumor),
                                   levels = design)
  fit1 <- lmFit(tmp_l,design)                 #拟合模型
  fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
  efit <- eBayes(fit2)
  tempOutput <- topTable(efit, coef=paste0(control,'-',tumor), n=Inf)
  degs <- na.omit(tempOutput)
  return(degs)
}

gsva_hallmark_diff<-.gsva_diff(gsva_hallmark,control = "Primary_Tumor",tumor = "Recurrent_Tumor")
gsva_kegg_diff<-.gsva_diff(gsva_kegg,control = "control",tumor = "T1")
gsva_go_diff<-.gsva_diff(gsva_go,control = "control",tumor = "T1")
gsva_reactome_diff<-.gsva_diff(gsva_reactome,control = "control",tumor = "T1")

gsva_hallmark_display<-gsva_hallmark_diff[1:20,c("logFC","P.Value")]
gsva_hallmark_display$name <- rownames(gsva_hallmark_display)
gsva_kegg_display<-gsva_kegg_diff[1:20,c("logFC","P.Value")]
gsva_kegg_display$name <- rownames(gsva_kegg_display)
gsva_go_display<-gsva_go_diff[1:20,c("logFC","P.Value")]
gsva_go_display$name <- rownames(gsva_go_display)
gsva_reactome_display<-gsva_reactome_diff[1:20,c("logFC","P.Value")]
gsva_reactome_display$name <- rownames(gsva_reactome_display)

pdf('./plots/rnaseq/1.degs/GSVA.N.vs.T.differential.pathways.barplots.pdf',width = 12)
p1<-ggplot(gsva_hallmark_display, aes(x= reorder(name,logFC), y=logFC,fill = -log10(P.Value))) +
  geom_bar(stat='identity', width=.7)  + scale_fill_gradient(low = "blue",high = "red")+
  labs(x="HALLMARK", y= "log2FC", title= "Differntial HALLMARK top 20")+
  theme_bw()+coord_flip()
p2<-ggplot(gsva_kegg_display, aes(x= reorder(name,logFC), y=logFC,fill = -log10(P.Value))) +
  geom_bar(stat='identity', width=.7)  + scale_fill_gradient(low = "blue",high = "red")+
  labs(x="KEGG", y= "log2FC", title= "Differntial KEGG top 20")+
  theme_bw()+coord_flip()
p3<-ggplot(gsva_go_display, aes(x= reorder(name,logFC), y=logFC,fill = -log10(P.Value))) +
  geom_bar(stat='identity', width=.7)  + scale_fill_gradient(low = "blue",high = "red")+
  labs(x="GO BP", y= "log2FC", title= "Differntial GO BP top 20")+
  theme_bw()+coord_flip()
p4<-ggplot(gsva_reactome_display, aes(x= reorder(name,logFC), y=logFC,str_wrap(name, width=2),fill = -log10(P.Value))) +
  geom_bar(stat='identity', width=.7)  + scale_fill_gradient(low = "blue",high = "red")+
  labs(x="REACTOME", y= "log2FC", title= "Differntial REACTOME top 20")+
  theme_bw()+coord_flip()+theme(axis.text.y = element_text(size=4))
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

#1.3.7 WGCNA analysis
.wgcna<-function(wgcna_data,control,tumor,gene_count){
  ## step 1 :
  tmp_l<-wgcna_data[,c(grep(paste("*_",control,sep=""),colnames(tpm_log)),
                       grep(paste("*_",tumor,sep=""),colnames(tpm_log)))]
  WGCNA_matrix = t(tmp_l[order(apply(tmp_l,1,mad), decreasing = T)[1:gene_count],])
  datExpr <- WGCNA_matrix
  ### step 2 :optimal beta value
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  cex1 = 0.9;
  pdf("./plots/rnaseq/1.degs/WGCNA.soft-thresholding power.Scatterplot.pdf")
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  if(is.na(sft$powerEstimate) == FALSE){
    #step3 Weight co-expression network
    net = blockwiseModules(
      datExpr,
      power = sft$powerEstimate,
      maxBlockSize = gene_count,
      TOMType = "unsigned", minModuleSize = 30,
      reassignThreshold = 0, mergeCutHeight = 0.25,
      numericLabels = TRUE, pamRespectsDendro = FALSE,
      saveTOMs = F,
      verbose = 3
    )
    mergedColors = labels2colors(net$colors)
    dynamicColors = labels2colors(net$unmergedColors)
    moduleColors=mergedColors
    # step4-Plot the dendrogram and the module colors underneath
    pdf("./plots/rnaseq/1.degs/WGCNA.dynamicColors+mergedColors.plot.pdf")
    par(mar = c(6, 8.5, 3, 3))
    plotDendroAndColors(net$dendrograms[[1]], cbind(dynamicColors[net$blockGenes[[1]]],mergedColors[net$blockGenes[[1]]]),
                        c("dynamic colors","merged colors"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()

    # step 5-Relationships between modules and traits
    datTraits<-data.frame(condition = str_remove(colnames(tmp_l),".+_"))
    datTraits$group<-rep(1,ncol(datTraits))
    datTraits$group[which(datTraits$condition == control)]<-0
    ##datTraits$nt<-datTraits$group
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    moduleColors <- labels2colors(net$colors)
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0);
    moduleTraitCor = cor(MEs, datTraits[,-1] , use = "p");
    colnames(moduleTraitCor) = "group"
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    pdf("./plots/rnaseq/1.degs/WGCNA.heatmap.plot.pdf",height = 10)
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(datTraits)[2],
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   font.lab.y = 0.5,
                   yColorWidth = strwidth("M"),
                   xColorWidth = strwidth("M"),
                   main = paste("Module-trait relationships"))
    dev.off()
    # step 6-Specific genetic analysis of modules for traits of interest
    group_nt = as.data.frame(datTraits$group);
    names(group_nt) = "conditions"
    module = str_remove(rownames(moduleTraitCor)[which(moduleTraitCor[,"group"]==max(moduleTraitCor[,"group"]))],"ME")
    modNames = substring(names(MEs), 3)
    geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    names(geneModuleMembership) = paste("MM", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    geneModuleMembership[1:4,1:4]

    geneTraitSignificance = as.data.frame(cor(datExpr, group_nt, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) = paste("GS.", names(group_nt), sep="");
    names(GSPvalue) = paste("p.GS.", names(group_nt), sep="");
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    pdf("./plots/rnaseq/1.degs/WGCNA.optimum.module.Scatterplot.pdf")
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for Luminal",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
    ## step 7
    # Select module genes
    probes = colnames(datExpr)
    inModule = (moduleColors==module);
    modProbes = probes[inModule];
    mod<-as.data.frame(modProbes)
    gene<-cbind(geneTraitSignificance,geneModuleMembership[,paste("MM",module,sep="")])
    gene1<-gene[which(rownames(gene)%in%mod[,1]=="TRUE"),]
    write.csv(gene1, './plots/rnaseq/1.degs/WGCNA.N.vs.T.optimal.gene.module.csv', row.names = T)
  }else{print("These data are not suitable for WGCNA analysis!")}
}

.wgcna(tpm_log,control="control",tumor="T1",gene_count=5000)

##output functional enrichment result
write.csv(ken@result, './plots/rnaseq/1.degs/KEGG.enrichment.all.deg.pvalue.table.csv')
write.csv(goCC@result, './plots/rnaseq/1.degs/GO.CC.enrichment.all.deg.pvalue.table.csv')
write.csv(goBP@result, './plots/rnaseq/1.degs/GO.BP.enrichment.all.deg.pvalue.table.csv')
write.csv(goMF@result, './plots/rnaseq/1.degs/GO.MF.enrichment.all.deg.pvalue.table.csv')
write.csv(gsea_hallmark, './plots/rnaseq/1.degs/HALLMARK.GSEA.enrichment.all.gene.pvalue.table.csv')
write.csv(gsea_kegg, './plots/rnaseq/1.degs/KEGG.GSEA.enrichment.all.gene.pvalue.table.csv')
write.csv(gsea_gobp, './plots/rnaseq/1.degs/GO.BP.GSEA.enrichment.all.gene.pvalue.table.csv')
write.csv(gsea_reactome, './plots/rnaseq/1.degs/REACTOME.GSEA.enrichment.all.gene.deg.pvalue.table.csv')
write.csv(gsva_hallmark, './plots/rnaseq/1.degs/HALLMARK.genesets.gsva.value.matrix.csv')
write.csv(gsva_kegg, './plots/rnaseq/1.degs/KEGG.genesets.gsva.value.matrix.csv')
write.csv(gsva_go, './plots/rnaseq/1.degs/GO.BP.genesets.gsva.value.matrix.csv')
write.csv(gsva_reactome, './plots/rnaseq/1.degs/REACTOME.genesets.gsva.value.matrix.csv')
write.csv(gsva_hallmark_diff, './plots/rnaseq/1.degs/HALLMARK.genesets.gsva.value.deg.pvalue.table.csv')
write.csv(gsva_kegg_diff, './plots/rnaseq/1.degs/KEGG.genesets.gsva.value.deg.pvalue.table.csv')
write.csv(gsva_go_diff, './plots/rnaseq/1.degs/GO.BP.genesets.gsva.value.deg.pvalue.table.csv')
write.csv(gsva_reactome_diff, './plots/rnaseq/1.degs/REACTOME.genesets.gsva.value.deg.pvalue.table.csv')

getwd()
save.image("/thinker/aid2/udata/all/gao/others/tcga/lgg/RNA_star_count_14_sample.RData")
load("/thinker/aid2/udata/all/gao/others/tcga/lgg/RNA_star_count.RData")




