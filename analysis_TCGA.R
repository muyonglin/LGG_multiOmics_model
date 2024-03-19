.libPaths(c("/thinker/aid2/udata/all/gao/R/R4.2.0/lib/R/library","/thinker/storage/udata/sunhr/R/R4.2.0/lib/R/library/", .libPaths()))
.libPaths(c("/thinker/aid2/udata/all/muyl/others/software/condaEnvs/tcgabiolinks/lib/R/library"),.libPaths())
.libPaths()
# .libPaths(c('/thinker/storage/biostaff/Biosoft/gao/R_lib161'
#             ,'/thinker/storage/udata/bing/RPackages'
#             ,'/thinker/storage/udata/bing/anaconda3/lib/R/library/'
#             ,'/thinker/storage/udata/muyl/software/R_lib161', .libPaths()))
library(stringr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
library(pheatmap)
library(maftools)
library(TCGAbiolinks)
library(pheatmap)
library(org.Hs.eg.db)
library(readr)
library(ggpubr)
library(BSgenome)
library(MutationalPatterns)
library(NMF)
library(tidyverse)
library(GenomicRanges)
library(maftools)

# mutation
# library(dplyr)
# library(DT)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

## 下载TCGA数据
# conda create -p /thinker/aid2/udata/all/muyl/others/software/condaEnvs/tcgabiolinks
# conda activate /thinker/aid2/udata/all/muyl/others/software/condaEnvs/tcgabiolinks
# conda install -c conda-forge r-base

## 0. set environment --------
## /thinker/aid2/udata/all/muyl/others/software/condaEnvs/tcgabiolinks
# ph <- '/thinker/aid2/udata/all/muyl/others/20.wujiaoxiang.fss.nanjing/'
# setwd(ph)
library(TCGAbiolinks)
# 1. download LGG files --------
# run in conda environment: /thinker/aid2/udata/all/muyl/others/software/condaEnvs/tcgabiolinks
# 1.1 download clinic information --------
clinic.index <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
write.csv(clinic.index, './tcga.lgg.clinic.information.table.csv', row.names = FALSE)
clinic.index <- read.csv("./tcga.lgg.clinic.information.table.csv")
clinic.index
# 1.2 download LGG mutation information --------
query.maf <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  legacy = FALSE,
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query.maf)
maf <- GDCprepare(query.maf)
library(stringr)
readr::write_tsv(maf, 'lgg.maf')

rm(list = ls())

lgg_maf_file ="/thinker/aid2/udata/all/gao/others/tcga/lgg/lgg.maf"
lgg_maf_file
lgg_maf = read.maf(lgg_maf_file)
getwd()

output_dir="/thinker/aid2/udata/all/gao/others/tcga/lgg/plots/"
dir.create(file.path(output_dir,'1.mutations/'), recursive = TRUE)
# pdf(file.path(output_dir, '1.mutations/wes.maf.summary.plot.pdf')
#     ,height = 10, width = 10)
pdf(file.path(output_dir, '1.mutations/wes.maf.summary.plot.pdf'))
plotmafSummary(lgg_maf,textSize = 0.6)
dev.off()
pdf(file.path(output_dir, '1.mutations/oncoplot.top.30.genes02.pdf'))
oncoplot(lgg_maf,
         sortByAnnotation=T,
         top=30,
         gene_mar = 6,
         borderCol=NULL,
         legend_height = 5)
         # , SampleNamefontSize=0.6
         # , sepwd_genes=1.5
         # , sepwd_samples=1.5
         # , sampleOrder = sampleOrder
         # , removeNonMutated = F
         # , gene_mar = 6)
         # , barcode_mar = 9)
dev.off()

# 3.1.3 driver gene oncoplot --------
# driver gene table
# ph <- '/home/ubuntu/efs/bed/drive_genes/'
# ph <- '/thinker/storage/biostaff/Biosoft/gao/drive_genes/'
driver_ph = '/thinker/storage/biostaff/Biosoft/gao/drive_genes/'
cell_driver_table <- fread(file.path(driver_ph, 'driverGene2018.csv'))
cell_driver_table <- as.data.frame(cell_driver_table)

nat_driver_table <- fread(file.path(driver_ph, 'Compendium_Cancer_Genes.tsv'))
nat_driver_table <- as.data.frame(nat_driver_table)

driver_gene_list <- unique(c(cell_driver_table$Gene, nat_driver_table$SYMBOL))
driver.maf <- subsetMaf(lgg_maf, genes=driver_gene_list)

pdf(file.path(output_dir, '1.mutations/oncoplot.top.30.driver.genes.pdf'))
oncoplot(driver.maf
         ,sortByAnnotation=T,
         top=30,
         gene_mar = 6,
         borderCol=NULL,
         legend_height = 5)

dev.off()
write.table(driver.maf@data,file.path(output_dir,"1.mutations/lgg_driver.maf_read_maf_data.xls"),sep = "\t",row.names = F)



write.table(lgg_maf@data,file.path(output_dir,"1.mutations/lgg_maf_read_maf_data.xls"),sep = "\t",row.names = F)
# 3.3 titv ----------
# 3.3.0 titv plot --------
dir.create(file.path(output_dir, '3.signature/')
           ,recursive=TRUE)
# 1.1 titv raw plot --------


laml.titv = titv(maf = lgg_maf, plot = FALSE, useSyn = TRUE)
pdf(file.path(output_dir, '3.signature/maf.titv.summary.plot.pdf'))
plotTiTv(res = laml.titv)
dev.off()
laml.titv$fraction.contribution
write.csv(laml.titv$fraction.contribution
          , file.path(output_dir, '3.signature/all.titv.fraction.contribution.csv'))
laml.titv$TiTv.fractions
write.csv(laml.titv$TiTv.fractions
          , file.path(output_dir, '3.signature/all.TiTv.fractions.csv'))
laml.titv$raw.counts
write.csv(laml.titv$raw.counts
          , file.path(output_dir, '3.signature/all.laml.titv$raw.counts.csv'))

# 3.3.1 signature analysis --------
lgg_maf@data$Tumor_Sample_Barcode
mut_mat <- trinucleotideMatrix(lgg_maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
mut_mat1 <- mut_mat$nmf_matrix
mut_mat1 <- t(mut_mat1)

mut_mat1 <- as.data.frame(mut_mat1)
mut_mat1 <- mut_mat1 + 0.0001
mut_mat1
# 96 profile plot
# plot_96_profile(mut_mat1, condensed = TRUE)
# 96 profiles
# pdf(file.path(output_dir, '3.signature/96.nucleotide.mutation.plots.all.samples.pdf'))
# plot_96_profile(mut_mat1, condensed = TRUE)
pdf(file.path(output_dir, '3.signature/96.nucleotide.mutation.plots.all.samples.pdf'))
plot_96_profile(mut_mat1)
dev.off()
write.csv(mut_mat1
          ,file.path(output_dir, '3.signature/96.nucleotide.mutation.plots.all.samples.csv'))

# signature
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/"
                , "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures
new_order = match(row.names(mut_mat1), cancer_signatures$Somatic.Mutation.Type)
new_order
cancer_signatures = cancer_signatures[as.vector(new_order),]
cancer_signatures
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type

cancer_signatures = as.matrix(cancer_signatures[,4:33])

fit_res <- fit_to_signatures(mut_mat1, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 10)

tmpdf <- fit_res$contribution[select,]
tmpdf <- as.data.frame(tmpdf)
tmpdf[,'rowsum'] <- apply(tmpdf, 1, sum)
tmpdf <- tmpdf[order(tmpdf[,'rowsum'], decreasing = T),]
tmpdf <- tmpdf[,1:(ncol(tmpdf) - 1)]
tmpdf <- as.matrix(tmpdf)

pdf(file.path(output_dir, '3.signature/cosmic.signature.relative.contribution.stacked.barplot.pdf'))
# plot_contribution(fit_res$contribution,
#                   coord_flip = T
# )
plot_contribution(fit_res$contribution
                  # coord_flip = T
)
dev.off()

fit_res$contribution
write.csv(fit_res$contribution
          ,file.path(output_dir, '3.signature/cosmic.signature.relative.contribution.stacked.barplot.csv'))

pdf(file.path(output_dir, '3.signature/cosmic.signature.relative.contribution.heatmap02.pdf'))
plot_contribution_heatmap(tmpdf,
                          cluster_samples = FALSE,
                          method = 'average',
                          show_rownames = F,
)

dev.off()

write.csv(tmpdf
          , file.path(output_dir, '3.signature/cosmic.signature.relative.contribution.heatmap.csv'))

# sample mutation profile similarity
cos_sim_samples <- cos_sim_matrix(mut_matrix1 = mut_mat1, mut_matrix2 = mut_mat1)
# plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)
# pdf(str_c(output_dir, '3.signature/sample.mutation.profile.pairwise.similarity.heatmap.pdf'))
# p <- pheatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)
# print(p)
# dev.off()

getwd()
pheatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA,show_rownames = F, show_colnames = F
         ,filename = file.path(output_dir, '3.signature/sample.mutation.profile.pairwise.similarity.heatmap.pdf'))


# KEGG --------------------------------------------------------------------
# GO kegg 富集通路
# c. kegg pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)
options(clusterProfiler.download.method = "wininet")
R.utils::setOption("clusterProfiler.download.method",'auto')
go_gene <- lgg_maf@data$Hugo_Symbol
go_gene
gene_L <- go_gene
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

# fn <- 'all.N.vs.T'
# /data/gao/K19109_WES
# 第一个为生物学途径（BP），表示分子功能的有序组合，以达到更广的生物功能，如有丝分裂或嘌呤代谢等；
# 第二个为细胞组份（CC），表示用于描述亚细胞结构、位置和大分子复合物，如核仁、端粒和识别起始的复合物等；
# 第三个为分子功能（MF），表示用于描述基因、基因产物的功能，如碳水化合物或ATP水解酶活性等。
dir.create(file.path(output_dir,'4.KEGG_GO/')
           ,recursive=TRUE)

pdf(file.path(output_dir,'4.KEGG_GO/enrichment.GO.KEGG.degs.pdf')
    ,height = 12, width = 12)
p <- dotplot(ken,title='KEGG enrich of mutation gene')
print(p)
p <- dotplot(goCC,title='GO CC enrich of mutation gene')
print(p)
p <- dotplot(goBP,title='GO BP enrich of mutation gene')
print(p)
p <- dotplot(goMF,title='GO MF enrich of mutation gene')
print(p)
dev.off()

# temp_ken_Infiltrat <- as.data.frame(ken)
# temp_ken_Infiltrat$geneID
# uniprot_ken_Infiltrat <- strsplit(temp_ken_Infiltrat$geneID,"/")
# symbol <- sapply(uniprot_ken_Infiltrat,function(x){
#   y=bitr(x, fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db")
#   #一对多，取第一个
#   y=y[!duplicated(y$SYMBOL),-1]
#   y=paste(y,collapse = "/")
# })
# symbol
# temp_ken_Infiltrat$geneID <- symbol
# temp_ken_Infiltrat
#
# write.table(temp_ken_Infiltrat,file="Figures/5.GO.KEGG20211015/kegg.cna_Infiltrat.maf.genes.xls",sep="\t",row.names=FALSE,quote =FALSE)
# write.table(goCC_Infiltrat,file="Figures/5.GO.KEGG20211015/goCC.cna_Infiltrat.maf.genes.xls",sep="\t",row.names=FALSE,quote =FALSE)
# write.table(goBP_Infiltrat,file="Figures/5.GO.KEGG20211015/goBP.cna_Infiltrat.maf.genes.xls",sep="\t",row.names=FALSE,quote =FALSE)
# write.table(goMF_Infiltrat,file="Figures/5.GO.KEGG20211015/goMF.cna_Infiltrat.maf.genes.xls",sep="\t",row.names=FALSE,quote =FALSE)

write.table(ken@result, file.path(output_dir,'4.KEGG_GO/kegg.enrichment.maf.genes.xls'),sep="\t",row.names=FALSE,quote =FALSE)
write.table(goCC@result, file.path(output_dir,'4.KEGG_GO/GO.CC.enrichment.all.mutated.maf.genes.xls'),sep="\t",row.names=FALSE,quote =FALSE)
write.table(goBP@result, file.path(output_dir,'4.KEGG_GO/GO.BP.enrichment.all.mutated.maf.genes.xls'),sep="\t",row.names=FALSE,quote =FALSE)
write.table(goMF@result, file.path(output_dir,'4.KEGG_GO/GO.MF.enrichment.all.mutated.maf.genes.xls'),sep="\t",row.names=FALSE,quote =FALSE)



# clinic data
lgg_clinical_index <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
lgg_clinical_index
# write.csv(clinical_index
#           , './ana_20211014/tcga.lihc.clinical.index.csv'
#           , row.names = FALSE)
# survival table
colnames(lgg_clinical_index)
surv_df <- lgg_clinical_index[, c('bcr_patient_barcode', 'days_to_death', 'days_to_last_follow_up', 'vital_status')]
surv_df

surv_df <- as.data.frame(surv_df)
getwd()
save.image("/thinker/aid2/udata/all/gao/others/tcga/lgg.RData")




