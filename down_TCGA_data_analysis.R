rm(list = ls());
# 0. setup environment  ---------------------------------------------------
.libPaths(c('/thinker/storage/biostaff/Biosoft/gao/R_lib161'
            ,'/thinker/storage/udata/bing/RPackages'
            ,'/thinker/storage/udata/bing/anaconda3/lib/R/library/'
            ,'/thinker/storage/udata/muyl/software/R_lib161', .libPaths()))
# .libPaths(c("/thinker/storage/udata/sunhr/R/R4.2.0/lib/R/library/", .libPaths()))

# 先安装devtools包
# install.packages("devtools")
# library(devtools)
# devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(stringr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
library(pheatmap)

library(TCGAbiolinks)
getwd()
# mutation
# maf_lihc <- GDCquery_Maf('LIHC', pipelines = 'mutect')
maf_lgg.mutect <- GDCquery_Maf("LGG", pipelines = "mutect2")
# maf_lgg.mutect <-
query.maf.hg19 <- GDCquery(project = "TCGA-LGG",
                           data.category = "Simple nucleotide variation",
                           data.type = "Simple somatic mutation",
                           # workflow.type = '',
                           access = "open",
                           legacy = TRUE
)

GDCdownload(query.maf.hg19)
maf_lgg <- GDCprepare(query.maf.hg19)

query.maf.hg19 <- GDCquery(project = "TCGA-LGG",
                           data.category = "Simple nucleotide variation",
                           data.type = "Simple somatic mutation",
                           # workflow.type = ''
                           access = "open",
                           legacy = TRUE
                           )
GDCdownload(query.maf.hg19)
maf_lgg <- GDCprepare(query.maf.hg19)
# prepare mutation manually
lihc.maf.df <- fread('./GDCdata/TCGA-LIHC/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/a630f0a0-39b3-4aab-8181-89c1dde8d3e2/TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf.gz')
# 用maftools内置的maf文件看看
lgg.maf = system.file('ext')
# 函数 ----------------------------------------------------------------------
getwd()
get_and_save_MAF <- function(project_name='LGG') {
  current.maf <- GDCquery_Maf(tumor=project_name, pipelines='mutect2')
  saveRDS(current.maf, file=str_interp('${project_name}.maf.RDS'))
}


# 批量下载 --------------------------------------------------------------------

get_and_save_MAF('LGG')




# mRNA   fpkm
query.star.count <- GDCquery(project = "TCGA-LGG",
                             experimental.strategy = "RNA-Seq",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "STAR - Counts"
                             # , legacy = TRUE
)
GDCdownload(query.star.count)

# query.htseq.count <- GDCquery(project = "TCGA-LGG",
#                              # experimental.strategy = "RNA-Seq",
#                              data.category = "Transcriptome Profiling",
#                              data.type = "Gene Expression Quantification",
#                              workflow.type = "HTSeq - Counts"
#                              # , legacy = TRUE
# )
# GDCdownload(query.htseq.count)

# lihc_Rnaseq <- GDCprepare(query, summarizedExperiment = FALSE)
# lihc_Rnaseq <- GDCprepare(query)
# prepare STAR results
file_table <- query.star.count$results[[1]]
file_table
# file_table[file_table$file_id == 'ff12abd3-0f45-4063-afa7-fa5cad973159',]
fn <- paste0('./GDCdata/TCGA-LGG/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/'
             , file_table[1, 'file_id'], '/', file_table[1, 'file_name'])
tmp_rc <- fread(fn)
tmp_rc <- tmp_rc[5:nrow(tmp_rc),]

lgg_star_count <- data.frame(gene_name=tmp_rc$gene_name, stringsAsFactors = FALSE)
lgg_star_count
for (i in 1:nrow(file_table)){
  fn <- paste0('./GDCdata/TCGA-LGG/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/'
               , file_table[i, 'file_id'], '/', file_table[i, 'file_name'])
  tmp_rc <- fread(fn)
  tmp_rc <- tmp_rc[5:nrow(tmp_rc),]

  lgg_star_count[, file_table[i, 'cases']] <- tmp_rc$unstranded
}

lihc_rna <- lihc_star_count
lihc_rna
# rownames(lihc_rna) <- lihc_rna$gene_name


# miRNA
query.mirna <- GDCquery(
  project = "TCGA-LGG",
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling",
  # barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
  data.type = "miRNA Expression Quantification"
)
GDCdownload(query.mirna)
lihc_mirna <- GDCprepare(query = query.mirna)

lihc_mirna_read_count <- lihc_mirna[,grepl('read_count',colnames(lihc_mirna))]
colnames(lihc_mirna_read_count) <- str_remove(colnames(lihc_mirna_read_count), 'read_count_')
rownames(lihc_mirna_read_count) <- lihc_mirna$miRNA_ID

# clinic data
lihc_clinical_index <- GDCquery_clinic(project = "TCGA-LGG", type = "clinical")
