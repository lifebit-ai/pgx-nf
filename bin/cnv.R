#!/usr/bin/env Rscript
library(PharmacoGxPrivate)
library(PharmacoGx)

options(stringsAsFactors=FALSE)
expprocessed= commandArgs(trailingOnly=TRUE)[1]
cnvinfo= commandArgs(trailingOnly=TRUE)[2]
cnvfeature= commandArgs(trailingOnly=TRUE)[3]

CNV_processed <- read.csv(expprocessed, sep="\t", check.names = FALSE)
CNV_processed$genename[12803] <- "PLPPR3"
CNV_processed$genename[16089] <- "SPHKAP"
CNV_processed$genename[25180] <- "USP17L2"
CNV_processed$genename[27091] <- "MED27"
rownames(CNV_processed) <- CNV_processed$genename
CNV_info <- read.csv(cnvinfo, row.names = 1)
colnames(CNV_processed[,6:82]) <- CNV_info$cellid
annot_gene <- read.csv(cnvfeature, row.names = 1)
tt <- annot_gene[match(CNV_processed$genename, annot_gene$gene_name), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(tt) <- CNV_processed$genename
tt$gene_name <- CNV_processed$genename
CNV_processed[,c("chrom","start","end","geneid","genename")] <- NULL
colnames(CNV_processed) <- rownames(CNV_info)
CNV_processed <- as.matrix(CNV_processed)
CNV_processed <- Biobase::ExpressionSet(CNV_processed)
fData(CNV_processed) <- tt
pData(CNV_processed) <- CNV_info
annotation(CNV_processed) <- "cnv"

save(CNV_processed, file="cnv_processed.RData")