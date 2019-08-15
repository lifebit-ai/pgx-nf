#!/usr/bin/env Rscript
library(PharmacoGxPrivate)
library(PharmacoGx)
library(readxl)
library(openxlsx)

options(stringsAsFactors=FALSE)
expprocessed= commandArgs(trailingOnly=TRUE)[1]
proteininfo= commandArgs(trailingOnly=TRUE)[2]
proteinfeature= commandArgs(trailingOnly=TRUE)[3]

RPPA_processed <- read.xlsx(expprocessed, sheet = 1, startRow = 3)
rownames(RPPA_processed) <- RPPA_processed$X1
RPPA_processed$X1 <- NULL
RPPA_processed <- t(RPPA_processed)
RPPA_processed <- as.matrix(RPPA_processed)
RPPA_processed <- Biobase::ExpressionSet(RPPA_processed)
protein_info <- read.csv(file=proteininfo, row.names = 1)
protein_feature <- read.csv(file=proteinfeature)
rownames(protein_feature) <- protein_feature$proteinid
colnames(RPPA_processed) <- protein_info$cellid
pData(RPPA_processed) <- protein_info
fData(RPPA_processed) <- protein_feature
annotation(RPPA_processed) <- "RPPA"


save(RPPA_processed, file="RPPA_processed.RData")