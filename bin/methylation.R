#!/usr/bin/env Rscript
library(PharmacoGxPrivate)
library(PharmacoGx)

options(stringsAsFactors=FALSE)

methylationexp = commandArgs(trailingOnly=TRUE)[1]
methylationinfo = commandArgs(trailingOnly=TRUE)[2]
methylationfeature = commandArgs(trailingOnly=TRUE)[3]

Methylation_processed <- read.csv(methylationexp,  sep="\t", skip = 8, check.names = FALSE)
rownames(Methylation_processed) <- Methylation_processed$TargetID
Methylation_processed <- Methylation_processed[, -grep("Signal", colnames(Methylation_processed))]
Methylation_processed <- Methylation_processed[, -grep("Detection", colnames(Methylation_processed))]
Methylation_info <- read.csv(methylationinfo, stringsAsFactors=FALSE, row.names=1)
Methylation_feature <- read.csv(methylationfeature, skip = 7)
Methylation_feature$Name <- NULL
colnames(Methylation_processed) <- sub(".AVG_Beta", "", colnames(Methylation_processed))
colnames(Methylation_processed) <- gsub(" ", "", colnames(Methylation_processed))
rownames(Methylation_feature) <- Methylation_feature$IlmnID
array_match <- Methylation_feature[match(Methylation_processed$TargetID, Methylation_feature$IlmnID),]
Methylation_processed$TargetID <- NULL
Methylation_processed$SYMBOL <- NULL
colnames(Methylation_processed) <- rownames(Methylation_info)
Methylation_processed  <- as.matrix(Methylation_processed)
Methylation_processed <-  Biobase::ExpressionSet(Methylation_processed)
pData(Methylation_processed) <- Methylation_info
fData(Methylation_processed) <- array_match
annotation(Methylation_processed) <- "Average Beta Values"


save(Methylation_processed, file="Methylation_processed.RData")