#!/usr/bin/env Rscript
library(PharmacoGxPrivate)
library(PharmacoGx)

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    #match to CurationCell
    
    matchToIDTableCELL <- function(ids,tbl, column) {
      sapply(ids, function(x) {
        myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating cell ids")
        }
        return(tbl[myx, "unique.cellid"])
      })
    }

options(stringsAsFactors=FALSE)
load(commandArgs(trailingOnly=TRUE)[1])
expmatrix = commandArgs(trailingOnly=TRUE)[2]
rnainfo = commandArgs(trailingOnly=TRUE)[3]
ensemblannot = commandArgs(trailingOnly=TRUE)[4]
expall = commandArgs(trailingOnly=TRUE)[5]
counts = commandArgs(trailingOnly=TRUE)[6]

RNA_exp_matrix_processed <- read.csv(expmatrix, sep = "\t", check.names = FALSE)
rownames(RNA_exp_matrix_processed) <- RNA_exp_matrix_processed$EnsEMBL_Gene_ID
RNA_info <- read.csv(file=rnainfo, row.names = 1)
annot_gene <- read.csv(ensemblannot, row.names = 1)
print(head(annot_gene))
annot_match <- annot_gene[match(rownames(RNA_exp_matrix_processed), rownames(annot_gene)),c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(annot_match) <- rownames(RNA_exp_matrix_processed)
annot_match$gene_id <- RNA_exp_matrix_processed$EnsEMBL_Gene_ID
RNA_exp_matrix_processed[,c("Gene_ID","FID","Seq_Name","EnsEMBL_Gene_ID")] <- NULL
colnames(RNA_exp_matrix_processed) <- rownames(RNA_info)
RNA_exp_matrix_processed <- as.matrix(RNA_exp_matrix_processed)
RNA_exp_matrix_processed <- Biobase::ExpressionSet(RNA_exp_matrix_processed)
pData(RNA_exp_matrix_processed) <- RNA_info
fData(RNA_exp_matrix_processed) <- annot_match
annotation(RNA_exp_matrix_processed) <- "RNAseq Exp Matrix"

RNA_exp_processed <- read.csv(file=expall, sep = "\t", check.names = FALSE)
rownames(RNA_exp_processed) <- RNA_exp_processed$EnsEMBL_Gene_ID
annot_match <- annot_gene[match(rownames(RNA_exp_processed), rownames(annot_gene)),c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(annot_match) <- rownames(RNA_exp_processed)
annot_match$gene_id <- RNA_exp_processed$EnsEMBL_Gene_ID
RNA_exp_processed[,c("Gene_ID","FID","Seq_Name","EnsEMBL_Gene_ID")] <- NULL
colnames(RNA_exp_processed) <- rownames(RNA_info)
RNA_exp_processed <- as.matrix(RNA_exp_matrix_processed)
RNA_exp_processed <- Biobase::ExpressionSet(RNA_exp_processed)
pData(RNA_exp_processed) <- RNA_info
fData(RNA_exp_processed) <- annot_match
annotation(RNA_exp_processed) <- "RNAseq Exp"

RNA_counts_processed <- read.csv(file=counts, sep = "\t", check.names = FALSE)
rownames(RNA_counts_processed) <- RNA_counts_processed$Gene
annot_match <- annot_gene[match(rownames(RNA_counts_processed), annot_gene$gene_name),c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(annot_match) <- rownames(RNA_counts_processed)
annot_match$gene_id <- RNA_counts_processed$Gene
RNA_counts_processed$Gene <- NULL
xx <- matchToIDTableCELL(colnames(RNA_counts_processed), curationCell, "GRAY.cellid")
colnames(RNA_counts_processed) <- xx
RNA_counts_processed <- as.matrix(RNA_counts_processed)
RNA_counts_processed <- Biobase::ExpressionSet(RNA_counts_processed)
pData(RNA_counts_processed) <- RNA_info
fData(RNA_counts_processed) <- annot_match
annotation(RNA_counts_processed) <- "RNAseq Counts"


save(RNA_exp_matrix_processed, RNA_exp_processed, RNA_counts_processed, file="RNAseq_processed.RData")