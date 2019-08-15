#!/usr/bin/env Rscript
library(PharmacoGxPrivate)
library(PharmacoGx)

options(stringsAsFactors=FALSE)
load(commandArgs(trailingOnly=TRUE)[1])

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
U133Aexp = commandArgs(trailingOnly=TRUE)[2]
U133Ainfo = commandArgs(trailingOnly=TRUE)[3]
U133Afeature = commandArgs(trailingOnly=TRUE)[4]

Exonexp = commandArgs(trailingOnly=TRUE)[5]
Exoninfo = commandArgs(trailingOnly=TRUE)[6]
Exonfeature = commandArgs(trailingOnly=TRUE)[7]


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

#RNA U133A Array

U133A_processed <- read.csv(U133Aexp, sep="\t")
U133A_sample_info <- read.csv(U133Ainfo, sep="\t")
rownames(U133A_sample_info) <- U133A_sample_info$Scan.Name
U133A_processed <- U133A_processed[-1,]
colnames(U133A_processed) <- gsub("^X", "",  colnames(U133A_processed))
array_design <- read.csv(U133Afeature, sep="\t", skip =17)
array_match <- array_design[match(U133A_processed$Scan.REF, array_design$Composite.Element.Name),]
rownames(array_match) <- array_match$Composite.Element.Name
rownames(U133A_processed) <- U133A_processed$Scan.REF
U133A_processed$Scan.REF <- NULL
U133A_sample_info <- U133A_sample_info[match(tolower(gsub(badchars, "",colnames(U133A_processed))), tolower(gsub(badchars, "",rownames(U133A_sample_info)))),]
rownames(U133A_sample_info) <- colnames(U133A_processed)
colnames(U133A_sample_info)[6] <- "cellid"
U133A_sample_info$cellid <- as.character(matchToIDTableCELL(U133A_sample_info$cellid, curationCell, "GRAY.cellid"))
rownames(U133A_sample_info) <- U133A_sample_info$cellid 
colnames(U133A_processed) <- U133A_sample_info$cellid
U133A_processed <- as.matrix(U133A_processed)
U133A_processed <-  Biobase::ExpressionSet(U133A_processed)
fData(U133A_processed) <- array_match
pData(U133A_processed) <- U133A_sample_info
annotation(U133A_processed) <- "U133A Array"


#RNA Exon Array

Exon_processed <- read.csv(Exonexp, sep="\t", check.names = FALSE)
Exon_processed <- Exon_processed[-1,]
Exon_processed <- Exon_processed[!is.na(Exon_processed$GeneID),]
rownames(Exon_processed) <- Exon_processed$GeneID
Exon_sample_info <- read.csv(Exoninfo, sep="\t")
rownames(Exon_sample_info) <- Exon_sample_info$Source.Name
annot_gene <- read.csv(Exonfeature, row.names = 1)
annot_match <- annot_gene[match(rownames(Exon_processed), annot_gene$gene_name), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(annot_match) <- Exon_processed$GeneID
annot_match$gene_name <- Exon_processed$GeneID
Exon_processed$GeneID <- NULL
colnames(Exon_sample_info)[2] <- "cellid"
Exon_sample_info <- Exon_sample_info[match((colnames(Exon_processed)), Exon_sample_info$cellid),]
Exon_sample_info$cellid <- as.character(matchToIDTableCELL(Exon_sample_info$cellid, curationCell, "GRAY.cellid"))
rownames(Exon_sample_info) <- Exon_sample_info$cellid
colnames(Exon_processed) <- Exon_sample_info$cellid
Exon_processed <- as.matrix(Exon_processed)
Exon_processed <-  Biobase::ExpressionSet(Exon_processed)
fData(Exon_processed) <- annot_match
pData(Exon_processed) <- Exon_sample_info
annotation(Exon_processed) <- "Exon Array"


save(U133A_processed, Exon_processed, file="RNA_processed.RData")

