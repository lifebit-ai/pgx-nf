#!/usr/bin/env Rscript
getDrugCuration <- function(verbose=FALSE){
  options(stringsAsFactors=FALSE)
  
  drug_all = commandArgs(trailingOnly=TRUE)[1]
  drug_all <- read.csv(file= drug_all, na.strings=c("", " ", "NA"))
  curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
  curationDrug <- curationDrug[ , c("unique.drugid", "GRAY.drugid")]
  rownames(curationDrug) <- curationDrug[ , "unique.drugid"]
  
  
  save(curationDrug, file="drug_cur.RData")
  
  return(curationDrug)

}

getDrugCuration(verbose = FALSE)