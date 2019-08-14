#!/usr/bin/env Rscript
getCellCuration <- function(verbose=FALSE){
options(stringsAsFactors=FALSE)
cell_all = commandArgs(trailingOnly=TRUE)[1]  
cell_all <- read.csv(file = cell_all, na.strings=c("", " ", "NA"))
curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
rownames(curationCell) <- curationCell[ , "unique.cellid"]

save(curationCell, file="cell_cur.RData")

return(curationCell)


}

getCellCuration(verbose = FALSE)