
getTissueCuration <- function(verbose=FALSE){
options(stringsAsFactors=FALSE)
cell_all = commandArgs(trailingOnly=TRUE)[1]
load(commandArgs(trailingOnly=TRUE)[2])  
cell_all <- read.csv(file = cell_all, na.strings=c("", " ", "NA"))
curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]

rownames(curationTissue) <- curationCell[ , "unique.cellid"]

save(curationTissue, file="tissue_cur.RData")

return(curationTissue)

}

getTissueCuration(verbose = FALSE)