#!/usr/bin/env Rscript
library(readxl)
library(openxlsx)
library(CoreGx)


getCelllineInfo <- function(verbose=FALSE){
  
options(stringsAsFactors=FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
load(commandArgs(trailingOnly=TRUE)[2])
print(head(curationCell))
cellline_info = commandArgs(trailingOnly=TRUE)[1]


#match to Cell Curation

matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}

cellineinfo <- read.xlsx(cellline_info, sheet = 1)
cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
rn <- cellineinfo[-1, 1]
cn <- t(cellineinfo[1, -1])
cn <- gsub(badchars, ".", cn)
cellineinfo <- cellineinfo[-1, -1]
dimnames(cellineinfo) <- list(rn, cn)
cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo[,1:10])
cellineinfo <- cellineinfo[which(!is.na(cellineinfo$Transcriptional.subtype)), ]
c1 <- matchToIDTableCELL(cellineinfo$cellid, curationCell, "GRAY.cellid")
c1 <- as.character(c1)
cellineinfo$cellid <- c1
cellineinfo$cellid[is.na(cellineinfo$cellid)]<-"NA"
rownames(cellineinfo) <-  cellineinfo$cellid
head(cellineinfo)

cellineinfo <- rbind(cellineinfo, c("SUM 190", "breast", "NA", "NA", "0","0","0","1","1","1","0","1"))
rownames(cellineinfo)[85] <- "SUM 190"

cellineinfo <- rbind(cellineinfo, c("HCC1500", "breast", "NA", "NA", "0","0","1","1","1","1","0","0"))
rownames(cellineinfo)[86] <- "HCC1500"

cellineinfo <- rbind(cellineinfo, c("DU-4475", "breast", "NA", "NA", "0","0","0","0","1","1","0","0"))
rownames(cellineinfo)[87] <- "DU-4475"

cellineinfo <- rbind(cellineinfo, c("HCC1007", "breast", "NA", "NA", "0","0","0","1","0","0","0","0"))
rownames(cellineinfo)[88] <- "HCC1007"

cellineinfo <- rbind(cellineinfo, c("MDA-MB-435", "breast", "NA", "NA", "0","0","0","1","0","0","0","0"))
rownames(cellineinfo)[89] <- "MDA-MB-435"

cellineinfo <- rbind(cellineinfo, c("HCC2157", "breast", "NA", "NA", "0","0","0","1","0","0","0","0"))
rownames(cellineinfo)[90] <- "HCC2157"

cellineinfo <- rbind(cellineinfo, c("184A1N4", "breast", "NA", "NA", "0","0","1","0","0","0","0","0"))
rownames(cellineinfo)[91] <- "184A1N4"



save(cellineinfo, file="cellline_info.RData")

return(cellineinfo)


}

getCelllineInfo(verbose = FALSE)