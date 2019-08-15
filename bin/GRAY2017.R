#!/usr/bin/env Rscript
getGRAYP <-
  function (
    verbose=FALSE,
    nthread=1){
    options(stringsAsFactors=FALSE)
    badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    z <- list()
    
    
    load(commandArgs(trailingOnly=TRUE)[1])
    load(commandArgs(trailingOnly=TRUE)[2])
    load(commandArgs(trailingOnly=TRUE)[3])
    
    #match to cell curation
    
    matchToIDTableCELL <- function(ids,tbl, column) {
      sapply(ids, function(x) {
        myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating cell ids")
        }
        return(tbl[myx, "unique.cellid"])
      })
    }
    
    
    #match to drug curation 
    
    matchToIDTableDRUG <- function(ids,tbl, column) {
      sapply(ids,function(x) {
        myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating drug ids")
        }
        return(tbl[myx, "unique.drugid"])
      })
    }
    
    
    removed_drugs <- c("Topotecan","L-779405","Crizotinib","PF-3084014","PF-3814735","PF-4691502","GSK2141795","GSK1059615","GSK1059868","GSK2119563","CGC-11144","Tykerb:IGF1R (1:1)","GSK-650394") 
    
    curationDrug <- curationDrug[which(!curationDrug$unique.drugid %in% removed_drugs),]
    
    ## cell information
    
    load(commandArgs(trailingOnly=TRUE)[4])

    
    
    
    ########################
    ####Drug Sensitivity####
    ########################
    
    print("starting sensitivity")
      
    load(commandArgs(trailingOnly=TRUE)[5])
    load(commandArgs(trailingOnly=TRUE)[6])
    
    sensitivity.profiles <- data.frame("auc_recomputed" = recomputed$AUC, "ic50_recomputed"=recomputed$IC50)
   
    #add published GR50, GEC50, GRmax, GRinf, GR Hill Coefficient, GR_AOC to sensitivity profiles
    badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]|[()]|[(]|[)]"
    gray.published <- gray.GRvalues
    gray.published <- gray.published[,c("Cell.Name","Small.Molecule.HMS.LINCS.ID.1","GR50","GEC50","GRmax","GRinf","GR.Hill.Coefficient", "GR_AOC")]
    remove.drug <- grep("Topotecan",gray.published$Small.Molecule.HMS.LINCS.ID.1)
    gray.published <-gray.published[-remove.drug,]
    
    curationCellTemp <- curationCell
    curationCellTemp$GRAY.cellid <- tolower(gsub(badchars, "",curationCellTemp$GRAY.cellid))
    
    curationDrugTemp <- curationDrug
    curationDrugTemp$GRAY.drugid <- tolower(gsub(badchars, "",curationDrugTemp$GRAY.drugid))
    
    published.cells.matched <- as.character(matchToIDTableCELL(tolower(gsub(badchars, "",gray.published$Cell.Name)), curationCellTemp, "GRAY.cellid"))
    published.drugs.matched <- as.character(matchToIDTableDRUG(tolower(gsub(badchars, "",gray.published$Small.Molecule.HMS.LINCS.ID.1)), curationDrugTemp, "GRAY.drugid"))
    
    gray.published$Cell.Name <- published.cells.matched
    gray.published$Small.Molecule.HMS.LINCS.ID.1 <- published.drugs.matched
    
    published_profiles <- sprintf("%s_%s",gray.published$Small.Molecule.HMS.LINCS.ID.1,gray.published$Cell.Name)
    ## handle replicates
    tt <- published_profiles
    for(i in 1:length(tt)) {
      xx <- which(tt == tt[i])
      if(length(xx) > 1) {
        for(j in 1:length(xx)) {
          tt[xx[j]] <- paste(tt[xx[j]], j, sep="_")
        }
      }
    }
    
    rownames(gray.published) <- tt
    sensitivity.profiles[,c("GR50_published","GEC50_published","GRmax_published","GRinf_published","GR.Hill.Coefficient_published","GR_AOC_published")] <- as.numeric(NA)
    
    sensitivity.profiles$GR50_published <- gray.published$GR50[match(rownames(sensitivity.profiles), rownames(gray.published))]
    sensitivity.profiles$GEC50_published <- gray.published$GEC50[match(rownames(sensitivity.profiles), rownames(gray.published))]
    sensitivity.profiles$GRmax_published <- gray.published$GRmax[match(rownames(sensitivity.profiles), rownames(gray.published))]
    sensitivity.profiles$GRinf_published <- gray.published$GRinf[match(rownames(sensitivity.profiles), rownames(gray.published))]
    sensitivity.profiles$GR.Hill.Coefficient_published <- gray.published$GR.Hill.Coefficient[match(rownames(sensitivity.profiles), rownames(gray.published))]
    sensitivity.profiles$GR_AOC_published <- gray.published$GR_AOC[match(rownames(sensitivity.profiles), rownames(gray.published))]
    
    
   
    slope <- NULL
    for(exp in rownames(raw.sensitivity)){
      slope <- c(slope, computeSlope(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"])) #computeSlope (returns normalized slope of drug response curve)
    }
    
    names(slope) <- rownames(raw.sensitivity)
    sensitivity.profiles <- cbind(sensitivity.profiles, "slope_recomputed"=slope)
    
    sensitivity.profiles[sensitivity.profiles == 'NULL'] <- NA
    
   
    druginfo <- data.frame("drugid"=curationDrug$unique.drugid)
    rownames(druginfo) <- druginfo$drugid
    
    print("finished sensitivity")
    
    #RNA-seq Processed Data
    
    load(commandArgs(trailingOnly=TRUE)[7])
    
    print("finished RNASeq")
    
    #RPPA Processed Data
    
    load(commandArgs(trailingOnly=TRUE)[8])
    
    print("finished RPPA")
    
    #RNA U133A Array & Exon Array Processed Data
    
    load(commandArgs(trailingOnly=TRUE)[9])
    
    print("finished RNA")
    
    #SNP Processed Data
    
    load(commandArgs(trailingOnly=TRUE)[10])
    
    print("finished SNP")
    
    #methylation processed data
    
    load(commandArgs(trailingOnly=TRUE)[11])
    
    print("finished methylation")
    
    z <- c(z,c(
      "rnaseq.exp.matrix"=RNA_exp_matrix_processed,
      "rnaseq.exp"=RNA_exp_processed,
      "rnaseq.counts"=RNA_counts_processed,
      "rppa"=RPPA_processed,
      "cnv"=CNV_processed,
      "rna.exon"=Exon_processed,
      "rna.u133a"=U133A_processed,
      "methylation"=Methylation_processed)
    )
    
    GRAY2017 <- PharmacoSet(molecularProfiles=z,
                        
                        name="GRAY", 
                        cell=cellineinfo, 
                        drug=druginfo, 
                        sensitivityInfo=sensitivity.info, 
                        sensitivityRaw=raw.sensitivity, 
                        sensitivityProfiles=sensitivity.profiles, 
                        sensitivityN=NULL,
                        curationCell=curationCell, 
                        curationDrug=curationDrug, 
                        curationTissue=curationTissue, 
                        datasetType="sensitivity")
    save(GRAY2017,file="GRAY_2017.RData")
    
    
    return (GRAY2017)
    
  }



library(PharmacoGxPrivate)
library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)




getGRAYP(verbose=FALSE,
         nthread=1)






