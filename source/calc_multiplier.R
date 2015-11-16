#calc_multiplier.R
#Calculates mutation rate multiplier
#Author: Boyang Zhao
#REQUIRES (datasets): 1) datasets_analyzed.RData (generated from analyze_basic.R)
#                     2) dProteinLengths (generated based on HGNC and UniProt FASTA files)
#GENERATES: dMAFcounts, dMultipliers

if(exists("config.avail") && config.avail){ #if configuration settings exist, keep these, clear the rest of workspace
  rm(list=ls()[grep("^[^config\\.].*$",ls())])
} else {
  rm(list=ls())
}

dset_dir <- "../datasets/"
dsetAnalyzed_dir <- paste(dset_dir, "analyzed/15.0301_p3_medtumor_all/", sep="")
dsetProteome_dir <- paste(dset_dir, "UniProt proteome/", sep="") #for dProteinLengths.RData
dsetCalc_dir <- paste(dsetAnalyzed_dir, "calc-series/null/", sep="") #output directory
source('util.R')

#settings
minRNAseqcutoff <- NULL #if NULL, then will not have any cutoffs and value of analyzeNotExpressed is irrelevant
analyzeNotExpressed <- TRUE
savePlots <- TRUE
datasetName <- "dMAFgenesMutExactSame"
dataset_matchName <- "dSNVentryMutExactSame"

#settings overwrite if configuration settings exist
if(exists("config.avail") && config.avail){
  dset_dir <- config.dset_dir
  dsetAnalyzed_dir <- config.dsetAnalyzed_dir
  dsetProteome_dir <- config.dsetProteome_dir
  dsetCalc_dir <- config.dsetCalc_dir
  
  minRNAseqcutoff <- config.calc.minRNAseqcutoff
  analyzeNotExpressed <- config.calc.analyzeNotExpressed
  savePlots <- config.calc.savePlots
  datasetName <- config.calc.datasetName
  dataset_matchName <- config.calc.dataset_matchName
}

#load datasets
load(paste(dsetAnalyzed_dir,"datasets_analyzed.RData",sep=""))
load(paste(dsetProteome_dir, "dProteinLengths.RData",sep="")) #load in dProteinLengths dataset

dataset <- eval(as.name(datasetName))
dataset_match <- eval(as.name(dataset_matchName))

#modify datasets based on cutoffs
if(!is.null(minRNAseqcutoff)){
  if(!analyzeNotExpressed){
    dataset <- dataset[(!is.na(dataset$mRNAvals_final) & dataset$mRNAvals_final >= minRNAseqcutoff) | 
                         is.na(dataset$mRNAvals_final),]
  } else {
    dataset <- dataset[dataset$mRNAvals_final < minRNAseqcutoff,]
  }
}

#quality controls
createDir(dsetCalc_dir) #create output dir if needed

###################################################
# Calculation methods
###################################################
calcCounts <- function(){
  #calculate mutation counts in atasets
  catn("Calculating mutation counts in dataset...")
  countfunc <- function(x){
    length(unique(x))
  }
  
  #count number of unique positions mutated per gene per tumor type
  posN_unique <- aggregate(dataset$aaPos, 
                           by=list(paste(dataset$geneSymbol,dataset$tumorType,sep="||")), 
                           FUN=countfunc)
  dMAFcounts <- data.frame(geneSymbol=sub("^(.*)\\|\\|.*$","\\1",posN_unique[[1]]),
                           tumorType=sub("^.*\\|\\|(.*)$","\\1",posN_unique[[1]]),
                           posN_unique=posN_unique[[2]],
                           stringsAsFactors=FALSE)
  protLengths_expanded <- dProteinLengths[match(dMAFcounts$geneSymbol, dProteinLengths$geneSymbol), 'length']
  
  dMAFcounts <- cbind(dMAFcounts, posNMatched_unique=rep(0,nrow(dMAFcounts)), protLength=protLengths_expanded, 
                      stringsAsFactors=FALSE)
  
  #loop through each tumor type, and get the exact match counts
  for(tumor in unique(dataset$tumorType)){
    id1 = paste(tumor, dataset_match$geneSymbol, dataset_match$aaPos, sep="_")
    id2 = paste(dataset$tumorType, dataset$geneSymbol, dataset$aaPos, sep="_")
    dataset_matchWorking <- dataset_match[id1 %in% id2,] #working dataset_match that is tumor specific
    
    if(nrow(dataset_matchWorking) > 0){
      #count number of unique positions matched per gene, and update dMAFcounts
      posN_unique <- aggregate(dataset_matchWorking$aaPos, 
                               by=list(dataset_matchWorking$geneSymbol), 
                               FUN=countfunc)
      lookup_idx <- match(dMAFcounts[dMAFcounts$tumorType == tumor,'geneSymbol'], posN_unique[[1]]) #indices
      lookup_NoNA <- !is.na(lookup_idx)  #true/false list
      dMAFcounts[dMAFcounts$tumorType == tumor,'posNMatched_unique'][lookup_NoNA] <- posN_unique[[2]][lookup_idx[lookup_NoNA]]
    }
  }

  #calculate mutation rates, match rates
  dMAFcounts <- cbind(dMAFcounts,
                      mutRate=dMAFcounts$posN_unique/dMAFcounts$protLength,
                      matchRate=dMAFcounts$posNMatched_unique/dMAFcounts$posN_unique,
                      stringsAsFactors=FALSE)
  
  #save data
  write.csv(dMAFcounts,file=paste(dsetCalc_dir, "dMAFcounts.csv",sep=""))
  
  return(dMAFcounts)

}

calcMultipliers <- function(dMAFcounts){
  #calculate mutation multipliers
  catn("Calculating mutation multipliers...")
  dMultipliers <- data.frame(tumorType=c(),
                             min_mutrate=c(),
                             max_mutrate=c(),
                             avg_mutrate=c(),
                             min_matrate=c(),
                             max_matrate=c(),
                             avg_matrate=c(),
                             multiplier=c(), stringsAsFactors=FALSE)
  for(tumor in unique(dMAFcounts$tumorType)){
    dMAFcounts_tumor <- dMAFcounts[dMAFcounts$tumorType==tumor,]
    mutrates <- dMAFcounts_tumor$mutRate
    matrates <- dMAFcounts_tumor$matchRate
    mlist <- data.frame(tumorType=tumor,
                        min_mutrate=min(mutrates, na.rm=TRUE),
                        max_mutrate=max(mutrates, na.rm=TRUE),
                        avg_mutrate=mean(mutrates, na.rm=TRUE),
                        min_matrate=min(matrates, na.rm=TRUE),
                        max_matrate=max(matrates, na.rm=TRUE),
                        avg_matrate=mean(matrates, na.rm=TRUE),
                        avg_multiplier=mean(matrates, na.rm=TRUE)/mean(mutrates, na.rm=TRUE),
                        stringsAsFactors=FALSE)
    dMultipliers <- rbind(dMultipliers, mlist)
  }
  
  #save data
  write.csv(dMultipliers,file=paste(dsetCalc_dir, "dMultipliers.csv",sep=""))
  
  return(dMultipliers)
}

summaryStats <- function(dMultipliers){
  catn("Generating summary statistics for multipliers...")
  plotfn <- 'hist_multipliers'
  openPDFdev(dsetCalc_dir, plotfn, savePlots)
  hist(dMultipliers$avg_multiplier,breaks=20, main="Distribution of multiplier values")
  saveFig(dsetCalc_dir, plotfn, savePlots)
  dev.off()
}

calcLMfit <- function(dMAFcounts){
  #NOT IN USE
  for(tumor in unique(dMAFcounts$tumorType)){
    dMAFcounts_tumor <- dMAFcounts[dMAFcounts$tumorType==tumor,]
    mutrates <- dMAFcounts_tumor$mutRate
    matrates <- dMAFcounts_tumor$matchRate
    
    plot(mutrates, matrates, xlab="mutation rates", ylab="match rates", main=tumor)
    lfit <- lm(matrates ~ mutrates)
    abline(lfit ,col="red")
    gfit1 <- glm(matrates ~ mutrates, family=gaussian)
    gfit2 <- glm(cbind(dMAFcounts_tumor$posNMatched_unique, 
                      dMAFcounts_tumor$posN_unique-dMAFcounts_tumor$posNMatched_unique) 
                ~ mutrates, family=binomial)
    
  }
}


###################################################
# Main
###################################################
dMAFcounts <- calcCounts()
dMultipliers <- calcMultipliers(dMAFcounts)
summaryStats(dMultipliers)
catn("Saving data...")
save(dMAFcounts, dMultipliers, file=paste(dsetCalc_dir, "dMultipliers.RData",sep=""))
