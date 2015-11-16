#analysis_scoring.R
#Analyze processed datasets: scoring functions
#Author: Boyang Zhao
#REQUIRES (datasets): 1) MutSigCV TCGA files (from TCGA, and ran through MutSigCV)
#                     2) datasets_analyzed.RData (generated from analyze_basic.R)
#                     3) dProteinLengths (protein lengths, generated based on HGNC and UniProt FASTA files)
#                     4) dMultipliers (mutation rate multipliers, generated from calc_multiplier.R)
#MODIFIES: dMAF to filter with RNA-seq expression data (but is not saved)
#GENERATES: sigList, sigListGene
#This R script accepts command line arguments (for specifying which tumor type to run),
#enabling parallelization on computing cluster

if(exists("config.avail") && config.avail){ #if configuration settings exist, keep these, clear the rest of workspace
  rm(list=ls()[grep("^[^config\\.].*$",ls())])
} else {
  rm(list=ls())
}

dset_dir <- "../datasets/"
dsetAnalyzed_dir <- paste(dset_dir, "analyzed/", sep="") #for datasets_analyzed.RData
dsetCalc_dir <- paste(dsetAnalyzed_dir, "calc-expr_null/", sep="") #for dMultipliers.RData
dsetTCGAMutSigCV_dir <- paste(dset_dir, "MutSigCV TCGA/", sep="") #for MutSigCV files
dsetProteome_dir <- paste(dset_dir, "UniProt proteome/", sep="") #for dProteinLengths.RData (if syntheticData=FALSE)
dsetScored_dir <- paste(dset_dir, "/", sep="") #output directory
dsetSyn_dir <- paste(dset_dir, "synthetic/", sep="") #for synthetic_datasets.RData (if syntheticData=TRUE)
source('util.R')
library(boot)

#settings
bootstrapN <- 1e3 #number of boostrap trials
permut_trialsN <- c(1e4,1e6) #number of permutation trials (for null distribution generation); a vector would mean [min max]
useMultiplier <- TRUE #if true, then overwrite mutRate_multiplier using dataset from dMultipliers.RData
mutRate_multiplier <- 10 #default mutation rate multiplier; unless overwritten (if useMultiplier is true)
minCasesNcutoff <- 0 #minimum number cases as cutoff (>=) to write to _sigListGene.csv
minRNAseqcutoff <- 1 #minimum normalized_RSEM values (upper quartile normalized) as cutoff (>=) for mRNA expression
analyzeNotExpressed <- NULL #if NULL, will partition dataset and analyze twice (once on expressed and once on not expressed)
syntheticData <- FALSE #using synthetic dataset instead
tumorTypesToAnalyze <- c() #if empty, analyze all tumor types, otherwise analyze the specified tumor type
workingDatasetName <- "dSNVentryMutExactSame"

#settings overwrite if configuration settings exist
if(exists("config.avail") && config.avail){
  dset_dir <- config.dset_dir
  dsetAnalyzed_dir <- config.dsetAnalyzed_dir
  dsetCalc_dir <- config.dsetCalc_dir
  dsetTCGAMutSigCV_dir <- config.dsetTCGAMutSigCV_dir
  dsetProteome_dir <- config.dsetProteome_dir
  dsetScored_dir <- config.dsetScored_dir
  dsetSyn_dir <- config.dsetSyn_dir
  
  bootstrapN <- config.scoring.bootstrapN
  permut_trialsN <- config.scoring.permut_trialsN
  useMultiplier <- config.scoring.useMultiplier
  mutRate_multiplier <- config.scoring.mutRate_multiplier
  minCasesNcutoff <- config.scoring.minCasesNcutoff
  minRNAseqcutoff <- config.scoring.minRNAseqcutoff
  analyzeNotExpressed <- config.scoring.analyzeNotExpressed
  syntheticData <- config.scoring.syntheticData
  tumorTypesToAnalyze <- config.scoring.tumorTypesToAnalyze
  workingDatasetName <- config.scoring.workingDatasetName
}

#initializing internal variables (DO NOT MODIFY) (prior to loading datasets)
workingDataset.parent <- NULL #this will get defined when reading in datasets below

#load datasets
if(!syntheticData){
  catn('Reading in datasets...')
  if(useMultiplier){load(paste(dsetCalc_dir,"dMultipliers.RData",sep=""))} #load in dMultipliers dataset
  load(paste(dsetProteome_dir, "dProteinLengths.RData",sep="")) #load in dProteinLengths dataset
  load(paste(dsetAnalyzed_dir,"datasets_analyzed.RData",sep="")) #dMAF and dSNV datasets (from analyze_basic.R)
} else {
  catn('Reading in synthetic datasets...')
  load(paste(dsetSyn_dir, "synthetic_dataset.RData",sep="")) #load in synthetic datasets: dMAF, dSNV, and dProteinLengths
}
workingDataset.parent <- eval(as.name(workingDatasetName))

#initializing internal variables (DO NOT MODIFY)
tumorTypes <- unique(dMAF$tumorType) #list of tumor types
tumorTypes <- tumorTypes[tumorTypes!='CELLLINE'] #exclude analysis of cell lines
if(length(tumorTypesToAnalyze) > 0){tumorTypes <- tumorTypes[tumorTypes %in% tumorTypesToAnalyze]}
preval <- 0 #for calculating progress
dMAFs <- c("dMAF") #variable name for dMAF (or a data slice of) to analyze
outputSuffices <- c("") #suffix for results filename
#dMAFS and outputSuffix are to enable the capability to run the analysis multiple times, based on different data slices of dMAF
#dMAFs define the variable name, and before analysis dMAF will be reassigned to the variable variable as specified in dMAFs

#if there is a command line value, then use value as the tumorType index, and only analyze this tumor type
#this provides capability to analyze on computing cluster with separate jobs for each tumor type
args <- commandArgs(TRUE)
if(length(args) > 0){
  catn('Using command line argument...')
  tumorTypes <- tumorTypes[as.numeric(args[1])]
}

#filter out genes not expressed
if(!is.null(minRNAseqcutoff)){
  if(!is.null(analyzeNotExpressed)){
    if(!analyzeNotExpressed){
      dMAF <- dMAF[(!is.na(dMAF$mRNAvals_final) & dMAF$mRNAvals_final >= minRNAseqcutoff) | is.na(dMAF$mRNAvals_final),]
    } else {
      dMAF <- dMAF[dMAF$mRNAvals_final < minRNAseqcutoff,]
      outputSuffices <- c("_notExpressed")
    }
  } else {
    catn("analyzeNotExpressed is set to NULL - will partition dataset to analyze separately RNA-seq expressed and not expressed data...")
    dMAFExpr <- dMAF[(!is.na(dMAF$mRNAvals_final) & dMAF$mRNAvals_final >= minRNAseqcutoff) | is.na(dMAF$mRNAvals_final),]
    dMAFNotEpr <- dMAF[dMAF$mRNAvals_final < minRNAseqcutoff,]
    
    dMAFs <- c("dMAFExpr", "dMAFNotEpr")
    outputSuffices <- c("", "_notExpressed")
  }
}
#NOTE: cannot just do dMAF[dMAF$mRNAvals_final >= 1,], this will add NAs to all columns for entries where mRNAvals_final are NAs
#Need to get the entries with and without NAs separately; requires guaranteed evaluation of operands left to right 

#quality controls
if(is.null(minRNAseqcutoff) && is.null(analyzeNotExpressed)){
  message("WARNING: both minRNAseqcutoff and analyzeNotExpressed are NULL. Since mininum RNAseq cutoff is not defined, cannot 
        partition dataset into expressed/not expressed; analysis will be done on entire dataset...")
}

createDir(dsetScored_dir) #create output dir if needed

###################################################
# Helper methods for statistical analyses
###################################################
calcPVal <- function(scores, refScore){
  #calculates one-sided p-val given distribution of scores, and reference score (refScore)
  pval <- sum(scores>=refScore)/length(scores)
  return(pval)
}

removeNotMatchedPosMut <- function(dataset, tumor){
  #remove any entries (from input dataset) where the residue positions + mutation are not found in dMAF of given tumor type
  id1 = paste(tumor, dataset$geneSymbol, dataset$aaPos, dataset$varAA, sep="_")
  id2 = paste(dMAF$tumorType, dMAF$geneSymbol, dMAF$aaPos, dMAF$varAA, sep="_")
  newdataset <- dataset[id1 %in% id2,]
  return(newdataset)
}

removeNotMatchedPos <- function(dataset, tumor){
  #remove any entries (from input dataset) where the residue positions are not found in dMAF of given tumor type
  id1 = paste(tumor, dataset$geneSymbol, dataset$aaPos, sep="_")
  id2 = paste(dMAF$tumorType, dMAF$geneSymbol, dMAF$aaPos, sep="_")
  newdataset <- dataset[id1 %in% id2,]
  return(newdataset)
}

removeMatchedPos <- function(dataset, refdataset){
  #remove any entries (from input dataset) with gene_pos also found in refdataset
  id1 = paste(tumor, dataset$geneSymbol, dataset$aaPos, sep="_")
  id2 = paste(refdataset$tumorType, refdataset$geneSymbol, refdataset$aaPos, sep="_")
  newdataset <- dataset[!(id1 %in% id2),]
  return(newdataset)
}

populateScores_perGene <- function(scores){
  #DEPRECATED
  #populate scores to all positions for each gene, given a gene-score pair
  #returns a list the same size as workingDataset
  #INPUT: scores is a list (scores and geneNames should be the same size)
  #DEPENDENCIES: workingDataset
  #REQUIRES: order need to be perserved between geneNames list (from unique method) and scores (from methods that 
  #          used the workingDatasetGene dataset as input; workingDatasetGene is generated using !duplicated approach)
  #          this was checked to be the case -> both duplicate and unique methods are stable and preserve order
  
  geneNames <- unique(workingDataset$geneSymbol)
  resultList <- sapply(1:nrow(workingDataset), function(idx){
    gene <- workingDataset[idx,'geneSymbol']
    scores[which(geneNames == gene)]
  })
  return(resultList)
}

populateScores_perGenePos <- function(scores){
  #DEPRECATED
  #populate scores to all positions for each gene, given a gene_pos-score pair
  #returns a list the same size as workingDataset
  idRef <- unique(paste(workingDataset$geneSymbol,workingDataset$aaPos,sep="_"))
  resultList <- sapply(1:nrow(workingDataset), function(idx){
    id <- paste(workingDataset[idx,'geneSymbol'], workingDataset[idx,'aaPos'], sep="_")
    scores[which(idRef == id)]
  })
  return(resultList)
}

getMutationsN <- function(caseID, tumor, gene){
  #get the number of mutations given tumor type, gene name, and case ID
  #DEPENDENCIES: dMAF
  mutationsN <- sum(dMAF$geneSymbol == gene & dMAF$tumorType == tumor & dMAF$caseID == caseID)
  return(mutationsN)
}

getCaseIDs <- function(tumor, gene, returnN = FALSE, asList = FALSE){
  #get case IDs for given gene and tumor type
  #if returnN is true, then return the number of cases, instead of list of case IDs
  #DEPENDENCIES: dMAF, tumor
  caseIDs <- unique(dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor,'caseID'])
  
  if(returnN){
    return(length(caseIDs))
  } else {
    if(asList){
      return(paste(caseIDs, collapse="|"))
    } else {
      return(caseIDs)
    }
  }
}

getCasesN <- function(idx, dataset){
  #DEPENDENCIES: tumor, preval
  calcProgress(idx, nrow(dataset), preval)
  gene <- dataset[idx, 'geneSymbol']
  casesN <- getCaseIDs(tumor, gene, returnN = TRUE)
  return(casesN)
}

getCasesList <- function(idx, dataset){
  #DEPENDENCIES: tumor, preval
  calcProgress(idx, nrow(dataset), preval)
  gene <- dataset[idx, 'geneSymbol']
  casesList <- getCaseIDs(tumor, gene, asList = TRUE)
  return(casesList)
}

getProtLength <- function(gene){
  #DEPENDENCIES: dProteinLengths
  if(sum(dProteinLengths$geneSymbol == gene) < 1){
    message(gettextf("WARNING: Cannot find protein length for %s", gene))
    return(NA)
  } else {
    return(dProteinLengths[dProteinLengths$geneSymbol == gene,'length'])
  }
}

###################################################
# Scoring methods
###################################################

###################################################
# Method based on permutation (per patient sample) to generate null distribution (DEPRECATED)
calcScore_permtPercase_single <- function(tumor, gene, pos){
  #DEPRECATED
  #calculate the fraction of samples where the position match
  #total <- length(unique(dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor, 'caseID']))
  total <- length(unique(dMAF[dMAF$geneSymbol == gene, 'caseID']))
  sameMatch <- length(unique(dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor & dMAF$aaPos == pos, 'caseID']))
  score <- sameMatch/total
  
  return(score)
}

permtPercaseNullDistr <- function(mutationsN_list, refDistr, refPos, trials_min, trials_max, refScore=NA, matchMultiple=FALSE){
  #DEPRECATED
  #generate the null distribution using permutation
  #refDistr: reference distribution, e.g. any positions on protein, or any positions of ones mutated in given tumor
  #refScore is NA if analyzing based on aggregated samples
  #refPos can be a vector of multiple positions, or just a single position
  scores <- c()
  for(idx in 1:trials_max){
    totalCasesN = length(mutationsN_list)
    sameMatch <- 0
    for(mutationsN in mutationsN_list){
      #for each case, generate a random mutation pattern
      if(protlength < mutationsN){
        message('WARNING: protlength < mutationsN')
        break
      }
      randPos <- sample(refDistr, mutationsN)
      
      refPos<-as.integer(refPos) #make sure all are integer type, before making the comparison
      matchesN <- sum(refPos %in% randPos)
      if(matchesN > 0){
        if(!matchMultiple){
          #if matchMultiple is false, then score is to just add 1
          sameMatch <- sameMatch + 1
        } else {
          #other, score is equal to the number of matches
          sameMatch <- matchesN
        }
      }
    }
    
    score <- sameMatch/totalCasesN
    scores <- c(scores, score)
    
    if(!is.na(refScore)){
      if(idx > trials_min && sum(refScore < scores) >= 10){
        #if the number of trial is greater than the minimum, and there are at least 10 scores from
        #randomly generated mutation profiles that have values equal to or greater than the observed reference score
        #then can the random sampling
        break
      }
    }
  }
  
  return(scores)
}

calcScore_permtPercase <- function(idx, dataset){
  #DEPRECATED
  #DEPENDENCIES: dProteinLengths, tumor, preval
  calcProgress(idx, nrow(dataset), preval)
  
  gene <- dataset[idx, 'geneSymbol']
  caseIDs <- getCaseIDs(tumor, gene)
  
  pval <- NA
  if(length(caseIDs) >= 2){ #need at least two cases
    
    #get the list of number of mutations, per case ID (each row a different case ID)
    mutationsN_caseID <- vapply(caseIDs, getMutationsN, 0, tumor=tumor, gene=gene)
    
    #per sample, get all the mutations at single position
    refScore <- calcScore_permtPercase_single(tumor, gene, dataset[idx,'aaPos'])
    posToCheck <- dataset[idx,'aaPos']
    
    #generate null distribution
    protLength <- getProtLength(gene)
    scores <- permtPercaseNullDistr(mutationsN_caseID, protLength, posToCheck, 1e4, 1e5, refScore)
  
    #calculate p-value
    pval <- calcPVal(scores, refScore)
  }
  
  return(pval)
}

###################################################
# Method based on probability of occurence from distribution of mutations observed
calcProbScore_byEntries <- function(refDistr, refPos){
  #calculates point probability score, given reference distribution
  #refDistr and refPos does not have to be unique
  distr_leftovers <- refDistr
  for(ref in refPos){
    idx <- which(distr_leftovers==ref)[1]
    distr_leftovers <- distr_leftovers[1:length(distr_leftovers) != idx]
  }
  countcomplement <- length(distr_leftovers)
  counts <- length(refDistr) - countcomplement
  prob <- counts/length(refDistr)
  
  return(prob)
}

calcProbScore_byPos <- function(refDistr, refPos){
  #calculates point probability score based on unique positions, given reference distribution
  refPos <- unique(refPos)
  counts <- sum(refDistr %in% refPos)
  prob <- counts/length(refDistr)
  
  return(prob)
}

probSD_boot_func <- function(distr, indices, refPos_unique){
  #calculates point probability given distribution and ref positions, using calcProbScore_byPos
  #DEPENDENCIES: called by boot in probSD_boot (method); uses calcProbScore_byPos (method)
  return(calcProbScore_byPos(distr[indices], refPos_unique))
}

probSD_boot <- function(distr, refPos_unique, N){
  #nonparametric bootstrap, and calculate standard deviation
  #DEPENDENCIES: called by mutSig_prob_single (method)
  bresults <- boot(distr, probSD_boot_func, R=N, refPos_unique=refPos_unique)
  probsd <- sd(bresults$t[,1])
  return(probsd)
}

calcProbScoreSD <- function(gene, refPos_unique, bootstrap=TRUE){
  #calculates probability score and score standard deviation (if boostrap is TRUE)
  #DEPENDENCIES: dMAF, tumor, bootstrapN, calcProbScore_byPos (method)
  
  #calculate point probability
  distr <- dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor,'aaPos']
  if(length(distr) < 1){message(gettextf("WARNING: no mutations found for %s, tumor type %s", gene, tumor))}
  prob <- calcProbScore_byPos(distr, refPos_unique)
  
  #calculate standard devation by bootstrapping
  probsd <- NA
  if(bootstrap){
    probsd <- probSD_boot(distr, refPos_unique, bootstrapN)
  }
  
  return(list(prob=prob,sd=probsd))
}

probNullDistr <- function(refDistr, refPos, trialsN, protLength){
  #generate the null distribution using permutation-based method
  #refDistr: reference distribution, e.g. any positions on protein, or any positions of ones mutated in given tumor
  #refPos can be a vector of multiple positions, or just a single position
  #DEPENDENCIES: tumor
  
  #quality controls
  if(length(refDistr) < 1 || length(refPos) < 1 || length(trialsN) < 1 || length(protLength) < 1){
    message('ERROR: one or more of inputs into mutSig_permtNulldistr is less than 1...')
    return(NULL)
  }
  
  #var definitions/changes
  if(exists("dMultipliers")){
    #update mutation rate multiplier is applicable
    mutRate_multiplier <- dMultipliers[dMultipliers$tumorType==tumor,"avg_multiplier"]
  }
  
  if(length(trialsN) == 2){
    #extract minimum and maximum number of trials for permutation
    trialsN_min <- trialsN[1]
    trialsN_max <- trialsN[2]
  } else {
    trialsN_max <- trialsN
    trialsN_min <- trialsN_max
  }
  
  #make sure all are integer type, before making the comparison
  refPos <- as.integer(refPos)
  refPos_unique <- unique(refPos)
  refDistr <- as.integer(refDistr)
  refDistr_unique <- unique(refDistr)
  mutationsN <- length(refPos)
  mutationsN_unique <- length(refPos_unique)
  
  mutRate <- length(refDistr_unique)/protLength
  
  #quality controls for mutRate_multiplier
  if(mutRate_multiplier*mutRate > 1){
    message(gettextf("WARNING: multiplier correction leads to mutation rate exceeding 1 for tumor %s, 
                   using corrected mutation rate of 1 instead...", tumor))
    mutRate <- 1
  } else if(mutRate_multiplier == 0){
    message(gettextf("WARNING: mutRate_multiplier for tumor %s is 0, using default multiplier value instead...", tumor))
  } else {
    mutRate <- mutRate_multiplier*mutRate
  }
  
  if(length(refDistr_unique) < mutationsN_unique){
    message('WARNING: length of reference mutations is less than number of matched mutations...')
  }
  
  refScore <- calcProbScore_byPos(refDistr, refPos_unique)
  scores <- c()
  
  #generate null distribution
  idx_firstPass <- -1 #used to keep track of the trial number when minimum sampling has been first met
  for(idx in 1:trialsN_max){
    #method 1a: sample positions based on uniform distribution from reference distribution
    #     if(length(refDistr) < 2){ #there is only one entry
    #       randPos <- refDistr
    #     } else {
    #       randPos <- sample(refDistr, mutationsN, replace=FALSE)
    #     }
    #method 1b: sample unique positions based on uniform distribution from reference distribution
    #     if(length(refDistr_unique) < 2){ #there is only one entry
    #       randPos <- refDistr_unique
    #     } else {
    #       randPos <- sample(refDistr_unique, mutationsN_unique, replace=FALSE)
    #     }
    
    #method 2: sample of match/nomatch based on binomial distribution (with prob ~ mutation rate)
    #then for each match, sample based on uniform distribution the unique positions

    matches <- sample(c(0,1), mutationsN_unique, prob = c(1-mutRate, mutRate), replace=TRUE)
    randPos <- rep(-1, mutationsN_unique)
    matchesN <- sum(matches)
    if(matchesN > 0){
      if(length(refDistr_unique) < 2){ #there is only one unique entry
        randPos[1:matchesN] <- rep(refDistr_unique, matchesN)
      } else {
        randPos[1:matchesN] <- sample(refDistr_unique, matchesN, replace=FALSE)
      }
    }
    
    #method 3: sample of match/nomatch based on binomial distribution (with prob ~ mutation rate)
    #then for each match, sample based on uniform distribution the positions
    #     matches <- sample(c(0,1), mutationsN, prob = c(1-mutRate, mutRate), replace=TRUE)
    #     randPos <- rep(-1, mutationsN)
    #     matchesN <- sum(matches)
    #     if(matchesN > 0){
    #       if(length(refDistr) < 2){ #there is only one entry
    #         randPos[1:matchesN] <- rep(refDistr, matchesN)
    #       } else {
    #         randPos[1:matchesN] <- sample(refDistr, matchesN, replace=FALSE)
    #       }
    #     }
    
    randPos_unique <- unique(randPos)
    score <- calcProbScore_byPos(refDistr, randPos_unique)
    scores <- c(scores, score)
    
    #after minimum number of trials have been met, check if this is sufficient sampling or more trials are needed
    if(idx > trialsN_min){
      pvalTmp <- calcPVal(scores, refScore)
      if(pvalTmp > 0){
        if(idx_firstPass == -1){
          idx_firstPass <- idx
        } else {
          if((idx-idx_firstPass) >= 1e4){
            #stop sampling if there has been at least 1e4 more trials since the criteria has been initially met
            break
          }
        }
      }
    }
  }
  
  pval <- calcPVal(scores, refScore)
  return(pval)
}

calcProbScore_pval <- function(gene, refPos_unique, perGene=FALSE){
  #calculates p-value
  #DEPENDENCIES: tumor, workingDataset, dMAF
  
  #get protein length
  protLength <- getProtLength(gene)
  
  #get reference disribution and positions
  #refDistr <- dProteinLengths[dProteinLengths$geneSymbol == gene, 'length'] #all positions across protein
  #refDistr <- unique(dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor,'aaPos']) #or unique mutated positions across protein
  refDistr <- dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor,'aaPos'] #or mutated positions across protein
  
  #expand refPos as a distribution
  refPos <- refDistr[refDistr %in% refPos_unique]
  
  #generate null distribution and calculate p-value
  pval <- probNullDistr(refDistr, refPos, permut_trialsN, protLength)
  
  return(pval)
}

calcProbStat <- function(idx, dataset, perGene=FALSE, calcPVal=FALSE){
  #main probability scoring method: calculates probability score, standard deviation, SNR, and p-value
  #DEPENDENCIES: tumor, preval, dMAF, workingDataset, workingDataset.parent
  calcProgress(idx, nrow(dataset), preval)
  
  gene <- dataset[idx, 'geneSymbol']
  refPos <- dataset[idx, 'aaPos']
  if(perGene){
    refPos <- workingDataset[workingDataset$geneSymbol == gene, 'aaPos']
  }
  refPos_unique <- unique(refPos)
  
  #calculate prob and probSD for given gene and refPos_unique
  val <- calcProbScoreSD(gene, refPos_unique)
  SNR <- val$prob/val$sd
  
  #generate null distribution and calculate p-value
  pval_prob <- NA
  if(calcPVal){
    pval_prob <- calcProbScore_pval(gene, refPos_unique, perGene)
  }
  
  #additional values
  uniquePosNdMAF <- length(unique(dMAF[dMAF$geneSymbol == gene, 'aaPos']))
  uniquePosNdMAFtumor <- length(unique(dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor,'aaPos']))
  uniquePosNExactMatch <- length(unique(workingDataset.parent[workingDataset.parent$geneSymbol == gene, 'aaPos']))
  uniquePosNExactMatchtumor <- length(unique(workingDataset[workingDataset$geneSymbol == gene, 'aaPos']))
  
  aaNotPos <- unique(dMAF[dMAF$geneSymbol == gene, 'aaPos'])[!(unique(dMAF[dMAF$geneSymbol == gene, 'aaPos']) %in% 
                  unique(dMAF[dMAF$geneSymbol == gene & dMAF$tumorType == tumor,'aaPos']))]
  exactMatchNonTumor <- aaNotPos %in% unique(workingDataset.parent[workingDataset.parent$geneSymbol == gene, 'aaPos'])
  
  contigMatrix <- c(uniquePosNExactMatchtumor, length(exactMatchNonTumor),
                    uniquePosNdMAFtumor-uniquePosNExactMatchtumor, length(aaNotPos)-length(exactMatchNonTumor))
  contigMatrix <- matrix(contigMatrix, nrow=2, ncol=2)
  FET <- fisher.test(contigMatrix)
  
  return(list(prob=val$prob,
              sd=val$sd,
              SNR=SNR,
              pval_prob=pval_prob,
              uniquePosNdMAF=uniquePosNdMAF,
              uniquePosNdMAFtumor=uniquePosNdMAFtumor,
              uniquePosNExactMatch=uniquePosNExactMatch,
              uniquePosNExactMatchtumor=uniquePosNExactMatchtumor,
              pval_FET=FET$p.value))
}

###################################################
# Main analysis method
###################################################
analyzeDataset <- function(workingDataset, workingDataset.parent){
  #analyze working dataset: not guaranteed to have unique gene and/or gene_pos
  uniqueGene_idx <- which(!duplicated(workingDataset$geneSymbol))
  uniqueGenePos_idx <- which(!duplicated(paste(workingDataset$geneSymbol,workingDataset$aaPos, sep="_")))
  
  #error checks
  if(nrow(workingDataset) < 1){
    message(gettextf("WARNING: working dataset for tumor %s is empty. Skipping analysis.", tumor))
    return(NULL)
  }
  
  #subsets of workingDatasetGene with unique gene, or unique gene_pos
  #subset dataset still contains all fields, but some are not relevant (e.g. positions in workingDatasetGene)
  workingDatasetGene <- workingDataset[uniqueGene_idx,]
  workingDatasetGenePos <- workingDataset[uniqueGenePos_idx,]
  
  #get indices to expand values (analyzed using workingDatasetGene or workingDatasetGenePos) to fill to workingDataset
  geneExpand_idx <- match(workingDataset$geneSymbol, workingDatasetGene$geneSymbol)
  genePosExpand_idx <- match(paste(workingDataset$geneSymbol,workingDataset$aaPos,sep="_"), 
                             paste(workingDatasetGenePos$geneSymbol,workingDatasetGenePos$aaPos,sep="_"))
  
  #per position analysis
  catn("Calculating probability score, SDs, and p-values for each position")
  preval <<- 0
  probResults <- lapply(1:nrow(workingDatasetGenePos), calcProbStat, 
                        dataset=workingDatasetGenePos, perGene=FALSE, calcPVal=FALSE)
  probResults <- unlist(probResults)
  prob <- probResults[names(probResults) == 'prob']
  probSD <- probResults[names(probResults) == 'sd']
  
  #per gene analysis
  catn("Calculating probability score, SDs, and p-values for each gene")
  preval <<- 0
  probGeneResults <- lapply(1:nrow(workingDatasetGene), calcProbStat, 
                            dataset=workingDatasetGene, perGene=TRUE, calcPVal=TRUE)
  probGeneResults <- unlist(probGeneResults)
  probGene <- probGeneResults[names(probGeneResults) == 'prob']
  probGeneSD <- probGeneResults[names(probGeneResults) == 'sd']
  probGeneSNR <- probGeneResults[names(probGeneResults) == 'SNR']
  pval_probGene <- probGeneResults[names(probGeneResults) == 'pval_prob']
  qval_probGene <- p.adjust(pval_probGene, method="BH")
  uniquePosNdMAF <- probGeneResults[names(probGeneResults) == 'uniquePosNdMAF']
  uniquePosNdMAFtumor <- probGeneResults[names(probGeneResults) == 'uniquePosNdMAFtumor']
  uniquePosNExactMatch <- probGeneResults[names(probGeneResults) == 'uniquePosNExactMatch']
  uniquePosNExactMatchtumor <- probGeneResults[names(probGeneResults) == 'uniquePosNExactMatchtumor']
  pval_FET <- probGeneResults[names(probGeneResults) == 'pval_FET']
  qval_FET <- p.adjust(pval_FET, method="BH")
  
  #get p_val from MutSigCV
  catn("Acquiring p_val from MutSigCV...")
  mutsigcv_fname <- paste(dsetTCGAMutSigCV_dir, tumor2filename(tumor),".maf_output.sig_genes.txt",sep="")
  #note: for COAD/READ, TCGA has separate datasets for each, but cBioPortal (and thus dMAF) only contains the grouped COAD/READ
  #the TCGA MutSigCV for COAD and READ was combined (and by taking the minimum q-value of the two), and was what was used here.
  if(file.exists(mutsigcv_fname)){
    dMAFMutSigCV <- read.table(mutsigcv_fname,header=TRUE,stringsAsFactors=FALSE)
    indices <- match(workingDatasetGene$geneSymbol, dMAFMutSigCV$gene)
    pval_MutSig <- dMAFMutSigCV[indices, 'p']
    qval_MutSig <- dMAFMutSigCV[indices, 'q']
  } else {
    message(gettextf("WARNING: MutsigCV file %s does not exist.", mutsigcv_fname))
    pval_MutSig <- rep(NA, nrow(workingDatasetGene))
    qval_MutSig <- pval_MutSig
  }
  
  #get number of cases
  catn("Acquiring number of cases...")
  preval <<- 0
  casesN <- vapply(1:nrow(workingDatasetGene), getCasesN, 0, dataset=workingDatasetGene)
  
  #catn("Getting case IDs list...")
  #preval <- 0
  #casesList <- vapply(1:nrow(workingDatasetGene), getCasesList, "", dataset=workingDatasetGene)
  #casesList <- populateScores_perGene(casesList)
  
  #saving results
  #per gene_position
  protLength <- vapply(workingDatasetGene$geneSymbol,getProtLength,0)
  protLength_e <- populateScores_perGene(protLength)
  
  sigList <- cbind(workingDataset,
                   prob=prob[genePosExpand_idx],
                   probSD=probSD[genePosExpand_idx],
                   probGene=probGene[geneExpand_idx],
                   probGeneSD=probGeneSD[geneExpand_idx],
                   SNR=probGeneSNR[geneExpand_idx],
                   pval_prob=pval_probGene[geneExpand_idx],
                   qval_prob=qval_probGene[geneExpand_idx],
                   pval_MutSig=pval_MutSig[geneExpand_idx],
                   qval_MutSig=qval_MutSig[geneExpand_idx],
                   protLength=protLength[geneExpand_idx],
                   uniquePosNdMAF=uniquePosNdMAF[geneExpand_idx],
                   uniquePosNdMAFtumor=uniquePosNdMAFtumor[geneExpand_idx],
                   uniquePosNExactMatch=uniquePosNExactMatch[geneExpand_idx],
                   uniquePosNExactMatchtumor=uniquePosNExactMatchtumor[geneExpand_idx],
                   pval_FET=pval_FET[geneExpand_idx],
                   qval_FET=qval_FET[geneExpand_idx],
                   casesN=casesN[geneExpand_idx])
  
  #per unique gene and remove genes with only one reported case
  sigListGene <- cbind(probGene=probGene,
                       probGeneSD=probGeneSD,
                       SNR=probGeneSNR,
                       pval_probGene=pval_probGene,
                       qval_probGene=qval_probGene,
                       pval_MutSig=pval_MutSig,
                       qval_MutSig=qval_MutSig,
                       protLength=protLength,
                       uniquePosNdMAF=uniquePosNdMAF,
                       uniquePosNdMAFtumor=uniquePosNdMAFtumor,
                       uniquePosNExactMatch=uniquePosNExactMatch,
                       uniquePosNExactMatchtumor=uniquePosNExactMatchtumor,
                       pval_FET=pval_FET,
                       qval_FET=qval_FET,
                       casesN=casesN)
  sigListGene <- as.data.frame(sigListGene, stringsAsFactors=FALSE)
  names(sigListGene) <- c('prob',
                          'probSD',
                          'SNR',
                          'pval_prob',
                          'qval_prob',
                          'pval_MutSig',
                          'qval_MutSig',
                          'protLength',
                          'uniquePosNdMAF',
                          'uniquePosNdMAFtumor',
                          'uniquePosNExactMatch',
                          'uniquePosNExactMatchtumor',
                          'pval_FET',
                          'qval_FET',
                          'casesN')
  row.names(sigListGene) <- workingDatasetGene$geneSymbol
  sigListGene <- sigListGene[sigListGene$casesN >= minCasesNcutoff,]
  
  catn("Writing results...")
  tryCatch({
    write.csv(sigList,file=paste(dsetScored_dir, tumor2filename(tumor), "_sigList",outputSuffix,".csv",sep=""))
    write.csv(sigListGene,file=paste(dsetScored_dir, tumor2filename(tumor), "_sigListGene",outputSuffix,".csv",sep=""))
    save(sigList, sigListGene, file=paste(dsetScored_dir, tumor2filename(tumor), "_sigList",outputSuffix,".RData",sep=""))
  }, warning = function(w){
    message(conditionMessage(w))
  }, error = function(err){
    message(conditionMessage(err))
  }, finally = {
  })
}

###################################################
# Analyze by tumor type
###################################################
for(idx in 1:length(dMAFs)){
  dMAF <- get(dMAFs[idx])
  outputSuffix <- outputSuffices[idx]
  
  for(tumor in tumorTypes){
    catn(gettextf("*** Analyzing %s %s***", tumor, outputSuffix))
    
    #define working dataset, and clean up dataset
    catn("Creating working dataset...")
    #create working dataset with exact gene_pos match between dMAF and dSNV
    workingDataset <- removeNotMatchedPos(workingDataset.parent, tumor)
    analyzeDataset(workingDataset, workingDataset.parent)
  }
}
