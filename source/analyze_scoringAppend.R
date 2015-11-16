#analysis_scoredAppend.R
#Append the scored list with significant MutSigCV hits that is not in the original overlay list
#Author: Boyang Zhao
#REQUIRES (datasets): 1) *_sigList.RData (generated from analyze_scoring.R)
#                     2) list of tumor types (from dMAF)
#GENERATES: sigList, sigListGene

if(exists("config.avail") && config.avail){ #if configuration settings exist, keep these, clear the rest of workspace
  rm(list=ls()[grep("^[^config\\.].*$",ls())])
} else {
  rm(list=ls())
}

dset_dir <- "../datasets/"
dsetTCGAMutSigCV_dir <- paste(dset_dir, "MutSigCV TCGA/", sep="") #for MutSigCV files
dsetAnalyzed_dir <- paste(dset_dir, "analyzed/", sep="") #for dMAF dataset
dsetScored_dir <- paste(dset_dir, "scored/", sep="") #both input and output directory
source("util.R")

#settings
MutSigCutoff <- 0.05 #minimum cutoff for MutSigCV qval (>=)
minRNAseqcutoff <- 1 #minimum normalized_RSEM values (upper quartile normalized) as cutoff (>=) for mRNA expression
analyzePairedExprData <- TRUE #for now always set to TRUE
tumorTypesToAnalyze <- c() #if empty, analyze all tumor types, otherwise analyze the specified tumor type

#quality controls
checkDirExists(dsetAnalyzed_dir)
checkDirExists(dsetScored_dir)

#load dMAF -> for getting the expression values
load(paste(dsetAnalyzed_dir,"datasets_analyzed.RData",sep=""))

#initializing internal variables (DO NOT MODIFY)
tumorTypes <- unique(dMAF$tumorType) #list of tumor types
tumorTypes <- tumorTypes[tumorTypes!='CELLLINE'] #exclude analysis of cell lines
if(length(tumorTypesToAnalyze) > 0){tumorTypes <- tumorTypes[tumorTypes %in% tumorTypesToAnalyze]}

###################################################
# Helper methods
###################################################
createAppendTable <- function(genesToAppend, dMAFMutSigCV){
  #given the dMAFMutSigCV dataset, and genes to append, pull out the create the data.frame that conforms with
  #format of sigListGene
  
  #expressed
  blankArray = rep(NA, length(genesToAppend))
  tableAppend <- data.frame(prob=blankArray,
                            probSD=blankArray,
                            SNR=blankArray,
                            pval_prob=blankArray,
                            qval_prob=blankArray,
                            pval_MutSig=dMAFMutSigCV[dMAFMutSigCV$gene %in% genesToAppend, 'p'],
                            qval_MutSig=dMAFMutSigCV[dMAFMutSigCV$gene %in% genesToAppend, 'q'],
                            protLength=blankArray,
                            uniquePosNdMAF=blankArray,
                            uniquePosNdMAFtumor=blankArray,
                            uniquePosNExactMatch=blankArray,
                            uniquePosNExactMatchtumor=blankArray,
                            pval_FET=blankArray,
                            qval_FET=blankArray,
                            casesN=blankArray)
  
  row.names(tableAppend) <- dMAFMutSigCV[dMAFMutSigCV$gene %in% genesToAppend, 'gene']
  
  return(tableAppend)
}

#get statistically significant MutSigCV genes
appendMutSigCVgenes <- function(tumor, sigListGene, criteriaExprPass){
  #Acquire MutSigCV genes that are significant and pass expression criteria
  #if criteriaExprPass is true, then filter as > minRNAseqcutoff, otherwise < minRNAseqcutoff

  #get p_val from MutSigCV
  mutsigcv_fname <- paste(dsetTCGAMutSigCV_dir, tumor2filename(tumor),".maf_output.sig_genes.txt",sep="")
  #note: for COAD/READ, TCGA has separate datasets for each, but cBioPortal (and thus dMAF) only contains the grouped COAD/READ
  #the TCGA MutSigCV for COAD and READ was combined (and by taking the minimum q-value of the two), and was what was used here.
  
  if(file.exists(mutsigcv_fname)){
    ######################
    # reading in datasets
    #read in MutSigCV
    dMAFMutSigCV <- read.table(mutsigcv_fname, header=TRUE, stringsAsFactors=FALSE)
    
    ######################
    # filters
    #filter for tumor type
    tumorIdx <- dMAF$tumorType == tumor

    #filter based on expression criteria
    if(criteriaExprPass){ #filter for expressed genes
      exprIdx <- (!is.na(dMAF$mRNAvals_final) & dMAF$mRNAvals_final >= minRNAseqcutoff) | is.na(dMAF$mRNAvals_final)
    } else { #filter for not expressed genes
      exprIdx <- dMAF$mRNAvals_final < minRNAseqcutoff
    }
    
    passExprTumor <- tumorIdx & exprIdx
    genes.passExprTumor <- dMAF[passExprTumor, 'geneSymbol']

    #filter MutSigCV
    genes.MutSigCV <- dMAFMutSigCV[dMAFMutSigCV$q <= MutSigCutoff, 'gene']
    
    #intersect, genes for tumor type that passed expression and MutSigCV significance
    MutSigCV.passedExpr <- genes.MutSigCV %in% genes.passExprTumor
    genes.MutSigCV.passedExpr <- genes.MutSigCV[MutSigCV.passedExpr]
    
    ######################
    # append
    genes.MutSigCV.add <- genes.MutSigCV.passedExpr[!(genes.MutSigCV.passedExpr %in% row.names(sigListGene))]
    
    t <- createAppendTable(genes.MutSigCV.add, dMAFMutSigCV)
    sigListGene.appended <- rbind(sigListGene, t)
    return(sigListGene.appended)
  } else {
    message(gettextf("WARNING: MutsigCV file %s does not exist.", mutsigcv_fname))
  }
  
  return(NULL)
}

###################################################
# Append by tumor type
###################################################
for(tumor in tumorTypes){
  catn(gettextf("*** Analyzing %s ***", tumor))
  
  #read in datasets (sigListGene in expressed and not expressed)
  if(analyzePairedExprData){
    #expressed
    outputSuffix <- ""
    sigList_filename <- paste(dsetScored_dir, tumor2filename(tumor), "_sigList.RData",sep="")
    if(file.exists(sigList_filename)){
      load(sigList_filename)
      sigListGene.appended <- appendMutSigCVgenes(tumor, sigListGene, TRUE) #note sigListGene.appended can be NULL if it is empty
      
      catn(gettextf("Expressed: original: %d; new %d",nrow(sigListGene),nrow(sigListGene.appended)))
      
      write.csv(sigListGene.appended, file=paste(dsetScored_dir, tumor2filename(tumor), "_sigListGeneAppended", outputSuffix, ".csv", sep=""))
      save(sigList, sigListGene, sigListGene.appended, file=paste(dsetScored_dir, tumor2filename(tumor), "_sigList",outputSuffix,".RData",sep=""))
    } else {
      message(gettextf("WARNING: sigList file %s does not exist.", sigList_filename))
    }
    
    #not expressed
    outputSuffix <- "_notExpressed"
    sigListNotExpr_filename <- paste(dsetScored_dir, tumor2filename(tumor), "_sigList_notExpressed.RData",sep="")
    if(file.exists(sigListNotExpr_filename)){
      load(sigListNotExpr_filename)
      sigListGene.appended <- appendMutSigCVgenes(tumor, sigListGene, FALSE)
      
      catn(gettextf("Not expressed: original: %d; new %d",nrow(sigListGene),nrow(sigListGene.appended)))
      
      write.csv(sigListGene.appended, file=paste(dsetScored_dir, tumor2filename(tumor), "_sigListGeneAppended", outputSuffix, ".csv", sep=""))
      save(sigList, sigListGene, sigListGene.appended, file=paste(dsetScored_dir, tumor2filename(tumor), "_sigList",outputSuffix,".RData",sep=""))
    } else {
      message(gettextf("WARNING: sigList file %s does not exist.", sigListNotExpr_filename))
    }
  }
  
}
