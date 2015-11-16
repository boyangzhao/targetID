#analysis_scored.R
#Analyze scored datasets (i.e. results from analyze_scoring.R)
#Author: Boyang Zhao
#REQUIRES (datasets): 1) *_sigList.RData (generated from analyze_scoring.R)
#                     2) list of tumor types (from dMAF)
#GENERATES: in 'dsetScoredSummary_dir' directory:
#               1) sigListGene.csv
#               2) 'siglist_all' in sigListGene.RData (aggregates and appends results with boolean pass filters)
#           in 'plots_dir' directory:
#               1) [plots pdf/png]
#MODIFIES: in dealing with NAs, Inf, etc in the dataset, the following are used in cleanups
#          1) for heatmaps and PCAs, a reduced dataset where any observation with at least NAs is removed
#          2) in preprocessing, for probSD is 0, this will yield a SNR of Inf. The SNR is recalculated using the min probSD
#             in the entire dataset (NOT IN USED RIGHT NOW)
#          3) for any entries with MutSigCV q-value (qval_MutSig) is NA. This is converted to 1.

if(exists("config.avail") && config.avail){ #if configuration settings exist, keep these, clear the rest of workspace
  rm(list=ls()[grep("^[^config\\.].*$",ls())])
} else {
  rm(list=ls())
}

dset_dir <- "../datasets/"
dsetAnalyzed_dir <- paste(dset_dir, "analyzed/", sep="") #for dMAF
dsetScored_dir <- paste(dset_dir, "scored/", sep="") #scored dataset
dsetScoredSummary_dir <- gettextf('%ssummary/summary_pval0.062_exprONLY_090415/', dsetScored_dir) #output directory
plots_dir <- gettextf('%splotsIndCancer_pval0.062_exprONLY_090415/', dsetScored_dir) #output directory
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
source("util.R")

#settings
savePlots <- TRUE
saveLog <- TRUE
minCasesNcutoff <- 2 #NULL means use all of given data; or impose additional cutoff on minimum number of cases (>=)
MutSigCutoff <- 0.05 #minimum cutoff for MutSigCV qval (>=)
probSigCutoff <- 0.06217 #default minimum adjusted p-value cutoff for probability score
SNRcutoff <- 2 #minimum cutoff for SNR value (>=)
usePairedExprData <- FALSE #if true, will analyze for both paired paritioned dataset results (expressed | not expressed)
                           #i.e. the final plots will contain both expressed and non-expressed genes.
tumorTypesToAnalyze <- c() #if empty, analyze all tumor types, otherwise analyze the specified tumor type

#settings overwrite if configuration settings exist
if(exists("config.avail") && config.avail){
  dset_dir <- config.dset_dir
  dsetAnalyzed_dir <- config.dsetAnalyzed_dir
  dsetScored_dir <- config.dsetScored_dir
  dsetScoredSummary_dir <- config.dsetScoredSummary_dir
  plots_dir <- config.plotsIndCancer_dir
  
  savePlots <- config.scored.savePlots
  saveLog <- config.scored.saveLog
  minCasesNcutoff <- config.scored.minCasesNcutoff
  MutSigCutoff <- config.scored.MutSigCutoff
  probSigCutoff <- config.scored.probSigCutoff
  SNRcutoff <- config.scored.SNRcutoff
  usePairedExprData <- config.scored.usePairedExprData
  tumorTypesToAnalyze <- config.scored.tumorTypesToAnalyze
}

#quality controls
checkDirExists(dsetAnalyzed_dir)
checkDirExists(dsetScored_dir)
createDir(dsetScoredSummary_dir)
if(savePlots){
  createDir(plots_dir)
}

#load datasets (to get access to dMAF); scored tumor specific datasets are loaded later
load(paste(dsetAnalyzed_dir,"datasets_analyzed.RData",sep=""))

#derived parameters
tumorTypes <- unique(dMAF$tumorType) #list of tumor types
tumorTypes <- tumorTypes[tumorTypes!='CELLLINE']
if(length(tumorTypesToAnalyze) > 0){tumorTypes <- tumorTypes[tumorTypes %in% tumorTypesToAnalyze]}
logFilename <- paste(dsetScoredSummary_dir,"analyze_scored.log",sep="")

#start diversion if applicable
if(saveLog){
  sink(logFilename, split=TRUE, append=FALSE)
}

###################################################
# Helper methods
###################################################
detProbSigCutoff <- function(vals){
  #determine pval/qval cutoff given a list of p-values/q-values
  #default method: use the minimum
  #DEPRECATED
  return(min(vals))
}

getIndxProbSig <- function(dataset){
  #get boolean indices of given dataset the significant hits based on SNR and significance cutoffs
  
  #use q-val and SNR
  return(!is.na(dataset$qval_prob) & !is.na(dataset$SNR) &
           dataset$qval_prob <= probSigCutoff & dataset$SNR >= SNRcutoff) 
}

getIndxMutSig <- function(dataset){
  #get boolean indices of given dataset the significant hits based MutSig significance cutoffs
  return(!is.na(dataset$qval_MutSig) & dataset$qval_MutSig <= MutSigCutoff)
}

dataframeNumeric <- function(dataset){
  #convert all columns are numeric, if not already
  rnames <- row.names(dataset)
  if(nrow(dataset) > 1){
    d <- as.data.frame(apply(dataset[,1:ncol(dataset)],2,function(x)as.numeric(x)))
  } else {
    cnames <- colnames(dataset)
    d <- as.data.frame(t(as.numeric(dataset)))
    colnames(d) <- cnames
  }
  
  row.names(d) <- rnames
  return(d)
}

SNRcorrect <- function(dataset){
  #if bootstrapped probSD is 0, set this to the smallest value in the dataset
  dataset <- cbind(dataset, SNRcorrected=dataset$SNR)
  SNRtoUpdate_idx <- !is.na(dataset$probSD) & dataset$probSD == 0
  minProbDB <- min(dataset$probSD[!SNRtoUpdate_idx], na.rm=TRUE)
  dataset$SNRcorrected[SNRtoUpdate_idx] <- dataset$prob[SNRtoUpdate_idx]/minProbDB
  return(dataset)
}

preprocessDataset <- function(dataset){
  #preprocess dataset
  
  #if minCasesNcutoff is defined, remove entries that is below the cutoff
  if(!is.null(minCasesNcutoff)){
    dataset <- dataset[(!is.na(dataset$casesN) & dataset$casesN >= minCasesNcutoff) | is.na(dataset$casesN),]
  }
  
  #dataset cleanups and extract properties
  idx_NAs <- is.na(dataset$qval_MutSig)
  if(sum(idx_NAs) > 0){
    message(gettextf("There are %.0f entries in MutSig qval with NAs, converting them to 1...", sum(idx_NAs)))
    dataset$qval_MutSig[idx_NAs] <- 1
  }
  
  return(dataset)
}

###################################################
# Main analysis of scored dataset
###################################################
analyzeScored <- function(tumor){
  #derived settings
  #return two lists: expr and nonexpr genes for given tumor
  
  plotsTumor_dir <- paste(plots_dir, tumor2filename(tumor), sep="")
  #create plots folder if they don't exist
  if(savePlots){
    createDir(plotsTumor_dir)
  }
  
  ##################################
  #reading in datasets
  notExprAvail <- FALSE
  datasetNotExpr <- NULL
  if(usePairedExprData){
    sigListNotExpr_filename <- paste(dsetScored_dir, tumor2filename(tumor), "_sigList_notExpressed.RData",sep="")
    if(file.exists(sigListNotExpr_filename)){
      load(sigListNotExpr_filename)
      sigListGeneNotExpr <- sigListGene
      sigListGeneNotExpr <- dataframeNumeric(sigListGeneNotExpr)
      datasetNotExpr <- preprocessDataset(SNRcorrect(sigListGeneNotExpr))
      
      if(exists('sigListGene.appended') & !is.null(sigListGene.appended)){
        sigListGeneNotExpr.appended <- sigListGene.appended
        sigListGeneNotExpr.appended <- dataframeNumeric(sigListGeneNotExpr.appended)
        datasetNotExpr.appended <- preprocessDataset(SNRcorrect(sigListGeneNotExpr.appended))
      }
      
      notExprAvail <- TRUE
    }
  }
  
  sigList_filename <- paste(dsetScored_dir, tumor2filename(tumor), "_sigList.RData",sep="")
  if(!file.exists(sigList_filename)){
    message(gettextf("WARNING: cannot open %s. Skipping analysis/plots generation...", sigList_filename))
    return(NULL)
  }
  load(sigList_filename) #load variables 'sigList' and 'sigListGene'
  sigListGene <- dataframeNumeric(sigListGene)
  datasetExpr <- preprocessDataset(SNRcorrect(sigListGene))
  
  if(exists('sigListGene.appended') & !is.null(sigListGene.appended)){
    sigListGene.appended <- dataframeNumeric(sigListGene.appended)
    datasetExpr.appended <- preprocessDataset(SNRcorrect(sigListGene.appended))
  }
  
  #default plot values
  colf <- rgb(0,0,0,70,maxColorValue=255)
  colb <- rgb(0,0,0,80,maxColorValue=255)
  colborder <- rgb(0,0,0,80,maxColorValue=255)
  colMutSig <- rgb(65,105,225,80,maxColorValue=255)
  colNovel <- rgb(178,34,34,80,maxColorValue=255)
  colNull <- rgb(150,150,150,70,maxColorValue=255)
  
  ##################################
  #Plots for combined dataset only
  #signal vs. significant
  if(notExprAvail){
    datasetCombined <- rbind(datasetExpr, datasetNotExpr)
  } else {
    datasetCombined <- datasetExpr
  }
  boolidx_probsig <- getIndxProbSig(datasetCombined)
  boolidx_mutsig <- getIndxMutSig(datasetCombined)
  colors <- rep(colNull,nrow(datasetCombined))
  colors[boolidx_probsig] <- colNovel
  colors[boolidx_mutsig] <- colMutSig
  xval <- datasetCombined$pval_prob
  xvalcutoff <- max(xval[boolidx_probsig])
  xvallab <- "Nominal p-value"
  textpos <- 4
  yscale <- c(0,ceiling(max(datasetCombined$SNRcorrected[xval!=0])))
  
  #log transform
  xval <- -log10(xval)
  xvalcutoff <- -log10(xvalcutoff)
  yscale[2] <- log2(yscale[2])
  xvallab <- "-log10[Nominal p-value]"
  textpos <- 2
  
  plotfn <- 'pvalvslog2SNR'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  plot(xval, log2(datasetCombined$SNRcorrected), xlab=xvallab, ylab="log2[SNR]", ylim=yscale,
       cex=2, col = colborder, bg = colors, pch=21)
  if((sum(boolidx_probsig)+sum(boolidx_mutsig)) > 0){
    text(xval[boolidx_probsig|boolidx_mutsig], log2(datasetCombined$SNRcorrected[boolidx_probsig|boolidx_mutsig]), 
         row.names(datasetCombined[boolidx_probsig|boolidx_mutsig,]), cex=0.75, pos=textpos, col="black")
  }
  abline(h=log2(SNRcutoff), col=colb, lty=2)
  abline(v=xvalcutoff, col=colb, lty=2)
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  #if there are any pval = 0, generate a separate plot for these hits
  if(sum(datasetCombined$pval_prob == 0) > 1){
    xval2 <- datasetCombined$pval_prob
    boolidx <- xval2==0
    xvallab <- "Nominal p-value"
    plotfn <- 'pvalvslog2SNR_zeros'
    openPDFdev(plotsTumor_dir, plotfn, savePlots)
    plot(xval2[boolidx], log2(datasetCombined[boolidx,]$SNRcorrected), 
         xlab=xvallab, ylab="log2[SNR]", ylim=yscale,
         cex=2, col = colborder, bg = colors[boolidx], pch=21)
    if((sum(boolidx_probsig)+sum(boolidx_mutsig)) > 0){
      text(xval2[boolidx&(boolidx_probsig|boolidx_mutsig)], 
           log2(datasetCombined$SNRcorrected[boolidx&(boolidx_probsig|boolidx_mutsig)]), 
           row.names(datasetCombined[boolidx&(boolidx_probsig|boolidx_mutsig),]), cex=0.75, pos=textpos, col="black")
    }
    saveFig(plotsTumor_dir, plotfn, savePlots)
    dev.off()
  }
  
  plotfn <- 'pvalvsprob'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  plot(xval, datasetCombined$prob, xlab=xvallab, ylab="prob",
       cex=2, col = colborder, bg = colors, pch=21)
  if((sum(boolidx_probsig)+sum(boolidx_mutsig)) > 0){
    text(xval[boolidx_probsig|boolidx_mutsig], datasetCombined$prob[boolidx_probsig|boolidx_mutsig],
         row.names(datasetCombined[boolidx_probsig|boolidx_mutsig,]), cex=0.75, pos=textpos, col="black")
  }
  #abline(h=0.2, col=colb, lty=2)
  #abline(v=xvalcutoff, col=colb, lty=2)
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  #cumulative distributions of p-values
  plotfn <- 'pval_ecdf'
  openPDFdev(plotsTumor_dir, plotfn, savePlots, height=4, width=8)
  g <- NULL
  if(notExprAvail){
    dToPlot <- melt(list(expressed=datasetExpr$pval_prob,notexpressed=datasetNotExpr$pval_prob))
    #gline <- geom_vline(xintercept=min(datasetNotExpr$pval_prob),linetype=2,size=0.8,color="#66676E")
    g <- ggplot(dToPlot, aes(value,group=L1,color=L1))
      
  }
  if(is.null(g)){
    dToPlot <- melt(list(expressed=datasetExpr$pval_prob))
    g <-ggplot(dToPlot, aes(value,group=L1,color=L1))
  }
  gh <- g + stat_ecdf(size=0.8) + labs(x="Nominal p-value (prob)", y="eCDF", color=gettextf("dataset (%s)",tumor)) + 
    ylim(0,1) + scale_x_log10() + theme(legend.position=c(0,1),legend.justification=c(0,1))
  print(gh)
  saveFig(plotsTumor_dir, plotfn, savePlots, height=500, width=900)
  dev.off()
  
  ##################################
  #Plots for expressed only
  dataset <- datasetExpr
  
  #extract values from dataset  
  casesNmax <- max(dataset$casesN)
  
  #reduced dataset with rows containing NAs removed
  rowIdxwithNAs <- apply(dataset,1,function(x){sum(is.na(x))})
  dataset_noNAs <- dataset[rowIdxwithNAs<1,]
  
  plotfn <- '2heatmap'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  cormatrix <- as.matrix(cor(na.omit(dataset[,c(1,2,4,6,8,9,10,11,12,13,15,16)]), method="spearman", use="pairwise"))
  if(sum(is.na(cormatrix))<1){
    heatmap.2(cormatrix, density.info="none", trace="none", col=brewer.pal(11,"RdBu"),
              cexRow=0.8,cexCol=0.8,margins=c(10,10))
    saveFig(plotsTumor_dir, plotfn, savePlots)
    dev.off()
  } else{
    message("WARNING: there are NAs in the correlation matrix")
  }
  
  #ranked SNR
  plotfn <- '3rankedSNR'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  labelsN = min(10, length(dataset$SNR))
  sorted_idx <- order(dataset$SNR)
  SNR.sorted <- dataset$SNR[sorted_idx]
  rownames.sorted <- rownames(dataset[sorted_idx,])
  plot(1:length(SNR.sorted), log2(SNR.sorted), xlab="Genes", ylab="log2[SNR]", cex=1.5, bg=colNull, col=colborder, pch=21)
  toLabel <- (length(SNR.sorted)-labelsN+1):length(SNR.sorted)
  text(toLabel, log2(SNR.sorted[toLabel]),
       rownames.sorted[toLabel], cex=0.75, pos=2, col="black")
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  #metrics vs. protein length
  plotfn <- '6protLengthvsSNR'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  smoothScatter(dataset$protLength, log2(dataset$SNR), xlab="Protein length (# of residues)", ylab="log2[SNR]")
  #plot(dataset$protLength, dataset$SNR, cex=2, bg=colb, col=colf, pch=21,
  #     xlab="Protein length (# of residues)", ylab="log2[SNR]")
  #text(dataset$protLength, dataset$SNR, row.names(dataset), cex=0.75, pos=4, col="black")
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  plotfn <- '7protLengthvsprob'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  smoothScatter(dataset$protLength, dataset$prob, xlab="Protein length (# of residues)", ylab="Prob")
  #plot(dataset$protLength, dataset$prob, cex=2, bg=colb, col=colf, pch=21,
  #     xlab="Protein length (# of residues)", ylab="Prob")
  #text(dataset$protLength, dataset$prob, row.names(dataset), cex=0.75, pos=4, col="black")
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  plotfn <- '8protLengthvspval'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  smoothScatter(dataset$protLength, dataset$pval_prob, xlab="Protein length (# of residues)", ylab="Nominal p-value (Prob)")
  #plot(dataset$protLength, dataset$pval_prob, cex=2, bg=colb, col=colf, pch=21,
  #     xlab="Protein length (# of residues)", ylab="Nominal p-value (Prob)")
  #text(dataset$protLength, dataset$qval_prob, row.names(dataset), cex=0.75, pos=4, col="black")
  saveFig(plotsTumor_dir, '8protLengthvspval', savePlots)
  dev.off()
  
  #number of unique positions
  plotfn <- '9pvalvsuniqPos'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  plot(dataset$uniquePosNdMAFtumor, dataset$pval_prob,
       xlab="Number of unique residues mutated in tumor",
       ylab="Nominal p-value (Prob)", cex=1.5, bg=colb, col=colf, pch=21)
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  plotfn <- '10pvalvsmutrate'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  plot(dataset$uniquePosNdMAFtumor/dataset$protLength, dataset$pval_prob, 
       xlab=gettextf("Number of unique residues mutated in tumor / protein length", tumor),
       ylab="Nominal p-value (Prob)", cex=2, bg=colb, col=colf, pch=21)
  lmfit <- lm(y ~ x, 
              data=data.frame(y=dataset$pval_prob,
                              x=dataset$uniquePosNdMAFtumor/dataset$protLength))
  #abline(a=coef(lmfit)[1], b=coef(lmfit)[2], col="red", lty=2)
  #text(dataset$uniquePosNdMAFtumor/dataset$protLength, dataset$pval_prob, row.names(dataset), cex=0.75, pos=4, col="black")
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  #plotfn <- '11uniqPosHist'
  #openPDFdev(plotsTumor_dir, plotfn, savePlots)
  #hist(dataset$uniquePosNdMAFtumor/dataset$uniquePosNdMAF,
  #     xlab=gettextf("Ratio of number of unique residues mutated in %s to all tumors", tumor),
  #     main="")
  #saveFig(plotsTumor_dir, plotfn, savePlots)
  #dev.off()
  
  plotfn <- '11SNRHist'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  hist(dataset$SNR,
      xlab=gettextf("SNR", tumor),
      main="")
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  #principal component analysis (using dataset_noNAs)
#   compvars <- c(1,2,4,6,8,9,10,11,12,13,15,16)
#   if(nrow(dataset_noNAs) >= length(compvars)){
#     plotfn <- '12PCA'
#     openPDFdev(plotsTumor_dir, plotfn, savePlots)
#     fit <- princomp(dataset_noNAs[,compvars], cor=FALSE, na.action=na.omit)
#     plot(fit$scores, cex=5*dataset_noNAs$casesN/casesNmax)
#     text(fit$scores, row.names(dataset_noNAs), cex=0.75, pos=4, col="black")
#     saveFig(plotsTumor_dir, plotfn, savePlots)
#     dev.off()
#     
#     plotfn <- '13PCAbiplot'
#     openPDFdev(plotsTumor_dir, plotfn, savePlots)
#     biplot(fit,cex=c(0.75,1))
#     saveFig(plotsTumor_dir, plotfn, savePlots)
#     dev.off()
#   } else {
#     message("Cannot perform PCA, number of observations (in reduced dataset) is too small.")
#   }
  
  #probSD plot
  plotfn <- '14probSDvsCases'
  openPDFdev(plotsTumor_dir, plotfn, savePlots)
  plot(dataset$casesN, dataset$probSD, xlab="Number of cases", ylab="probSD (bootstrap)", cex=2, bg=colb, col=colf, pch=21)
  saveFig(plotsTumor_dir, plotfn, savePlots)
  dev.off()
  
  ##################################
  #Return list of significant genes that pass threshold
  siglist_out <- NULL
  siglistNotExpr_out <- NULL
  
  dataset <- if(exists('datasetExpr.appended')) datasetExpr.appended else datasetExpr
  siglist_out <- cbind(dataset,
                       MutSigPass=rep(0, nrow(dataset)),
                       probPass=rep(0, nrow(dataset)),
                       tumorType=rep(tumor, nrow(dataset)), stringsAsFactors=FALSE)
  siglist_out$MutSigPass[getIndxMutSig(dataset)] <- 1
  siglist_out$probPass[getIndxProbSig(dataset)] <- 1
  
  dataset <- if(exists('datasetNotExpr.appended')) datasetNotExpr.appended else datasetNotExpr
  if(!is.null(dataset)){
    siglistNotExpr_out <- cbind(dataset,
                                tumorType=rep(tumor, nrow(dataset)), stringsAsFactors=FALSE)
  }
  
  return(list(siglist=siglist_out, notExpr=siglistNotExpr_out))
}

###################################################
# Analyze by tumor type
###################################################
catn(gettextf("*** Analyzing %s tumor types...", length(tumorTypes)))

#analysis
siglist_all <- c()
notExpr_all <- c()
for(tumor in tumorTypes){
  catn(gettextf("Analyzing and generating plots for %s...", tumor))
  
  #analyze data, generate plots
  r <- analyzeScored(tumor)
  siglist_tumor <- r$siglist
  notExpr_tumor <- r$notExpr
  
  #combine results
  if(!is.null(siglist_tumor)){
    #stop using gene names as row.names, so siglist_tumor can be merged
    siglist_tumor <- cbind(geneSymbol=row.names(siglist_tumor), siglist_tumor, stringsAsFactors=FALSE)
    row.names(siglist_tumor) <- NULL
    
    #merge siglist_tumor (across tumor types)
    siglist_all <- rbind(siglist_all, siglist_tumor)
  }
  
  if(!is.null(notExpr_tumor)){
    #stop using gene names as row.names, so siglist_tumor can be merged
    notExpr_tumor <- cbind(geneSymbol=row.names(notExpr_tumor), notExpr_tumor, stringsAsFactors=FALSE)
    row.names(notExpr_tumor) <- NULL
    
    #merge siglist_tumor (across tumor types)
    notExpr_all <- rbind(notExpr_all, notExpr_tumor)
  }
}

graphics.off() #close all devices
catn("Writing results...")
write.csv(siglist_all,file=paste(dsetScoredSummary_dir, "scored_sigListGene.csv",sep=""))
if(!is.null(notExpr_all)){
  write.csv(notExpr_all,file=paste(dsetScoredSummary_dir, "scored_notExprListGene.csv",sep=""))
}
save(siglist_all,notExpr_all,file=paste(dsetScoredSummary_dir,"sigListGene.RData",sep=""))
if(saveLog){
  #close diversion
  sink()
}
