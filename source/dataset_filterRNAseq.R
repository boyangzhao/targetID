#dataset_filterRNAseq.R
#Author: Boyang Zhao
#Read in processed datasets and add RNA-seq data to dMAF
#Note: requires calling cBioPortal to get RNA-Seq data
#REQUIRES: datasets_processed_p1.RData generated from dataset_process.R
#MODIFIES: dMAF
#GENERATES: dHUMSAVAR, dMAF, dSNV

rm(list=ls())
library(cgdsr)
source("util.R")
dsetProcessed_dir <- "../datasets/processed/"
load(paste(dsetProcessed_dir,"datasets_processed_p1.RData",sep=""))

###################################################
# Functions for RNA-seq acquisition
###################################################
#Acquire mRNA RNA-Seq V2 (normalized_count) from TCGA through cBioPortal
acqRNAseq <- function(dataset){
  #given dataset, query cBioPortal to get the corresponding RNA-seq (using caseID)
  #returns a list of RNA-seq values, the same size (i.e. number of rows) as given dataset
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  preval<-0
  mRNAvals <- c()
  
  for(idx in 1:nrow(dataset)){
    calcProgress(idx, nrow(dataset), preval)
    gene <- dataset$geneSymbol[idx]
    tumor <- dataset$tumorType[idx]
    tumor[tumor == 'COAD/READ'] <- 'coadread'
    tumor <- tolower(tumor)
    case <- dataset$caseID[idx]
    datasource <- sub("^([A-Za-z0-9]*)-.*$","\\1",dataset$caseID[idx])
    datasource <- tolower(datasource)
    cancerstudy <- paste(tumor, '_', datasource, sep="")
    geneticprofile <- paste(cancerstudy, '_rna_seq_v2_mrna', sep="") #correspond to rsem.genes.normalized_results
    
    if(datasource == 'tcga'){
      #only focus on TCGA
      
      mRNAval <- getProfileData(mycgds, gene, geneticprofile, cases=case)
      if(nrow(mRNAval) > 0){
        mRNAval <- mRNAval[1,1]
        if(is.na(as.numeric(mRNAval))){
          caselist <- getCaseLists(mycgds,cancerstudy)[getCaseLists(mycgds,cancerstudy)$case_list_id==geneticprofile,'case_ids']
          caselist <- unlist(strsplit(caselist,' '))
          if(!(case %in% caselist)){
            mRNAval = "CASENOTAVAILABLE"
          }
        }
      } else {
        mRNAval = "PROFILEEMPTY"
      }
    } else {
      mRNAval = "NOTTCGA"
    }
    mRNAvals <- c(mRNAvals, mRNAval)
  }
  
  return(mRNAvals)
}

###################################################
# Incorporate RNA-seq data into dMAF
###################################################
mRNAvals <- acqRNAseq(dMAF)

if(nrow(dMAF) != length(mRNAvals)){
  stop('Lengths of dMAF and mRNAvals do not match.')
}

#Append two columns to dMAF, one with mRNA values directly from cBioPortal query above,
#and one with the final mRNA values to use, which fill in non-numeric mRNAvals with median of mRNAvals of given gene
dMAF <- cbind(dMAF, mRNAvals=mRNAvals, mRNAvals_final=rep(NA,nrow(dMAF)), stringsAsFactors=FALSE)

#Get median mRNA values for each gene, and generate a gene_mRNA_pair lookup table consist of gene:mRNA value pairs
dataset_withmRNA_idx <- !is.na(as.numeric(dMAF$mRNAvals))
dataset_withmRNA <- dMAF[dataset_withmRNA_idx,]
dataset_withmRNA_genes <- unique(dataset_withmRNA$geneSymbol)

genesN_noexpression <- 0
if(length(dataset_withmRNA_genes) < length(unique((dMAF$geneSymbol)))){
  genesN_noexpression <- length(unique((dMAF$geneSymbol))) - length(dataset_withmRNA_genes)
  print(gettextf("There are %s genes without any mRNA expression data",genesN_noexpression))
}

getmRNA_medvals <- function(gene){
  median(as.numeric(dataset_withmRNA[dataset_withmRNA$geneSymbol == gene, 'mRNAvals']))
}
mRNA_medvals <- vapply(dataset_withmRNA_genes, getmRNA_medvals, 0)
gene_mRNA_pair <- data.frame(geneSymbol=dataset_withmRNA_genes,
                             mRNA_medvals=mRNA_medvals, stringsAsFactors=FALSE)

#Fill in entries in mRNAvals_final column of dataset
#lookup corresponding index in gene_mRNA_pair for each gene in dataset
expand_idx <- match(dMAF$geneSymbol, gene_mRNA_pair$geneSymbol)

#for entries with no mRNA value, fill in median mRNA value for that gene
dMAF[!dataset_withmRNA_idx,'mRNAvals_final'] <- gene_mRNA_pair[expand_idx[!dataset_withmRNA_idx],'mRNA_medvals']

#for entries with mRNA value, grab the value from mRNAvals column and copy over to mRNAvals_final column
dMAF[dataset_withmRNA_idx,'mRNAvals_final'] <- dMAF[dataset_withmRNA_idx,'mRNAvals']

#convert mRNAvals_final column to type numeric
dMAF$mRNAvals_final<-as.numeric(dMAF$mRNAvals_final)

#quality controls
NAsN <- sum(is.na(dMAF$mRNAvals_final))
if(NAsN > 0){
  genesNAs <- unique(dMAF[is.na(dMAF$mRNAvals_final),'geneSymbol'])
  print(gettextf("Note: There are %s NAs in dMAF mRNAvals_final column, corresponding to %s genes.", 
                 NAsN, length(genesNAs)))
}

if(genesN_noexpression < 1 && NAsN > 0){
  print("ERROR: Something is wrong, all genes have at least one piece of mRNA expression info but there are NAs in 
        final dataset")
}

write.csv(mRNAvals,file=paste(dsetProcessed_dir,"mRNAvals.csv",sep=""))
write.csv(gene_mRNA_pair,file=paste(dsetProcessed_dir,"mRNAvals_medpair.csv",sep=""))
write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p2.csv",sep=""))
save(mRNAvals,gene_mRNA_pair,file=paste(dsetProcessed_dir,"datasets_processed_p2_supp.RData",sep=""))
save(dHUMSAVAR,dMAF,dSNV,file=paste(dsetProcessed_dir,"datasets_processed_p2.RData",sep=""))
