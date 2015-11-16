#dataset_filterRNAseq2.R
#Read in processed dataset (from dataset_filterRNAseq.R) and segment dataset
#based on different ways of dealing with missing RNA-seq data, and TCGA vs. non TCGA
#Author: Boyang Zhao
#REQUIRES: dMAF (from datasets_processed_p2.RData, generated from dataset_filterRNAseq.R)
#GENERATES: datasets_processed_p3_medall_all.RData, 
#           datasets_processed_p3_medtumor_all.RData,
#           datasets_processed_p3_medtumor_TCGAonly.RData

rm(list=ls())
dsetProcessed_dir <- "../datasets/processed/" #input/output directory
source('util.R')

#settings
#derived parameters

#load datasets
load(paste(dsetProcessed_dir,"datasets_processed_p2.RData",sep=""))
dataset <- dMAF

#adjustments
#rename header mRNAvals_final to mRNAvals.adj.all
names(dataset)[names(dataset) == "mRNAvals_final"] <- "mRNAvals.adj.all"
dataset <- cbind(dataset,
                 mRNAvals.adj.perTumor=rep(NA,nrow(dataset)),
                 mRNAvals_final=rep(NA,nrow(dataset)),
                 stringsAsFactors=FALSE)

#dataset slices
dataset.genetumor <- paste(dataset$geneSymbol,'_',dataset$tumorType,sep="")
dataset.genetumor.uniq <- unique(dataset.genetumor)
dataset.withmRNA.idx <- !is.na(as.numeric(dataset$mRNAvals))
dataset.withmRNA <- dataset[dataset.withmRNA.idx,]
dataset.withmRNA.genetumor <- paste(dataset.withmRNA$geneSymbol,'_',dataset.withmRNA$tumorType,sep="")
dataset.withmRNA.genetumor.uniq <- unique(dataset.withmRNA.genetumor)
dataset.withmRNA.N <- nrow(dataset.withmRNA)
dataset.N <- nrow(dataset)

###############################
#general stat on entries with and without expression data
catn(gettextf("There are %s out of %s (%s%%) entries without expression data",
              dataset.N-dataset.withmRNA.N,
              dataset.N,
              round((dataset.N-dataset.withmRNA.N)*100/dataset.N)))

genetumor.noexpression.N <- 0
if(length(dataset.withmRNA.genetumor.uniq) < length(dataset.genetumor.uniq)){
  genetumor.noexpression <- dataset.genetumor.uniq[!(dataset.genetumor.uniq %in% dataset.withmRNA.genetumor.uniq)]
  genetumor.noexpression.N <- length(genetumor.noexpression)
  geneNames <- unique(sub("^(.*)_.*$","\\1",genetumor.noexpression))
  tumorTypes <- unique(sub("^.*_(.*)$","\\1",genetumor.noexpression))
  
  catn(gettextf("There are %s out of %s (%s%%) unique gene_tumor pairs without any mRNA expression data",
                genetumor.noexpression.N,
                length(dataset.genetumor.uniq),
                round(genetumor.noexpression.N*100/length(dataset.genetumor.uniq))))
  catn(gettextf("Corresponding to %s genes and %s tumor types", length(geneNames), length(tumorTypes)))
}

###############################
#generate reference table of expression values for gene_tumor
catn("Generating mRNA (per tumor) reference table...")
getmRNA_medvals <- function(genetumor.idx){
  calcProgress(genetumor.idx, length(dataset.withmRNA.genetumor.uniq), preval)
  genetumor <- dataset.withmRNA.genetumor.uniq[genetumor.idx]
  median(as.numeric(dataset.withmRNA[dataset.withmRNA.genetumor == genetumor, 'mRNAvals']))
}
preval<-0
mRNA_medvals <- vapply(1:length(dataset.withmRNA.genetumor.uniq), getmRNA_medvals, 0)
genetumor.mRNAvals <- data.frame(genetumor=dataset.withmRNA.genetumor.uniq,
                                 mRNA_medvals=mRNA_medvals, stringsAsFactors=FALSE)

###############################
#Fill in entries in mRNAvals.adj.perTumor column of dataset
catn("Filling in mRNA values for missing entries...")

#lookup corresponding index in genetumor.mRNAvals for each gene in dataset
expand_idx <- match(dataset.genetumor, genetumor.mRNAvals$genetumor)

#for entries with no mRNA value, fill in median mRNA value for that gene (tumor-specific)
dataset[!dataset.withmRNA.idx,'mRNAvals.adj.perTumor'] <- genetumor.mRNAvals[expand_idx[!dataset.withmRNA.idx],'mRNA_medvals']

#for entries with mRNA value, grab the value from mRNAvals column and copy over to mRNAvals_final column
dataset[dataset.withmRNA.idx,'mRNAvals.adj.perTumor'] <- dataset[dataset.withmRNA.idx,'mRNAvals']

###############################
#Fill in entries in mRNAvals.adj.final column of dataset
#fill first based on dealing missing mRNA with median value per gene per tumor
#then based on dealing missing mRNA with median value per gene for all tumors for any remaining entries
NAidx <- is.na(dataset$mRNAvals.adj.perTumor)
dataset$mRNAvals_final[!NAidx] <- dataset$mRNAvals.adj.perTumor[!NAidx]
dataset$mRNAvals_final[NAidx] <- dataset$mRNAvals.adj.all[NAidx]

#convert mRNAvals.adj.final column to type numeric
dataset$mRNAvals_final <- as.numeric(dataset$mRNAvals_final)

#quality controls
NAsN <- sum(is.na(dataset$mRNAvals_final))
if(NAsN > 0){
  genesNAs <- unique(dataset[is.na(dataset$mRNAvals_final),'geneSymbol'])
  message(gettextf("Note: There are %s NAs in dMAF mRNAvals_final column, corresponding to %s genes.", 
                   NAsN, length(genesNAs)))
}

###############################
#slice data with TCGA-only data
dataset.TCGAonly <- dataset[dataset$mRNAvals != "NOTTCGA",]


###############################
#save datasets
catn("Saving datasets...")
write.csv(mRNA_medvals,file=paste(dsetProcessed_dir,"mRNAvals_medpair_genetumor.csv",sep=""))
save(mRNA_medvals,file=paste(dsetProcessed_dir,"datasets_processed_p3_supp.RData",sep=""))

#with dMAF, missing mRNA value for gene based on pan-cancer median mRNA value for given gene
#this is the untouched version of dMAF from p2
write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p3_medall_all.csv",sep=""))
save(dHUMSAVAR,dMAF,dSNV,file=paste(dsetProcessed_dir,"datasets_processed_p3_medall_all.RData",sep=""))

#with dMAF, missing mRNA value for gene based on tumor-specific median mRNA value for given gene, and 
#if not applicable, then take pan-cancer median value
dMAF <- dataset
write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p3_medtumor_all.csv",sep=""))
save(dHUMSAVAR,dMAF,dSNV,file=paste(dsetProcessed_dir,"datasets_processed_p3_medtumor_all.RData",sep=""))

#with dMAF, filtered for only TCGA entries
#missing mRNA value for gene based on tumor-specific median value for given gene, and if not applicable,
#then take pan-cancer median value
dMAF <- dataset.TCGAonly
write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p3_medtumor_TCGAonly.csv",sep=""))
save(dHUMSAVAR,dMAF,dSNV,file=paste(dsetProcessed_dir,"datasets_processed_p3_medtumor_TCGAonly.RData",sep=""))
