#dataset_exclusionFilters.R
#Author: Boyang Zhao
#Filter variants out based on PanCancer12

rm(list=ls())
dsetProcessed_dir <- "../datasets/processed/" #input/output directory
dsetMAF_dir <- "../datasets/pancancer12/" #input/output directory
source('util.R')

#load datasets
load(paste(dsetProcessed_dir,"datasets_processed_p3_medtumor_TCGAonly.RData",sep=""))

dMAF.filtered <- data.frame()

#loop through MAFs
tumors <- c('blca', 'brca', 'cesc', 'coadread', 'gbm', 'hnsc', 'kirc', 'kirp', 'laml',
            'lgg', 'luad', 'lusc', 'ov', 'paad', 'prad', 'skcm', 'stad', 'thca', 'ucec')

for(tumor in tumors){
  print(gettextf('Reading %s...',tumor))
  #d12 <- read.table(gettextf('./syn1729383_maf/%s_cleaned_filtered.maf',tumor), sep='\t', header=T, stringsAsFactors=F)
  d12 <- read.csv(gettextf('%ssyn1729383_maf/%s_cleaned_filtered.csv', dsetMAF_dir, tumor), header=T, stringsAsFactors=F)
  if(tumor == 'coadread'){
    tumorName = 'COAD/READ'
  } else {
    tumorName = toupper(tumor)
  }
  
  #find intersection between dMAF and pancancer12 filtered dataset, and
  #filter out any entries for given tumor type where variant is not found in pancancer12
  d <- dMAF[dMAF$tumorType == tumorName,]
  d12.posvar <- paste(d12$Chromosome, d12$Start_Position, d12$reference, d12$variant, sep='_')
  d.posvar <- paste(d$chr, d$start, d$refAllele, d$varAllele, sep='_')
  intersect.posvar <- intersect(d12.posvar, d.posvar)
  
  #append
  dMAF.filtered <- rbind(dMAF.filtered, d[d.posvar %in% intersect.posvar,])
}

#save new dMAF
dMAF <- dMAF.filtered #overview dMAF and save
write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p3_medtumor_TCGAonly_filtered.csv",sep=""))
save(dHUMSAVAR,dMAF,dSNV,file=paste(dsetProcessed_dir,"datasets_processed_p3_medtumor_TCGAonly_filtered.RData",sep=""))
