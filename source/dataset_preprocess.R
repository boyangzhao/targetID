#dataset_preprocess.R
#Author: Boyang Zhao
#Read in and combine raw datasets, and save proprocessed datasets to processed folder.
#REQUIRES: dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw generated from dataset_acq.R
#GENERATES: dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw

rm(list=ls())
dsetRaw_dir <- "../datasets/local/"
dsetProcessed_dir <- "../datasets/processed/" #output directory
source('util.R')

#quality controls
createDir(dsetProcessed_dir) #create output dir if needed

###########################################
#Preprocess dMAFraw dataset
print("Preprocessing dMAFraw...")
#combine dMAFraw csv files
filenames = list.files(path=dsetRaw_dir,pattern="^dMAFraw[0-9]*_[0-9]*\\.csv")
idx<-order(as.integer(sub("^dMAFraw([0-9]*)_[0-9]*\\.csv","\\1",filenames)))
filenames <- filenames[idx] #sorted

dMAFraw <- data.frame(matrix(nrow=0,ncol=22))
for(filename in filenames){
  fileToOpen = paste(dsetRaw_dir,filename,sep="")
  print(paste("Reading ",fileToOpen,"...",sep=""))
  d <- read.csv(fileToOpen,header=TRUE,stringsAsFactors=FALSE)
  dMAFraw <- rbind(dMAFraw,d)
}

rm(d)

#save dMAFraw to dsetProcessed_dir folder
write.csv(dMAFraw,file=paste(dsetProcessed_dir,"dMAFraw.csv",sep=""))

###########################################
#Preprocess dSNVraw dataset
print("Preprocessing dSNVraw...")
#retrive dSNVraw from dsetRaw_dir folder
dSNVraw <- read.csv(paste(dsetRaw_dir,"dSNVraw.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
#save dSNVraw to dsetProcessed_dir folder
write.csv(dSNVraw,file=paste(dsetProcessed_dir,"dSNVraw.csv",sep=""))

###########################################
#Preprocess dHUMSAVARraw dataset
print("Preprocessing dHUMSAVARraw...")
#retrive dHUMSAVARraw from dsetRaw_dir folder
dHUMSAVARraw <- read.csv(paste(dsetRaw_dir,"dHUMSAVARraw.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
#save dSNVraw to dsetProcessed_dir folder
write.csv(dHUMSAVARraw,file=paste(dsetProcessed_dir,"dHUMSAVARraw.csv",sep=""))

###########################################
#Read in dRefSeqraw dataset
print("Preprocessing UCSC RefSeq genes data...")
dRefSeqraw <- read.table(paste(dsetRaw_dir,"refseq_genes.txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
write.csv(dRefSeqraw,file=paste(dsetProcessed_dir,"refseq_genes.txt",sep=""))

###########################################
#Save datasets
save(dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw, file=paste(dsetProcessed_dir,"datasets_preprocessed.RData",sep=""))
