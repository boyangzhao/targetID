#dataset_preprocess.R
#Author: Boyang Zhao
#Read in and combine raw datasets, and save proprocessed datasets to processed folder.
#REQUIRES: dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw generated from dataset_acq.R
#GENERATES: dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw

rm(list=ls())
dsetRaw_dir <- "./"
dsetProcessed_dir <- "./" #output directory
source('util.R')

#quality controls
createDir(dsetProcessed_dir) #create output dir if needed

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
#Save datasets
save(dSNVraw, dHUMSAVARraw, file=paste(dsetProcessed_dir,"datasets_preprocessed.RData",sep=""))
