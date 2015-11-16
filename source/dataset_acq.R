#dataset_acq.R
#Author: Boyang Zhao
#Read in and preprocess datasets
#Note: requires calling cBioPortal to get mutation data
#REQUIRES: ClinVar dataset, HUMSAVAR dataset, UCSC RefSeq, [cBioPortal API]
#Generates: dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw (saved in [dataset dir]/local/)

rm(list=ls())
dset_dir <- "../datasets/"
library(cgdsr)

###########################################
#Read in Mendelian allele variants
print("Reading in SNV data...")
dSNVraw <- read.csv(paste(dset_dir,"ClinVar/variant_summary.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)

#get list of unique gene names in dSNP
geneList = unique(dSNVraw$GeneSymbol)

write.csv(dSNVraw,file=paste(dset_dir,"local/","dSNVraw.csv",sep=""))

###########################################
# #Read in somatic mutations in tumor samples, from TCGA MAF files
# #loop through each file in TCGA file directory, and read in data
# filenames <- list.files(path=paste(dset_dir,"TCGA/",sep=""), pattern="*.maf",full.names=T,recursive=FALSE)
# 
# dMAFraw <- data.frame(matrix(nrow=0,ncol=14))
# 
# #Parse MAF data
# for(filename in filenames){
#   #called per MAF file
#   print(paste("reading",filename))
#   tcgaSM <- read.table(paste(dset_dir,"TCGA/", filename,sep=""), sep="\t",header=TRUE,stringsAsFactors=FALSE)
#   
#   #extract TCGA tumor type name, import TCGA data to dMAFraw
#   tumortypeName <- substr(filename,start=9,stop=regexpr("_",filename)-1)
#   
#   #list layout: tumor, gene name, gene id, variant class, variant type, ref allele, 
#   #             tumor allele1, tumor allele2, chromosome #, start, end,
#   #             mutation status, tumor sample barcode, matched sample barcode
#   info_idx <- c(1,2,9,10,11,12,13,5,6,7,26,16,17)
#   dMAFraw <- rbind(dMAFraw, cbind(tumortypeName,tcgaSM[info_idx]))
#   
#   rm(tcgaSM)
# }
# 
# #rename some of the headers
# names(dMAFraw)[c(1,2,9,10,11)] <- c('TumorType','GeneSymbol','Chromosome','Start','Stop')

#Read in somatic mutations, using API from cBioPortal
print("Reading in TCGA data... Using cBioPortal API...")
cgdsh <- CGDS("http://www.cbioportal.org/public-portal/")
cancerstudies <- getCancerStudies(cgdsh)[,1]

dMAFraw <- data.frame(matrix(nrow=0,ncol=22))
c = 0
total = length(geneList)
for(geneName in geneList){
  if((c*100/total) %% 2 == 0){
    print(paste(c*100/total,"% complete...",sep=""))
  }
  
  for(cancerstudy in cancerstudies){
    mutdata <- getMutationData(cgdsh, paste(cancerstudy,"_all",sep=""), paste(cancerstudy,"_mutations",sep=""), geneName)
    if(nrow(mutdata)>0){
      dMAFraw = rbind(dMAFraw,mutdata)
    }
  }
  c = c + 1
}

rm(mutdata)

write.csv(dMAFraw,file=paste(dset_dir,"local/","dMAFraw.csv",sep=""))

###########################################
#Read in dHUMSAVAR dataset
print("Reading in HUMSVAR data...")
dHUMSAVARraw <- read.csv(paste(dset_dir,"humsavar/humsavar082914.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
write.csv(dHUMSAVARraw,file=paste(dset_dir,"local/","dHUMSAVARraw.csv",sep=""))

###########################################
#Read in dRefSeqraw dataset
print("Reading in UCSC RefSeq genes data...")
dRefSeqraw <- read.table(paste(dset_dir,"UCSC/refseq_genes083114.txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
write.csv(dRefSeqraw,file=paste(dset_dir,"local/","refseq_genes.txt",sep=""))

###########################################
save(dSNVraw, dMAFraw, dHUMSAVARraw, dRefSeqraw, file=paste(dset_dir,"local/","datasets.RData",sep=""))
