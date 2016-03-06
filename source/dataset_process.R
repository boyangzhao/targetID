#dataset_process.R
#Author: Boyang Zhao and Justin Pritchard
#Read in preprocessed datasets, apply filters, and standardize field values
#Note: requires mutationassessor API to get aa position mutation information
#REQUIRES: datasets_preprocessed.RData generated from dataset_preprocess.R (specifically dSNVraw, dMAFraw, dHUMSAVARraw)
#GENERATES: dHUMSAVAR, dMAF, dSNV

rm(list=ls())
dsetRaw_dir <- "../datasets/local/"
dsetProcessed_dir <- "../datasets/processed/"
source("util.R")
require(RCurl)

#load preprocessed datasets
load(paste(dsetProcessed_dir,"datasets_preprocessed.RData",sep=""))

###################################################
# Helper methods
###################################################
aa3to1 <- function(inputData){
  #convert 3-letter AA to 1-letter AA
  AA3 = c("Ala","Arg","Asn","Asp","Cys","Glu","Gln","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Ter")
  AA1 = c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")
  
  for(i in c(1:length(AA3))){
    inputData <- gsub(AA3[i],AA1[i],inputData,ignore.case=FALSE)
  }
  
  return(inputData)
}


aaMutSplit <- function(AAmutData){
  #split AA mutation data into three separate lists
  #e.g. L162V, slit into L 162 V
  #applicable mutations: AA->*, *->AA, AA->AA, AA->del, AA->fs, AA->splice)
  
  rowsTotal = length(AAmutData)
  
  #initialize lists
  aaPos <- rep("",rowsTotal)
  refAA <- rep("",rowsTotal)
  varAA <- rep("",rowsTotal)
  
  #get indices of applicable mutation AA entries to change
  idxToRec <- grep("^[A-Z\\*][0-9]*([A-Z\\*]|del|fs|_splice)$",AAmutData,ignore.case=TRUE) #find entries with applicable mutations to record
  idxNotToRec <- which(!(1:rowsTotal %in% idxToRec)) #complement of idxToRec
  
  #delimit with "-"
  toRec <- sub("^([A-Z\\*])([0-9]*)([A-Z\\*]|del|fs|_splice)$","\\1-\\2-\\3",AAmutData,ignore.case=TRUE)
  
  #mark empty or not applicable entries with "NA-NA-NA", so the split below works
  toRec[idxNotToRec] <- "NA-NA-NA"
  
  toRec_split <- strsplit(toRec,"-") #break into triplets: refAA, aaPos, varAA
  toRec_unlisted <- unlist(toRec_split) #collapse list of list to a single list
  toRec_unlisted_idx <- (idxToRec-1)*3 #adjust the index, since now each position is broken into three parts (items)
  
  refAA[idxToRec] <- toRec_unlisted[toRec_unlisted_idx+1]
  aaPos[idxToRec] <- toRec_unlisted[toRec_unlisted_idx+2]
  varAA[idxToRec] <- toRec_unlisted[toRec_unlisted_idx+3]
  
  #collect all results into list
  aaMutInfo <- list(refAA,aaPos,varAA)
  names(aaMutInfo) <- c("ref","pos","var")
  
  return(aaMutInfo)
}

domainSplit <- function(domainData){
  #parse domain@position from mutation assessor, split into three separate lists:
  #   name, start pos of domain, and end pos of domain
  #INPUT: string (format: domain id // domain name // start // end)
  #       e.g. "PF00129 // Class I Histocompatibility antigen, domains alpha 1 and 2 // 26 // 202"
  
  domainSplit <- list(name="",
                      start="", 
                      end="")
  
  if(domainData != ""){
    domainSplit_t <- strsplit(domainData, "//")
    domainSplit_t <- domainSplit_t[[1]]
    domainSplit <- list(name=domainSplit_t[1],
                        start=domainSplit_t[3], 
                        end=domainSplit_t[4])
  }
  
  return(domainSplit)
}

###########################################
# Process dMAF dataset
###########################################
print("Processing dMAF dataset...")

#fields (same between dMAF and dSNV): ID, geneSymbol, chr, start, end, refAllele, varAllele, aaPos, refAA, varAA
#fields (different between dMAF and dSNV): origin, mutationType, aaChange
#fields (unique to dMAF): tumorType, caseID

#ID, geneSymbol, chr, start, end
idx1 <- c(1,3,14,15,16)

#refAllele, varAllele
refAllele <- dMAFraw$reference_allele
varAllele <- dMAFraw$variant_allele

#some reference_allele or variant_allele has values "TRUE", this should be "T", the nucleotide; correct this:
idxToChange <- which(refAllele == "TRUE")
idxNoNA <- which(dMAFraw$xvar_link[idxToChange] != "NA")
idxToCheck <- idxToChange[idxNoNA]
nt <- sub(".*=hg19,[A-Z0-9]*,[0-9]*,([A-Z]),[A-Z]&.*","\\1",dMAFraw$xvar_link[idxToCheck])
if(length(which(nt != "T")) == 0){
  print("FYI: all reference_allele=\"TRUE\" are indeed meant to be nucleotide T. (checked for ones where xvar_link != NA)")
}
refAllele[idxToChange] = "T" #change TRUE to T

idxToChange <- which(varAllele=="TRUE")
idxNoNA <- which(dMAFraw$xvar_link[idxToChange] != "NA")
idxToCheck <- idxToChange[idxNoNA]
nt <- sub(".*=hg19,[A-Z0-9]*,[0-9]*,[A-Z],([A-Z])&.*","\\1",dMAFraw$xvar_link[idxToCheck])
if(length(which(nt != "T")) == 0){
  print("FYI: all variant_allele=\"TRUE\" are indeed meant to be nucleotide T. (checked for ones where xvar_link != NA)")
}
varAllele[idxToChange] = "T" #change TRUE to T

#aaChange
idx_aaChange <- c(9)

#aaPos, refAA, varAA (for mutations AA->*, *->AA, AA->AA, AA->del, AA->fs, AA->splice)
aaMutInfo <- aaMutSplit(dMAFraw$amino_acid_change)

#origin, mutationType
idx2 <- c(6,7)

#tumorType
tumortype <- substr(dMAFraw$genetic_profile_id,1,regexpr("_",dMAFraw$genetic_profile_id)-1) 
tumortype <- toupper(tumortype)
tumortype[tumortype=="COADREAD"] <- "COAD/READ"

#caseID
idx3 <- c(4)

nameLabels = c("ID", "geneSymbol", "chr", "start", "end", "refAllele", "varAllele", 
               "aaPos", "refAA", "varAA", "origin", "mutationType", "aaChange", "tumorType", "caseID")

dMAF = data.frame(dMAFraw[,idx1],
                  refAllele, varAllele,
                  aaMutInfo$pos, aaMutInfo$ref, aaMutInfo$var,
                  dMAFraw[,idx2],
                  dMAFraw[,idx_aaChange], tumortype,
                  dMAFraw[,idx3],stringsAsFactors=FALSE)
names(dMAF) <- nameLabels

write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p0.csv",sep=""))

###########################################
# Process dHUMSAVARraw dataset
###########################################
print("Processing dHUMSAVAR dataset...")

#fields (unique to dHUMSAVAR): ID
#fields (same between dMAF and dHUMSAVAR): geneSymbol, aaPos, refAA, varAA
#fields (unique to dHUMSAVAR): aachange, SPID, FTID, dsSNP, diseaseName

#focus only on pathogenic SNPs
dHUMSAVARraw <- dHUMSAVARraw[dHUMSAVARraw$varType == "Disease",]

#geneSymbol
idx1 <- c(1,2)

#aaChange, aaPos, refAA, varAA
aaChange <- sub("^[p]\\.(.*)$","\\1",dHUMSAVARraw$aaChange)
aaChange <- aa3to1(aaChange) #convert 3-letter AA to 1-letter AA

aaMutInfo <- aaMutSplit(aaChange)

#SPID (Swiss-Prot ID), FTID (Swiss-Prot feature ID), dsSNP, diseaseName
idx2 <- c(3,4,7,8)

nameLabels = c("ID","geneSymbol", "aaPos", "refAA", "varAA", "aaChange", 
               "SPID", "FTID", "dsSNP", "diseaseName")
dHUMSAVAR = data.frame(dHUMSAVARraw[,idx1],
                       aaMutInfo$pos, aaMutInfo$ref, aaMutInfo$var, aaChange,
                       dHUMSAVARraw[,idx2], stringsAsFactors=FALSE)
names(dHUMSAVAR) <- nameLabels

write.csv(dHUMSAVAR,file=paste(dsetProcessed_dir,"dHUMSAVAR_p0.csv",sep=""))


###########################################
# Process dSNV dataset
###########################################
print("Processing dSNV dataset...")

#fields (same between dMAF and dSNV): ID, geneSymbol, chr, start, end, refAllele, varAllele, aaPos, refAA, varAA
#fields (different between dMAF and dSNV): origin, mutationType, aaChange,
#fields (unique to dSNV): ntChange, domain, domainStart, domainEnd, clinicalSignificance, cytogenetic, otherIDs, mappingIssues

#focus on assembly 37 (same as what TCGA data are built on)
dSNVworking <- dSNVraw[dSNVraw$Assembly == "GRCh37",]
#focus only on pathogenic SNPs
dSNVworking <- dSNVworking[dSNVworking$ClinicalSignificance != "Benign",]
#remove somatic SNPs
dSNVworking <- dSNVworking[dSNVworking$Origin != "somatic",]
dSNVworking <- dSNVworking[dSNVworking$Origin != "somatic;unknown",]
dSNVworking <- dSNVworking[dSNVworking$Origin != "somatic;uncertain",]
#focus only on entries with known NM_ identifier
dSNVworking <- dSNVworking[grep("^NM_.*",dSNVworking$Name),]

rowsTotal <- nrow(dSNVworking)

#ID
idx1 <- c(1)
#chr, start, end
idx2 <- c(15,16,17)

#ntChange
ntChange <- sub(".*:[cgn]\\.(.*)\\(.*","\\1",dSNVworking$Name) #get index of records that contains .c, .g, or .n info
idxToRec1 <- grep(".*:[cgn]\\..*\\(.*",dSNVworking$Name)
ntChange <- sub(".*:[cgn]\\.(.*)$","\\1",ntChange)
idxToRec2 <- grep(".*:[cgn]\\..*$",ntChange)
idxToRec = union(idxToRec1,idxToRec2)
idxNotToRec <- which(!(1:rowsTotal %in% idxToRec)) #complement of idxToRec

ntChange[idxNotToRec] <- "" #remove records that does not contain .c, .g, or .n info

#refAllele, varAllele (for point mutations)
refAllele <- rep("",rowsTotal) #initialize lists
varAllele <- rep("",rowsTotal) #initialize lists

idxToRec <- grep(".*[A-Z]>[A-Z].*",ntChange)
refAllele[idxToRec] = sub(".*([A-Z])>[A-Z].*","\\1",ntChange)[idxToRec]
varAllele[idxToRec] = sub(".*[A-Z]>([A-Z]).*","\\1",ntChange)[idxToRec]

# #aaChange, aaPos, refAA, varAA
# deprecated -> used for parsing through dSNVraw dataset; will get these info from mutation assessor instead
# idxToRec <- grep(".*\\([p]\\.(.*)\\).*",dSNVworking$Name) #get index of records that contains .p info
# idxNotToRec <- which(!(1:rowsTotal %in% idxToRec)) #complement of idxToRec
# aaChange <- sub(".*\\([p]\\.(.*)\\).*","\\1",dSNVworking$Name)
# aaChange[idxNotToRec] <- "" #remove records that does not contain .p info
# aaChange <- aa3to1(aaChange) #convert 3-letter AA to 1-letter AA
# 
# aaMutInfo <- aaMutSplit(aaChange)

aaMutInfo <- list(ref=c(),pos=c(),var=c())
domainInfo <- list(name=c(),start=c(),end=c())
aaChange <- c()
geneSymbol <- c()
mappingissuesN = 0
mappingissues = c()
preval <- 0
for(idx in 1:rowsTotal){
  val <- floor(idx*100/rowsTotal)
  if(val %% 2 == 0 && preval != val){
    print(paste(val,"% complete...",sep=""))
    preval <- val
  }
  
  url <- "http://mutationassessor.org/?cm=var&var=hg19,%s,%s,%s,%s&fts=all&frm=txt "
  url <- gettextf(url,dSNVworking[idx,15],dSNVworking[idx,16],refAllele[idx],varAllele[idx])
  rsp <- getURL(url)
  
  entries <- strsplit(rsp,"\n",fixed=TRUE)
  entries <- entries[[1]]
  entries <- strsplit(entries,"\\t",fixed=TRUE)
  
  data<-data.frame(matrix(unlist(entries[[-1]]),ncol=length(entries[[1]]),byrow=TRUE),stringsAsFactors=FALSE)
  names(data)<-entries[[1]]
  
  #most of the mapping issues are irrelevant here -> e.g. issued with indel, nonsense, etc.
  mappingissues <- c(mappingissues, data[["Mapping issue"]])
  if(data[["Mapping issue"]] != ""){
    #print(paste("Mapping issue detected... for ",url,sep=""))
    mappingissuesN = mappingissuesN + 1
  }
  
  aaMutInfoEntry <- aaMutSplit(data[["AA variant"]])
  domainInfoEntry <- domainSplit(data[["domain@position"]])
  
  geneSymbol <- c(geneSymbol, data[["Gene"]])
  aaChange <- c(aaChange, data[["AA variant"]])
  aaMutInfo$ref <- c(aaMutInfo$ref, aaMutInfoEntry$ref)
  aaMutInfo$pos <- c(aaMutInfo$pos, aaMutInfoEntry$pos)
  aaMutInfo$var <- c(aaMutInfo$var, aaMutInfoEntry$var)
  domainInfo$name <- c(domainInfo$name, domainInfoEntry$name)
  domainInfo$start <- c(domainInfo$start, domainInfoEntry$start)
  domainInfo$end <- c(domainInfo$end, domainInfoEntry$end)
  
}

idx3 <- c(13,3) #origin, mutationType
idx4 <- c(7,18,25) #clinicalSignificance, cytogenetic, otherIDs

nameLabels = c("ID", "geneSymbol", "chr", "start", "end", "refAllele", "varAllele",
               "aaPos", "refAA", "varAA", "origin", "mutationType", "aaChange", "ntChange", 
               "domain", "domainStart", "domainEnd",
               "clinicalSignificance", "cytogenetic", "otherIDs", "mappingIssues")

dSNV = data.frame(dSNVworking[,idx1],
                  geneSymbol,
                  dSNVworking[,idx2],
                  refAllele, varAllele,
                  aaMutInfo$pos, aaMutInfo$ref, aaMutInfo$var,
                  dSNVworking[,idx3],
                  aaChange, ntChange, domainInfo$name, domainInfo$start, domainInfo$end,
                  dSNVworking[,idx4],
                  mappingissues, stringsAsFactors=FALSE)
names(dSNV) <- nameLabels

write.csv(dSNV,file=paste(dsetProcessed_dir,"dSNV_p0.csv",sep=""))

save(dSNV, dMAF, dHUMSAVAR, file=paste(dsetProcessed_dir,"datasets_processed_p0.RData",sep=""))

###########################################
# Additional cleanups, and processing by comparing across datasets (generates datasets_processed_p1)
###########################################
print("Additional processing (if applicable) across datasets")
rm(list=ls())
dsetProcessed_dir <- "../datasets/processed/"
load(paste(dsetProcessed_dir,"datasets_processed_p0.RData",sep=""))
require(stringr)
source("util_misc.R")

###################################################
#Remove silent mutations
removeSilentMut <- function(dataset1, dName){
  silentMutN <- sum(dataset1$refAA == dataset1$varAA & dataset1$refAA != "")
  if(silentMutN > 0){
    print(gettextf("Removing %.0f silent mutations from %s dataset",silentMutN,dName))
    dataset1 <- dataset1[dataset1$refAA != dataset1$varAA,]
  } else {
    print(gettextf("No silent mutations found in %s",dName))
  }
  return(dataset1)
}

dSNV <- removeSilentMut(dSNV, "dSNV")
dMAF <- removeSilentMut(dMAF, "dMAF")
dHUMSAVAR <- removeSilentMut(dHUMSAVAR, "dSNV")

###################################################
#Remove duplicate entries for same patient (i.e. duplicate gene_aaChange_caseID) (JR edit)
aachange_ID <- paste(dMAF$geneSymbol,dMAF$aaChange,dMAF$caseID,sep="_")
duplicated <- duplicated(aachange_ID)
print(gettextf("Removing %.0f duplicate entries from dMAF...",sum(duplicated)))
dMAF <- dMAF[!duplicated,]

###################################################
#Remove entries where refAA doesn't match within each dataset
dH_idxToCompr <- dHUMSAVAR$aaPos != "" & dHUMSAVAR$refAA != "-" & nchar(dHUMSAVAR$refAA) == 1
dHid <- paste(dHUMSAVAR$geneSymbol,dHUMSAVAR$aaPos,sep="_")
dHid2 <- paste(dHUMSAVAR$geneSymbol,dHUMSAVAR$aaPos,dHUMSAVAR$refAA,sep="_")
if(length(unique(dHid)) != length(unique(dHid2))){
  print("There are multiple entries in dHUMSAVAR with same geneSymbol_aaPos_refAA...")
  dHid_dup <- duplicated(dHid)
  dHid2_dup <- duplicated(dHid2)
  idxDiff <- xor(dHid_dup,dHid2_dup)
  difflist <- dHUMSAVAR[idxDiff,]
  difflist_u <- unique(paste(difflist$geneSymbol,difflist$aaPos,sep="_"))
  
  print(gettextf("Removing %.0f positions from dHUMSAVAR dataset...",length(difflist_u)))
  dHUMSAVAR <- dHUMSAVAR[!(dHid %in% difflist_u) & dH_idxToCompr,]
} else {
  print("There are no multiple entries in dHUMSAVAR with same geneSymbol_aaPos_refAA. Continuing...")
}

dM_idxToCompr <- dMAF$aaPos != "" & dMAF$refAA != "-" & nchar(dMAF$refAA) == 1
dMid <- paste(dMAF$geneSymbol,dMAF$aaPos,sep="_")
dMid2 <- paste(dMAF$geneSymbol,dMAF$aaPos,dMAF$refAA,sep="_")
if(length(unique(dMid)) != length(unique(dMid2))){
  print("There are multiple entries in dMAF with same geneSymbol_aaPos_refAA...")
  dMid_dup <- duplicated(dMid)
  dMid2_dup <- duplicated(dMid2)
  idxDiff <- xor(dMid_dup,dMid2_dup)
  difflist <- dMAF[idxDiff,]
  difflist_u <- unique(paste(difflist$geneSymbol,difflist$aaPos,sep="_"))
  
  print(gettextf("Removing %.0f positions from dMAF dataset...",length(difflist_u)))
  dMAF <- dMAF[!(dMid %in% difflist_u) & dM_idxToCompr,]
} else {
  print("There are no multiple entries in dMAF with same geneSymbol_aaPos_refAA. Continuing...")
}

###################################################
#Specific corrections
#Before cross-check datasets for refAA mismatch, correct entries (if possible) if genes get mapped to different isoforms
print("Correcting discrepancies in protein mapping between dMAF and dHUMSAVAR...")

#PTPN11 was mapped to isoform 2 in cBioPortal but isoform 1 in HUMSAVAR
#here, isoform naming is based on UniProt database; need to correct aaPos, refAA, and aaChange so both datasets refer to
#consistent protein (in this case, isoform 1)
dMAF <- correctPTPN11(dMAF)

#other types of specific corrections
#gene "TRUE" is supposed to be "T"
dMAF[dMAF$geneSymbol=="TRUE","geneSymbol"] <- "T"

###################################################
#Remove entries where refAA doesn't match between dMAF and dHUMSAVAR
dH_idxToCompr <- dHUMSAVAR$aaPos != "" & dHUMSAVAR$refAA != "-" & nchar(dHUMSAVAR$refAA) == 1
dM_idxToCompr <- dMAF$aaPos != "" & dMAF$refAA != "-" & nchar(dMAF$refAA) == 1

dHid <- paste(dHUMSAVAR$geneSymbol,dHUMSAVAR$aaPos,sep="_")
dMid <- paste(dMAF$geneSymbol,dMAF$aaPos,sep="_")
mutSame_dH <- dHid %in% dMid & dH_idxToCompr #index in relation to dHUMSAVAR
mutSame_dM <- dMid %in% dHid & dM_idxToCompr #index in relation to dMAF

dHid2 <- paste(dHUMSAVAR$geneSymbol,dHUMSAVAR$aaPos,dHUMSAVAR$refAA,sep="_")
dMid2 <- paste(dMAF$geneSymbol,dMAF$aaPos,dMAF$refAA,sep="_")
mutSame2_dH <- dHid2 %in% dMid2 & dH_idxToCompr #index in relation to dHUMSAVAR
mutSame2_dM <- dMid2 %in% dHid2 & dM_idxToCompr #index in relation to dMAF

diff <- xor(mutSame_dH,mutSame2_dH)
diffN <- sum(diff)
if(diffN > 0){
  #print(gettextf("Removing %.0f entries from dHUMSAVAR where refAllele doesn't match across dHUMSAVAR and dMAF",diffN))
  #dHUMSAVAR <- dHUMSAVAR[!diff,] %this only removes entries; need to be more strigent, remove the entire gene instead.
  
  genesToRemove <- unique(dHUMSAVAR[diff,"geneSymbol"])
  dHUMSAVAR <- dHUMSAVAR[!(dHUMSAVAR$geneSymbol %in% genesToRemove),]
  print(gettextf("Removing %.0f genes from dHUMSAVAR where refAllele doesn't match across dHUMSAVAR and dMAF",length(genesToRemove)))
  
} else {
  print("All entries from dHUMSAVAR have refAllele matched across dHUMSAVAR and dMAF")
}

diff <- xor(mutSame_dM,mutSame2_dM)
diffN <- sum(diff)
if(diffN > 0){
  #print(gettextf("Removing %.0f entries from dMAF where refAllele doesn't match across dHUMSAVAR and dMAF",diffN))
  #dMAF <- dMAF[!diff,] %this only removes entries; need to be more strigent, remove the entire gene instead.
  
  genesToRemove <- unique(dMAF[diff,"geneSymbol"])
  dMAF <- dMAF[!(dMAF$geneSymbol %in% genesToRemove),]
  print(gettextf("Removing %.0f genes from dMAF where refAllele doesn't match across dHUMSAVAR and dMAF",length(genesToRemove)))
  
} else {
  print("All entries from dMAF have refAllele matched across dHUMSAVAR and dMAF")
}

#save data
write.csv(dMAF,file=paste(dsetProcessed_dir,"dMAF_p1.csv",sep=""))
write.csv(dHUMSAVAR,file=paste(dsetProcessed_dir,"dHUMSAVAR_p1.csv",sep=""))
save(dSNV, dMAF, dHUMSAVAR, file=paste(dsetProcessed_dir,"datasets_processed_p1.RData",sep=""))
