#analysis_basic.R
#Analyze processed datasets: basic statistics
#Author: Boyang Zhao and Justin Pritchard
#REQUIRES: dHUMSAVAR and dMAF (from datasets_processed_p3_*.RData, generated from dataset_filterRNAseq2.R)
#GENERATES: datasets_analyzed.RData

if(exists("config.avail") && config.avail){ #if configuration settings exist, keep these, clear the rest of workspace
  rm(list=ls()[grep("^[^config\\.].*$",ls())])
} else {
  rm(list=ls())
}

dset_dir <- "../datasets/"
dsetProcessed_dir <- paste(dset_dir, "processed/", sep="")
dsetAnalyzed_dir <- paste(dset_dir, "analyzed/", sep="") #output directory
source('util.R')
require(stringr)

#settings
savePlots <- TRUE
processedDatasetName <- "datasets_processed_p3_medtumor_all.RData"

#settings overwrite if configuration settings exist
if(exists("config.avail") && config.avail){
  dset_dir <- config.dset_dir
  dsetProcessed_dir <- config.dsetProcessed_dir 
  dsetAnalyzed_dir <- config.dsetAnalyzed_dir
  savePlots <- config.basic.savePlots
  processedDatasetName <- config.basic.processedDatasetName
}

#derived parameters
plots_dir <- paste(dsetAnalyzed_dir,"plots/",sep="")

#load datasets
load(paste(dsetProcessed_dir, processedDatasetName,sep=""))

#additional cleanups
#remove from dSNV, entries with blank genesymbols, generally these refer to indels, etc.
#dSNV <- dSNV[dSNV$geneSymbol != "",]
dSNV <- dHUMSAVAR
geneList = unique(dSNV$geneSymbol)

#quality controls
createDir(dsetAnalyzed_dir) #create output dir if needed
if(savePlots){createDir(plots_dir)} #create output plots dir if needed

######################################
# Basic statistics
######################################
catn("Comparing each gene across datasets...")
genesCompr <- list(inCommon=c(),dSNVOnly=c(),dMAFOnly=c())
genesComprCorr <- list(dSNV=c(),dMAF=c()) #count of number of mutations per gene (each row = a gene)

for(geneName in geneList){
  #loop through each gene and compare across datasets
  in_dSNV_N <- sum(dSNV$geneSymbol==geneName)
  in_dMAF_N <- sum(dMAF$geneSymbol==geneName)
  
  genesComprCorr$dSNV <- c(genesComprCorr$dSNV,in_dSNV_N)
  genesComprCorr$dMAF <- c(genesComprCorr$dMAF,in_dMAF_N)
  
  if (in_dSNV_N > 0 && in_dMAF_N > 0)
    genesCompr$inCommon <- c(genesCompr$inCommon, geneName)
  else if (in_dSNV_N > 0)
    genesCompr$dSNVOnly <- c(genesCompr$dSNVOnly, geneName)
  else if (in_dMAF_N > 0)
    genesCompr$dMAFOnly <- c(genesCompr$dMAFOnly, geneName)
}

#summarize genesCompr
genesComprN <- lapply(genesCompr,length)
genesComprN <- unlist(genesComprN)
names(genesComprN) <- c("In Common", "In dSNV Only", "In dMAF Only")

#slice of data.frame for only genes in common
dSNVgenesComOnly <- dSNV[dSNV$geneSymbol %in% genesCompr$inCommon,]
dMAFgenesComOnly <- dMAF[dMAF$geneSymbol %in% genesCompr$inCommon,]

write.csv(dSNVgenesComOnly,file=paste(dsetAnalyzed_dir,"1.dSNVgenesComOnly.csv",sep=""))
write.csv(dMAFgenesComOnly,file=paste(dsetAnalyzed_dir,"1.dMAFgenesComOnly.csv",sep=""))

#generate summary plots
plotfn <- "basic_stat1"
openPDFdev(plots_dir, plotfn, savePlots)
barplot(genesComprN, ylab="Frequency (Genes)") #unique list of genes
saveFig(plots_dir, plotfn, savePlots)
dev.off()

plotfn <- "basic_stat2"
openPDFdev(plots_dir, plotfn, savePlots)
plot(genesComprCorr$dSNV, genesComprCorr$dMAF, log="xy",
     main=gettextf("r = %.2f",cor(genesComprCorr$dSNV, genesComprCorr$dMAF)),
     xlab="Freq of mutations for each gene (Mendelian)", ylab="Freq of mutations for each gene (TCGA)")
saveFig(plots_dir, plotfn, savePlots)
dev.off()

#plot per tumor type (note: time intensive)
# tumortypes <- unique(dMAF$tumorType)
# dimN = round(sqrt(length(tumortypes)))
# par(mfrow=c(dimN,dimN)) 
# genesComprCorrPerTumor = list()
# for(tumortype in tumortypes[1:3]){
#   for(geneName in geneList){
#     dMAFforGene <- dMAF[dMAF$geneSymbol==geneName,]
#     in_dMAF_N_perTumor <- sum(dMAFforGene$tumorType == tumortype)
#     genesComprCorrPerTumor[[eval(tumortype)]] <- c(genesComprCorrPerTumor[[eval(tumortype)]],in_dMAF_N_perTumor)
#   }
#   
#   plot(genesComprCorr$dSNV, genesComprCorrPerTumor[[eval(tumortype)]], log="xy",
#        main=gettextf("%s, r = %.2f",tumortype,cor(genesComprCorr$dSNV, genesComprCorr$dMAF)),
#        xlab="Freq of mutations for each gene (Mendelian)", ylab="Freq of mutations for each gene (TCGA)")
# }
# par() #reset par settings

######################################
catn("Examining mutations per base/residue position")
mutCompare <- function(dataset1, dataset2, field1, pos, ref, title){
  #Compare between dSNV and dMAF, using field1_pos as identifier for comparison
  #REQUIRES: dSNV, dMAF
  #          field1 and pos are fields in both datasets
  #          pos should be a position value; will ignore any entries where pos in dataset1 is empty
  #          order of dataset1 and dataset2 matters; mutCompare checks each entry of dataset1 in dataset2
  #OUTPUT: list of 1) genes common to both 2) aa position
  
  #create identifier, for comparison
  d1 <- paste(dataset1[[eval(field1)]],dataset1[[eval(pos)]],sep="_")
  d2 <- paste(dataset2[[eval(field1)]],dataset2[[eval(pos)]],sep="_")
  #compare the two datasets, for the same field1_field2 (e.g. gene_aaPos)
  mutSame <- d1 %in% d2 & dataset1[[eval(pos)]] != ""
  mutSameN = sum(mutSame)
  
  #now also create identifier with the ref position nt/residue
  d1b <- paste(dataset1[[eval(field1)]],dataset1[[eval(pos)]],dataset1[[eval(ref)]],sep="_")
  d2b <- paste(dataset2[[eval(field1)]],dataset2[[eval(pos)]],dataset2[[eval(ref)]],sep="_")
  #compare the two datasets, for the same field1_pos_ref (e.g. gene_aaPos_refAA)
  mutSameb <- d1b %in% d2b & dataset1[[eval(pos)]] != ""
  mutSamebN = sum(mutSameb)
  
  #a TRUE value here suggests that there is discrepancy between ref nt/AA between the two datasets
  error <- xor(mutSame,mutSameb)
  
  errorN <- sum(error)
  if(errorN > 0){
    message(gettextf("There are %.0f entries where ref nt/AA does not match and are removed from analysis", errorN))
  }
  
  genesInCommonN <- length(unique(dataset1[mutSameb & !error,"geneSymbol"]))
  genesDiffN <- length(unique(dataset1[!mutSameb & !error,"geneSymbol"]))
  data <- c(mutSamebN, nrow(dataset1[!error,])-mutSamebN)
  
  plotfn <- "basic_stat3"
  openPDFdev(plots_dir, plotfn, savePlots)
  pie(data,
      labels=c(gettextf("Mutations at same position\n(%.0f, %.1f%%) (%0.f genes)",
                        data[1], data[1]*100/sum(data), genesInCommonN),
               gettextf("Mutations at different positions\n(%.0f, %.1f%%) (%0.f genes)",
                        data[2],data[2]*100/sum(data), genesDiffN)
               ),
      main=title)
  saveFig(plots_dir, plotfn, savePlots)
  dev.off()
  
  return(list(geneSymbol=dataset1[mutSameb & !error,"geneSymbol"],
              pos=dataset1[mutSameb & !error, pos],
              refAA=dataset1[mutSameb & !error, ref]))
}

##############
#For genes in common in dMAF and dSNV
#same base position
# ntMutSameComOnly <- mutCompare(dSNVgenesComOnly, dMAFgenesComOnly, "chr","start","refAllele",
#                                "DNA, same position\n (for only genes mutated in both)")
#same residue position
aaMutSameComOnly <- mutCompare(dSNVgenesComOnly, dMAFgenesComOnly, "geneSymbol","aaPos","refAA",
                               "Amino acid, same position\n (for only genes mutated in both)")

#slice of data.frame for only genes (at least one mutated position is shared across datasets)
dSNVgenesMutSame <- dSNV[dSNV$geneSymbol %in% aaMutSameComOnly$geneSymbol,]
dMAFgenesMutSame <- dMAF[dMAF$geneSymbol %in% aaMutSameComOnly$geneSymbol,]

#slice of data.frame for only entries (with same gene and same position mutated in both datasets)
dSNVentryMutSame <- dSNV[paste(dSNV$geneSymbol,dSNV$aaPos,dSNV$refAA,sep="_") %in% 
                           paste(aaMutSameComOnly$geneSymbol,aaMutSameComOnly$pos,aaMutSameComOnly$ref,sep="_"),]
dMAFentryMutSame <- dMAF[paste(dMAF$geneSymbol,dMAF$aaPos,dMAF$refAA,sep="_") %in%
                           paste(aaMutSameComOnly$geneSymbol,aaMutSameComOnly$pos,aaMutSameComOnly$ref,sep="_"),]

write.csv(dSNVgenesMutSame,file=paste(dsetAnalyzed_dir,"2.dSNVgenesMutSame.csv",sep=""))
write.csv(dMAFgenesMutSame,file=paste(dsetAnalyzed_dir,"2.dMAFgenesMutSame.csv",sep=""))
write.csv(dSNVentryMutSame,file=paste(dsetAnalyzed_dir,"2.dSNVentryMutSame.csv",sep=""))
write.csv(dMAFentryMutSame,file=paste(dsetAnalyzed_dir,"2.dMAFentryMutSame.csv",sep=""))

##############
#For only genes mutated on the same residue in both dMAF and dSNV
#bases with similar residue properties
catn("Examining mutations per base/residue position, with same/similar AA properties...")

#define aa groups with similar properties
aa_acidic <- c("D","E")
aa_basic <- c("H","K","R")
aa_nonpolar <- c("A","C","G","I","L","M","F","P","W","V")
aa_polar <- c("N","Q","S","T","Y")
aa_stop <- c("*")
aaGroups <- list(aa_acidic, aa_basic, aa_nonpolar, aa_polar, aa_stop)

mutGroupCompr <- function(geneMutItem){
  #for each given entry in SNV, find the same gene and aa mutated position in dMAF
  #check to see if the variant mutated is in the same aa_group between given entry (from dSNV) and dMAF entry
  #INPUT: geneMutItem is an entry from dSNVentryMutSame, called from apply function
  #         geneMutItem is coming as single vector, since degenerate dimension is dropped;
  #         treat geneMutItem as list instead of as data.frame
  #OUTPUT: vector of boolean 1) exact match 2) shared properties 3) error

  aaMutExactSame <- FALSE #mutation is exactly the same between dSNV and dMAF
  aaMutSameGroup <- FALSE #mutation has similar aa properties between dSNV and dMAF
  error <- FALSE
  
  #slice of dataframe with same gene and aa mutated in dMAF, as that given in geneMutItem
  #when using apply for data frame, the input geneMutItem seems to get leading spaces for aaPos, needs to trim
  #the input values before concatenation
  dMAFentryMutSameItem <- dMAFentryMutSame[paste(dMAFentryMutSame$geneSymbol,dMAFentryMutSame$aaPos,sep="_") == 
                                             paste(str_trim(geneMutItem["geneSymbol"]),str_trim(geneMutItem["aaPos"]),sep="_"),]
  
  if(nrow(dMAFentryMutSameItem) < 1){
    error <- TRUE
  } else {
  
    #refSame <- geneMutItem['refAA'] == dMAFentryMutSameItem$refAA #this is already checked before in mutCompare
    varSame <- geneMutItem['varAA'] == dMAFentryMutSameItem$varAA
    
    if(sum(varSame) > 0){ #at least one of the variant mutations in dMAF is also found in dSNV
      aaMutExactSame <- TRUE      
    } else {
      for(aaGroup in aaGroups){
        #check if variant residue from both datasets are part of the same aa group
        #for dMAFentryMutSameItem$varAA, dMAFentryMutSameItem could have multiple entries,
        #  require just at least one to be in group to return TRUE
        boolN <- sum(dMAFentryMutSameItem$varAA %in% aaGroup)
        if(boolN > 0 && geneMutItem['varAA'] %in% aaGroup){
            aaMutSameGroup <- TRUE
            break
        }
      }
    }
    
  }
  
  return(c(aaMutExactSame, aaMutSameGroup, error))
}

#go through each entry in dSNVentryMutSame, and check for the mutated aa property to that in dMAFentryMutSame
aaMutSameGroup_bool <- apply(dSNVentryMutSame, 1, mutGroupCompr)
if(sum(aaMutSameGroup_bool[3,]) > 0){
  message("Something is wrong, one or more of the dSNVentryMutSame entries are not found in dMAFentryMutSame.")
}

#keep a list of gene, and position for entries with the properties is the same
aaMutExactSame <- dSNVentryMutSame[aaMutSameGroup_bool[1,],c("geneSymbol","aaPos","varAA")]
aaMutSameGroup <- dSNVentryMutSame[aaMutSameGroup_bool[2,],c("geneSymbol","aaPos","varAA")]
aaMutDiff <- dSNVentryMutSame[!(aaMutSameGroup_bool[1,] | aaMutSameGroup_bool[2,]),
                              c("geneSymbol","aaPos","varAA")]

#slice of data.frame for only genes (with at least one same mutation across both datasets)
dSNVgenesMutExactSame <- dSNV[dSNV$geneSymbol %in% aaMutExactSame$geneSymbol,]
dMAFgenesMutExactSame <- dMAF[dMAF$geneSymbol %in% aaMutExactSame$geneSymbol,]

#slice of data.frame for only genes (with at least one mutation with shared aa properties in both datasets)
dSNVgenesMutSameGroup <- dSNV[dSNV$geneSymbol %in% aaMutSameGroup$geneSymbol,]
dMAFgenesMutSameGroup <- dMAF[dMAF$geneSymbol %in% aaMutSameGroup$geneSymbol,]

#slice of data.frame for only entries (with same gene, same mutation with shared aa properties in both datasets)
dSNVentryMutExactSame <- dSNV[paste(dSNV$geneSymbol,dSNV$aaPos,dSNV$varAA,sep="_") %in%
                                paste(aaMutExactSame$geneSymbol,aaMutExactSame$aaPos,aaMutExactSame$varAA,sep="_"),]
dMAFentryMutExactSame <- dMAF[paste(dMAF$geneSymbol,dMAF$aaPos,dMAF$varAA,sep="_") %in%
                                paste(aaMutExactSame$geneSymbol,aaMutExactSame$aaPos,aaMutExactSame$varAA,sep="_"),]

#slice of data.frame for only entries (with same gene, same mutation with shared aa properties in both datasets)
dSNVentryMutSameGroup <- dSNV[paste(dSNV$geneSymbol,dSNV$aaPos,sep="_") %in%
                                paste(aaMutSameGroup$geneSymbol,aaMutSameGroup$aaPos,sep="_"),]
dMAFentryMutSameGroup <- dMAF[paste(dMAF$geneSymbol,dMAF$aaPos,sep="_") %in%
                                paste(aaMutSameGroup$geneSymbol,aaMutSameGroup$aaPos,sep="_"),]


write.csv(dSNVgenesMutExactSame,file=paste(dsetAnalyzed_dir,"3.dSNVgenesMutExactSame.csv",sep=""))
write.csv(dMAFgenesMutExactSame,file=paste(dsetAnalyzed_dir,"3.dMAFgenesMutExactSame.csv",sep=""))
write.csv(dSNVentryMutExactSame,file=paste(dsetAnalyzed_dir,"3.dSNVentryMutExactSame.csv",sep=""))
write.csv(dMAFentryMutExactSame,file=paste(dsetAnalyzed_dir,"3.dMAFentryMutExactSame.csv",sep=""))

write.csv(dSNVgenesMutSameGroup,file=paste(dsetAnalyzed_dir,"4.dSNVgenesMutSameGroup.csv",sep=""))
write.csv(dMAFgenesMutSameGroup,file=paste(dsetAnalyzed_dir,"4.dMAFgenesMutSameGroup.csv",sep=""))
write.csv(dSNVentryMutSameGroup,file=paste(dsetAnalyzed_dir,"4.dSNVentryMutSameGroup.csv",sep=""))
write.csv(dMAFentryMutSameGroup,file=paste(dsetAnalyzed_dir,"4.dMAFentryMutSameGroup.csv",sep=""))

mutExactSameN <- sum(aaMutSameGroup_bool[1,])
mutSameGroupN <- sum(aaMutSameGroup_bool[2,])

data <- c(mutExactSameN, mutSameGroupN, 
          length(aaMutSameGroup_bool[1,])-(mutExactSameN+mutSameGroupN))

dataGenes <- c(length(unique(dSNVgenesMutExactSame$geneSymbol)),
               length(unique(dSNVgenesMutSameGroup$geneSymbol)),
               length(unique(dSNV[dSNV$geneSymbol %in% aaMutDiff$geneSymbol,"geneSymbol"])))

plotfn <- "basic_stat4"
openPDFdev(plots_dir, plotfn, savePlots)
pie(data,
    labels=c(gettextf("Mutations with same AA\n(%.0f, %.1f%%) (%.0f genes)",
                      data[1], data[1]*100/sum(data), dataGenes[1]),
             gettextf("Mutations with similar AA properties\n(%.0f, %.1f%%) (%.0f genes)",
                      data[2], data[2]*100/sum(data), dataGenes[2]),
             gettextf("Mutations with different AA properties\n(%.0f, %.1f%%) (%.0f genes)",
                      data[3], data[3]*100/sum(data), dataGenes[3])
    ),
    main="Amino acid, similar properties\n (for only positions mutated in both SNV and MAF)")
saveFig(plots_dir, plotfn, savePlots)
dev.off()

#Save datasets
catn("Saving datasets...")
save(geneList, dMAF, dSNV, 
     dSNVgenesComOnly, dMAFgenesComOnly,
     dSNVgenesMutSame, dMAFgenesMutSame,
     dSNVentryMutSame, dMAFentryMutSame,
     dSNVgenesMutExactSame, dMAFgenesMutExactSame,
     dSNVentryMutExactSame, dMAFentryMutExactSame,
     dSNVgenesMutSameGroup, dMAFgenesMutSameGroup,
     dSNVentryMutSameGroup, dMAFentryMutSameGroup,
     file=paste(dsetAnalyzed_dir,"datasets_analyzed.RData",sep=""))
