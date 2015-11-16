#calc_synData.R
#Generates synthetic dataset for in silico analyses
#Author: Boyang Zhao

rm(list=ls())
dset_dir <- "../datasets/"
synData_dir <- gettextf('%ssynthetic/', dset_dir) #output directory
source("util.R")

#settings
savePlots <- TRUE

#initialize internal variables
preval <- 0

gcd <- function(val){
  #determine GCD using euclidean method for a pair of values
  #val is a two-element vector, with first element >= second element
  r <- val[1]%%val[2]
  if(r){
    gcd(c(val[2],r))
  } else {
    return(val[2])
  }
}

gcdSet <- function(val){
  #determine GCD for a group of numbers
  #performs pairwise gcd calculation, and determine the minimum gcd
  #val is a vector of values
  
  gcdVals <- c()
  n <- length(val)
  pairwiseCompr.indices <- combn(1:n,2)
  for(idx in 1:ncol(pairwiseCompr.indices)){
    val.indices <- pairwiseCompr.indices[,idx]
    gcdVals <- c(gcdVals, gcd(val[val.indices])) #append the gcd determined for this pair to 'gcdVals' list
  }
  
  return(min(gcdVals))
}

adjSetbyGCD <- function(valSet, multiFactor){
  #adjust set of values using gcd
  valSet.adj <- ceiling(valSet*multiFactor) #scale by factor, and will ceil all values (so tolerance need to be defined by multiFactor)
  gcd <- gcdSet(valSet.adj)
  valSet.adj <- valSet.adj/gcd #change set to integers by dividing by gcd
  
  return(valSet.adj)
}

generateEvents <- function(protLen,
                           matchedN.unique,
                           probScore,
                           matchDistr, #default uniform distribution
                           bgNoise,
                           totalEventsMultiplier,
                           nonmatchDistr = 'uniform',
                           mutN.unique.max = NULL,
                           sampling=FALSE){
  #Generates the following variables
  #INPUTS:
  # protLen: protein length
  # matchedN.unique: number of unique positions that are matches
  # probScore: probability score
  # matchDistr: distribution of the matched positions; lenght of this vector need to be equal to value of 'matchedN.unique'
  #             if 'matchDistr' is empty, will equal to an uniform distribution
  # bgNoise: background noise, ranges [0,1]; any values falling outside of this range will default to the boundary values
  # totalEventsMultiplier: amplify the signal by boosting all the events uniformly by this value
  # nonmatchDistr: nonmatch distribution, takes the following takes: 'uniform'
  # mutN.unique.max: maximum number of mutations; if NULL, then limit is only constrained by other parameters
  #                  if value here is smaller than matchedN.unique, default value will be matchedN.unique + 1
  # sampling: if TRUE, then entries will be stochastically determined with sampling drawn from distribution
  #OUTPUTS:
    #matchedN
    #matchedN.unique
    #nonmatchedN
    #nonmatchedN.unique
    #
    #pos.matched
    #pos.matched.unique
    #pos.nonmatched
    #pos.nonmatched.unique
  
  #quality controls
  if(totalEventsMultiplier < 1){ totalEventsMultiplier <- 1 }
  if(bgNoise > 1){ bgNoise <- 1}
  if(bgNoise < 0){ bgNoise <- 0}
  if(length(matchDistr) < 1){
    matchDistr = rep(1/matchedN.unique, matchedN.unique) #default to uniform distribution
  }
  if(matchedN.unique != length(matchDistr)){
    message("matchedN.unique value needs to match length of matchDistr")
  }
  if(protLen < 1){ protLen <- 100}
  if(!is.null(mutN.unique.max) && mutN.unique.max <= matchedN.unique){
    mutN.unique.max <- matchedN.unique + 1
  }
  
  #####
  ## Per unit calculation
  #####
  # Determine entries for matched positions
  posToDrawFrom <- 1:protLen
  
  #####
  # Calculate event counts, and adjust counts accordingly to ensure all values are integers
  if(matchedN.unique > 1){
    matchDistr.adj <- adjSetbyGCD(matchDistr, 100)
    matchedN <- sum(matchDistr.adj)
  } else {
    matchDistr.adj <- c(1)
    matchedN <- 1
  }

  nonmatchedN <- (matchedN/probScore) - matchedN #total number of entries as not matches given probScore constraint
  
  if(nonmatchedN != 0){
    if(nonmatchedN < 1){
      multifactor <- 1/nonmatchedN #smallest factor needed to make sure there is at least one event for not matched
      multifactor <- 10^ceiling(log10(1/nonmatchedN)) #turn this into multiplies of 10
    } else {
      multifactor <- 10 #will take care of any decimals in nonmatchedN just in case; this effectively also defines a tolerance
    }
    distrSet.adj <- adjSetbyGCD(c(matchDistr.adj, nonmatchedN), multifactor)
  
    len <- length(distrSet.adj)
    matchDistr.adj <- distrSet.adj[1:len-1] #updated matchDistr.adj after consideration of nonmatchedN
    nonmatchedN <- distrSet.adj[len] #updated nonmatchedN
    matchedN <- sum(matchDistr.adj) #updated matchedN
  }

  #####
  # Determine matched positions
  pos.matched.unique <- sample(posToDrawFrom, matchedN.unique, replace=FALSE)
  if(matchedN>1){
    if(sampling && length(pos.matched.unique) > 1){
      pos.matched <- sample(pos.matched.unique, matchedN, prob=matchDistr.adj, replace=TRUE) #stochastic
    } else {
      pos.matched <- rep(pos.matched.unique, matchDistr.adj) #deterministic
    }
  } else {
    #matchedN=1, i.e. there is just one match total
    pos.matched <- pos.matched.unique
  }

  #####
  # Determine entries for not matched positions
  nonmatchedN.unique <- 0
  pos.nonmatched <- c()
  pos.nonmatched.unique <- c()
  if(nonmatchedN != 0){
    posToDrawFrom <- posToDrawFrom[!(posToDrawFrom %in% pos.matched.unique)]
    nonmatchedN.unique <- round(bgNoise * nonmatchedN)
    
    #checking bounds for nonmatchedN.unique
    if(nonmatchedN.unique<1){ nonmatchedN.unique <- 1} #lower bound
    if(nonmatchedN.unique>length(posToDrawFrom)){ nonmatchedN.unique <- length(posToDrawFrom)} #upper bound
    if(!is.null(mutN.unique.max)){
      nonmatchedN.unique.allowable <- max(1,mutN.unique.max - matchedN.unique)
      if(nonmatchedN.unique > nonmatchedN.unique.allowable){
        #based on the maximum number of allowable mutations set, and the already existing number of unique positions
        #mutated that is a match, the number of unique poisitions that is not a match need to be readjusted here
        #to satisfy constraint
        nonmatchedN.unique <- nonmatchedN.unique.allowable
      }
    }
    
    #####
    # Determine non-matched positions
    pos.nonmatched.unique <- sample(posToDrawFrom, nonmatchedN.unique, replace=FALSE)
    
    switch(nonmatchDistr,
           uniform={
             if(nonmatchedN>1){
               if(sampling && length(pos.nonmatched.unique) > 1){
                 pos.nonmatched <- sample(pos.nonmatched.unique, nonmatchedN, replace=TRUE) #stochastic
               } else {
                 pos.nonmatched <- rep(pos.nonmatched.unique, nonmatchedN/nonmatchedN.unique) #deterministic
               }
             } else {
               #nonmatchedN=1, i.e. there is just one non-match total
               pos.nonmatched <- pos.nonmatched.unique
             }
           }
    )
  }

  #replicate synthetic entries by 'totalEventsMultiplier' to boost total number of events
  pos.matched <- rep(pos.matched, totalEventsMultiplier)
  pos.nonmatched <- rep(pos.nonmatched, totalEventsMultiplier)
  matchedN <- matchedN*totalEventsMultiplier
  nonmatchedN <- nonmatchedN*totalEventsMultiplier
  
  events <- list(protLen=protLen,
                 matchedN = matchedN,
                 matchedN.unique = matchedN.unique,
                 nonmatchedN = nonmatchedN,
                 nonmatchedN.unique = nonmatchedN.unique,
                 pos.matched = pos.matched,
                 pos.matched.unique = pos.matched.unique,
                 pos.nonmatched = pos.nonmatched,
                 pos.nonmatched.unique = pos.nonmatched.unique)
  
  return(events)
}

convertEventsToDataset <- function(events, geneSymbol, datasets=c()){
  #INPUT: events from generateEvents, geneSymbol
  #       if datasets is given, will append newly generated dataset from events to this
  #OUTPUT datasets
  
  aaPos.dMAF <- c(events$pos.matched, events$pos.nonmatched)
  dMAF.len <- length(aaPos.dMAF)
  aaPos.dSNV <- events$pos.matched
  dSNV.len <- length(aaPos.dSNV)
  
  dMAF <- data.frame("geneSymbol"=rep(geneSymbol, dMAF.len),
                     "chr"=rep('-', dMAF.len),
                     "start"=rep('-', dMAF.len),
                     "end"=rep('-', dMAF.len),
                     "refAllele"=rep('-', dMAF.len),
                     "varAllele"=rep('-', dMAF.len),
                     "aaPos"=aaPos.dMAF,
                     "refAA"=rep('A', dMAF.len),
                     "varAA"=rep('B', dMAF.len),
                     "origin"=rep('-', dMAF.len),
                     "mutationType"=rep('Missense_Mutation', dMAF.len),
                     "aaChange"=rep('-', dMAF.len),
                     "tumorType"=rep('-', dMAF.len),
                     "caseID"=rep('-', dMAF.len),
                     "mRNAvals"=rep(1, dMAF.len),
                     "mRNAvals_final"=rep(1, dMAF.len),
                     stringsAsFactors=FALSE)

  dSNV <- data.frame("geneSymbol"=rep(geneSymbol, dSNV.len),
                     "aaPos"=aaPos.dSNV,
                     "refAA"=rep('A', dSNV.len),
                     "varAA"=rep('B', dSNV.len),
                     "aaChange"=rep('-', dSNV.len),
                     "SPID"=rep('-', dSNV.len),
                     "FTID"=rep('-', dSNV.len),
                     "dsSNP"=rep('-', dSNV.len),
                     "diseaseName"=rep('-', dSNV.len),
                     stringsAsFactors=FALSE)

  dProteinLengths <- data.frame("geneSymbol"=geneSymbol,
                                "UniProtID"="-",
                                "refSeqID"="-",
                                "length"=events$protLen,
                                stringsAsFactors=FALSE)
  
  if(length(datasets) < 1){
    datasets <- list(dMAF=c(),dSNV=c(),dProteinLengths=c())
  }
  
  datasets$dMAF <- rbind(datasets$dMAF, dMAF)
  datasets$dSNV <- rbind(datasets$dSNV, dSNV)
  datasets$dProteinLengths <- rbind(datasets$dProteinLengths, dProteinLengths)
  
  return(datasets)
}

generateSynData_cases <- function(){
  #generates synthetic dMAF and dSNV datasets
  catn("Generating synthetic datasets (based on select cases)...")
  
  #format: generateEvents(protLen, matchedN.unique, probScore, matchDistr, bgNoise, totalEventsMultiplier)
  
  #considerations
  #base case
  events <- generateEvents(100, 2, 0.8, c(0.5,0.5), 0, 1)
  datasets <- convertEventsToDataset(events, "Gene1")
  
  #same prob score, but different distributions of match, low background
  events <- generateEvents(100, 2, 0.8, c(0.1,0.9), 0, 1)
  datasets <- convertEventsToDataset(events, "Gene2", datasets)
  
  #same prob score, different protein size
  events <- generateEvents(500, 2, 0.8, c(0.5,0.5), 0, 1)
  datasets <- convertEventsToDataset(events, "Gene3", datasets)
  
  #same prob score, but different distributions of match, mid background
  events <- generateEvents(100, 2, 0.8, c(0.5,0.5), 0.5, 1)
  datasets <- convertEventsToDataset(events, "Gene4", datasets)
  
  #same prob score, but different distributions of match, high background
  events <- generateEvents(100, 2, 0.8, c(0.5,0.5), 1, 1)
  datasets <- convertEventsToDataset(events, "Gene5", datasets)
  
  #same prob score, high SNR
  events <- generateEvents(100, 2, 0.8, c(0.5,0.5), 0, 10)
  datasets <- convertEventsToDataset(events, "Gene6", datasets)
  
  #one out of one match
  events <- generateEvents(100, 1, 1, c(), 0, 1)
  datasets <- convertEventsToDataset(events, "Gene7", datasets)
  
  dMAF <- datasets$dMAF
  dSNV <- datasets$dSNV
  dProteinLengths <- datasets$dProteinLengths
  
  return(list(dMAF=dMAF, dSNV=dSNV, dProteinLengths=dProteinLengths))
}


generateSynData_sampling <- function(N=1000){
  #generates synthetic dataset through latin hypercube sampling
  #INPUTS:
  # N: number of samples
  
  catn("Generating synthetic datasets (based on sampling)...")
  
  #format: generateEvents(protLen, matchedN.unique, probScore, matchDistr, bgNoise, nonmatchDistr, totalEventsMultiplier)
  
  library(lhs)
  paramSet <- maximinLHS(N, 5)
  protLen.min <- 10
  protLen.max <- 8000
  matchedN.unique.min <- 1 # unit: number of positions
  matchedN.unique.max <- 100 #unit: number of positions
  matchedN.unique.maxPercent <- 0.8 #unit: percent; % of maximum possible (i.e. protein length)
  probScore.min <- 0.01
  probScore.max <- 0.99
  bgNoise.min <- 0
  bgNoise.max <- 1
  totalEventsMultiplier.min <- 1
  totalEventsMultiplier.max <- 50
  mutN.unique.max.val <- 100
  
  paramSet[,1] <- round(paramSet[,1]*(protLen.max-protLen.min)+protLen.min)
  paramSet[,2] <- round(paramSet[,2]*(
    min(paramSet[,1]*matchedN.unique.maxPercent,matchedN.unique.max) - matchedN.unique.min)+matchedN.unique.min
    )
  paramSet[,3] <- paramSet[,3]*(probScore.max-probScore.min)+probScore.min
  paramSet[,4] <- paramSet[,4]*(bgNoise.max-bgNoise.min)+bgNoise.min
  paramSet[,5] <- round(paramSet[,5]*(totalEventsMultiplier.max-totalEventsMultiplier.min)+totalEventsMultiplier.min)
  
  datasets <- c()
  for(idx in 1:N){
    calcProgress(idx, N, preval)
    events <- generateEvents(paramSet[idx,1], #protLen
                             paramSet[idx,2], #matchedN.unique
                             paramSet[idx,3], #probScore
                             c(), #matchDistr
                             paramSet[idx,4], #bgNoise
                             paramSet[idx,5], #totalEventsMultiplier
                             mutN.unique.max = mutN.unique.max.val,
                             sampling = TRUE)
    datasets <- convertEventsToDataset(events, gettextf("Gene%s",idx), datasets)
  }
  
  dMAF <- datasets$dMAF
  dSNV <- datasets$dSNV
  dProteinLengths <- datasets$dProteinLengths
  
  return(list(dMAF=dMAF, dSNV=dSNV, dProteinLengths=dProteinLengths))
}

###################################################
# Main
###################################################
#generate synthetic data (simple)
r<-generateSynData_cases()
dMAF <- r$dMAF
dSNV <- r$dSNV
dProteinLengths <- r$dProteinLengths
createDir(paste(synData_dir, "cases/",sep=""))

plots_dir <- gettextf('%scases/plots/', synData_dir)
plotMutationsInDatasets(dMAF, dSNV, dProteinLengths, plots_dir, savePlots)

catn("Saving synthetic datasets...")
write.csv(dMAF,file=paste(synData_dir, "cases/", "dMAF_cases.csv",sep=""))
write.csv(dSNV,file=paste(synData_dir, "cases/", "dSNV_cases.csv",sep=""))
write.csv(dProteinLengths,file=paste(synData_dir, "cases/", "dProteinLengths_cases.csv",sep=""))
save(dMAF, dSNV, dProteinLengths, file=paste(synData_dir, "cases/", "synthetic_dataset.RData",sep=""))

#generate synthetic data (sampling)
r<-generateSynData_sampling(1000)
dMAF <- r$dMAF
dSNV <- r$dSNV
dProteinLengths <- r$dProteinLengths
createDir(paste(synData_dir, "sampling/",sep=""))

plots_dir <- gettextf('%ssampling/plots/', synData_dir)
plotMutationsInDatasets(dMAF, dSNV, dProteinLengths, plots_dir, savePlots)

catn("Saving synthetic datasets...")
write.csv(dMAF,file=paste(synData_dir, "sampling/", "dMAF_sampling.csv",sep=""))
write.csv(dSNV,file=paste(synData_dir, "sampling/", "dSNV_sampling.csv",sep=""))
write.csv(dProteinLengths,file=paste(synData_dir, "sampling/", "dProteinLengths_sampling.csv",sep=""))
save(dMAF, dSNV, dProteinLengths, file=paste(synData_dir, "sampling/", "synthetic_dataset.RData",sep=""))
