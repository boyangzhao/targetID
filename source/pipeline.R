#pipeline.R
rm(list=ls())
config.avail = TRUE
lockBinding("config.avail", globalenv())


#######################################
# Configurations
#directory locations
config.dset_dir <- "../datasets/"
config.dsetProcessed_dir <- paste(config.dset_dir, "pipeline_test/processed/", sep="")
config.dsetAnalyzed_dir <- paste(config.dset_dir, "pipeline_test/analyzed/", sep="")
config.dsetScored_dir <- paste(config.dset_dir, "pipeline_test/scored/", sep="") 
config.dsetSyn_dir <- paste(config.dset_dir, "pipeline_test/synthetic/", sep="") #for synthetic_datasets.RData (if syntheticData=TRUE)

config.dsetCalc_dir <- paste(config.dsetAnalyzed_dir, "calc/", sep="") #for dMultipliers.RData
config.dsetTCGAMutSigCV_dir <- paste(config.dset_dir, "MutSigCV TCGA/", sep="") #for MutSigCV files
config.dsetProteome_dir <- paste(config.dset_dir, "UniProt proteome/", sep="") #for dProteinLengths.RData

config.dsetScoredSummary_dir <- gettextf('%ssummary/', config.dsetScored_dir)
config.dsetSurvivalSummary_dir <- gettextf('%ssummary/stratatumorTypeTRUE_mut_age_gender/', config.dsetScored_dir)

config.plotsIndCancer_dir <- gettextf('%splotsIndCancer/', config.dsetScored_dir)
config.plotsPanCancer_dir <- gettextf('%splotsPanCancer/', config.dsetScored_dir)
config.plotsSurvival_dir <- gettextf('%splotsSurvival/stratatumorTypeTRUE_mut_age_gender/', config.dsetScored_dir)

#for analyze_basic.R
config.basic.savePlots <- TRUE
config.basic.processedDatasetName <- "pipeline_test.RData"

#for calc_multiplier.R
config.calc.minRNAseqcutoff <- NULL #if NULL, then will not have any cutoffs and value of analyzeNotExpressed is irrelevant
config.calc.analyzeNotExpressed <- TRUE
config.calc.savePlots <- TRUE
config.calc.datasetName <- "dMAFgenesMutExactSame"
config.calc.dataset_matchName <- "dSNVentryMutExactSame"

#for analyze_scoring.R
config.scoring.bootstrapN <- 1e3 #number of boostrap trials
config.scoring.permut_trialsN <- c(1e4,1e6) #number of permutation trials (for null distribution generation); a vector would mean [min max]
config.scoring.useMultiplier <- TRUE #if true, then overwrite mutRate_multiplier using dataset from dMultipliers.RData
config.scoring.mutRate_multiplier <- 10 #default mutation rate multiplier; unless overwritten (if useMultiplier is true)
config.scoring.minCasesNcutoff <- 0 #minimum number cases as cutoff (>=) to write to _sigListGene.csv
config.scoring.minRNAseqcutoff <- 5 #minimum normalized_RSEM values (upper quartile normalized) as cutoff (>=) for mRNA expression
config.scoring.analyzeNotExpressed <- NULL #if NULL, will partition dataset and analyze twice (once on expressed and once on not expressed)
config.scoring.syntheticData <- FALSE #using synthetic dataset instead
config.scoring.tumorTypesToAnalyze <- c() #if empty, analyze all tumor types, otherwise analyze the specified tumor type
config.scoring.workingDatasetName <- "dSNVentryMutSame"

#for analyze_scored.R
config.scored.savePlots <- TRUE
config.scored.saveLog <- TRUE
config.scored.minCasesNcutoff <- 2 #NULL means use all of given data; or impose additional cutoff on minimum number of cases (>=)
config.scored.MutSigCutoff <- 0.05 #minimum cutoff for MutSigCV qval (>=)
config.scored.probSigCutoff <- 0.1 #default minimum cutoff for probability score pval/qval (>=)
config.scored.SNRcutoff <- 2 #minimum cutoff for SNR value (>=)
config.scored.usePairedExprData <- TRUE #if true, will search for paired paritioned dataset results (expressed | not expressed),
                                 #and use not expressed results to set the cuttoff pval/qval for prob (probSigCutoff)
config.scored.tumorTypesToAnalyze <- c() #if empty, analyze all tumor types, otherwise analyze the specified tumor type

#for analyze_scoredPanCancer.R
config.pancancer.savePlots <- TRUE
config.pancancer.saveLog <- FALSE
config.pancancer.tumorTypesExclude <- c("CELLLINE", "ESCA", "PAAD", "LAML") #tumor types to exclude

#for analyze_scoredSurvival.R
config.surv.saveLog <- TRUE
config.surv.savePlots <- TRUE #for survivalAnalysis_indv
config.surv.strata_tumorType <- TRUE #for survivalAnalysis_indv
config.surv.coxmodel <- "mut_age_gender" #for survivalAnalysis_indv; possible values 'mut', 'mut_age', 'mut_age_gender'
config.surv.verboseOutput <- TRUE #for survivalAnalysis_indv
config.surv.tumorTypesExclude <- c("CELLLINE", "ESCA", "PAAD") #tumor types to exclude

#######################################

message("[1/6] Running basic summary statistics...")
source("analyze_basic.R")
message("[2/6] Calculating proportionality coefficients...")
source("calc_multiplier.R")
message("[3/6] Calculating scores...")
source("analyze_scoring.R")
message("[4/6] Analyzing scored results...")
source("analyze_scored.R")
message("[5/6] Performing pan-cancer analyses...")
source("analyze_scoredPanCancer.R")
message("[6/6] Performing survival analyses...")
source("analyze_scoredSurvival.R")

