#util.R
#shared utilities
#Author: Boyang Zhao

calcProgress <- function(idx, total, preval){
  #Calculates progress
  #DEPENDENCIES: preval needs to be initialized as global variable
  
  if(total < 1){
    message("Warning: in progress calculation, total is zero")
  } else if(idx > total){
    message("Warning: in progress calculation, index is greater than total")
  } else {
    val <- floor(idx*100/total)
    if(val %% 2 == 0 && preval != val){
      catn(paste(val,"% complete...",sep=""))
      preval <<- val
    }
  }
}

## File handling methods
createDir <- function(dirname){
  #create directory if it doesn't exist
  
  #need to remove the trailing "/" for dir before checking for dir existence
  if(!file.exists(sub("^(.*)/$","\\1",dirname))){
    dir.create(dirname)
  }
}

checkDirExists <- function(dirname){
  #check if direcotry exists, if not throws error message
  if(!file.exists(sub("^(.*)/$","\\1",dirname))){
    stop(gettextf("Directory '%s' does not exist...\nStopping execution...", dirname))
  }
}

tumor2filename <- function(tumor){
  #convert given tumor value to a filename-friendly value, if needed
  if(tumor == "COAD/READ"){
    return("COADREAD")
  }
  
  return(tumor)
}

catn <- function(txt){
  #append line break to text and call cat
  cat(paste(txt,"\n",sep=""))
}

## Plotting methods
saveFig <- function(figdir, figname, savePlots, height=800, width=800, res=150){
  #save figure given directory and filename
  if(savePlots){
    dev.print(png, height=height, width=width, units="px", res=res, filename=gettextf('%s/%s.png', figdir, figname))
  }
}

openPDFdev <- function(figdir, figname, savePlots, width=7, height=7){
  #opens PDF device
  if(savePlots){
    pdf(gettextf('%s/%s.pdf', figdir, figname), width=width, height=height)
    dev.control("enable") #turn on recording so dev can be copied
  }
}

plotMutationsInDatasets <- function(dRef, dMut, dProteinLengths, plots_dir = "", savePlots = FALSE){
  #plots mutations in given datasets
  #dRef: reference dataset (e.g. dMAF)
  #dMut: annoatation dataset, for match (e.g. dSNV)
  #dProteinLengths: protein lengths dataset
  
  library(mutationsplotter)
  #if mutationsplotter is not installed, uncomment lines below
  #library(devtools)
  #install_github('boyangzhao/mutationsplotter')
  
  catn("Plotting mutations in given datasets...")
  
  #quality controls
  createDir(plots_dir)
  
  genes <- unique(dRef$geneSymbol)
  
  if(length(genes) < 1){
    message("There is nothing to plot. Reference dataset is empty.")
  } else {
    for(geneSymbol in genes){
      dataslice <- dRef[dRef$geneSymbol==geneSymbol,]
      tally <- aggregate(dataslice$aaPos, by=list(dataslice$aaPos), FUN=length)
      mutPos <- tally[[1]]
      posCount <- tally[[2]]
      matchPos <- unique(dMut[dMut$geneSymbol==geneSymbol,'aaPos'])
      protLen <- dProteinLengths[dProteinLengths$geneSymbol == geneSymbol, 'length']
        
      plotfn <- paste("mut", geneSymbol, sep=".")
      openPDFdev(plots_dir, plotfn, savePlots, height=4, width=8)
      plotMutations(mutPos, posCount, protLen=protLen, 
                    annotatePos=matchPos, annotateSymbol=NULL, annotateCol='default', 
                    cex=1.5)
      saveFig(plots_dir, plotfn, savePlots, height=700, width=1500, res=200)
      dev.off()
    }
  }
}

