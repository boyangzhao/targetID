#parseProteome.R
#Merge uniprot_proteome (fasta file) with HGNC
#Generates a list of protein lengths based on hugo symbol
#Procedure: use HGNC to look up the uniprotID for each hugo symbol, and then use uniprotID to lookup the protein
#sequence in fasta file

library("seqinr")
prot <- read.fasta(file="homo_sapiens_proteome.fasta")
hgnc <- read.csv("HGNC.csv",header=TRUE,stringsAsFactors=FALSE)


#parse fasta to create a table with uniprotID and protein seq length
parsefunc <- function(x){
  id <- attr(x,'name')
  uniprotid <- sub("^.*\\|(.*)\\|.*_.*$", "\\1", id)
  protlength <- length(x)
  
  list(id=uniprotid, length=protlength)
}

protlist_tmp <- vapply(prot, parsefunc, list('a',0))
protlist <- data.frame(UniProtID=unlist(protlist_tmp[1,]), length=unlist(protlist_tmp[2,]))

#for each HGNC gene, lookup the protein seq length
hgnc_approved <- hgnc[hgnc$Status == 'Approved',]
hgnc_approved <- hgnc_approved[!duplicated(hgnc_approved$geneSymbol),]
hgnc_approved <- hgnc_approved[hgnc_approved$UniProtID != "",]

getProtLength <- function(x){
  protlist[protlist$UniProtID==x,'length']
}

dataset.merged<- merge(hgnc_approved, protlist, by="UniProtID")

dProteinLengths <- data.frame(geneSymbol=dataset.merged$geneSymbol,
                              UniProtID=dataset.merged$UniProtID,
                              refSeqID=dataset.merged$RefSeqID,
                              length=dataset.merged$length)
write.csv(dProteinLengths,file="dProteinLengths.csv")
save(dProteinLengths, file="dProteinLengths.RData")