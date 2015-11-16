#util_misc.R
#miscellaneous utilities
#Author: Boyang Zhao

correctPTPN11 <- function(dMAF){
  #PTPN11 was mapped to isoform 2 in cBioPortal but isoform 1 in HUMSAVAR
  #here, isoform naming is based on UniProt database; need to correct aaPos,and aaChange so both datasets refer to
  #consistent protein (in this case, isoform 1)
  
  PTPN11_seq = "MTSRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLNCADPTSERWFHGHLSGKEAEKLLTEKGKHGSFLVRESQSHPGDFVLSVRTGDDKGESNDGKSKVTHVMIRCQELKYDVGGGERFDSLTDLVEHYKKNPMVETLGTVLQLKQPLNTTRINAAEIESRVRELSKLAETTDKVKQGFWEEFETLQQQECKLLYSRKEGQRQENKNKNRYKNILPFDHTRVVLHDGDPNEPVSDYINANIIMPEFETKCNNSKPKKSYIATQGCLQNTVNDFWRMVFQENSRVIVMTTKEVERGKSKCVKYWPDEYALKEYGVMRVRNVKESAAHDYTLRELKLSKVGQALLQGNTERTVWQYHFRTWPDHGVPSDPGGVLDFLEEVHHKQESIMDAGPVVVHCSAGIGRTGTFIVIDILIDIIREKGVDCDIDVPKTIQMVRSQRSGMVQTEAQYRFIYMAVQHYIETLQRRIEEEQKSKRKGHEYTNIKYSLADQTSGDQSPLPPCTPTPPCAEMREDSARVYENVGLMQQQKSFR"
  
  #408-411 is missing in isoform 2
  #for positions 408 and above in isoform 2, it is +4 in reference to isoform 1
  #therefore, add +4 to all positions greater than 407 for PTPN11 in dMAF
  idxToCorrect <- dMAF$geneSymbol=="PTPN11" & as.numeric(dMAF$aaPos)>407
  dMAF[idxToCorrect,"aaPos"] <- as.character(as.numeric(dMAF[idxToCorrect,"aaPos"])+4)
  
  dMAF[idxToCorrect,"aaChange"] <- paste(dMAF[idxToCorrect,"refAA"],
                                         dMAF[idxToCorrect,"aaPos"],
                                         dMAF[idxToCorrect,"varAA"],
                                         sep="")
  return(dMAF)
}
