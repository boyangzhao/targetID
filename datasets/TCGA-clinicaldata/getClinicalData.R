#getClinicalData.R
#Author: Boyang Zhao
#Retrieve clinical data from cBioPortal

rm(list=ls())
library(cgdsr)

mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
cancerstudies <- getCancerStudies(mycgds)

dClinicalDataRaw <- c()
for(cancerstudy in cancerstudies[,1]){
  #only focus on tcga
  datasource <- sub("^[a-z0-9]*_([a-z0-9_]*)$", "\\1", cancerstudy)
  
  if(datasource == "tcga" || datasource == "tcga_pub"){
    print(gettextf("Retrieving clinical data for %s...", cancerstudy))
    
    tumorType <- toupper(sub("^([a-z0-9]*)_[a-z0-9_]*$", "\\1", cancerstudy))
    
    #retrieving clinical data
    clinicaldata <- getClinicalData(mycgds, paste(cancerstudy,'_all',sep=""))
    
    len <- nrow(clinicaldata)
    if(len > 0){
      #collect and add to dClinicalDataRaw data structure
      write.csv(clinicaldata,file=gettextf("individualTumors/%s.csv",cancerstudy))
      save(clinicaldata, file=gettextf("individualTumors/%s.RData",cancerstudy))
      
      dataToAdd <- data.frame(caseID = row.names(clinicaldata),
                              age = ifelse(rep(!is.null(clinicaldata$AGE),len), clinicaldata$AGE, rep(NA,len)),
                              gender = ifelse(rep(!is.null(clinicaldata$GENDER),len), clinicaldata$GENDER, rep(NA,len)),
                              race = ifelse(rep(!is.null(clinicaldata$RACE),len), clinicaldata$RACE, rep(NA,len)),
                              dfs_months = ifelse(rep(!is.null(clinicaldata$DFS_MONTHS),len), clinicaldata$DFS_MONTHS, rep(NA,len)),
                              dfs_status = ifelse(rep(!is.null(clinicaldata$DFS_STATUS),len), clinicaldata$DFS_STATUS, rep(NA,len)),
                              os_months = ifelse(rep(!is.null(clinicaldata$OS_MONTHS),len), clinicaldata$OS_MONTHS, rep(NA,len)),
                              os_status = ifelse(rep(!is.null(clinicaldata$OS_STATUS),len), clinicaldata$OS_STATUS, rep(NA,len)),
                              tumorType = tumorType,
                              cancerStudy=cancerstudy,
                              stringsAsFactors = FALSE)
      dClinicalDataRaw <- rbind(dClinicalDataRaw, dataToAdd)
    }
  }
}

#format change for caseID, change "." to "-"
#format is as follows: TCGA.TSS.participantID.sampleID
#drop the sampleID
dClinicalData <- dClinicalDataRaw
dClinicalData$caseID <- sub("^([A-Za-z0-9]*)\\.([A-Za-z0-9]*)\\.([A-Za-z0-9]*)\\.([A-Za-z0-9]*)$", "\\1-\\2-\\3", dClinicalData$caseID)

#remove duplicates
dClinicalData_unique <- dClinicalData[!duplicated(dClinicalData$caseID),1:ncol(dClinicalData)-1]

#save data
write.csv(dClinicalDataRaw,file="dClinicalDataRaw.csv")
write.csv(dClinicalData,file="dClinicalData.csv")
write.csv(dClinicalData_unique,file="dClinicalData_unique.csv")
save(dClinicalDataRaw, dClinicalData,dClinicalData_unique,file="datasets_clinicalData.RData")
