suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyquant))  # Loads tidyverse, tidquant, financial pkgs, xts/zoo
suppressPackageStartupMessages(library(dataRetrieval))
suppressPackageStartupMessages(library(gmailr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))

gageInfo <- read_csv('data/USGSgageStationInformation.csv')


## Functions ##

# Pull NWIS data
NWISpull <- function(gageNo,start,end){
  allUnitData <- readNWISuv(siteNumbers=gageNo,
                            parameterCd=c("00010", "00095", "00300", "00400", "63680"),
                            startDate=as.Date(start,"%Y-%m-%d"),
                            endDate=as.Date(end,"%Y-%m-%d"),
                            tz='America/New_York')
  allUnitData <- renameNWISColumns(allUnitData)
}






for(i in 1:24){
  print(i)
  # Pull all gage data available since it began recording
  gageData <- NWISpull(paste(0,gageInfo[i,6],sep=""), "2017-08-01",Sys.Date())%>%
    dplyr::select(agency_cd,site_no,dateTime,'Turb_Inst') %>%
    tq_mutate(
      select     = Turb_Inst,
      mutate_fun = rollapply, 
      # rollapply args
      width      = 7,
      align      = "right",
      by.column  = FALSE,
      FUN        = max,
      # FUN args
      na.rm      = T,
      # tq_mutate args
      col_rename = "rw30_max")
  
  gageInfo$Turbidity_90th[i] <- quantile(gageData$Turb_Inst,prob=0.9,na.rm=T)
  gageInfo$Turbidity_95th[i] <- quantile(gageData$Turb_Inst,prob=0.95,na.rm=T)
  gageInfo$Turbidity_99th[i] <- quantile(gageData$Turb_Inst,prob=0.99,na.rm=T)
  gageInfo$rw30_90th[i] <- quantile(gageData$rw30_max,prob=0.90,na.rm=T)
  gageInfo$rw30_95th[i] <- quantile(gageData$rw30_max,prob=0.95,na.rm=T)
  gageInfo$rw30_99th[i] <- quantile(gageData$rw30_max,prob=0.99,na.rm=T)

  
}

#write.csv(gageInfo,'data/gageInfoEVJ.csv')

