suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyquant))  # Loads tidyverse, tidquant, financial pkgs, xts/zoo
suppressPackageStartupMessages(library(dataRetrieval))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(twitteR))

setwd ("~/PipelineMonitoring_Server")
gageInfo <- read_csv('data/gageInfoEVJ.csv')
twitterAuth <- readRDS('data/twitterAuth.RDS')
#source('rmd2r.R')
#rmd2rscript('parameterFunctions.Rmd')
#source('parameterFunctions.R')
source('parameterFunctions[rmd2r].R')

# Set up oauth for server manually
setup_twitter_oauth(consumer_key= twitterAuth$consumer_key,
                    consumer_secret = twitterAuth$consumer_secret,
                    access_token = twitterAuth$access_token,
                    access_secret = twitterAuth$access_secret)


# Inputs to functions and loop
last2Hours <- format(Sys.time()-7200,format="%Y-%m-%d %H:00:00") # Need to filter by the last 2 hours bc exceedance could start mid-hour and be lost if constantly resetting on hour
start_time <- Sys.time()
a <- 1:(nrow(gageInfo)-2)
n <- length(a)
allGageResults <- list() # make an empty list to store all gage data dataframes

for(i in a[seq(1,n,2)]){ # make a sequence 1:26 with only odd numbers bc paired gages are listed below upstream gage
  print(i)               # only 26 bc the last two records of gageInfo don't have actual gages
  upstreamDataForTurbidity <- NWISpull(paste(0,gageInfo$`USGS Station ID`[i],sep=''), Sys.Date()-1,Sys.Date())#format(Sys.time()-21600,format="%Y-%m-%d %H:%M:%S"), format(Sys.time(),format="%Y-%m-%d %H:%M:%S"))
  upstreamData <- dplyr::filter(upstreamDataForTurbidity,dateTime >= last2Hours)
  # Special step to get gage height from Lafayette gage for Roanoke gage pair only
  if(gageInfo$`USGS Station ID`[i]=='205450393'){
    upstreamDataLafayetteForTurbidity <- NWISpull("02054500", Sys.Date()-1,Sys.Date())%>%#format(Sys.time()-21600,format="%Y-%m-%d %H:%M:%S"), format(Sys.time(),format="%Y-%m-%d %H:%M:%S"))
      select(dateTime,GH_Inst,GH_Inst_cd)
    upstreamDataLafayette  <- dplyr::filter(upstreamDataLafayetteForTurbidity,dateTime >= last2Hours) 
    upstreamData <- full_join(upstreamData,upstreamDataLafayette,by='dateTime') %>% 
      dplyr::select(agency_cd,site_no,dateTime,Wtemp_Inst,Wtemp_Inst_cd,GH_Inst,GH_Inst_cd,everything())
    upstreamDataForTurbidity <- full_join(upstreamDataForTurbidity,upstreamDataLafayetteForTurbidity,by='dateTime') %>% 
      dplyr::select(agency_cd,site_no,dateTime,Wtemp_Inst,Wtemp_Inst_cd,GH_Inst,GH_Inst_cd,everything())
    }
  downstreamDataForTurbidity <- NWISpull(paste(0,gageInfo$`USGS Station ID`[i+1],sep=''),Sys.Date()-1,Sys.Date())
  downstreamData <- dplyr::filter(downstreamDataForTurbidity,dateTime >= last2Hours)
  
  gageResults <- dataScan(upstreamData, downstreamData, 
                          WQclassGage1 = gageInfo$WQS_Class[i], 
                          WQclassGage2 = gageInfo$WQS_Class[i+1], 
                          pHspecialStandards = gageInfo$pH_SpecialStandards[i],
                          pHrangeAllowance = gageInfo$pH_RangeAllowance[i], 
                          SpCond_Designation = gageInfo$SpCond_Designation[i], 
                          turbidityBaseline = gageInfo$Turbidity_Baseline[i], 
                          turbidity99th1 = gageInfo$Turbidity_99th[i],
                          turbidity99th2 = gageInfo$Turbidity_99th[i+1])
  
  # Send Notification if no data could be compared upstream/downstream for pull
  noGageDataTweet(gageResults,i,last2Hours)
  
  # Send Notification if any 5 minute flags are blown
  anyExceedance(gageResults,upstreamData,downstreamData)
  
  # Send Notification if any Temperature hourly change flags are blown
  tempTimeTweet(gageResults,upstreamData,downstreamData) 
  
  # Send notification if any turbidity change flags are blown
  turbTimeTweet(gageResults,upstreamData,downstreamData,turbidity99th1=gageInfo$Turbidity_99th[i],
                turbidity99th2=gageInfo$Turbidity_99th[i+1],upstreamDataForTurbidity,
                downstreamDataForTurbidity)
  # Save data, if in interactive testing mode
  #allGageResults[[i]] <- gageResults 
}
# Running all sites takes how long?
Sys.time()-start_time

#test <- allGageResults[[1]]
#write.csv(test,'data/test.csv')