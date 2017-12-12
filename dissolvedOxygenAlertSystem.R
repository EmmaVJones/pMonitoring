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

# Find Threshold Exceedances,  ################################################ programmed as > not >= !!!!!!!!!
numericThreshold <- function(x, threshold){
  suppressWarnings(
    if(min(x,na.rm=T)==Inf){NA
    }else{ifelse(min(x,na.rm=T)>threshold,1,0) })}
percentThreshold <- function(x,threshold){  
  suppressWarnings(
    if(min(x,na.rm=T)==Inf){NA
    }else{ifelse(min(x,na.rm=T)>threshold,1,0)  })}



# Dissolved Oxygen Analysis
dissolvedOxygenByWQSclass <- function(x,minDO){
  together <- mutate(x,upstreamMinDOviolation=ifelse(upstream<minDO,1,0),
                     downstreamMinDOviolation=ifelse(downstream<minDO,1,0))
  
  return(together)
  
}


# Temperature Analysis
temperatureByWQSclass <- function(x,maxT,natTchange,maxHourlyTchange){
  together <- mutate(x,upstreamMaxTviolation=ifelse(upstream>maxT,1,0),# 9VAC25-260-50. Numerical Criteria for Maximum Temperature
                     downstreamMaxTviolation=ifelse(downstream>maxT,1,0), # 9VAC25-260-50. Numerical Criteria for Maximum Temperature
                     riseAboveNaturalTviolation=ifelse(numericDiff>natTchange,1,0))# 9VAC25-260-60. Rise Above Natural Temperature
  # 9VAC25-260-70. Maximum Hourly Temperature Change
  tidyverse_diff_rollstats <- together %>%
    tq_mutate(
      select     = numericDiff,
      mutate_fun = rollapply, 
      # rollapply args
      width      = 12,
      align      = "right",
      by.column  = FALSE,
      FUN        = numericThreshold,
      # FUN args
      threshold  = maxHourlyTchange,
      # tq_mutate args
      col_rename = "hourlyTchangeViolation")
  return(tidyverse_diff_rollstats)
}


temperature <- function(upstreamData, downstreamData, parameter, WQclass){
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,Turb_Inst)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100)
  
  # Max Temperature Violations
  if(WQclass==6){dat <- temperatureByWQSclass(together,maxT = 20, natTchange = 1, maxHourlyTchange = 0.5)%>%
    mutate(WQSapplied=WQclass)}
  if(WQclass==5){dat <- temperatureByWQSclass(together,maxT = 21, natTchange = 3, maxHourlyTchange = 2)%>%
    mutate(WQSapplied=WQclass)}
  if(WQclass==4){dat <- temperatureByWQSclass(together,maxT = 31, natTchange = 3, maxHourlyTchange = 2)%>%
    mutate(WQSapplied=WQclass)}
  if(WQclass==3){dat <- temperatureByWQSclass(together,maxT = 32, natTchange = 3, maxHourlyTchange = 2)%>%
    mutate(WQSapplied=WQclass)}
  
  return(dat)
}
#test <- temperature(upstreamData,downstreamData,'Wtemp_Inst',6)











## Apply threshold exceedance functions across time series, multitple parameters
dataScan <- function(upstreamData, downstreamData, WQclass, turbidityBaseline, turbidityUnitIncrease, turbidityPercentIncrease){
  
  ## Temperature Analysis ##
  tempResults <- temperature(upstreamData,downstreamData,'Wtemp_Inst',WQclass) # will have to preprogram WQS class depending on site
  
  ## Turbidity Analysis ##
  # Choose turbidity method based on baseline parameter value
  turbidityResults <- turbidity(upstreamData, downstreamData, 'Turb_Inst', turbidityBaseline, turbidityUnitIncrease, turbidityPercentIncrease)
  
  return(tempResults)#turbidityResults)
  # Eventually will need to build way for notifications to be combined across parameters to send to
  # email function to sift through and trigger email if violations occur
  
}

#test <- dataScan(upstreamData, downstreamData, WQclass=6, turbidityBaseline=40, turbidityUnitIncrease=6, turbidityPercentIncrease=15)








## Test Analyses ##
last2Hours <- format(Sys.time()-7200,format="%Y-%m-%d %H:00:00") # Need to filter by the last 2 hours bc exceedance could start mid-hour and be lost if constantly resetting on hour

upstreamData <- NWISpull('03171597', Sys.Date()-1,Sys.Date())%>%#format(Sys.time()-21600,format="%Y-%m-%d %H:%M:%S"), format(Sys.time(),format="%Y-%m-%d %H:%M:%S"))
  dplyr::filter(dateTime >= last2Hours)
downstreamData <- NWISpull('0317159760',Sys.Date()-1,Sys.Date())%>%
  dplyr::filter(dateTime >= last2Hours)

#test <- dataScan(upstreamData, downstreamData, WQclass=6, turbidityBaseline=40, turbidityUnitIncrease=6, turbidityPercentIncrease=15)
test <- dataScan(upstreamData, downstreamData, WQclass=6, turbidityBaseline=40, turbidityUnitIncrease=6, turbidityPercentIncrease=15)


