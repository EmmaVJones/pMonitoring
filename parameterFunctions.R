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

NArecords <- function(x){
  sum(is.na(x))
}

turbidityNumericThreshold <- function(x,threshold){
  exceed <- numericThreshold(x,threshold)
  NAs <- NArecords(x)
  z <- c(turbidity_Exceedance=exceed,turbidity_NAs=NAs)
  
  return(z)
}

turbidityPercentThreshold <- function(x,threshold){
  exceed <- numericThreshold(x,threshold)
  NAs <- NArecords(x)
  z <- c(turbidity_Exceedance=exceed,turbidity_NAs=NAs)
  
  return(z)
}


temperatureHourlyChange <- function(x,threshold){
  Tchange <- max(x, na.rm=T)-min(x, na.rm=T)
  NAs <- NArecords(x)
  violation <- ifelse(Tchange>threshold,1,0)
  z <- c(T_maxHourlyChange=Tchange, T_hourlyChangeViolation=violation,T_NAs=NAs)
  return(z)
}



# Temperature Analysis
temperatureByWQSclass <- function(x,maxT,natTchange,maxHourlyTchange){
  if(nrow(x)<13){
    together <- mutate(x,T_upstreamMaxViolation=ifelse(upstream>maxT,1,0),# 9VAC25-260-50. Numerical Criteria for Maximum Temperature
                       T_downstreamMaxViolation=ifelse(downstream>maxT,1,0), # 9VAC25-260-50. Numerical Criteria for Maximum Temperature
                       T_riseAboveNaturalViolation=ifelse(numericDiff>natTchange,1,0), # 9VAC25-260-60. Rise Above Natural Temperature
                       # since there is less than an hour's worth of data to analyze in this scenario, report back NA's
                       T_maxHourlyChange= as.numeric(NA),
                       T_upstreamHourlyChangeViolation= as.numeric(NA),
                       T_upstreamNAs= as.numeric(NA),
                       T_maxHourlyChange.1= as.numeric(NA),
                       T_downstreamHourlyChangeViolation= as.numeric(NA),
                       T_downstreamNAs= as.numeric(NA))
     }else{
      together <- mutate(x,T_upstreamMaxViolation=ifelse(upstream>maxT,1,0),# 9VAC25-260-50. Numerical Criteria for Maximum Temperature
                         T_downstreamMaxViolation=ifelse(downstream>maxT,1,0), # 9VAC25-260-50. Numerical Criteria for Maximum Temperature
                         T_riseAboveNaturalViolation=ifelse(numericDiff>natTchange,1,0)) %>% # 9VAC25-260-60. Rise Above Natural Temperature
        tq_mutate(  # 9VAC25-260-70. Maximum Hourly Temperature Change
          select     = upstream,
          mutate_fun = rollapply, 
          # rollapply args
          width      = 13,
          align      = "right",
          by.column  = FALSE,
          FUN        = temperatureHourlyChange,
          # FUN args
          threshold  = maxHourlyTchange) %>%
        tq_mutate( # 9VAC25-260-70. Maximum Hourly Temperature Change
          select     = downstream,
          mutate_fun = rollapply, 
          # rollapply args
          width      = 13,
          align      = "right",
          by.column  = FALSE,
          FUN        = temperatureHourlyChange,
          # FUN args
          threshold  = maxHourlyTchange) %>%
        rename(T_upstreamHourlyChangeViolation= T_hourlyChangeViolation, T_upstreamNAs= T_NAs,
               T_downstreamHourlyChangeViolation= T_hourlyChangeViolation.1, T_downstreamNAs=T_NAs.1)
     }
      
      return(together)
      
}

temperature <- function(upstreamData, downstreamData, parameter, WQclass){
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    filter(upstream >= 0 & upstream < 40 | is.na(upstream)  & downstream >= 0 & downstream < 40  | is.na(downstream)) %>% # Ensure data is realistic, filter out wonky probe readings
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100,
           window= lag(dateTime,12),
           T_valid1hrWindow=as.double(dateTime-lag(dateTime,12),units='mins')) # lag 11 bc lag already grabs 1 row above
  
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




# Dissolved Oxygen Analysis
dissolvedOxygenByWQSclass <- function(x,minDO){
  together <- mutate(x,DO_upstreamMinViolation=ifelse(upstream<minDO,1,0), # 9VAC25-260-50. Numerical Criteria for Dissolved Oxygen
                     DO_downstreamMinViolation=ifelse(downstream<minDO,1,0),  # 9VAC25-260-50. Numerical Criteria for Dissolved Oxygen
                     DO_upDownDifferenceViolation = ifelse(abs(numericDiff)>1,1,0)) # Flag if upstream/downstream difference > 1 mg/L
  return(together)
}

dissolvedOxygen <- function(upstreamData, downstreamData, parameter, WQclass){
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime')) %>%
    filter(upstream >= 0 & upstream < 25  | is.na(upstream) & downstream >= 0 & downstream < 25 | is.na(downstream)) %>% # Ensure data is realistic, filter out wonky probe readings
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100)
  
  # Max Temperature Violations
  if(WQclass==6){dat <- dissolvedOxygenByWQSclass(together,minDO = 8)%>%mutate(WQSapplied=WQclass)}
  if(WQclass==5){dat <- dissolvedOxygenByWQSclass(together,minDO = 7)%>%mutate(WQSapplied=WQclass)}
  if(WQclass==4){dat <- dissolvedOxygenByWQSclass(together,minDO = 7)%>%mutate(WQSapplied=WQclass)}
  if(WQclass==3){dat <- dissolvedOxygenByWQSclass(together,minDO = 7)%>%mutate(WQSapplied=WQclass)}
  
  return(dat)
}




# pH Analysis
pHbyWQS <- function(x,pHspecialStandards,pHrangeAllowance){
  # Adjust standards if WQS special standard for site
  min_pH <- ifelse(pHspecialStandards=="Y",6.5,6)
  max_pH <- ifelse(pHspecialStandards=="Y",9.5,9)
  
  together <- mutate(x,pH_upstreamViolation=ifelse(upstream > min_pH & upstream < max_pH,0,1), # 9VAC25-260-50. Numerical Criteria for pH
                     pH_downstreamViolation=ifelse(downstream > min_pH & upstream < max_pH,0,1),  # 9VAC25-260-50. Numerical Criteria for pH
                     pH_upDownDifferenceViolation = ifelse(abs(numericDiff)>pHrangeAllowance,1,0)) # Flag if upstream/downstream difference > allowed range
  return(together)
}

pH <- function(upstreamData, downstreamData, parameter, pHspecialStandards, pHrangeAllowance){
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime')) %>%
    filter(upstream >= 0 & upstream <= 14 | is.na(upstream) & downstream >= 0 & downstream <= 14 | is.na(downstream)) %>% # Ensure data is realistic, filter out wonky probe readings
    mutate(parameter=parameter,
           numericDiff=downstream-upstream, 
           pctDiff=(numericDiff/downstream)*100)
  
  pHbyWQS(together, pHspecialStandards, pHrangeAllowance)
}




#Specific Conductivity Analysis
SpCondByDesignation <- function(x,SpCond_Designation){
  upDownDifferenceThreshold <- ifelse(SpCond_Designation==500,100,50)
  
  mutate(x,spCond_upstreamViolation=ifelse(upstream > SpCond_Designation,1,0), 
         spCond_downstreamViolation=ifelse(downstream > SpCond_Designation,1,0),
         spCond_upDownNumericDifferenceViolation = ifelse(numericDiff>upDownDifferenceThreshold,1,0)) %>% # Flag if upstream/downstream numeric difference > allowed range
    rowwise()%>% # Make sure the upDownDifferenceThreshold p
    mutate(spCond_upDownPercentDifferenceViolation = ifelse(upDownDifferenceThreshold==100,
                                                            ifelse(pctDiff>20,1,0), "Not Applicable")) # Flag if upstream/downstream percent difference > allowed range, only for high conductivity sites
}


SpCond <- function(upstreamData, downstreamData, parameter, SpCond_Designation){
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime')) %>%
    filter(upstream >= 0 & upstream < 2000 | is.na(upstream) & downstream >= 0 & downstream < 2000 | is.na(downstream)) %>% # Ensure data is realistic, filter out wonky probe readings
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100)
  
  SpCondByDesignation(together, SpCond_Designation)
}




# Turbidity analysis
turbidity <- function(upstreamData, downstreamData, parameter, turbidityBaseline, turbidity95th){
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    filter(upstream >= 0 & upstream < 4000 | is.na(upstream) & downstream >= 0 & downstream < 4000 | is.na(downstream) ) %>% # Ensure data is realistic, filter out wonky probe readings, as of 12/19/2017 the max(turbidity at all sites)= 1170, so 4000 is very inclusive
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100,
           window= lag(dateTime,6),
           turbidity_valid30minuteWindow=dateTime-lag(dateTime,6)) # lag 6 bc lag already grabs 1 row above
  
  # if < 7 records the rolling functions will fail, report out NA's if not enough data to run full metrics
  if(nrow(together)<7){
    tidyverse_diff_rollstats <- mutate(together,turbidity_Exceedance=as.numeric(NA),turbidity_NAs=as.numeric(NA),
                                       turbidity_ExceedanceType=as.numeric(NA),turbidity_upstreamExceed95th=as.numeric(NA),
                                       turbidity_downstreamExceed95th=as.numeric(NA))
                       
  }else{
    if(turbidityBaseline <= 40){
      tidyverse_diff_rollstats <- together %>%
        tq_mutate(
          select     = numericDiff,
          mutate_fun = rollapply, 
          # rollapply args
          width      = 7,
          align      = "right",
          by.column  = FALSE,
          FUN        = turbidityNumericThreshold,
          # FUN args
          threshold  = 5.999999999) %>% mutate(turbidity_ExceedanceType='Numeric')
    }else{
      tidyverse_diff_rollstats <- together %>%
        tq_mutate(
          select     = pctDiff,
          mutate_fun = rollapply, 
          # rollapply args
          width      = 7,
          align      = "right",
          by.column  = FALSE,
          FUN        = turbidityPercentThreshold, 
          # FUN args
          threshold  = 14.99999999) %>% mutate(turbidity_ExceedanceType='Percent')
    }
    
    # Exceed 95th percentile upstream or downstream
    tidyverse_diff_rollstats <- tidyverse_diff_rollstats %>%
      tq_mutate(
        select     = upstream,
        mutate_fun = rollapply, 
        # rollapply args
        width      = 7,
        align      = "right",
        by.column  = FALSE,
        FUN        = turbidityNumericThreshold,
        # FUN args
        threshold  = turbidity95th) %>%
      rename(turbidity_upstreamExceed95th=turbidity_Exceedance.1,turbidity_upstreamNAs=turbidity_NAs.1) %>%
      tq_mutate(
        select     = downstream,
        mutate_fun = rollapply, 
        # rollapply args
        width      = 7,
        align      = "right",
        by.column  = FALSE,
        FUN        = turbidityNumericThreshold,
        # FUN args
        threshold  = turbidity95th) %>%
      rename(turbidity_downstreamExceed95th=turbidity_Exceedance.1,turbidity_downstreamNAs=turbidity_NAs.1) 
    # Only count violations if they in fact exceed 30 minutes, go back and double check flags correspond with 6 continuous rows of data (no time gaps)
    # mutate(window= lag(dateTime,6),
    #    valid30minuteWindow=dateTime-lag(dateTime,6)) # lag 6 bc lag already grabs 1 row above
  }
  
  return(tidyverse_diff_rollstats)
}


## Apply threshold exceedance functions across time series, multitple parameters
dataScan <- function(upstreamData, downstreamData, WQclass, pHspecialStandards, pHrangeAllowance, SpCond_Designation, turbidityBaseline, turbidity95th){
  
  ## Temperature Analysis ##
  T_Results <- temperature(upstreamData,downstreamData,'Wtemp_Inst',WQclass) %>% # will have to preprogram WQS class depending on site
    select(agency_cd,dateTime,site_no.x,site_no.y,WQSapplied,T_valid1hrWindow,T_upstreamMaxViolation,T_downstreamMaxViolation,
           T_upstreamHourlyChangeViolation,T_upstreamNAs,T_downstreamHourlyChangeViolation,T_downstreamNAs,T_riseAboveNaturalViolation)
  
  ## Dissolved Oxygen Analysis ##
  DO_Results <- dissolvedOxygen(upstreamData,downstreamData,'DO_Inst',WQclass) %>% # will have to preprogram WQS class depending on site
    select(agency_cd,dateTime,site_no.x,site_no.y,DO_upstreamMinViolation,DO_downstreamMinViolation,DO_upDownDifferenceViolation)
  
  ## pH Analysis ##
  pH_Results <- pH(upstreamData,downstreamData,'pH_Inst',pHspecialStandards,pHrangeAllowance) %>% # will have to preprogram WQS special standards depending on site
    select(agency_cd,dateTime,site_no.x,site_no.y,pH_upstreamViolation,pH_downstreamViolation,pH_upDownDifferenceViolation)
  
  ## Specific Conductivity Analysis ##
  spCond_Results <- SpCond(upstreamData,downstreamData,'SpecCond_Inst', SpCond_Designation) %>% 
    select(agency_cd,dateTime,site_no.x,site_no.y,spCond_upstreamViolation,spCond_downstreamViolation,spCond_upDownNumericDifferenceViolation,
           spCond_upDownPercentDifferenceViolation)
  
  ## Turbidity Analysis ##
  # Choose turbidity method based on baseline parameter value
  turbidity_Results <- turbidity(upstreamData, downstreamData, 'Turb_Inst', turbidityBaseline, turbidity95th) %>% 
    select(agency_cd,dateTime,site_no.x,site_no.y,turbidity_valid30minuteWindow,turbidity_NAs,turbidity_ExceedanceType,turbidity_Exceedance,
           turbidity_upstreamExceed95th,turbidity_upstreamNAs,turbidity_downstreamExceed95th,turbidity_downstreamNAs)
  
  allParameters <- plyr::join_all(list(T_Results,DO_Results, pH_Results, spCond_Results, turbidity_Results),by=c('agency_cd','dateTime','site_no.x','site_no.y'))
  
  return(allParameters)
}

#------------------------## Notification Functions ##---------------------------------------#
# Based on column name, output upstream or downstream gage number
updownDifference <-function(colName,gageResults){
  if(length(grep('upstream',colName))){return(gageResults$site_no.x[1])}
  if(length(grep('downstream',colName))){return(gageResults$site_no.y[1])}}  

# Based on column name, output parameter code (for weblink)
parameterCodeDecipher <- function(colName){
  if(length(grep('T',colName,ignore.case = F))){return('00010')}
  if(length(grep('DO',colName,ignore.case = F))){return('00300')}
  if(length(grep('pH',colName,ignore.case = F))){return('00400')}
  if(length(grep('spCond',colName,ignore.case = F))){return('00095')}
  if(length(grep('turbidity',colName,ignore.case = F))){return('63680')}}

# Based on gage number, output stream name
gageNumberToStream <- function(gageNumber){
  gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]}

# Send tweet if any metrics that require 5 min exceedance for notification
anyExceedance <- function(gageResults,upstreamData,downstreamData){
  minmaxData <- select(gageResults,dateTime,site_no.x,site_no.y,T_upstreamMaxViolation,
                       T_downstreamMaxViolation,T_riseAboveNaturalViolation,DO_upstreamMinViolation,
                       DO_downstreamMinViolation,DO_upDownDifferenceViolation,
                       pH_upstreamViolation,pH_downstreamViolation,pH_upDownDifferenceViolation,
                       spCond_upstreamViolation,spCond_downstreamViolation,
                       spCond_upDownNumericDifferenceViolation,spCond_upDownPercentDifferenceViolation)
  
  for(i in 4:length(minmaxData)){# start at 4 bc first columns are gage numbers and date
    columnName <- colnames(minmaxData)[i]
    PCcolName <- gsub("Violation","Notification",columnName) # PC name for tweet
    date1 <- as.Date(minmaxData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
    date2 <- as.Date(minmaxData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
    parameterCode <- parameterCodeDecipher(columnName) # change parameter from column name to USGS parameter code for link
    uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
    if(length(grep('T',columnName,ignore.case = F))){selection <- c("Wtemp_Inst","Wtemp_Inst_cd")}
    if(length(grep('DO',columnName,ignore.case = F))){selection <- c("DO_Inst","DO_Inst_cd")}
    if(length(grep('pH',columnName,ignore.case = F))){selection <- c("pH_Inst","pH_Inst_cd")}
    if(length(grep('spCond',columnName,ignore.case = F))){selection <- c("SpecCond_Inst","SpecCond_Inst_cd")}
    if(length(grep('turbidity',columnName,ignore.case = F))){selection <- c("Turb_Inst","Turb_Inst_cd")}
    
    
    if(length(grep("riseAboveNatural",columnName)) > 0 | length(grep('upDown',columnName)) > 0){
      gageNumber1 <- minmaxData$site_no.x[1]
      gageNumber2 <- minmaxData$site_no.y[1]
      streamName <- gageInfo[grep(as.numeric(gageNumber1),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_',parameterCode,'=on&site_no=',gageNumber1,'%2C',gageNumber2,'&format=gif_mult_sites',sep='')
      datToSave <- minmaxData[,c(1:3,i)]
      upstreamRAW <- select(upstreamData,site_no,dateTime,selection)
      names(upstreamRAW)[c(1,3:4)] <- paste(names(upstreamRAW)[c(1,3:4)],".x",sep="")
      downstreamRAW <- select(downstreamData,site_no,dateTime,selection)
      names(downstreamRAW)[c(1,3:4)] <- paste(names(downstreamRAW)[c(1,3:4)],".y",sep="")
      datToSave <- plyr::join_all(list(datToSave,upstreamRAW,downstreamRAW),by='dateTime')
     }else{
      gageNumber <- updownDifference(columnName,minmaxData)
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_',parameterCode,'=on&site_no=',gageNumber,sep = "")
      #Pull together data to save for archive
      if(length(grep('upstream',columnName)) > 0){
        otherRAWdata  <- select(upstreamData,dateTime,site_no,selection)
        datToSave <- minmaxData[,c(1:2,i)] %>%
          full_join(otherRAWdata,by='dateTime')}
      if(length(grep('downstream',columnName)) > 0){
        otherRAWdata  <- select(downstreamData,dateTime,site_no,selection)
        datToSave <- minmaxData[,c(1,3,i)] %>%
          full_join(otherRAWdata,by='dateTime')}
      }
    if(1 %in% minmaxData[,i]){ # If there is a single 5 min violation for a metric, send tweet and save analysis
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)
    }}
}

# Tweet if any hourly change violations (based on stream WQS)
tempTimeTweet <- function(gageResults,upstreamData,downstreamData){
  tempTimeData <- select(gageResults,dateTime,site_no.x,site_no.y,WQSapplied,T_valid1hrWindow,
                         T_upstreamHourlyChangeViolation,T_upstreamNAs,
                         T_downstreamHourlyChangeViolation,T_downstreamNAs)
  
  # Only use data from correct 1hr window
  validData <- filter(tempTimeData,T_valid1hrWindow==60)
  if(nrow(validData)>0){
    # Upstream Analysis
    validUP <- filter(validData,T_upstreamHourlyChangeViolation>0, # Limit dataset to only violations
                      T_upstreamNAs <=4) # Only allow up to 4 missing Temperature reading per hour for any hour to count toward max change rate violation
    # Only send tweet if at least 1 valid row (1 hr violation) & <= 4 missing Temperature readings within hour
    if(nrow(validUP)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'T_upstreamHourlyChangeNotification'
      gageNumber <- unique(validData$site_no.x)[1]
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_00010=on&site_no=',gageNumber,sep = "")
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- tempTimeData[,c(1,2,4:7)]
      upstreamRAW <- select(upstreamData,site_no,dateTime,Wtemp_Inst,Wtemp_Inst_cd) %>%
        dplyr::rename(site_no.x=site_no,Wtemp_Inst.x=Wtemp_Inst,Wtemp_Inst_cd.x=Wtemp_Inst_cd)
      datToSave <- plyr::join_all(list(datToSave,upstreamRAW),by='dateTime')
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)
    }
    # Downstream Analysis
    validDOWN <- filter(validData,T_downstreamHourlyChangeViolation>0, # Limit dataset to only violations
                        T_downstreamNAs <=4) # Only allow up to 4 missing Temperature reading per hour for any hour to count toward max change rate violation
    # Only send tweet if at least 1 valid row (1 hr violation) & <= 4 missing Temperature readings within hour
    if(nrow(validDOWN)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'T_downstreamHourlyChangeNotification'
      gageNumber <- unique(validData$site_no.y)[1]
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_00010=on&site_no=',gageNumber,sep = "")
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- tempTimeData[,c(1,3:5,8,9)]
      downstreamRAW <- select(downstreamData,site_no,dateTime,Wtemp_Inst,Wtemp_Inst_cd) %>%
        dplyr::rename(site_no.y=site_no,Wtemp_Inst.y=Wtemp_Inst,Wtemp_Inst_cd.y=Wtemp_Inst_cd)
      datToSave <- plyr::join_all(list(datToSave,downstreamRAW),by='dateTime')
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)}
  }
}

# Send Notification if no data could be compared upstream/downstream for pull
noGageDataTweet <- function(gageResults){
  if(nrow(filter(gageResults,!is.na(site_no.x))) == 0 | nrow(filter(gageResults,!is.na(site_no.y))) == 0){
    uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
    # find Culprit
    if(nrow(filter(gageResults,!is.na(site_no.x))) == 0){columnName <- 'upstreamNoData'}
    if(nrow(filter(gageResults,!is.na(site_no.y))) == 0){columnName <- 'downstreamNoData'}
    if(nrow(filter(gageResults,!is.na(site_no.x))) == 0 &
       nrow(filter(gageResults,!is.na(site_no.y))) == 0){columnName <- 'NoData'}
    
    gageNumberExists <- as.character(
      filter(gageInfo,`USGS Station ID` %in% 
               c(substring(unique(gageResults$site_no.x),2),substring(unique(gageResults$site_no.y),2) ))%>%
        filter(!is.na(`USGS Station ID`)) %>% select(`USGS Station ID`))
    streamName <- gageInfo[grep(as.numeric(gageNumberExists),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
    gageNumberMissing <- paste(0,as.character(
      filter(gageInfo,`Stream Name` %in% streamName) %>%
        filter(`USGS Station ID` != gageNumberExists) %>% select(`USGS Station ID`)),sep='')
    date1 <- as.Date(gageResults$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
    date2 <- as.Date(gageResults$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
    parameterCode <- '&cb_00010=on&cb_00095=on&cb_00300=on&cb_00400=on&cb_63680=on'
    weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                     parameterCode,'&site_no=',gageNumberMissing,sep = "")
    tweet(paste(streamName,columnName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
    # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
    write.csv(gageResults,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)}
}

# Turbidity tweet
turbTimeTweet <- function(gageResults,upstreamData,downstreamData){
  turbTimeData <- select(gageResults,dateTime,site_no.x,site_no.y,WQSapplied,turbidity_valid30minuteWindow,
                         turbidity_NAs,turbidity_ExceedanceType,turbidity_Exceedance,
                         turbidity_upstreamExceed95th,turbidity_upstreamNAs,
                         turbidity_downstreamExceed95th,turbidity_downstreamNAs)
  
  # Only use data from correct 1hr window
  validData <- filter(turbTimeData,turbidity_valid30minuteWindow==30)
  
  # upstream downstream change comparison
  if(nrow(validData)>0){
    validComparison <- filter(validData,turbidity_Exceedance>0, # Limit dataset to only violations
                              turbidity_NAs <= 3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    if(nrow(validComparison)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'turbidity_upstreamDownstreamNotification'
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      gageNumber1 <- unique(validData$site_no.x)[1]
      gageNumber2 <- unique(validData$site_no.y)[1]
      streamName <- gageInfo[grep(as.numeric(gageNumber1),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_63680=on&site_no=',gageNumber1,'%2C',gageNumber2,'&format=gif_mult_sites',sep='')
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- turbTimeData[,c(1:3,5:8)]
      upstreamRAW <- select(upstreamData,site_no,dateTime,Turb_Inst,Turb_Inst_cd) %>%
        dplyr::rename(site_no.x=site_no,Turb_Inst.x=Turb_Inst,Turb_Inst_cd.x=Turb_Inst_cd)
      downstreamRAW <- select(downstreamData,site_no,dateTime,Turb_Inst,Turb_Inst_cd) %>%
        dplyr::rename(site_no.y=site_no,Turb_Inst.y=Turb_Inst,Turb_Inst_cd.y=Turb_Inst_cd)
      datToSave <- plyr::join_all(list(datToSave,upstreamRAW,downstreamRAW),by='dateTime')
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)
    }
    
    # Upstream absolute threshold Analysis
    validUP <- filter(validData,turbidity_upstreamExceed95th>0, # Limit dataset to only violations
                      turbidity_upstreamNAs <=3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    # Only send tweet if at least 1 valid row (1 hr violation) & <= 3 missing turbidity readings within hour
    if(nrow(validUP)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'turbidity_upstreamExceed95th'
      gageNumber <- unique(validData$site_no.x)[1]
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_63680=on&site_no=',gageNumber,sep = "")
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- turbTimeData[,c(1,2,5,9,10)]
      upstreamRAW <- select(upstreamData,site_no,dateTime,Turb_Inst,Turb_Inst_cd) %>%
        dplyr::rename(site_no.x=site_no,Turb_Inst.x=Turb_Inst,Turb_Inst_cd.x=Turb_Inst_cd)
      datToSave <- plyr::join_all(list(datToSave,upstreamRAW),by='dateTime')
      
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)
    }
    # Downstream Analysis
    validDOWN <- filter(validData,turbidity_downstreamExceed95th>0, # Limit dataset to only violations
                        turbidity_downstreamNAs <=3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    # Only send tweet if at least 1 valid row (1 hr violation) & <= 4 missing turbidity readings within hour
    if(nrow(validDOWN)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'turbidity_downstreamExceed95th'
      gageNumber <- unique(validData$site_no.x)[1]
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_63680=on&site_no=',gageNumber,sep = "")
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- turbTimeData[,c(1,3,5,11,12)]
      downstreamRAW <- select(downstreamData,site_no,dateTime,Turb_Inst,Turb_Inst_cd) %>%
        dplyr::rename(site_no.y=site_no,Turb_Inst.y=Turb_Inst,Turb_Inst_cd.y=Turb_Inst_cd)
      datToSave <- plyr::join_all(list(datToSave,downstreamRAW),by='dateTime')
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)}
  }
}



#-----------------------------------------------------------------------------------------------
