# Pull NWIS data
NWISpull <- function(gageNo,start,end){
  allUnitData <- readNWISuv(siteNumbers=gageNo,
                            parameterCd=c("00010", "00095", "00300", "00400", "63680","00065"),
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
temperatureByWQSclass <- function(x, maxHourlyTchange1, maxHourlyTchange2){
  # 9VAC25-260-70. Maximum Hourly Temperature Change, UPSTREAM
  if(nrow(filter(x,!is.na(upstream))) < 13){
    x <- mutate(x,T_maxHourlyChange=NA,T_hourlyChangeViolation=NA,T_NAs=NA)
  }else{x <- tq_mutate(x,
                       select     = upstream,
                       mutate_fun = rollapply, 
                       # rollapply args
                       width      = 13,
                       align      = "right",
                       by.column  = FALSE,
                       FUN        = temperatureHourlyChange,
                       # FUN args
                       threshold  = maxHourlyTchange1)}
  # 9VAC25-260-70. Maximum Hourly Temperature Change, DOWNSTREAM  
  if(nrow(filter(x,!is.na(downstream))) < 13){
    x <- mutate(x,T_maxHourlyChange.1=NA,T_hourlyChangeViolation.1=NA,T_NAs.1=NA)
  }else{x <- tq_mutate(x,
                       select     = downstream,
                       mutate_fun = rollapply, 
                       # rollapply args
                       width      = 13,
                       align      = "right",
                       by.column  = FALSE,
                       FUN        = temperatureHourlyChange,
                       # FUN args
                       threshold  = maxHourlyTchange2)}
  x <- rename(x, T_upstreamHourlyChangeViolation= T_hourlyChangeViolation, T_upstreamNAs= T_NAs,
              T_downstreamHourlyChangeViolation= T_hourlyChangeViolation.1, T_downstreamNAs=T_NAs.1)
  return(x)
}


temperature <- function(upstreamData, downstreamData, parameter, WQclassGage1, WQclassGage2){
  # Add Wtemp fields if no Wtemp data retrieved from probes, populate with NA 
  if(unique(c('Wtemp_Inst',"Wtemp_Inst_cd") %in% names(upstreamData))){
    UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  if(unique(c('Wtemp_Inst',"Wtemp_Inst_cd") %in% names(downstreamData))){
    DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  # Establish thresholds based on WQS for max Temperature Violations
  if(WQclassGage1==6){maxT1 <- 20; natTchange <- 1; maxHourlyTchange1 <- 0.5}
  if(WQclassGage2==6){maxT2 <- 20; natTchange <- 1; maxHourlyTchange2 <- 0.5}
  if(WQclassGage1==5){maxT1 <- 21; natTchange <- 3; maxHourlyTchange1 <- 2}
  if(WQclassGage2==5){maxT2 <- 21; natTchange <- 3; maxHourlyTchange2 <- 2}
  if(WQclassGage1==4){maxT1 <- 31; natTchange <- 3; maxHourlyTchange1 <- 2}
  if(WQclassGage2==4){maxT2 <- 31; natTchange <- 3; maxHourlyTchange2 <- 2}
  if(WQclassGage1==3){maxT1 <- 32; natTchange <- 3; maxHourlyTchange1 <- 2}
  if(WQclassGage2==3){maxT2 <- 32; natTchange <- 3; maxHourlyTchange2 <- 2}
  # natTchange called same parameter name bc want the larger # to override smaller 
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    # Ensure data is realistic, replace wonky probe readings with NA so they are skipped
    # Filtering out bad data is not a good idea because you will lose entire rows of data and subsequently
    #  throw off other analyses down the line
    mutate(upstream = replace(upstream, upstream < 0 | upstream > 40, NA ),
                downstream = replace(downstream, downstream < 0 | downstream > 40, NA )) %>%
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100,
           window= lag(dateTime,12),
           T_valid1hrWindow=as.double(dateTime-lag(dateTime,12),units='mins'), # lag 12 bc lag already grabs 1 row above
           T_upstreamMaxViolation=ifelse(upstream>maxT1,1,0),# 9VAC25-260-50. Numerical Criteria for Maximum Temperature
           T_downstreamMaxViolation=ifelse(downstream>maxT2,1,0), # 9VAC25-260-50. Numerical Criteria for Maximum Temperature
           T_riseAboveNaturalViolation=ifelse(numericDiff>natTchange,1,0)) # 9VAC25-260-60. Rise Above Natural Temperature
    
  temperatureByWQSclass(together, maxHourlyTchange1, maxHourlyTchange2)%>%
    mutate(WQSapplied=ifelse(WQclassGage1==WQclassGage2,WQclassGage1,paste(WQclassGage1,";",WQclassGage2,sep="")))
}




# Dissolved Oxygen Analysis
dissolvedOxygenByWQSclass <- function(x,minDO1,minDO2){
  together <- mutate(x,DO_upstreamMinViolation=ifelse(upstream<minDO1,1,0), # 9VAC25-260-50. Numerical Criteria for Dissolved Oxygen
                     DO_downstreamMinViolation=ifelse(downstream<minDO2,1,0),  # 9VAC25-260-50. Numerical Criteria for Dissolved Oxygen
                     DO_upDownDifferenceViolation = ifelse(abs(numericDiff)>1,1,0)) # Flag if upstream/downstream difference > 1 mg/L
  return(together)
}

dissolvedOxygen <- function(upstreamData, downstreamData, parameter, WQclassGage1, WQclassGage2){
  # Add DO_Inst fields if no DO_Inst data retrieved from probes, populate with NA 
  if(unique(c('DO_Inst',"DO_Inst_cd") %in% names(upstreamData))){
    UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  if(unique(c('DO_Inst',"DO_Inst_cd") %in% names(downstreamData))){
    DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  # Establish thresholds based on WQS for min DO violations
  if(WQclassGage1==6){minDO1 <- 8}
  if(WQclassGage2==6){minDO2 <- 8}
  if(WQclassGage1 %in% c(3,4,5)){minDO1 <- 7}
  if(WQclassGage2 %in% c(3,4,5)){minDO2 <- 7}
  
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime')) %>%
    # Ensure data is realistic, replace wonky probe readings with NA so they are skipped
    # Filtering out bad data is not a good idea because you will lose entire rows of data and subsequently
    #  throw off other analyses down the line
    mutate(upstream = replace(upstream, upstream < 0 | upstream > 25, NA ),
           downstream = replace(downstream, downstream < 0 | downstream > 25, NA )) %>%
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100)
  
  dissolvedOxygenByWQSclass(together,minDO1,minDO2)%>%
    mutate(WQSapplied=ifelse(WQclassGage1==WQclassGage2,WQclassGage1,paste(WQclassGage1,";",WQclassGage2,sep="")))
}




# pH Analysis
pHbyWQS <- function(x,pHspecialStandards,pHrangeAllowance){
  # Adjust standards if WQS special standard for site
  min_pH <- ifelse(pHspecialStandards=="Y",6.5,6)
  max_pH <- ifelse(pHspecialStandards=="Y",9.5,9)
  
  together <- mutate(x,pH_upstreamViolation=ifelse(upstream > min_pH & upstream < max_pH,0,1), # 9VAC25-260-50. Numerical Criteria for pH
                     pH_downstreamViolation=ifelse(downstream > min_pH & downstream < max_pH,0,1),  # 9VAC25-260-50. Numerical Criteria for pH
                     pH_upDownDifferenceViolation = ifelse(abs(numericDiff)>pHrangeAllowance,1,0)) # Flag if upstream/downstream difference > allowed range
  return(together)
}

pH <- function(upstreamData, downstreamData, parameter, pHspecialStandards, pHrangeAllowance){
  # Add pH_Inst fields if no pH_Inst data retrieved from probes, populate with NA 
  if(unique(c('pH_Inst',"pH_Inst_cd") %in% names(upstreamData))){
    UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  if(unique(c('pH_Inst',"pH_Inst_cd") %in% names(downstreamData))){
    DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime')) %>%
    # Ensure data is realistic, replace wonky probe readings with NA so they are skipped
    # Filtering out bad data is not a good idea because you will lose entire rows of data and subsequently
    #  throw off other analyses down the line
    mutate(upstream = replace(upstream, upstream < 0 | upstream > 14, NA ),
           downstream = replace(downstream, downstream < 0 | downstream > 14, NA )) %>%
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
  # Add SpecCond_Inst fields if no SpecCond_Inst data retrieved from probes, populate with NA 
  if(unique(c('SpecCond_Inst',"SpecCond_Inst_cd") %in% names(upstreamData))){
    UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  if(unique(c('SpecCond_Inst',"SpecCond_Inst_cd") %in% names(downstreamData))){
    DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime')) %>%
    # Ensure data is realistic, replace wonky probe readings with NA so they are skipped
    # Filtering out bad data is not a good idea because you will lose entire rows of data and subsequently
    #  throw off other analyses down the line
    mutate(upstream = replace(upstream, upstream < 0, NA ),
           downstream = replace(downstream, downstream < 0, NA )) %>%
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100)
  
  SpCondByDesignation(together, SpCond_Designation)
}




# Turbidity analysis
turbidity <- function(upstreamData, downstreamData, parameter, turbidityBaseline, turbidity99th1, turbidity99th2){
  # only run function if Turbidity data came from both datasets
  if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(upstreamData))){
    UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
    # If gage height data available at gage then grab that, too
    if("GH_Inst" %in% names(upstreamData)){
      if(unique(upstreamData$site_no)[!is.na(unique(upstreamData$site_no))]=="0205450393"){
        gh <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,GH_Inst)%>%
          mutate(site_noGH='02054500')%>% dplyr::select(agency_cd,site_noGH,everything(),-site_no)
      }else{
        gh <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,GH_Inst)%>%
          rename(site_noGH=!!names(.[2]))}
      
    }else{gh <- dplyr::select(upstreamData,agency_cd,dateTime)%>%
        mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}
  }else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(downstreamData))){
    DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
    # If gage height data available at gage then grab that, too
    if("GH_Inst" %in% names(downstreamData)){
      gh <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,GH_Inst) %>%
        rename(site_noGH=!!names(.[2]))
    }else{
      if(unique(upstreamData$site_no)[!is.na(unique(upstreamData$site_no))]=="0205450393"){
        gh <- gh
        }else{
          gh <- dplyr::select(downstreamData,agency_cd,dateTime)%>%
            mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}}
      
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    full_join(gh,by=c('agency_cd','dateTime')) %>%
    # Ensure data is realistic, replace wonky probe readings with NA so they are skipped
    # Filtering out bad data is not a good idea because you will lose entire rows of data and subsequently
    #  throw off other analyses down the line
    mutate(upstream = replace(upstream, upstream < 0 | upstream > 1500, NA ),
           downstream = replace(downstream, downstream < 0 | downstream > 1500, NA )) %>%
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100,
           window= lag(dateTime,6),
           turbidity_valid30minuteWindow=dateTime-lag(dateTime,6)) # lag 6 bc lag already grabs 1 row above
  
  # if < 7 records the rolling functions will fail, report out NA's if not enough data to run full metrics
  if(nrow(together)<7){
    tidyverse_diff_rollstats <- mutate(together,turbidity_Exceedance=as.numeric(NA),turbidity_NAs=as.numeric(NA),
                                       turbidity_ExceedanceType=as.numeric(NA),turbidity_upstreamExceed99th=as.numeric(NA),
                                       turbidity_downstreamExceed99th=as.numeric(NA))
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
      }}
  #  Prevent tq_mutate from bombing out if entire dataset fed to it is NA's
  if(nrow(filter(together,!is.na(upstream))) < 7){
    tidyverse_diff_rollstats <- mutate(tidyverse_diff_rollstats,turbidity_upstreamExceed99th=NA,turbidity_upstreamNAs=NA)
  }else{
    # Exceed 99th percentile upstream or downstream
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
        threshold  = turbidity99th1) %>%
      rename(turbidity_upstreamExceed99th=turbidity_Exceedance.1,turbidity_upstreamNAs=turbidity_NAs.1)
  }
  #  Prevent tq_mutate from bombing out if entire dataset fed to it is NA's
  if(nrow(filter(together,!is.na(downstream))) < 7){
    tidyverse_diff_rollstats <- mutate(tidyverse_diff_rollstats,turbidity_downstreamExceed99th=NA,turbidity_downstreamNAs=NA)
  }else{
    tidyverse_diff_rollstats <- tidyverse_diff_rollstats %>%
      tq_mutate(
        select     = downstream,
        mutate_fun = rollapply, 
        # rollapply args
        width      = 7,
        align      = "right",
        by.column  = FALSE,
        FUN        = turbidityNumericThreshold,
        # FUN args
        threshold  = turbidity99th2) %>%
      rename(turbidity_downstreamExceed99th=turbidity_Exceedance.1,turbidity_downstreamNAs=turbidity_NAs.1) 
  }
  return(tidyverse_diff_rollstats)
}


## Apply threshold exceedance functions across time series, multitple parameters
dataScan <- function(upstreamData, downstreamData, WQclassGage1, WQclassGage2, pHspecialStandards, pHrangeAllowance, SpCond_Designation, turbidityBaseline, turbidity99th1, turbidity99th2){
  
  ## Temperature Analysis ##
  T_Results <- temperature(upstreamData,downstreamData,'Wtemp_Inst',WQclassGage1,WQclassGage2) %>% # will have to preprogram WQS class depending on site
    select(agency_cd,dateTime,site_no.x,site_no.y,WQSapplied,T_valid1hrWindow,T_upstreamMaxViolation,T_downstreamMaxViolation,
           T_upstreamHourlyChangeViolation,T_upstreamNAs,T_downstreamHourlyChangeViolation,T_downstreamNAs,T_riseAboveNaturalViolation)
  
  ## Dissolved Oxygen Analysis ##
  DO_Results <- dissolvedOxygen(upstreamData,downstreamData,'DO_Inst',WQclassGage1,WQclassGage2) %>% # will have to preprogram WQS class depending on site
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
  turbidity_Results <- turbidity(upstreamData, downstreamData, 'Turb_Inst', turbidityBaseline, turbidity99th1, turbidity99th2) %>% 
    select(agency_cd,dateTime,site_no.x,site_no.y,site_noGH,GH_Inst,turbidity_valid30minuteWindow,turbidity_NAs,turbidity_ExceedanceType,turbidity_Exceedance,
           turbidity_upstreamExceed99th,turbidity_upstreamNAs,turbidity_downstreamExceed99th,turbidity_downstreamNAs)
  
  allParameters <- plyr::join_all(list(T_Results,DO_Results, pH_Results, spCond_Results, turbidity_Results),by=c('agency_cd','dateTime','site_no.x','site_no.y'))
  
  return(allParameters)
}

#------------------------## Notification Functions ##---------------------------------------#
# Based on column name, output upstream or downstream gage number
updownDifference <-function(colName,gageResults){
  if(length(grep('upstream',colName))){
    return(unique(gageResults$site_no.x)[!is.na(unique(gageResults$site_no.x))])}
  if(length(grep('downstream',colName))){
    return(unique(gageResults$site_no.y)[!is.na(unique(gageResults$site_no.y))])}}  

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
    PCcolName <- gsub("Violation","DataReview",columnName) # PC name for tweet
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
      gageNumber1 <- unique(minmaxData$site_no.x)[!is.na(unique(minmaxData$site_no.x))]
      gageNumber2 <- unique(minmaxData$site_no.y)[!is.na(unique(minmaxData$site_no.y))]
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
      PCcolName <- 'T_upstreamHourlyChangeDataReview'
      gageNumber <- unique(validData$site_no.x)[!is.na(unique(validData$site_no.x))]
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
      PCcolName <- 'T_downstreamHourlyChangeDataReview'
      gageNumber <- unique(validData$site_no.y)[!is.na(unique(validData$site_no.y))]
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
noGageDataTweet <- function(gageResults,i,last2Hours){
  if(nrow(filter(gageResults,!is.na(site_no.x))) == 0 | nrow(filter(gageResults,!is.na(site_no.y))) == 0){
    uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
    date1 <- as.Date(last2Hours)-3 # Need to hardwire date into link so get same results well after original notification
    date2 <- as.Date(last2Hours) # Need to hardwire date into link so get same results well after original notification
    parameterCode <- '&cb_00010=on&cb_00095=on&cb_00300=on&cb_00400=on&cb_63680=on'
    
    # find Culprit
    if(nrow(filter(gageResults,!is.na(site_no.x))) == 0){
      columnName <- 'upstreamNoData'
      gageNumberMissing <- paste(0,gageInfo$`USGS Station ID`[i],sep = "")}
    if(nrow(filter(gageResults,!is.na(site_no.y))) == 0){
      columnName <- 'downstreamNoData'
      gageNumberMissing <- paste(0,gageInfo$`USGS Station ID`[i+1],sep = "")}
    if(nrow(filter(gageResults,!is.na(site_no.x))) == 0 &
       nrow(filter(gageResults,!is.na(site_no.y))) == 0){columnName <- 'NoData'}
    # If both gages are missing, give special tweet
    if(columnName=='NoData'){
      streamName <- gageInfo$`Stream Name`[i]
      gageNumberMissing1 <- paste(0,gageInfo$`USGS Station ID`[i],sep="")
      gageNumberMissing2 <- paste(0,gageInfo$`USGS Station ID`[i+1],sep = "")
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       parameterCode,'&site_no=',gageNumberMissing1,sep = "")
      weblink2 <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                        parameterCode,'&site_no=',gageNumberMissing2,sep = "")
      tweet(paste(streamName,columnName,weblink," & ",weblink2,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
      # fix gageResults so we know which gages blank data associated with
      gageResults$site_no.x <- gageNumberMissing1
      gageResults$site_no.y <- gageNumberMissing2
      }else{weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       parameterCode,'&site_no=',gageNumberMissing,sep = "")}
    
    tweet(paste(streamName,columnName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T)
    
    # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
    write.csv(gageResults,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)}
}

# turbidity Plots
pairedTurbidityPlot <- function(turbidityData){
  #Figure out which gage the GH data comes from
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  ggplot(turbidityData,aes(x = dateTime)) + 
    geom_point(aes(y=upstream, colour='Upstream Gage'))  + 
    geom_point(data=turbidityData, aes(dateTime,downstream,colour='Downstream Gage')) +
    geom_line(data=turbidityData, aes(dateTime,GH_Inst,colour=paste(updown,' Gage Height',sep="")),linetype=2) +
    scale_y_continuous(sec.axis = sec_axis(trans=~.+0,name = 'Gage Height (ft)')) +
    scale_color_manual(values=c('blue','red','gray'))+
    labs(y="Turbidity (NTU)",x="Date",colour='Parameter')
}

upstreamTurbidityPlot <- function(turbidityData,turbidity99th1){
  #Figure out which gage the GH data comes from
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  ggplot(turbidityData,aes(x = dateTime)) + 
    geom_point(aes(y=upstream, colour='Upstream Gage'))  + 
    geom_line(data=turbidityData,aes(dateTime,turbidity99th1,colour='turbidity99th')) +
    geom_line(data=turbidityData, aes(dateTime,GH_Inst,colour=paste(updown,' Gage Height',sep="")),linetype=2) +
    scale_y_continuous(sec.axis = sec_axis(trans=~.+0,name = 'Gage Height (ft)')) +
    scale_color_manual(values=c('gray','red','blue'))+
    labs(y="Turbidity (NTU)",x="Date",colour='Parameter')
}

downstreamTurbidityPlot <- function(turbidityData,turbidity99th2){
  #Figure out which gage the GH data comes from
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  ggplot(turbidityData,aes(x = dateTime)) + 
    geom_point(data=turbidityData, aes(dateTime,downstream,colour='Downstream Gage')) +
    geom_line(data=turbidityData,aes(dateTime,turbidity99th2,colour='turbidity99th')) +
    geom_line(data=turbidityData, aes(dateTime,GH_Inst,colour=paste(updown,' Gage Height',sep="")),linetype=2) +
    scale_y_continuous(sec.axis = sec_axis(trans=~.+0,name = 'Gage Height (ft)')) +
    scale_color_manual(values=c('blue','gray','red'))+
    labs(y="Turbidity (NTU)",x="Date",colour='Parameter')
}

# Turbidity tweet
turbTimeTweet <- function(gageResults,upstreamData,downstreamData,turbidity99th1,turbidity99th2){
  turbTimeData <- select(gageResults,dateTime,site_no.x,site_no.y,site_noGH,GH_Inst,WQSapplied,turbidity_valid30minuteWindow,
                         turbidity_NAs,turbidity_ExceedanceType,turbidity_Exceedance,
                         turbidity_upstreamExceed99th,turbidity_upstreamNAs,
                         turbidity_downstreamExceed99th,turbidity_downstreamNAs)
  
  # Only use data from correct 1hr window
  validData <- filter(turbTimeData,turbidity_valid30minuteWindow==30)
  
  
  if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(upstreamData))){
    UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,'Turb_Inst')%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
    # If gage height data available at gage then grab that, too
    if("GH_Inst" %in% names(upstreamData)){
      if(unique(upstreamData$site_no)[!is.na(unique(upstreamData$site_no))]=="0205450393"){
        gh <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,GH_Inst)%>%
          mutate(site_noGH='02054500')%>% dplyr::select(agency_cd,site_noGH,everything(),-site_no)
      }else{
        gh <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,GH_Inst)%>%
          rename(site_noGH=!!names(.[2]))}
      
    }else{gh <- dplyr::select(upstreamData,agency_cd,dateTime)%>%
      mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}
  }else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(downstreamData))){
    DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,'Turb_Inst')%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
    # If gage height data available at gage then grab that, too
    if("GH_Inst" %in% names(downstreamData)){
      gh <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,GH_Inst) %>%
        rename(site_noGH=!!names(.[2]))
    }else{
      if(unique(upstreamData$site_no)[!is.na(unique(upstreamData$site_no))]=="0205450393"){
        gh <- gh
      }else{
        gh <- dplyr::select(downstreamData,agency_cd,dateTime)%>%
          mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}}
    
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  
  ## organize Turbidity data for ggplots
  #if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(upstreamData))){
  #  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,Turb_Inst)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  #  # If gage height data available at gage then grab that, too
  #  if("GH_Inst" %in% names(upstreamData)){
  #    if(unique(upstreamData$site_no)[!is.na(unique(upstreamData$site_no))]=="0205450393"){
  #      gh <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,GH_Inst)%>%
  #        mutate(site_noGH='02054500')%>% dplyr::select(agency_cd,site_noGH,everything(),-site_no)
  #    }else{
  #      gh <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,GH_Inst)%>%
  #        rename(site_noGH=!!names(.[2]))}
  #    
  #    }else{gh <- dplyr::select(upstreamData,agency_cd,dateTime)%>%
  #        mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}
  #}else{UP <- select(upstreamData,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
  #if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(downstreamData))){
  #  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,Turb_Inst)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  #  # If gage height data available at gage then grab that, too
  #  if("GH_Inst" %in% names(downstreamData)){
  #    gh <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,GH_Inst) %>%
  #      rename(site_noGH=!!names(.[2]))
  #  }else{gh <- dplyr::select(downstreamData,agency_cd,dateTime)%>%
  #    mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}
  #}else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    full_join(gh,by=c('agency_cd','dateTime')) %>%
    mutate(turbidity99th1=turbidity99th1,turbidity99th2=turbidity99th2)
  
  
  # upstream downstream change comparison
  if(nrow(validData)>0){
    validComparison <- filter(validData,turbidity_Exceedance>0, # Limit dataset to only violations
                              turbidity_NAs <= 3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    if(nrow(validComparison)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'turbidity_upstreamDownstreamDataReview'
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      gageNumber1 <- unique(validData$site_no.x)[!is.na(unique(validData$site_no.x))]
      gageNumber2 <- unique(validData$site_no.y)[!is.na(unique(validData$site_no.y))]
      streamName <- gageInfo[grep(as.numeric(gageNumber1),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_63680=on&site_no=',gageNumber1,'%2C',gageNumber2,'&format=gif_mult_sites',sep='')
      file_name <- paste('notifications/images/',uniqueTweetIdentifier,'.jpg',sep='')
      jpeg(file_name)
      print(pairedTurbidityPlot(together))
      dev.off()
      
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T,
            mediaPath=paste('notifications/images/',uniqueTweetIdentifier,'.jpg',sep=''))
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      write.csv(together,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)
    }
    
    # Upstream absolute threshold Analysis
    validUP <- filter(validData,turbidity_upstreamExceed99th>0, # Limit dataset to only violations
                      turbidity_upstreamNAs <=3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    # Only send tweet if at least 1 valid row (1 hr violation) & <= 3 missing turbidity readings within hour
    if(nrow(validUP)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'turbidity_upstreamExceed99th'
      gageNumber <- unique(validData$site_no.x)[!is.na(unique(validData$site_no.x))]
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_63680=on&site_no=',gageNumber,sep = "")
      file_name <- paste('notifications/images/',uniqueTweetIdentifier,'.jpg',sep='')
      jpeg(file_name)
      print(upstreamTurbidityPlot(together))
      dev.off()
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T,
            mediaPath=paste('notifications/images/',uniqueTweetIdentifier,'.jpg',sep=''))
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- turbTimeData[,c(1,2,7,11,12)]
      upstreamRAW <- dplyr::select(together,-c(site_no.y,downstream))
      datToSave <- plyr::join_all(list(datToSave,upstreamRAW),by='dateTime')
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)
    }
    # Downstream Analysis
    validDOWN <- filter(validData,turbidity_downstreamExceed99th>0, # Limit dataset to only violations
                        turbidity_downstreamNAs <=3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    # Only send tweet if at least 1 valid row (1 hr violation) & <= 4 missing turbidity readings within hour
    if(nrow(validDOWN)>0){
      uniqueTweetIdentifier <- format(as.numeric(Sys.time()),nsmall=4) # Tweet function will not work if the status is the same as a previous status
      PCcolName <- 'turbidity_downstreamExceed99th'
      gageNumber <- unique(validData$site_no.x)[!is.na(unique(validData$site_no.x))]
      streamName <- gageInfo[grep(as.numeric(gageNumber),gageInfo$`USGS Station ID`),]$`Stream Name`[1]
      date1 <- as.Date(validData$dateTime[1])-3 # Need to hardwire date into link so get same results well after original notification
      date2 <- as.Date(validData$dateTime[1]) # Need to hardwire date into link so get same results well after original notification
      weblink <- paste('waterdata.usgs.gov/va/nwis/uv?period=&begin_date=',date1,'&end_date=',date2,
                       '&cb_63680=on&site_no=',gageNumber,sep = "")
      file_name <- paste('notifications/images/',uniqueTweetIdentifier,'.jpg',sep='')
      jpeg(file_name)
      print(downstreamTurbidityPlot(together))
      dev.off()
      tweet(paste(streamName,PCcolName,weblink,'UID:',uniqueTweetIdentifier,sep=" "),bypassCharLimit=T,
            mediaPath=paste('notifications/images/',uniqueTweetIdentifier,'.jpg',sep=''))
      # Save a record of the data that caused the tweet, indexed by uniqueTweetIdentifier
      datToSave <- turbTimeData[,c(1,3,7,13,14)]
      downstreamRAW <- dplyr::select(together,-c(site_no.x,upstream))
      datToSave <- plyr::join_all(list(datToSave,downstreamRAW),by='dateTime')
      write.csv(datToSave,paste('notifications/',uniqueTweetIdentifier,'.csv',sep=""),row.names = F)}
  }
}



#-----------------------------------------------------------------------------------------------
