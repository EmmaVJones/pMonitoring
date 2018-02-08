1=5

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
together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
  full_join(gh,by=c('agency_cd','dateTime')) %>%
  mutate(turbidity99th1=turbidity99th1,turbidity99th2=turbidity99th2)


turbidityData <- together


#pairedTurbidityPlot <- function(turbidityData){
  #Figure out which gage the GH data comes from
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  up <- dplyr::select(turbidityData,dateTime,upstream)%>%
    gather(gage,measure,upstream)%>%mutate(Turbidity99=turbidity99th1,panel=2)#'Upstream Turbidity (NTU)')
  down <- dplyr::select(turbidityData,dateTime,downstream)%>%
    gather(gage,measure,downstream)%>%mutate(Turbidity99=turbidity99th2,panel=3)#'Downstream Turbidity (NTU)')
  gh <- dplyr::select(turbidityData,dateTime,GH_Inst)%>%
    rename(measure=GH_Inst)%>% mutate(gage=updown,Turbidity99=NA,panel=1)%>% select(dateTime,gage,Turbidity99,measure,panel)#paste(updown,'Gage Height (ft)',sep=' ')) 
  dat <- rbind(up,down,gh)
 correctName <- c(`1`=paste(updown,'Gage Height (ft)',sep=' '),
                   `2`='Upstream Turbidity (NTU)',`3`='Downstream Turbidity (NTU)')
                   
# 3 plots
plot3_last2hours <- ggplot(dat,mapping = aes(x=dateTime,y=measure))+
  facet_grid(panel~.,scale='free',labeller = as_labeller(correctName))+
  geom_point(data=up,stat = 'identity',colour='coral')+
  geom_line(data=up,aes(dateTime,Turbidity99),colour='red')+
  geom_point(data=down,stat = 'identity',colour='blueviolet')+
  geom_line(data=down,aes(dateTime,Turbidity99),colour='blue')+
  geom_point(data=gh,stat = 'identity',colour='gray')
    
}
  
    
up <- dplyr::select(turbidityData,dateTime,upstream)%>%
  gather(gage,measure,upstream)%>%mutate(Turbidity99=turbidity99th1,panel='Turbidity (NTU)')
down <- dplyr::select(turbidityData,dateTime,downstream)%>%
  gather(gage,measure,downstream)%>%mutate(Turbidity99=turbidity99th2,panel='Turbidity (NTU)')
gh <- dplyr::select(turbidityData,dateTime,GH_Inst)%>%
  rename(measure=GH_Inst)%>% mutate(gage=updown,Turbidity99=NA,panel=paste(updown,'Gage Height (ft)',sep=' ')) %>% select(dateTime,gage,Turbidity99,measure,panel)
dat <- rbind(up,down,gh)


# 2 plots
plot2 <- ggplot(dat,mapping = aes(x=dateTime,y=measure))+
  facet_grid(panel~.,scale='free')+
  
  
  geom_point(data=up,stat = 'identity',colour='coral')+
  geom_line(data=up,aes(dateTime,Turbidity99),colour='red')+
  
  geom_point(data=down,stat = 'identity',colour='blueviolet')+
  geom_line(data=down,aes(dateTime,Turbidity99),colour='blue')+
  geom_point(data=gh,stat = 'identity',colour='gray')
  
  


# 2 plots with entire record & last 2 hours highlighted
upstreamData1 <- upstreamDataForTurbidity
downstreamData1 <- downstreamDataForTurbidity


if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(upstreamData1))){
  UP <- dplyr::select(upstreamData1,agency_cd,site_no,dateTime,'Turb_Inst')%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  # If gage height data available at gage then grab that, too
  if("GH_Inst" %in% names(upstreamData1)){
    if(unique(upstreamData1$site_no)[!is.na(unique(upstreamData1$site_no))]=="0205450393"){
      gh <- dplyr::select(upstreamData1,agency_cd,site_no,dateTime,GH_Inst)%>%
        mutate(site_noGH='02054500')%>% dplyr::select(agency_cd,site_noGH,everything(),-site_no)
    }else{
      gh <- dplyr::select(upstreamData1,agency_cd,site_no,dateTime,GH_Inst)%>%
        rename(site_noGH=!!names(.[2]))}
    
  }else{gh <- dplyr::select(upstreamData1,agency_cd,dateTime)%>%
    mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}
}else{UP <- select(upstreamData1,agency_cd,site_no,dateTime)%>%mutate(upstream=NA)}
if(unique(c('Turb_Inst',"Turb_Inst_cd") %in% names(downstreamData1))){
  DOWN <- dplyr::select(downstreamData1,agency_cd,site_no,dateTime,'Turb_Inst')%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  # If gage height data available at gage then grab that, too
  if("GH_Inst" %in% names(downstreamData1)){
    gh <- dplyr::select(downstreamData1,agency_cd,site_no,dateTime,GH_Inst) %>%
      rename(site_noGH=!!names(.[2]))
  }else{
    if(unique(upstreamData1$site_no)[!is.na(unique(upstreamData1$site_no))]=="0205450393"){
      gh <- gh
    }else{
      gh <- dplyr::select(downstreamData1,agency_cd,dateTime)%>%
        mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}}
  
}else{DOWN <- select(downstreamData1,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}


together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
  full_join(gh,by=c('agency_cd','dateTime')) %>%
  mutate(turbidity99th1=turbidity99th1,turbidity99th2=turbidity99th2)
turbidityData <- together

pairedTurbidityPlot_2days <- function(turbidityData){
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  up <- dplyr::select(turbidityData,dateTime,upstream)%>%
    gather(gage,measure,upstream)%>%mutate(Turbidity99=turbidity99th1,panel=2)#'Upstream Turbidity (NTU)')
  down <- dplyr::select(turbidityData,dateTime,downstream)%>%
    gather(gage,measure,downstream)%>%mutate(Turbidity99=turbidity99th2,panel=3)#'Downstream Turbidity (NTU)')
  gh <- dplyr::select(turbidityData,dateTime,GH_Inst)%>%
    rename(measure=GH_Inst)%>% mutate(gage=updown,Turbidity99=NA,panel=1)%>% select(dateTime,gage,Turbidity99,measure,panel)#paste(updown,'Gage Height (ft)',sep=' ')) 
  dat <- rbind(up,down,gh)
  dat2 <- dplyr::filter(dat,dateTime >= last2Hours & panel==2) %>%
    mutate(a=min(dateTime),b=max(dateTime),c=min(measure,na.rm=T),
           d=max(max(measure,na.rm=T),max(Turbidity99)))
  dat3 <- dplyr::filter(dat,dateTime >= last2Hours & panel==3) %>%
    mutate(a=min(dateTime),b=max(dateTime),c=min(measure,na.rm=T),
           d=max(max(measure,na.rm=T),max(Turbidity99)))
  
  #dat2 <- mutate(dat,c=min(measure,na.rm=T),d=max(measure,na.rm=T))%>%
  #  dplyr::filter(dateTime >= last2Hours) %>%
  #  mutate(a=min(dateTime),b=max(dateTime))
  #dat3 <- dat2 %>% mutate(panel=2)
  #dat2 <- dat2 %>% mutate(panel=3)
  
  correctName <- c(`1`=paste(updown,'Gage Height (ft)',sep=' '),
                   `2`='Upstream Turbidity (NTU)',`3`='Downstream Turbidity (NTU)')
  
  
  plot3_alldata2hr <- ggplot(dat,mapping = aes(x=dateTime,y=measure))+
    facet_grid(panel~.,scale='free',labeller = as_labeller(correctName))+
    geom_rect(data=dat2,aes(xmin=a,xmax=b,ymin=c,ymax=d,fill=TRUE),
              alpha=0.2)+
    geom_rect(data=dat3,aes(xmin=a,xmax=b,ymin=c,ymax=d,fill=TRUE),
              alpha=0.2)+
    
    geom_point(data=up,stat = 'identity',colour='coral')+
    geom_line(data=up,aes(dateTime,Turbidity99),colour='red')+
    
    geom_point(data=down,stat = 'identity',colour='blueviolet')+
    geom_line(data=down,aes(dateTime,Turbidity99),colour='blue')+
    geom_point(data=gh,stat = 'identity',colour='gray') +
    
    
    scale_fill_manual(values='gray')+guides(fill=F)
  return(plot3_alldata2hr)
}



file_name <- paste('figures/BCU3plots_2hrs.jpg',sep='')
jpeg(file_name)
print(plot3_last2hours)
dev.off()


file_name <- paste('figures/BCUplot.jpg',sep='')
jpeg(file_name)
print(plot2)
dev.off()


file_name <- paste('figures/BCU3plots_alldata.jpg',sep='')
jpeg(file_name)
print(plot3_alldata2hr )
dev.off()

























# Changes that would need to be made to incorporate 2 days of data in turb plots

# all parameterAlertSystem
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



# Send notification if any turbidity change flags are blown
turbTimeTweet(gageResults,upstreamData,downstreamData,turbidity99th1=gageInfo$Turbidity_99th[i],
              turbidity99th2=gageInfo$Turbidity_99th[i+1],upstreamDataForTurbidity,
              downstreamDataForTurbidity)






# parameter Functions

dataManipulationForTurbTweets_2hr <- function(upstreamData,downstreamData,turbidity99th1,turbidity99th2){
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
        gh <- gh}
      if(exists('gh')){gh <- gh}else{
        gh <- dplyr::select(downstreamData,agency_cd,dateTime)%>%
          mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}}
    
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    full_join(gh,by=c('agency_cd','dateTime')) %>%
    mutate(turbidity99th1=turbidity99th1,turbidity99th2=turbidity99th2)
  return(together)
}


dataManipulationForTurbTweets_2days <- function(upstreamData,downstreamData,turbidity99th1,turbidity99th2){
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
        gh <- gh}
      if(exists('gh')){gh <- gh}else{
        gh <- dplyr::select(downstreamData,agency_cd,dateTime)%>%
          mutate(site_noGH=NA,GH_Inst=NA)%>%dplyr::select(agency_cd,site_noGH,dateTime,GH_Inst)}}
    
  }else{DOWN <- select(downstreamData,agency_cd,site_no,dateTime)%>%mutate(downstream=NA)}
  
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    full_join(gh,by=c('agency_cd','dateTime')) %>%
    mutate(turbidity99th1=turbidity99th1,turbidity99th2=turbidity99th2)
  return(together)
}

pairedTurbidityPlot_2days <- function(turbidityData,turbidity99th1,turbidity99th2){
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  up <- dplyr::select(turbidityData,dateTime,upstream)%>%
    gather(gage,measure,upstream)%>%mutate(Turbidity99=turbidity99th1,panel=2)#'Upstream Turbidity (NTU)')
  down <- dplyr::select(turbidityData,dateTime,downstream)%>%
    gather(gage,measure,downstream)%>%mutate(Turbidity99=turbidity99th2,panel=3)#'Downstream Turbidity (NTU)')
  gh <- dplyr::select(turbidityData,dateTime,GH_Inst)%>%
    rename(measure=GH_Inst)%>% mutate(gage=updown,Turbidity99=NA,panel=1)%>% select(dateTime,gage,Turbidity99,measure,panel)#paste(updown,'Gage Height (ft)',sep=' ')) 
  dat <- rbind(up,down,gh)
  dat2 <- dplyr::filter(dat,dateTime >= last2Hours & panel==2) %>%
    mutate(a=min(dateTime),b=max(dateTime),c=min(measure,na.rm=T),
           d=max(max(measure,na.rm=T),max(Turbidity99)))
  dat3 <- dplyr::filter(dat,dateTime >= last2Hours & panel==3) %>%
    mutate(a=min(dateTime),b=max(dateTime),c=min(measure,na.rm=T),
           d=max(max(measure,na.rm=T),max(Turbidity99)))
  correctName <- c(`1`=paste(updown,'Gage Height (ft)',sep=' '),
                   `2`='Upstream Turbidity (NTU)',`3`='Downstream Turbidity (NTU)')
  
  
  plot3_alldata2hr <- ggplot(dat,mapping = aes(x=dateTime,y=measure))+
    facet_grid(panel~.,scale='free',labeller = as_labeller(correctName))+
    geom_rect(data=dat2,aes(xmin=a,xmax=b,ymin=c,ymax=d,fill=TRUE),
              alpha=0.2)+
    geom_rect(data=dat3,aes(xmin=a,xmax=b,ymin=c,ymax=d,fill=TRUE),
              alpha=0.2)+
    
    geom_point(data=up,stat = 'identity',colour='coral')+
    geom_line(data=up,aes(dateTime,Turbidity99),colour='red')+
    
    geom_point(data=down,stat = 'identity',colour='blueviolet')+
    geom_line(data=down,aes(dateTime,Turbidity99),colour='blue')+
    geom_point(data=gh,stat = 'identity',colour='gray') +
    
    
    scale_fill_manual(values='gray')+guides(fill=F)
  return(plot3_alldata2hr)
}


upstreamTurbidityPlot2day <- function(turbidityData,turbidity99th1){
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  up <- dplyr::select(turbidityData,dateTime,upstream)%>%
    gather(gage,measure,upstream)%>%mutate(Turbidity99=turbidity99th1,panel=2)#'Upstream Turbidity (NTU)')
  gh <- dplyr::select(turbidityData,dateTime,GH_Inst)%>%
    rename(measure=GH_Inst)%>% mutate(gage=updown,Turbidity99=NA,panel=1)%>% select(dateTime,gage,Turbidity99,measure,panel)#paste(updown,'Gage Height (ft)',sep=' ')) 
  dat <- rbind(up,gh)
  dat2 <- dplyr::filter(dat,dateTime >= last2Hours & panel==2) %>%
    mutate(a=min(dateTime),b=max(dateTime),c=min(measure,na.rm=T),
           d=max(max(measure,na.rm=T),max(Turbidity99)))
  correctName <- c(`1`=paste(updown,'Gage Height (ft)',sep=' '),
                   `2`='Upstream Turbidity (NTU)')
  
  upstream2dayplot <- ggplot(dat,mapping = aes(x=dateTime,y=measure))+
    facet_grid(panel~.,scale='free',labeller = as_labeller(correctName))+
    geom_rect(data=dat2,aes(xmin=a,xmax=b,ymin=c,ymax=d,fill=TRUE),
              alpha=0.2)+
    geom_point(data=up,stat = 'identity',colour='coral')+
    geom_line(data=up,aes(dateTime,Turbidity99),colour='red')+
    geom_point(data=gh,stat = 'identity',colour='gray') +
    scale_fill_manual(values='gray')+guides(fill=F)
  return(upstream2dayplot)
}


downstreamTurbidityPlot2day <- function(turbidityData,turbidity99th2){
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.x)[!is.na(unique(turbidityData$site_no.x))]){updown <- "Upstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] ==
     unique(turbidityData$site_no.y)[!is.na(unique(turbidityData$site_no.y))]){updown <- "Downstream"}
  if(unique(turbidityData$site_noGH)[!is.na(unique(turbidityData$site_noGH))] == "02054500"){
    updown <- 'Lafayette'}
  
  down <- dplyr::select(turbidityData,dateTime,downstream)%>%
    gather(gage,measure,downstream)%>%mutate(Turbidity99=turbidity99th2,panel=3)#'Downstream Turbidity (NTU)')
  gh <- dplyr::select(turbidityData,dateTime,GH_Inst)%>%
    rename(measure=GH_Inst)%>% mutate(gage=updown,Turbidity99=NA,panel=1)%>% select(dateTime,gage,Turbidity99,measure,panel)#paste(updown,'Gage Height (ft)',sep=' ')) 
  dat <- rbind(down,gh)
  dat3 <- dplyr::filter(dat,dateTime >= last2Hours & panel==3) %>%
    mutate(a=min(dateTime),b=max(dateTime),c=min(measure,na.rm=T),
           d=max(max(measure,na.rm=T),max(Turbidity99)))
  correctName <- c(`1`=paste(updown,'Gage Height (ft)',sep=' '),`3`='Downstream Turbidity (NTU)')
  
  
  downstream2dayplot <- ggplot(dat,mapping = aes(x=dateTime,y=measure))+
    facet_grid(panel~.,scale='free',labeller = as_labeller(correctName))+
    geom_rect(data=dat3,aes(xmin=a,xmax=b,ymin=c,ymax=d,fill=TRUE),
              alpha=0.2)+
    geom_point(data=down,stat = 'identity',colour='blueviolet')+
    geom_line(data=down,aes(dateTime,Turbidity99),colour='blue')+
    geom_point(data=gh,stat = 'identity',colour='gray') +
    scale_fill_manual(values='gray')+guides(fill=F)
  return(downstream2dayplot)
}

# Turbidity tweet
turbTimeTweet <- function(gageResults,upstreamData,downstreamData,turbidity99th1,turbidity99th2,
                          upstreamDataForTurbidity,downstreamDataForTurbidity){
  turbTimeData <- select(gageResults,dateTime,site_no.x,site_no.y,site_noGH,GH_Inst,WQSapplied,turbidity_valid30minuteWindow,
                         turbidity_NAs,turbidity_ExceedanceType,turbidity_Exceedance,
                         turbidity_upstreamExceed99th,turbidity_upstreamNAs,
                         turbidity_downstreamExceed99th,turbidity_downstreamNAs)
  
  # Only use data from correct 1hr window
  validData <- filter(turbTimeData,turbidity_valid30minuteWindow==30)
  
  # upstream downstream change comparison
  if(nrow(validData)>0){
    validComparison <- filter(validData,turbidity_Exceedance>0, # Limit dataset to only violations
                              turbidity_NAs <= 3) # Only allow up to 3 missing turbidity reading per hour for any hour to count toward max change rate violation
    
    # Data manipulation for plots
    together <- dataManipulationForTurbTweets_2hr(upstreamData,downstreamData,turbidity99th1,turbidity99th2)
    together2days <- dataManipulationForTurbTweets_2days(upstreamDataForTurbidity,downstreamDataForTurbidity,turbidity99th1,turbidity99th2)
    
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
      print(pairedTurbidityPlot_2days(together2days,turbidity99th1,turbidity99th2))
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
      print(upstreamTurbidityPlot2day(together2days,turbidity99th1))
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
      print(downstreamTurbidityPlot2day(together2days,turbidity99th2))
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









