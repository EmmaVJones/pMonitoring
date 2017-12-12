suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyquant))  # Loads tidyverse, tidquant, financial pkgs, xts/zoo
suppressPackageStartupMessages(library(dataRetrieval))
suppressPackageStartupMessages(library(gmailr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))

# Email authorization
gmail_auth("compose",secret_file="pipelinemonitoring.json")


## Functions required ##

# Pull NWIS data
NWISpull <- function(gageNo,start,end){
  allUnitData <- readNWISuv(siteNumbers=gageNo,
                            parameterCd=c("00010", "00095", "00300", "00400", "63680"),
                            startDate=as.Date(start,"%Y-%m-%d"),
                            endDate=as.Date(end,"%Y-%m-%d"),
                            tz='America/New_York')
  allUnitData <- renameNWISColumns(allUnitData)
}

# Find Threshold Exceedances
numericThreshold <- function(x, threshold){ ifelse(min(x,na.rm=T)>threshold,1,0) }
percentThreshold <- function(x,threshold){ ifelse(min(x,na.rm=T)>=threshold,1,0) }

# Apply threshold exceedance functions across time series
dataScan <- function(upstreamData, downstreamData, parameter, baseline, unitIncrease, percentIncrease){
  
  UP <- dplyr::select(upstreamData,agency_cd,site_no,dateTime,parameter)%>%rename(upstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  DOWN <- dplyr::select(downstreamData,agency_cd,site_no,dateTime,Turb_Inst)%>%rename(downstream=!!names(.[4])) # change parameter to general name to make further manipulations easier
  
  together <- full_join(UP,DOWN,by=c('agency_cd','dateTime'))%>%
    mutate(parameter=parameter,
           numericDiff=downstream-upstream,
           pctDiff=(numericDiff/downstream)*100)
  
  # Choose method based on baseline parameter value
  if(min(UP$upstream, na.rm = T) >= baseline && min(DOWN$downstream, na.rm = T) >= baseline){
    tidyverse_diff_rollstats <- together %>%
      tq_mutate(
        select     = pctDiff,
        mutate_fun = rollapply, 
        # rollapply args
        width      = 12,
        align      = "right",
        by.column  = FALSE,
        FUN        = percentThreshold,
        # FUN args
        threshold  = percentIncrease,
        # tq_mutate args
        col_rename = "exceedPCTdiff"
      )
  }else{
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
        threshold  = unitIncrease,
        # tq_mutate args
        col_rename = "exceedNUMERICdiff"
      )
  }
  return(tidyverse_diff_rollstats)
}

# Send Email Alerts
sendEmail_function <- function(x){
  exceedances <- x%>%dplyr::filter(.[[10]]==1) # Only get exceedance records
  
  # Record important gage/parameter information
  gage <- exceedances[1,5]
  param <- exceedances[1,7]
  
  # If any exceedances are noted, send email alert
  if(nrow(exceedances)>0){ 
    addresses <- suppressMessages(read_csv("NotificationDocumentation/addresses_TEST.csv")) # email address list
    datetime <- format(Sys.time(), "%Y-%m-%d_%H%M%S") # keep %H%M%S same for the .RDS file as well to enable file matching
    datetimeTEXT <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    emailSubject <- paste(param,"at USGS Gage:",gage,sep=" ")
    email_sender <- 'Monitoring Server <MonitoringSystem2017@gmail.com>' # your Gmail address
    #optional_bcc <- 'Anonymous <emma.v.jones@icloud.com>'     # for me, TA address
    body <- paste("Hi, %s.\n\nPlease visit link to dashboard to review ",emailSubject,".\n\n",
                  "Notification Time:",datetimeTEXT,
                  "\n\nCheers,\n",
                  "VDEQ Pipeline Monitoring Team")
    
    emailData <- addresses %>%
      mutate(
        To = sprintf('%s <%s>', name, email),
        #Bcc = optional_bcc,
        From = email_sender,
        Subject = sprintf('Notification for %s', emailSubject),
        body = sprintf(body, name, emailSubject)) %>%
      select(To, #Bcc, 
             From, Subject, body)
    
    # Keep a record of email composition
    write_csv(emailData, paste("NotificationDocumentation/composedEmails/Notification",datetime,".csv",sep=""))
    
    emails <- emailData %>%
      pmap(mime)
    #use_secret_file("pipelinemonitoring.json") # moved to top of script
    
    safe_send_message <- safely(send_message) # Safely send message (i.e. don't bomb out rest of list if one email isn't correct)
    
    # Record of what emails were sent correctly
    sent_mail <- emails %>%
      map(safe_send_message)
    
    # Keep a record of what was send to whom and when
    saveRDS(sent_mail,paste("NotificationDocumentation/sendReceipts/sent-emails_",datetime,".RDS", sep = ""))
  }
}



## Run Analyses ##
last2Hours <- format(Sys.time()-7200,format="%Y-%m-%d %H:00:00") # Need to filter by the last 2 hours bc exceedance could start mid-hour and be lost if constantly resetting on hour

upstreamData <- NWISpull('0205450393', Sys.Date()-1,Sys.Date())%>%#format(Sys.time()-21600,format="%Y-%m-%d %H:%M:%S"), format(Sys.time(),format="%Y-%m-%d %H:%M:%S"))
  dplyr::filter(dateTime >= last2Hours)
downstreamData <- NWISpull('0205450495',Sys.Date()-1,Sys.Date())%>%
  dplyr::filter(dateTime >= last2Hours)


test <- dataScan(upstreamData, downstreamData, 'Turb_Inst', 40, -0.6, 15)
sendEmail_function(test)