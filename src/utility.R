
## removes random effects from an epimodel object
rm_re <- function(object) {
  len <- NROW(as.matrix(object))
  locs <- grep("b\\[", object$orig_names)
  nchains <- length(object$stanfit@sim$samples)
  for (loc in locs) {
    for (chain in 1:nchains) {
      object$stanfit@sim$samples[[chain]][[loc]] <- rep(0., len)
    }
  }
  return(object)
}

rm_re_rw <- function(object) {
  len <- NROW(as.matrix(object))
  locs <- c(grep("b\\[", object$orig_names),grep("ac_", object$orig_names))
  nchains <- length(object$stanfit@sim$samples)
  for (loc in locs) {
    for (chain in 1:nchains) {
      object$stanfit@sim$samples[[chain]][[loc]] <- rep(0., len)
    }
  }
  return(object)
}



## add additional data for 3 weeks
get_newdata <- function(res,areaname="ltla"){
    df <- res$fit$data
    case_date <- df$date[which(!is.na(df$Cases_week))[1]]
    death_date <- df$date[which(!is.na(df$Deaths_week))[1]]
    if (is.element("admissions_week",colnames(df)))
        admissions_date <- df$date[which(!is.na(df$admissions_week))[1]]

    a <- as.character(df[[areaname]][1])
    df <- subset(df[df[[areaname]] == a,],select=-group)
#    c("ltla", "date", "Cases_week", "Deaths_week","averageMobility","averageRt","logitRt")]
    N <- 29
    dfnew <- data.frame(a, date=max(df$date) + 1:N,
                        Cases_week=NA, Deaths_week=NA)
    if (is.element("admissions_week",colnames(df)))
        dfnew$admissions_week <- NA
    colnames(dfnew)[1] <- areaname
    for (j in colnames(df)){
        if (!is.element(j,colnames(dfnew))) ## carry last value forward
            dfnew[,j] <- df[,j][df$date==max(df$date)]
    }
#                            averageMobility=df$averageMobility[df$date==max(df$date)],
#                        averageRt=df$averageRt[df$date==max(df$date)],
#                        logitRt=df$logitRt[df$date==max(df$date)],

    df <- rbind(df, dfnew)
    w <- (df$date> res$last_obs) &(as.integer(df$date - case_date) %% 7 == 0)
    df$Cases_week[w] <- -1
    w <- (df$date> res$last_obs) & (as.integer(df$date - death_date) %% 7 == 0)
    df$Deaths_week[w] <- -1
    if (is.element("admissions_week",colnames(df))){
        w <- (df$date> res$last_obs) & (as.integer(df$date - admissions_date) %% 7 == 0)
        df$admissions_week[w] <- -1
    }
    
    df
}


## add additional data for 3 weeks
get_newdata_past <- function(res){
    df <- res$fit$data
    case_date <- df$date[which(!is.na(df$Cases_week))[1]]
    death_date <- df$date[which(!is.na(df$Deaths_week))[1]]
    a <- as.character(df$ltla[1])
    df <- subset(df[df$ltla == a,],select=-group)
#    df <- df[df$ltla == a,c("ltla", "date", "Cases_week", "Deaths_week","averageMobility","averageRt","logitRt")]
    N <- 29
    dfnew <- data.frame(ltla=a, date=max(df$date) + 1:N,
                        Cases_week=NA, Deaths_week=NA)
    for (j in colnames(df)){
        if (!is.element(j,colnames(dfnew))) ## carry last value forward
            dfnew[,j] <- df[,j][df$date==max(df$date)]
    }
    df <- rbind(df, dfnew)
    Cases_start <- df$date[min(which(!is.na(df$Cases_week)))]
    w <- (df$date> res$last_obs|df$date<Cases_start) &(as.integer(df$date - case_date) %% 7 == 0)
    df$Cases_week[w] <- -1
    w <- (df$date> res$last_obs) & (as.integer(df$date - death_date) %% 7 == 0)
    df$Deaths_week[w] <- -1
    df
}


getcaseend <- function(data){
    nobs <- table(data$date)
    Caseend <-  max(as.Date(names(nobs)[nobs>=0.7*max(nobs)]))-2 ## remove last 2 days - to avoid problems through reporting delays - that usually is 3 days including the day on which data was downloaded

## remove one more date if there is a strong drop in cases
    if (sum(data$Cases[data$date==Caseend-1],na.rm=TRUE)*0.6>
        sum(data$Cases[data$date==Caseend],na.rm=TRUE)) 
        Caseend <- Caseend-1
    Caseend
}


### IFR IAR i2o_rates utility functions
add_i2o_rates <- function(data,i2o_rates,region="England"){
    i2o_rates <- i2o_rates %>% filter(region==!!region)
    IFR <- i2o_rates%>%filter(type=="IFR") %>% select(-type) %>%
        rename(ifr=value)
#    if (any(IFR$region!="England")) stop("Currently only England regions allowed")
    IFR <- IFR %>% select(-region)
    if (length(unique(IFR$ifr))==1) IFR$ifr <- jitter(IFR$ifr)
    data <- data %>% left_join(IFR,by=c("date")) %>% tidyr::fill(ifr,.direction="downup")
    
    
    IAR <- i2o_rates%>%filter(type=="IAR") %>% select(-type) %>%
        rename(iar=value)
 #   if (any(IAR$region!="England")) stop("Currently only England regions allowed")
    IAR <- IAR %>% select(-region)
    if (length(unique(IAR$iar))==1) IAR$iar <- jitter(IAR$iar)
    data <- data %>% left_join(IAR,by=c("date")) %>% tidyr::fill(iar,.direction="downup")
    
}

get_IFR_sd <- function(i2o_rates){
    IFR_sd <- i2o_rates%>%filter(type=="IFR_sd")
    if (dim(IFR_sd)[1]!=1) stop("expected exactly one IFR_sd")
    IFR_sd$value[1]
}

get_IAR_sd <- function(i2o_rates){
    IAR_sd <- i2o_rates%>%filter(type=="IAR_sd")
    if (dim(IAR_sd)[1]!=1) stop("expected exactly one IAR_sd")
    IAR_sd <- IAR_sd$value[1]
}

add_i2o_rates_option <- function(option_list,default="data/i2o_rates.rds"){
    append(option_list,make_option(c("--i2o_rates"),action="store",default=default,help="i2o rates file [default \"%default\"]"))

}


gatherresults <- function(j,areaname="ltla"){
    tryCatch({
        cat(j," ")
        load(paste(basedir,j,sep=""))
        case_date <- lastFullWeek 
        death_date <- lastFullWeek 

        newdata <- get_newdata(res,areaname=areaname)
        if (is.element("admissions_week",colnames(newdata)))
            admissions_date <- lastFullWeek
        
        a <- as.character(newdata[[areaname]][1])
#        df <- df[df$ltla == a,c("ltla", "date","averageMobility","averageRt")]
#        N <- 21
        ## dfnew <- data.frame(ltla=a,
        ##                     date=max(df$date) + 1:N,
        ##                     averageMobility=df$averageMobility[df$date==max(df$date)],
        ##                     averageRt=df$averageRt[df$date==max(df$date)])
        ## df <- rbind(df, dfnew)
        ## df$week <- ceiling((df$date - res$last_obs_date)/7)
        ## newdata <- df
        newdata[,c("Cases_week", "Deaths_week")] <- NA ##to align with observations
        w <- as.integer(newdata$date - case_date) %% 7 == 0
        newdata$Cases_week[w] <- -1
        w <- as.integer(newdata$date - death_date) %% 7 == 0
        newdata$Deaths_week[w] <- -1
        if (is.element("admissions_week",colnames(newdata))){
            newdata[,"admissions_week"] <- NA 
            w <- as.integer(newdata$date - admissions_date) %% 7 == 0
            newdata$admissions_week[w] <- -1
        }
        
        mat <- posterior_predict(res$fit, newdata=newdata, types="Cases_week")
        cat(a,"\n")
        dates <- mat$time
        mat <- t(mat$draws)

        mat_d<- posterior_predict(res$fit, newdata=newdata, types="Deaths_week")
        dates_d<- mat_d$time
        mat_d<- t(mat_d$draws)


        if (is.element("admissions_week",colnames(newdata))){
            mat_a <- posterior_predict(res$fit, newdata=newdata, types="admissions_week")
            dates_a <- mat_a$time
            mat_a <- t(mat_a$draws)
        }
        
        newdata[,c("Cases_week", "Deaths_week")] <- -1
        if (is.element("admissions_week",colnames(newdata))){
            newdata[,"admissions_week"] <- -1
        }
        rt_post <- posterior_rt(rm_re(res$fit), newdata=newdata)
        dates_rt <- rt_post$time
        rt <- t(rt_post$draws)

        infections <- posterior_infections(res$fit, newdata=newdata)
        dates_inf <- infections$time
        mat_inf <- t(infections$draws)
        
        resadd <- list()
        for (k in names(reportinfo)){
            if (reportinfo[[k]]$type=="infections"){
                w <- (dates_inf>=reportinfo[[k]]$dates[1]&
                      dates_inf<=reportinfo[[k]]$dates[2])
                future_cases <- colSums(mat_inf[w,,drop=F], na.rm=T)
                if (reportinfo[[k]]$"by100k"){
                    future_cases <- future_cases/popltla[a]*1e5
                }
                
                resadd[[k]] <- quantile(future_cases,
                                        c(0.5,(1-reportinfo[[k]]$CI)/2,
                                          1-(1-reportinfo[[k]]$CI)/2))
                names(resadd[[k]]) <- c("Median","CILow","CIUp")
                next
            }
            if (reportinfo[[k]]$type=="Pgx"&&reportinfo[[k]]$which=="Infections"){
                w <- (dates_inf>=reportinfo[[k]]$dates[1]&
                      dates_inf<=reportinfo[[k]]$dates[2])
                future_cases <- colSums(mat_inf[w,,drop=F], na.rm=T)
                if (reportinfo[[k]]$"by100k"){
                    future_cases <- future_cases/popltla[a]*1e5
                }
                
                resadd[[k]] <- mean(future_cases>=reportinfo[[k]]$x)
                names(resadd[[k]]) <- NULL
                next
            }
            if (reportinfo[[k]]$type=="obs"){
                df <- res$fit$data
                w <- (df$date>=reportinfo[[k]]$dates[1]&
                          df$date<=reportinfo[[k]]$dates[2])
                resadd[[k]] <-  sum(df[w, reportinfo[[k]]$which], na.rm=T)
                if (reportinfo[[k]]$"by100k"){
                    resadd[[k]] <- resadd[[k]]/popltla[a]*1e5
                }
                names(resadd[[k]]) <- NULL
            }
            if (reportinfo[[k]]$type=="pred"){
                if (reportinfo[[k]]$which=="Cases_week"){
                    w <- (dates>=reportinfo[[k]]$dates[1]&
                          dates<=reportinfo[[k]]$dates[2])
                    future_cases <- colSums(mat[w,,drop=F], na.rm=T)
                } else if (reportinfo[[k]]$which=="Deaths_week"){
                    w <- (dates_d>=reportinfo[[k]]$dates[1]&
                          dates_d<=reportinfo[[k]]$dates[2])
                    future_cases <- colSums(mat_d[w,,drop=F], na.rm=T)
                } else if (reportinfo[[k]]$which=="admissions_week"){
                    w <- (dates_a>=reportinfo[[k]]$dates[1]&
                          dates_a<=reportinfo[[k]]$dates[2])
                 
                    future_cases <- colSums(mat_a[w,,drop=F], na.rm=T)
                }  else
                    stop(paste("Type ", reportinfo[[k]]$type, "not allowed"))

                if (reportinfo[[k]]$"by100k"){
                    future_cases <- future_cases/popltla[a]*1e5
                }
                
                resadd[[k]] <- quantile(future_cases,
                                        c(0.5,(1-reportinfo[[k]]$CI)/2,
                                          1-(1-reportinfo[[k]]$CI)/2))
                names(resadd[[k]]) <- c("Median","CILow","CIUp")
                
            }
            if (reportinfo[[k]]$which=="Ravg") {
                w <- which(dates_rt>=reportinfo[[k]]$dates[1]&
                               dates_rt<=reportinfo[[k]]$dates[2])
                rt_atdate <- colMeans(rt[w,,drop=F])
                if (reportinfo[[k]]$type=="Pgx"){
                    resadd[[k]] <- mean(rt_atdate>reportinfo[[k]]$x)
                    names(resadd[[k]]) <- NULL
                    
                }else{
                    resadd[[k]] <- quantile(rt_atdate,
                                            c(0.5,(1-reportinfo[[k]]$CI)/2,
                                              1-(1-reportinfo[[k]]$CI)/2))
                    names(resadd[[k]]) <- c("Median","CILow","CIUp")
                }
            }else if (reportinfo[[k]]$type=="Pgx"){
                w <- (dates>=reportinfo[[k]]$dates[1]&
                          dates<=reportinfo[[k]]$dates[2])
                future_cases <- colSums(mat[w,,drop=F], na.rm=T)
                if (reportinfo[[k]]$"by100k"){
                    future_cases <- future_cases/popltla[a]*1e5
                }
                
                resadd[[k]] <- mean(future_cases>=reportinfo[[k]]$x)
                names(resadd[[k]]) <- NULL
                
            }
            if (reportinfo[[k]]$which=="R") {
                w <- which(dates_rt==reportinfo[[k]]$date)
                rt_atdate <- rt[w,,drop=F]
                resadd[[k]] <- quantile(rt_atdate,
                                        c(0.5,(1-reportinfo[[k]]$CI)/2,
                                          1-(1-reportinfo[[k]]$CI)/2))
                names(resadd[[k]]) <- c("Median","CILow","CIUp")
                
            }
        }
        #        if (is.na(startpred)){ ##determine start of prediction period
        #            startpred <- max(dates)-predlen+1
        #        }else{
        #            if (startpred!=max(dates)-predlen+1)
        #                stop("date period not consistent")
        #        }
        
        #        future_cases <- colSums(mat[dates >=startpred,])
        #        qtl <- quantile(future_cases, c(0.5,.05,0.95))
        #        qtl <- signif(qtl,3)
        list(resadd=resadd,a=a)
    },error=function(e) NA)
}
