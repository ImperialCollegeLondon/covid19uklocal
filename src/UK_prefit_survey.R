## for epidemia with the new observation format and infections_normal branch

## script is currently tested for England
library(optparse)
library(here)
library(epidemia)
library(dplyr)
library(zoo)
library(gridExtra)
library(tidyverse)
library(zoo)
options(mc.cores=parallel::detectCores())
source(here('UK_Linelist/ifriar.R'))

option_list <- list(
    make_option(c("--input"),action="store", default=paste(here(),"data/uk-public-ltla-combined.rds",sep=""),help="Input file [default \"%default\"]"),
    make_option(c("--input_survey"),action="store", default=paste(here(),"data/uk-public-survey.rds",sep=""),help="Input file [default \"%default\"]"),
    make_option(c("--nchains"),action="store", type="integer", default=4,help="Number of Chains [default \"%default\"]"),
    make_option(c("--iter"),action="store", type="integer", default=6000,help="Number of iterations [default \"%default\"]"),
    make_option(c("--region"),action="store", default="England",help="Region [default \"%default\"]"),
    make_option(c("--thin"),action="store", type="integer", default=2,help="Amount of thinning of results [default \"%default\"]"),
    make_option(c("--output"),action="store", default=here("data"),
                help="Output directory [default \"%default\"]"),
    make_option(c("--areaname"),action="store", default="Area_name",
                help="Column name for area in input file [default \"%default\"]")

)
opt <- parse_args(OptionParser(option_list=option_list))



dataorig <- readRDS(opt$input)
dataorig[[opt$areaname]] <- gsub("_"," ", dataorig[[opt$areaname]])
pops <- read.csv2(paste(here(),"data/modified_population.csv",sep="/"),sep=",") # doesn't work on hpc
pops <- data.frame(ltla=gsub("_"," ",as.character(pops$AREA)),pop=pops$Y2018,stringsAsFactors=FALSE)
pop <- pops$pop
names(pop) <- pops$ltla
dataorig$pop <- pop[dataorig[[opt$areaname]]]

data <- ungroup(dataorig)
data <- select(data, !!opt$areaname, Cdate, Cases, Deaths, pop) # , Tests
data <- rename(data, ltla=!!opt$areaname, date=Cdate, Cases=Cases, Deaths_week=Deaths)

source(here("UK_Linelist/public/utility.R"))
Caseend <- getcaseend(data)

cat("Case End=",as.character(Caseend),"\n")

data <- group_by(data, ltla)
data <- arrange(data, date, .by_group=T)

a <- gsub("_", " ",opt$region)

data <- data[is.element(data$ltla,a),]

data <- data[data$date<=Caseend,] ## cut off data at end of observation
## reading survey data
data_survey <- readRDS(opt$input_survey)
data_survey <- 
    data_survey %>%
    rename(date = Cdate, ltla = Region_name)
## reading ons, infections (right now this is the only info but can change)
data_ons <- 
    data_survey %>%
    filter(Type == 'Infections' & Survey == 'ONS') %>%
    select(ltla,date,Value, LowerInterval, UpperInterval) %>%
    rename(Infections=Value, ILowerInterval = LowerInterval, IUpperInterval = UpperInterval)

## joining data 
data <- left_join(data, data_ons)

## reading react data of total infections
data_react <- 
    data_survey %>%
    filter(Type == 'Total_Infections' & Survey == 'REACT') %>%
    select(ltla,date,Value, LowerInterval, UpperInterval) %>%
    rename(TotalInfections=Value, TILowerInterval = LowerInterval, TIUpperInterval = UpperInterval)
data <- left_join(data, data_react)
## weeks ends with the last case
data$week <- ceiling((data$date - Caseend)/7) ## to allign with Caseend
data$week[data$week==0] <- -1 ##to stabilise the end.

data$Cases_week <- NA

#data$Deaths_week <- NA

#data$Tests_week <- NA

for (ltla in a){
    w <- data$ltla==ltla
    
    o <- order(data$date[w]) ##ensure dates are ordered
    data[w,] <- data[w,][o,]
    
    ##compute cumulative obs per week
    data$Cases_week[w] <- sapply(1:sum(w),
                                 function(i)
                                     if (as.integer((data$date[w][i]-Caseend))%% 7==0){
                                         sum(data$Cases[w][i:(max(1,i-6))])   
                                     }else {
                                         NA
                                     }
    )
    
    
    
    start <- which(cumsum(ifelse(is.na(data$Deaths_week[w]),0,data$Deaths_week[w]))>=10)
    if (length(start)>=1){
        ostart <- min(data$date[w][start])
    }else{
        ostart <- as.Date("2030-01-01")
    }
    epistart <- ostart -30
    data$week[w] <- pmax(data$week[w],data$week[w][which(data$date[w]==ostart)-14])# start random walk a bit later than the epidemic
    
    Casestart <- as.Date("2020-06-01") ### cases start in June in the model
    data$Cases[w][data$date[w]<Casestart] <- NA
    data$Cases[w][data$date[w]>Caseend] <- NA
    
    data$Cases_week[w][data$date[w]<Casestart] <- NA
    data$Cases_week[w][data$date[w]>Caseend] <- NA
    data$Deaths_week[w][data$date[w]>Caseend] <- NA
    data$Deaths_week[w][data$date[w]<ostart] <- NA
    data$ltla[w][data$date[w]<epistart] <- NA
    data <- data[!is.na(data$ltla),]
}

obs <- list()
obs$Cases <- epiobs(formula=Cases(ltla,date) ~ 1,
                    i2o=c(0,0,0,rep(1/10,10))*.21)  ##iar .3 added to "i2o"

i2o2week <- function(i2o)
    rowSums(sapply(0:6, function (k) c(rep(0,k),i2o,rep(0,6-k))))

# code for deaths formula
deaths_last_date <- tail(data$date[which(!is.na(data$Deaths_week))], n=1)
ifr_dates = seq.Date(as.Date("2020-06-01"), deaths_last_date-1, by="2 week")
ifr_dates_formula <- lapply(ifr_dates, function(x) paste0("I(date>'", x, "')"))
ifr_dates_formula <- formula(paste0("Deaths_week(ltla, date) ~ 1 + ", 
                                    paste(ifr_dates_formula, collapse=" + ")))
## Poisson model
obs$Deaths_week <- epiobs(formula=ifr_dates_formula,
                          link="identity",
                          family="neg_binom",
                          prior = rstanarm::normal(0, .1),
                          prior_aux = rstanarm::normal(location=5,2),
                          prior_intercept=rstanarm::normal(1,.05),
                          i2o=i2o2week(EuropeCovid$obs$deaths$i2o*.009)) ##react ifr

## formula for cases
cases_last_date <- tail(data$date[which(!is.na(data$Cases_week))], n=1)
iar_dates = seq.Date(as.Date("2020-06-01"), cases_last_date-1, by="2 week")
iar_dates_formula <- lapply(iar_dates, function(x) paste0("I(date>'", x, "')"))
iar_dates_formula <- formula(paste0("Cases_week(ltla, date) ~ 1 + ", 
                                    paste(iar_dates_formula, collapse=" + ")))
obs$Cases_week <- epiobs(formula=iar_dates_formula,
                         link="identity",
                         family="neg_binom",
                         prior_aux = rstanarm::normal(location=5,2),
                         prior_intercept=rstanarm::normal(1,.25),
                         prior = rstanarm::normal(0, .3),
                         i2o=i2o2week(obs$Cases$i2o))
## ONS infections as observations
sd <- (mean(data$IUpperInterval - data$ILowerInterval, na.rm=TRUE)/4) * sqrt(14) 
obs$infections <- epiobs(formula = Infections(ltla, date) ~1,
                         i2o = 1, # direct linking of the infections
                         prior_intercept = rstanarm::normal(1, .1),
                         prior_aux = rstanarm::normal(sd, sd/5),
                         family = "normal",
                         link = "identity",)

## REACT total infections as observations
sd <- (mean(data$TIUpperInterval - data$TILowerInterval, na.rm=TRUE)/4)
obs$total_infections <- epiobs(formula = TotalInfections(ltla, date) ~1,
                         i2o = rep(1,nrow(data)), # direct linking of the infections
                         prior_intercept = rstanarm::normal(1, .1),
                         prior_aux = rstanarm::normal(sd, sd/5),
                         family = "normal",
                         link = "identity",)

args <- list()
args$data <- data
args$obs <- list(Cases_week=obs$Cases_week,Deaths_week=obs$Deaths_week, infections=obs$infections, total_infections=obs$total_infections)
#args$obs <- list(Deaths=obs$Deaths,Cases=obs$Cases)
args$rt <- epirt(formula=R(ltla,date) ~ rw(time = week, prior_scale = 0.2)+I(date>="2020-03-23"))

args$pops <- pops
args$si <- EuropeCovid$si
args$group_subset <- a #"Leicester"
args$algorithm <- "sampling"
args$sampling_args <- list(iter=opt$iter, chains=opt$nchains,thin=opt$thin,control=list(adapt_delta=0.93,max_treedepth=13))
args$data$ltla <- factor(args$data$ltla) ##temporary fix - country needs to be categorical

args$init_run <- list(iter=8000, chains=1)#, control=list(adapt_delta=0.95,max_treedepth=15))

args$alpha <- 1
args$beta <- 10
args$pop_adjust <- FALSE
#args$init_run <- TRUE

fit <- do.call(epim, args)
grid.arrange(plot_rt(fit),
             plot_infections(fit), 
             plot_obs(fit, type="Deaths_week"),
             plot_obs(fit, type="Cases_week"),
             plot_obs(fit, type="Infections"),
             nrow=5)

res <- list(fit=fit,
            model=fit$formula,
            last_obs_date=fit$data$date[max(which(!is.na(fit$data$Cases)))],
            today=max(fit$data$date),
            ltla=a,
            time=time)

rt <- posterior_rt(fit)
res$meanRt <- tibble(date=rt$time,averageRt=colMeans(rt$draws))

ifr <- posterior_ifriar(fit, "Deaths_week") * sum(fit$obs$Deaths_week$i2o) / 7 #Weekly observations
iar <- posterior_ifriar(fit, "Cases_week") * sum(fit$obs$Cases_week$i2o) / 7 #Weekly observations

# Use these objects for ifr and iar
ifr_quantiles <- get_ifriar_quantiles(fit, ifr, c(90))
iar_quantiles <- get_ifriar_quantiles(fit, iar, c(90))

start_date <- min(args$data$date)
first_ifr_date <- min(as.Date(ifr_quantiles$date))
first_iar_date <- min(as.Date(iar_quantiles$date))
end_date <- max(args$data$date)

ifr_all <- 
    ifr_quantiles %>%
    select(date, mean) %>%
    mutate(date = as.Date(date)) %>%
    complete(date = seq.Date(first_ifr_date, end_date, by="day"), fill = list(mean = NA)) %>%
    mutate(mean = na.locf(mean)) %>%
    complete(date = seq.Date(start_date, first_ifr_date, by="day"), fill = list(mean = NA)) %>%
    mutate(mean = na.locf(mean, fromLast = TRUE)) %>%
    mutate('region' = a, type = 'IFR') %>%
    rename(value = mean) %>%
    add_row(region=a,type="IFR_sd",date=NA,value= (mean(ifr_quantiles$upper - ifr_quantiles$lower)/4)/mean(ifr_quantiles$mean))
iar_all <- 
    iar_quantiles %>%
    select(date, mean) %>%
    mutate(date = as.Date(date)) %>%
    complete(date = seq.Date(first_iar_date, end_date, by="day"), fill = list(mean = NA, lower = NA, upper = NA, tag = 'CI', level = 90)) %>%
    mutate(mean = na.locf(mean)) %>%
    complete(date = seq.Date(start_date, first_iar_date, by="day"), fill = list(mean = NA)) %>%
    mutate(mean = na.locf(mean, fromLast = TRUE)) %>%
    mutate('region' = a, type = 'IAR') %>%
    rename(value = mean) %>%
    add_row(region=a,type="IAR_sd",date=NA,value= (mean(iar_quantiles$upper - iar_quantiles$lower)/4)/mean(iar_quantiles$mean))

res$i2o_rates <- bind_rows(ifr_all, iar_all)
# save result to file
save(res, file=paste(opt$output,"/fm-", gsub(" ","_",a), ".rdata",sep=""))
saveRDS(res$i2o_rates, file =paste(opt$output,"/i2o_rates.rds", sep=""))
