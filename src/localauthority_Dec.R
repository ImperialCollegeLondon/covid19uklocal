## for epidemia with the new observation format

## script runs for a single ltla
library(optparse)
library(here)
library(epidemia)
library(dplyr)
library(zoo)
options(mc.cores=parallel::detectCores())

source(here("UK_Linelist/public/utility.R"))


option_list <- list(
    make_option(c("--ltla"),action="store", type="integer", default=-1,help="Which LTLA to run this for [default \"%default\"]"),
    make_option(c("--ltlaname"),action="store", default="Leicester",help="Which LTLA to run this for [default \"%default\"]"),
    make_option(c("--input"),action="store", default=paste(here(),"/data/uk-public-ltla-combined.rds",sep=""),help="Input file [default \"%default\"]"),
    make_option(c("--nchains"),action="store", type="integer", default=4,help="Number of Chains [default \"%default\"]"),
    make_option(c("--iter"),action="store", type="integer", default=1000,help="Number of iterations [default \"%default\"]"),
    make_option(c("--thin"),action="store", type="integer", default=1,help="Amount of thinning of results [default \"%default\"]"),
    make_option(c("--UKfits"),action="store", default=here("fits"),  help="fits for UK model [default \"%default\"]"),
    make_option(c("--fit_params"),action="store", default=here("fits/fit_params.csv"),  help="population file [default \"%default\"]"),
    make_option(c("--output"),action="store", default=here("fits"),
                help="Output directory [default \"%default\"]")
)
option_list <- add_i2o_rates_option(option_list)
opt <- parse_args(OptionParser(option_list=option_list))

dataorig <- readRDS(opt$input)

pops <- read.csv2(here("data/modified_population.csv"),sep=",") 
pops <- data.frame(ltla=as.character(pops$AREA),pop=pops$Y2018,stringsAsFactors=FALSE)
pop <- pops$pop
names(pop) <- pops$ltla
dataorig$pop <- pop[dataorig$Area_name]

data <- ungroup(dataorig)
data <- select(data, Area_name, Cdate, Cases, Deaths, pop)
data <- rename(data, ltla=Area_name, date=Cdate, Cases=Cases, Deaths_week=Deaths)

Caseend <- getcaseend(data)

cat("Case End=",as.character(Caseend),"\n")


data <- group_by(data, ltla)
data <- arrange(data, date, .by_group=T)
if (opt$ltla==-1){
    a <- opt$ltlaname
    a <- gsub("_"," ",a)
}else{
    i <- opt$ltla
    a <- unique(data$ltla)[i]
}

data <- data[is.element(data$ltla,a),]

Caseend <- min(max(data$date)-2,Caseend)

data <- data[data$date<=Caseend,] 

## weeks ends with the last case
data$week <- ceiling((data$date - (Caseend-3))/7) ## to allign with Caseend
data$week[data$week>0] <- 0 ##to stabilise the end

data$Cases_week <- NA

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
        ostart <- min(min(data$date[w][start]),as.Date("2020-06-01"))
    }else{
        ostart <- as.Date("2020-06-01")
    }
    ostart <- as.Date("2020-12-15")
    cat("Observation start=",as.character(ostart),"\n")
    epistart <- ostart -30
    cat("Start of Epidemic=",as.character(epistart),"\n")

    Casestart <- as.Date("2020-12-05") ### cases start in June in the model
    data$Cases[w][data$date[w]<Casestart] <- NA
    data$Cases[w][data$date[w]>Caseend] <- NA

    data$Cases_week[w][data$date[w]<Casestart] <- NA
    data$Cases_week[w][data$date[w]>Caseend] <- NA
    data$Deaths_week[w][data$date[w]>Caseend] <- NA
    data$Deaths_week[w][data$date[w]<ostart] <- NA

    data$ltla[w][data$date[w]<epistart] <- NA
    data <- data[!is.na(data$ltla),]
}


## Load Regional fit
region <- dataorig$Region_name[min(which(dataorig$Area_name==a))]
load(paste(opt$UKfits,"/fm-",gsub(" ","_",region),".rds",sep=""))
data <- left_join(data, res$meanRt)
firstentryaRt <- min(which(!is.na(data$averageRt)))
data$averageRt[1:firstentryaRt] <- data$averageRt[firstentryaRt] ##fill NAs at beginning with first value
wNA <- which(is.na(data$averageRt))
wNA <- wNA[wNA<10]
data$averageRt[wNA] <- data$averageRt[max(wNA)+1]
wNA <- which(is.na(data$averageRt))
wNA <- wNA[wNA>10]
if (length(wNA)>0){
    data$averageRt[wNA] <- data$averageRt[min(wNA)-1]
}
data$logitRt <- log(data$averageRt/3.28/2/(1-data$averageRt/3.28/2))
data$averageRt <- data$averageRt-3.28
##stop regional trend a month before the end of observation to get regional trend.
data$logitRt[data$date>(Caseend-45)] <- data$logitRt[data$date==Caseend-45]

#### Parsing i2o_rates into data
i2o_rates <- readRDS(opt$i2o_rates)
data <- data %>% add_i2o_rates(i2o_rates)
IFR_sd <- get_IFR_sd(i2o_rates)
IAR_sd <- get_IAR_sd(i2o_rates)

obs <- list()
obs$Cases <- epiobs(formula=Cases(ltla,date) ~ 1,
                    i2o=c(0,0,0,rep(1/10,10)))

i2o2week <- function(i2o)
    rowSums(sapply(0:6, function (k) c(rep(0,k),i2o,rep(0,6-k))))


obs$Deaths_week <- epiobs(formula=Deaths_week(ltla,date) ~ 0+ifr,
                          link="identity",
                          family="quasi_poisson",
                          prior_aux = rstanarm::normal(location=3,2),
                          prior=rstanarm::normal(1,IFR_sd,autoscale=FALSE),
                          i2o=i2o2week(EuropeCovid$obs$deaths$i2o)) 
obs$Cases_week <- epiobs(formula=Cases_week(ltla,date) ~ 0 + iar,
                        link="identity",
                         family="quasi_poisson",
                            prior_aux = rstanarm::normal(location=3,2),
                          prior=rstanarm::normal(1,IAR_sd,autoscale=FALSE),
                         i2o=i2o2week(obs$Cases$i2o))  

obs$logit_Deaths_week <- epiobs(formula=Deaths_week(ltla,date) ~ 0 + I(log(ifr/(1-ifr))),
                                link="logit",
                                family="quasi_poisson",
                                prior_aux = rstanarm::normal(location=3,2),
                                prior=rstanarm::normal(1,IAR_sd,autoscale=FALSE),
                                i2o=2*i2o2week(EuropeCovid$obs$deaths$i2o)) 
obs$logit_Cases_week <- epiobs(formula=Cases_week(ltla,date) ~ 0 + I(log(iar/(1-iar))),
                               link="logit",
                               family="quasi_poisson",
                               prior_aux = rstanarm::normal(location=3,2),
                               prior=rstanarm::normal(1,IFR_sd,autoscale=FALSE),
                               i2o=2*i2o2week(obs$Cases$i2o))  


args <- list()
args$data <- data
args$obs <- list(Cases_week=obs$Cases_week,Deaths_week=obs$Deaths_week)
args$prior_covariance <- rstanarm::decov(shape=2,scale=0.15)
args$rt <- epirt(formula=R(ltla,date) ~  rw(time=week,prior_scale=.15)+logitRt,prior=rstanarm::normal(1,.05),prior_intercep=rstanarm::normal(0,.1))

args$pops <- pops
args$si <- c(EuropeCovid$si[1:33],sum(EuropeCovid$si[34:length(EuropeCovid$si)]))
args$group_subset <- a 
args$algorithm <- "sampling"
adapt_delta <- 0.92
args$alpha <- 1
args$beta <- 10

args$data$ltla <- factor(args$data$ltla) ##temporary fix - country needs to be categorical

args$init_run <- list(iter=1000, chains=1)

fitparams <- read.csv(opt$fit_params,sep="&")
fitparams$Area_name <- gsub("_"," ",fitparams$Area_name)
fitparams <- fitparams %>% filter(Area_name==gsub("_"," ",a))
if (dim(fitparams)[1]>0){
    opt$iter <- max(fitparams$iter)
    opt$thin <- max(fitparams$thin)
    adapt_delta <- max(fitparams$adapt_delta)
}

args$sampling_args <- list(iter=opt$iter, chains=opt$nchains,thin=opt$thin,control=list(adapt_delta=adapt_delta,max_treedepth=11))


cat("iter=",opt$iter, "thin= ", opt$thin, " adapt_delta=",adapt_delta,"\n")
time <- system.time({fit <- do.call("epim", args)})

sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
Rhat <- max(rstan::summary(fit$stanfit)$summary[,"Rhat"])
divergent <- sum(sapply(sampler_params, function(x) x[,"divergent__"]))
cat("Rhat=",Rhat," divergent steps=", divergent,"\n")


res <- list(fit=fit,
            model=fit$formula,
            last_obs_date=fit$data$date[max(which(!is.na(fit$data$Cases)))],
            today=max(fit$data$date),
            ltla=a,
            time=time)


if (Rhat>=1.2||divergent>=opt$nchains*opt$iter/2.*.005){
    if (divergent>=opt$nchains*opt$iter/2.*.005) {
        adapt_delta <- (adapt_delta+1.)/2.
    } else if (Rhat>=1.2 && opt$iter<6000){
        opt$iter <- opt$iter*2
        opt$thin <- opt$thin*2
    }
    write.table(data.frame(Area_name=a,
                           iter=opt$iter,
                           thin=opt$thin,
                           adapt_delta=adapt_delta),
                sep="&",append=TRUE,col.names=FALSE,row.names=FALSE,
                file=opt$fit_params)
    dir.create(file.path(opt$output,"failedruns"), showWarnings = FALSE)
    
# save result to file
    save(res, file=paste(opt$output,"/failedruns/fm-", gsub(" ","_",a), Sys.time(), ".rds",sep=""))

    stop("Sampling not successful; parameters adjusted")
}


warnings()



# save result to file
save(res, file=paste(opt$output,"/fm-", gsub(" ","_",a), ".rdata",sep=""))

