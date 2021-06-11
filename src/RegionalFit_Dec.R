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
    make_option(c("--input"),action="store", default=paste(here(),"/data/uk-public-ltla-combined.rds",sep=""),help="Input file [default \"%default\"]"),
    make_option(c("--nchains"),action="store", type="integer", default=4,help="Number of Chains [default \"%default\"]"),
    make_option(c("--iter"),action="store", type="integer", default=2000,help="Number of iterations [default \"%default\"]"),
    make_option(c("--region"),action="store", default="UK",help="Region [default \"%default\"]"),
    make_option(c("--thin"),action="store", type="integer", default=2,help="Amount of thinning of results [default \"%default\"]"),
    make_option(c("--output"),action="store", default=here("fits"),
                help="Output directory [default \"%default\"]")
)
option_list <- add_i2o_rates_option(option_list,here("data/i2o_rates.rds"))
opt <- parse_args(OptionParser(option_list=option_list))



dataorig <- readRDS(opt$input)
dataorig$Area_name <- gsub("_"," ", dataorig$Area_name)
pops <- read.csv2(paste(here(),"data/modified_population.csv",sep="/"),sep=",") # doesn't work on hpc
pops <- data.frame(ltla=gsub("_"," ",as.character(pops$AREA)),pop=pops$Y2018,stringsAsFactors=FALSE)
pop <- pops$pop
names(pop) <- pops$ltla
dataorig$pop <- pop[dataorig$Area_name]

data <- ungroup(dataorig)
data <- select(data, Area_name, Cdate, Cases, Deaths, pop) # , Tests
data <- rename(data, ltla=Area_name, date=Cdate, Cases=Cases, Deaths_week=Deaths)

Caseend <- getcaseend(data)

cat("Case End=",as.character(Caseend),"\n")

data <- group_by(data, ltla)
data <- arrange(data, date, .by_group=T)

a <- gsub("_", " ",opt$region)

data <- data[is.element(data$ltla,a),]

data <- data[data$date<=Caseend,] ## cut off data at end of observation

## weeks ends with the last case
data$week <- ceiling((data$date - (Caseend-3))/7) ## to allign with Caseend
data$week[data$week>0] <- 0 ##to stabilise the end


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
    ostart <- as.Date("2020-12-15")
    epistart <- ostart -30
    data$week[w] <- pmax(data$week[w],data$week[w][which(data$date[w]==ostart)-14])# start random walk a bit later than the epidemic
    
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




data <- 
    data %>% 
    mutate( changePoint = as.numeric(format(date, '%m'))-8) %>% 
    mutate( changePoint = if_else(changePoint <0, 0, 1))



#### Parsing i2o_rates into data
i2o_rates <- readRDS(opt$i2o_rates)
data <- data %>% add_i2o_rates(i2o_rates)
IFR_sd <- get_IFR_sd(i2o_rates)
IAR_sd <- get_IAR_sd(i2o_rates)


obs <- list()
## obs$Deaths <- epiobs(formula=Deaths(ltla,date) ~ 1,
##                     i2o=EuropeCovid$obs$deaths$i2o/100) ##ifr .01 added to i2o
obs$Cases <- epiobs(formula=Cases(ltla,date) ~ 1,
                    i2o=c(0,0,0,rep(1/10,10)))  ##iar .3 added to "i2o"

i2o2week <- function(i2o)
    rowSums(sapply(0:6, function (k) c(rep(0,k),i2o,rep(0,6-k))))


## Poisson model
obs$Deaths_week <- epiobs(formula=Deaths_week(ltla,date) ~ 0+ifr,# + factor(changePoint),
                          link="identity",
                          # link="logit",
                          #                          family="poisson",
                          family="neg_binom",
                          prior_aux = rstanarm::normal(location=5,2),
                          prior=rstanarm::normal(1,IFR_sd,autoscale=FALSE),
                          i2o=i2o2week(EuropeCovid$obs$deaths$i2o)) ##ifr .01 added to i2o
# i2o=i2o2week(EuropeCovid$obs$deaths$i2o*.0136*2))
# obs$Cases_week <- epiobs(formula=Cases_week(ltla,date) ~ 1 + I(Tests_week/pop * 1e3),
# obs$Cases_week <- epiobs(formula=Cases_week(ltla,date) ~ 1 + I(Cases_week/Tests_week * 100),
obs$Cases_week <- epiobs(formula=Cases_week(ltla,date) ~ 0+iar , #+ I(Tests_week/pop * 1e3),
                         link="identity",
                         #                         family="poisson",
                         family="neg_binom",
                         prior_aux = rstanarm::normal(location=5,2),
                          prior=rstanarm::normal(1,IAR_sd,autoscale=FALSE),
                         i2o=i2o2week(obs$Cases$i2o))  ##iar .19 added to "i2o"



args <- list()
args$data <- data
args$obs <- list(Cases_week=obs$Cases_week,Deaths_week=obs$Deaths_week)
#args$obs <- list(Deaths=obs$Deaths,Cases=obs$Cases)
args$rt <- epirt(formula=R(ltla,date) ~ rw(time = week, prior_scale = 0.2)+I(date>="2020-03-23"))

args$pops <- pops
args$si <- EuropeCovid$si
args$group_subset <- a #"Leicester"
args$algorithm <- "sampling"
args$sampling_args <- list(iter=opt$iter, chains=opt$nchains,thin=opt$thin,control=list(adapt_delta=0.96,max_treedepth=12))
args$data$ltla <- factor(args$data$ltla) ##temporary fix - country needs to be categorical

args$init_run <- list(iter=2000, chains=1)#, control=list(adapt_delta=0.95,max_treedepth=15))

args$alpha <- 1
args$beta <- 10
args$pop_adjust <- TRUE

for (j in 1:10){
    
    time <- system.time({fit <- do.call("epim", args)})
    args$sampling_args$iter <- min(args$sampling_args$iter*2,8000)
                               
    sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
    Rhat <- max(rstan::summary(fit$stanfit)$summary[,"Rhat"])
    divergent <- sum(sapply(sampler_params, function(x) x[,"divergent__"]))
    cat("Rhat=",Rhat," divergent steps=", divergent,"\n")
    if (Rhat<1.2&&divergent<40) break;
}


warnings()

donotrun <- function(){
    print(fit,2)
    library(gridExtra)
    
    grid.arrange(plot_rt(fit),plot_obs(fit,"Cases_week"),plot_obs(fit,"Deaths_week"),ncol=1)
}

res <- list(fit=fit,
            model=fit$formula,
            last_obs_date=fit$data$date[max(which(!is.na(fit$data$Cases)))],
            today=max(fit$data$date),
            ltla=a,
            time=time)

rt <- posterior_rt(fit)
res$meanRt <- tibble(date=rt$time,averageRt=colMeans(rt$draws))

save(res,file=paste(opt$output,"/fm-",opt$region,".rdata",sep=""))

