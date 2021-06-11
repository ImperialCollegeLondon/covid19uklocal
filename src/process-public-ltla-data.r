library(tidyr)
library(dplyr)
library(readr)
library(janitor)
library(lubridate)
library(magrittr)
library(ggpubr)
library(tidyr)
library(readxl)
library(tools)
library(ggplot2)
library(flextable)
library(stringr)
library(utils)
library(httr)
library(here)
library(rvest)
library(aweek)

# scraplinks code taken from https://gist.github.com/paulrougieux/e1ee769577b40cd9ed9db7f75e9a2cc2
scraplinks <- function(url){
  # Create an html document from the url
  webpage <- xml2::read_html(url)
  # Extract the URLs
  url_ <- webpage %>%
    rvest::html_nodes("a") %>%
    rvest::html_attr("href")
  # Extract the link text
  link_ <- webpage %>%
    rvest::html_nodes("a") %>%
    rvest::html_text()
  return(tibble(link = link_, url = url_))
}

get_data_web <- function(url, url_page, file_to_disk, retry_link = NULL){
  tryCatch({
    r <- RETRY("GET", url, 
               write_disk(file_to_disk, overwrite=TRUE))
    
    if (http_error(r)) {
      
      if (is.null(retry_link)){
        stop("Error downloading file")
      } else{
        get_data_web(retry_link, url_page, file_to_disk)
      }
    }
  },
  error = function(e) {
    stop(sprintf("Error downloading file '%s': %s, please check %s",
                 url, e$message, url_page))
  })
}

get_regional_data <- function(df_tibble = NULL, region = NULL){
  if (is.null(df_tibble) || is.null(region)){
    print('Need a dataframe and region to subset')
  } else{
    df_subset <-  
      df_tibble %>%
      filter(Region_name == region) %>%
      group_by(Cdate) %>%   
      summarise_at(vars(Cases,Deaths), sum) %>% 
      mutate(Area_name = region, Region_name = 'UK') %>%
      ungroup()
  }
}

process_survey_data <- function(){
  # population of England by age bands we need
  pop_england <- 
    tibble('Region_name'='England', '010'=6264662, '1116'=3896599, '1724'=5346009, '2534'=7609363, '3549'=10863751, '5069'=13486687, '70above'=7556976) %>% 
    pivot_longer(-Region_name, names_to = 'subgroup', values_to = 'Population')
  # using web scraping to get latest link everytime, instead of relaying on date formats
  base_url_ons <- "https://www.ons.gov.uk"
  url_page_ons <- "https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata"
    ons_links <- scraplinks(url_page_ons)
    url_ons <- paste0(base_url_ons,filter(ons_links, grepl('xlsx',link))$url)
    file_ons <- here('data/ons_survey.xlsx')
    if (length(url_ons)>1){
        tmp <- lapply(url_ons,function(u){
            get_data_web(url = u,
                         url_page = url_page_ons,
                         file_to_disk = file_ons
                         )
                                        # reading data and processing it to get what we want
            df_ons <- 
                readxl::read_excel(file_ons, sheet = '2b', skip = 9, guess_max = 20) %>%
                select(c(1,6:8)) %>%
                filter_all(any_vars(!is.na(.)))
            names(df_ons) <- c("Cdate", "Value", "LowerInterval", "UpperInterval")
            df_ons}
            )
        df_ons <- do.call(bind_rows,tmp)         
    }else{
        get_data_web(url = url_ons, 
                     url_page = url_page_ons,
                     file_to_disk = file_ons
                     )
  # reading data and processing it to get what we want
        df_ons <- 
            readxl::read_excel(file_ons, sheet = '2b', skip = 6, guess_max = 20) %>%
            select(c(1,6:8)) %>%
            filter_all(any_vars(!is.na(.)))
        names(df_ons) <- c("Cdate", "Value", "LowerInterval", "UpperInterval")
    }
  df_ons <- 
    df_ons %>%
    mutate(Value = as.double(Value),
           LowerInterval = as.double(LowerInterval),
           UpperInterval = as.double(UpperInterval)
    ) %>%
    mutate(Cdate=as.Date(Cdate)) %>%
    filter(Cdate >= as.Date('2020-08-14')) %>%
    mutate(Type='Infections', 
           Survey = 'ONS', 
           Region_name = 'England',
           Bi_week = as.integer((Cdate-min(Cdate)))%% 14) %>%
    filter(Bi_week == 0) %>%
    select(-Bi_week)
  df_react <- 
    tibble('Cdate' = as.Date('2020-06-30'), 
           "Value" = 3360000, 
           "LowerInterval" = 3210000, 
           "UpperInterval" = 3510000, 
           Type='Total_Infections', 
           Survey = 'REACT', 
           Region_name = 'England') %>%
    bind_rows(tibble('Cdate' = as.Date('2020-09-26'), 
              "Value" = 45000, 
              "LowerInterval" = 41000, 
              "UpperInterval" = 53000, 
              Type='Infections', 
              Survey = 'REACT', 
              Region_name = 'England'),
              tibble('Cdate' = as.Date('2020-10-21'), 
                     "Value" = 96000, 
                     "LowerInterval" = 86000, 
                     "UpperInterval" = 105000, 
                     Type='Infections', 
                     Survey = 'REACT', 
                     Region_name = 'England'),
              tibble('Cdate' = as.Date('2020-11-18'), 
                     "Value" = 72000, 
                     "LowerInterval" = 58000, 
                     "UpperInterval" = 78000, 
                     Type='Infections', 
                     Survey = 'REACT', 
                     Region_name = 'England')
              )
  df_survey <- full_join(df_ons, df_react)
  saveRDS(df_survey, file = paste(here(),'data/uk-public-survey.rds', sep= '/'))
  df_ons_Over25 <-
    readxl::read_excel(file_ons, sheet = '1g', skip = 6, guess_max = 20)
  names(df_ons_Over25) <- c('Cdate', "Positivity_010_mean", "Positivity_010_lower", "Positivity_010_upper", 
                            "Positivity_1116_mean", "Positivity_1116_lower", "Positivity_1116_upper", 
                            "Positivity_1724_mean", "Positivity_1724_lower", "Positivity_1724_upper",
                            "Positivity_2534_mean", "Positivity_2534_lower", "Positivity_2534_upper", 
                            "Positivity_3549_mean", "Positivity_3549_lower", "Positivity_3549_upper", 
                            "Positivity_5069_mean", "Positivity_5069_lower", "Positivity_5069_upper",
                            "Positivity_70above_mean", "Positivity_70above_lower", "Positivity_70above_upper")
  df_ons_age <- 
    df_ons_Over25 %>%
    filter_all(any_vars(!is.na(.))) %>%
    pivot_longer( cols = Positivity_010_mean:Positivity_70above_upper,
                  names_to = c("subgroup", "type"),
                  names_pattern = "Positivity_?(.*)_(.*)",
                  values_to = "Value") %>%
    pivot_wider(names_from = type, values_from = Value) %>%
    left_join(pop_england) %>%
    mutate(mean = mean * Population, lower = lower * Population, upper = upper * Population, Cdate = as.Date(Cdate), 
           Bi_week = as.integer((Cdate-min(Cdate)))%% 14) %>%
    filter(Bi_week == 0 | (Cdate %in% df_survey$Cdate[df_survey$Survey=='ONS'])) %>%
    select(-c(lower, upper, Population, Bi_week)) %>%
    mutate(Cdate = case_when(Cdate==as.Date('2020-09-12') ~ as.Date('2020-09-06'), TRUE ~ Cdate)) %>% #remove this temp fix
    filter((Cdate %in% df_survey$Cdate[df_survey$Survey=='ONS']))
  
  # create aggregated data
  df_ons_age_agg <-
    df_ons_age %>%
    group_by(Cdate, Region_name) %>%
    summarise_at(vars(mean), sum) %>%
    rename(total = mean)
  # joining last two frames
  
  df_ons_Over25 <-
    df_ons_age %>%
    left_join(df_ons_age_agg) %>%
    mutate(mean = mean / total) %>%
    filter(!subgroup %in% c('010', '1116', '1724')) %>%
    group_by(Cdate, Region_name) %>%
    summarise_at(vars(mean), sum) %>%
    rename(factor = mean) %>%
    mutate(subgroup = 'Over25', Survey = 'ONS') %>%
    bind_rows(tibble('Cdate' = unique(df_ons_age$Cdate), 'factor' = 1, 'Region_name' = 'England', 'subgroup' = 'All', Survey = 'ONS')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-06-30'), 'factor' = 0.7, 'Region_name' = 'England', 'subgroup' = 'Over25', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-06-30'), 'factor' = 1, 'Region_name' = 'England', 'subgroup' = 'All', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-09-26'), 'factor' = 0.7, 'Region_name' = 'England', 'subgroup' = 'Over25', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-09-26'), 'factor' = 1, 'Region_name' = 'England', 'subgroup' = 'All', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-10-21'), 'factor' = 0.7, 'Region_name' = 'England', 'subgroup' = 'Over25', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-10-21'), 'factor' = 1, 'Region_name' = 'England', 'subgroup' = 'All', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-11-18'), 'factor' = 0.7, 'Region_name' = 'England', 'subgroup' = 'Over25', Survey = 'REACT')) %>%
    bind_rows(tibble('Cdate' = as.Date('2020-11-18'), 'factor' = 1, 'Region_name' = 'England', 'subgroup' = 'All', Survey = 'REACT'))
  saveRDS(df_ons_Over25, file = paste(here(),'data/uk-public-survey-factors.rds', sep= '/'))
  file.remove(file_ons)
}
process_ltla_public_data <- function(){
  # getting data for local areas
  local_areas_data <- read_csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&metric=newDeaths28DaysByDeathDate&format=csv")
  # getting data for regions
  regions_data <- read_csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newCasesBySpecimenDate&metric=newDeaths28DaysByDeathDate&format=csv")
  # getting nations data
  nations_data <- read_csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&metric=newCasesBySpecimenDate&metric=newDeaths28DaysByDeathDate&format=csv")
  # getting UK data
  uk_data <- read_csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesBySpecimenDate&metric=newDeaths28DaysByDeathDate&format=csv")
  
  # changing name Comhairle nan Eilean Siar to match earlier data for local areas
  local_areas_data <- mutate(local_areas_data, areaName = case_when(areaName=="Comhairle nan Eilean Siar" ~ "Na h-Eileanan Siar", TRUE ~ areaName))
  # Buckinghamshire doesn't exist in data (adding it, it's a mystery area it gets created and removed so many times)
  buck_areas <- c("South Bucks", "Chiltern", "Wycombe", "Aylesbury Vale")
  buck_data <- 
    filter(local_areas_data, areaName%in%buck_areas) %>%
    group_by(date) %>%
    summarise_at(vars(newCasesBySpecimenDate,newDeaths28DaysByDeathDate), sum, na.rm=TRUE) %>%
    ungroup() %>%
    mutate(areaCode = "EN0700000", areaType="ltla", areaName = "Buckinghamshire")
  # changing name yorkshire and the humber to match earlier data for regions
  regions_data <- mutate(regions_data, areaName = case_when(areaName=="Yorkshire and The Humber" ~ "Yorkshire and Humber", TRUE ~ areaName))
  # changing name in United Kingdom to UK in overall data
  uk_data <- mutate(uk_data, areaName = "UK")
  
  # stacking all data above each other and creating weekly deaths as thats what we had in previous version
  df_data <- 
    bind_rows(local_areas_data, regions_data, nations_data,uk_data,buck_data) %>%
    rename(Cdate = date, Area_name = areaName, Cases = newCasesBySpecimenDate, Deaths = newDeaths28DaysByDeathDate) %>%
    select(Cdate, Area_name, Cases, Deaths) %>%
    inner_join(read_csv((here("data/la_to_regions_website.csv")))) %>%
    mutate(week = isoweek(Cdate),
           year = case_when((year(Cdate)==2021  & isoweek(Cdate) == 53) ~ 2020, TRUE ~  year(Cdate))) %>%
    group_by(Area_name, Region_name, week, year) %>%
    arrange(Area_name, Cdate) %>%
    mutate(Deaths = c(rep(NA, length(Deaths)-1), 1) * cumsum(Deaths)) %>%
    ungroup() %>%
    select(-c(week, year))
  # somehow local areas of Wales do not have any death data on UK site so for now I am just taking old data from ONS
  deaths_file = paste(here(),"data/ltla_weekly_deaths.xlsx", sep='/')
  base_url_ons <- "https://www.ons.gov.uk"
  url_page_deaths <- "https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/causesofdeath/datasets/deathregistrationsandoccurrencesbylocalauthorityandhealthboard"
  ons_deaths_links <- scraplinks(url_page_deaths)
  url_deaths<- paste0(base_url_ons,filter(ons_deaths_links, grepl('xlsx',link))$url)
  if (length(url_deaths)>1){
    get_data_web(url = url_deaths[grep("2020/",url_deaths)], url_page = url_page_deaths, file_to_disk = deaths_file)
    df_ons <- read_excel(deaths_file, sheet = 4, skip = 3)%>%mutate(year=2020)
    get_data_web(url = url_deaths[grep("2021/",url_deaths)], url_page = url_page_deaths, file_to_disk = deaths_file)
    df_ons <- df_ons %>% bind_rows(read_excel(deaths_file, sheet = 4, skip = 3)%>%mutate(year=2021))
  }else{
    get_data_web(url = url_deaths, url_page = url_page_deaths, file_to_disk = deaths_file)
    df_ons <- read_excel(deaths_file, sheet = 4, skip = 3)%>%mutate(year=2020)
  }
  file.remove(deaths_file)
  # get all local areas in wales
  wales_data <- 
    df_data %>%
    filter(Region_name=="Wales") %>%
    select(-Deaths)
  df_wales_deaths <- 
    df_ons %>%
    filter(`Geography type`=='Local Authority' & `Cause of death`=='COVID 19') %>%
    rename(Area_name = `Area name`, Week_number = `Week number`, Deaths = `Number of deaths`) %>%
    filter(Area_name %in% unique(wales_data$Area_name)) %>%
    group_by(Area_name, year, Week_number) %>%
    summarise_at(vars(Deaths), sum) %>%
    ungroup() %>%
    mutate(Cdate = get_date(week=Week_number, year = year, day = 7)) %>%
    select(Area_name,Cdate, Deaths)
  
  # joining cases and deaths for wales
  wales_data <- left_join(wales_data,df_wales_deaths)
  # first removing wales from original and then binding rows
  df_final <- 
    df_data %>%
    filter(!Region_name=="Wales") %>%
    bind_rows(wales_data)
  saveRDS(df_final, here("data/uk-public-ltla-combined.rds"))
}
process_ltla_public_data()
#process_survey_data()
