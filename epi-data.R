###
###   INGESTION OF CLINICAL DATA
###   (not all data sets can be shared, hence this code will not run)
### 


library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr) 
library(patchwork)


source('hosp-ped.R')
source('Rt.R')


plot_flu_history <- function(cases) {
  g = cases %>% 
    filter(type=='pos', virus %in% c('H1N1','H3N2', 'B')) %>%
    group_by(virus, date) %>% 
    summarise(y = sum(v, na.rm=TRUE), .groups='drop') %>% 
    ggplot(aes(x=date, y=y, color=virus)) + 
    geom_line(size=1) + 
    theme(panel.grid.minor = element_blank())+
    scale_x_date(date_breaks = '1 year', date_labels = '%Y')+
    scale_color_brewer(palette = 'Set1')+
    labs(title='influenza in Canada', x='', y='num. postive tests')
  g
}

#' Influenza cases by subtype
get_cases <- function(subtype) {

  message('Retrieving cases...', appendLF = F)
    
  cases = readRDS('data/RVDSS-influenza-2022-02-25.rds') %>%
    mutate(yss = as.numeric(substr(season,1,4))) %>% 
    mutate(virus = case_when(
      (name=='PosFLUA2') ~ 'H3N2',
      (name=='PosFLUA3') ~ 'unknown',
      (name=='PosFLUAP') ~ 'H1N1',
      (name=='PosFLUB')  ~ 'B',
      # add a "N/A" label for tests
      (name=='TestFLUAB')  ~ 'N/A',
    )) %>% 
    mutate(type = ifelse(grepl('^Pos', name),'pos','test'))
  
  # Just an opportunity to plot the time series...
  pdf('plots/plot-flu-history.pdf', width=12)
  plot_flu_history(cases)
  dev.off()
  
  # Aggregate at the National level
  cases.can = cases  %>%
    group_by(yss, date, type, virus) %>% 
    summarise(n = sum(v), .groups='drop')
  
  
  # Subset of tests and positive tests of subtype of interest
  # Drop "virus" column
  # Transform to wide format of aggregated cases
  cases.can.wide = cases.can %>%
    filter(virus %in% c(subtype, "N/A")) %>%
    select(yss, date, type, n) %>%
    pivot_wider(names_from = type, values_from = n) %>%
    mutate(positivity = pos/test)
  
  message('done.')
  return(list(
    longformat = cases.can,
    wideformat = cases.can.wide
  ))
}


get_vax_cvg <- function(variables) {
  
  message('Retrieving vaccine coverage...', appendLF = F)
  
  # Import and format vaccine coverage datasets
  
  sivc_coverage <- read_excel("data/sivc_vaccine_coverage_2015-2021.xlsx") %>%
    pivot_longer(!age.group, names_to = "season", values_to = "percent") %>%
    mutate(yss = as.numeric(substring(season, 8, 11)),
           date = ymd(paste0(yss+1,'-02-01'))) %>%
    rename(vax.cvg.sivc = percent) %>%
    filter(grepl("^All adults", age.group)) %>%
    select(yss, date, vax.cvg.sivc)
  
  statcan_coverage <- read.csv("data/statcan-cchs_vaccine-coverage_2015-2020.csv") %>%
    rename(Year = REF_DATE) %>%
    filter(GEO == "Canada (excluding territories)",
           Age.group == "Total, 12 years and over",
           Sex == "Both sexes",
           Indicators == "Influenza immunization in the past 12 months",
           UOM == "Percent") %>%
    mutate(Characteristics = case_when(
      (Characteristics == "Low 95% confidence interval, percent") ~ "Low.95.CI",
      (Characteristics == "High 95% confidence interval, percent") ~ "High.95.CI",
      (Characteristics == "Percent") ~ "Coverage")) %>%
    select(Year, Characteristics, VALUE) %>%
    pivot_wider(names_from = Characteristics, values_from = VALUE) %>%
    mutate(Date = ymd(paste0(Year,'-02-01')),
           vax.cvg.sc    = Coverage/100,
           vax.cvg.sc.lo = Low.95.CI/100,
           vax.cvg.sc.hi = High.95.CI/100,
           yss = Year - 1)
  message('done.')
  return(list(
    sivc    = sivc_coverage,
    statcan = statcan_coverage
  ))
  
}

get_vax_eff <- function(subtype) {
  
  message('Retrieving vaccine efficacy...', appendLF = F)
  
  # Import and format vaccine effectiveness dataset
  spsn_eff <- read_excel("data/spsn_vaccine-eff_2004-05_2019-20.xlsx", na = "-") %>%
    mutate(yss         = as.numeric(substring(season, 1, 4))) %>%
    mutate(
           date        = ymd(paste0(yss+1,'-02-01')),
           eff         = eff/100,
           vax.eff.lo  = low.95.ci/100,
           vax.eff.hi  = high.95.ci/100) %>%
    rename(vax.eff = eff) %>%
    filter(type == subtype)
  message('done.')
  return(spsn_eff)
}


get_deaths <- function() {
  
  message('Retrieving deaths...', appendLF = F)
  
  # Import influenza and pneumonia mortality data
  
  statcan_deaths <- read.csv("data/statcan_leading-causes-death_2000-2020.csv") %>%
    # Influenza deaths occurring in a calendar year N
    # are assumed to be associated with the season
    # that started in calendar year N-1:
    mutate(yss = ref_date - 1) %>%  
    filter(grepl("^Influenza", Leading.causes.of.death..ICD.10.) | 
             grepl("^Total", Leading.causes.of.death..ICD.10.)) %>%
    mutate(Age.at.time.of.death = gsub(" ", ".", 
                                       trimws(str_to_sentence(str_extract(Age.at.time.of.death, "[^,]+$")))),
           Characteristics = case_when(
             Characteristics == "Rank of leading causes of death" ~ "Rank.cause.death",
             Characteristics == "Number of deaths" ~ "Number.deaths",
             Characteristics == "Percentage of deaths" ~ "Percentage.deaths",
             Characteristics == "Age-specific mortality rate per 100,000 population" ~ "Mortality.rate")) %>%
    mutate(Characteristics2 = paste0(substr(Leading.causes.of.death..ICD.10., 1, 3),
                                     "_", Characteristics, "_", Age.at.time.of.death)) %>%
    filter(Sex == "Both sexes") %>%
    select(yss, Characteristics2, VALUE) %>%
    rename(death.rate = VALUE) %>%
    #
    # TODO: DOUBLE CHECK THIS IS WHAT WE WANT:
    filter(grepl('^Inf_Mortality.+All.ages', Characteristics2)) %>% 
    select(-Characteristics2)
  
  message('done.')  
  return(statcan_deaths)
}


categorize_positivity <- function(cases.can.wide) {
message('Categorize positivity...', appendLF = FALSE)
  # Peak positivity (using data frame from above)
  peak.pos.by.season <- cases.can.wide %>%
    group_by(yss) %>%
    summarise(peak.pos = max(positivity, na.rm = TRUE)) %>%
    mutate(peak.pos.sev = case_when(
      peak.pos >= quantile(peak.pos, 0.75) ~ "Severe",
      peak.pos >= quantile(peak.pos, 0.50) ~ "Moderate",
      peak.pos >= quantile(peak.pos, 0.25) ~ "Mild",
      TRUE ~ "Low",
    ))
  message('done.')
  return(peak.pos.by.season)
}



get_epi_data <- function(subtype, Rt.all.seasons = FALSE) {
  
  cases = get_cases(subtype)
  pkpos = categorize_positivity(cases$wideformat) %>% 
    mutate(subtype=subtype)
  
  h = hosp_ped_norm_pop(subtype) 
  
  R0 = calc_R_multi(subtype, cases, all.seasons = Rt.all.seasons)
  
  return(list(
    subtype = subtype,
    cases   = cases,
    pkpos   = pkpos,
    R0      = R0,
    vax.cvg = get_vax_cvg(),
    vax.eff = get_vax_eff(subtype), 
    deaths  = get_deaths(),
    hosp    = h$can %>% rename(yss = year)
  ))
}


