###
###   INGEST PEDIATRIC HOSPITALIZATION DATA
###   (data cannot be publicly shared, so code will not run)
###


library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)

library(haven)


get_hosp_ped <- function() {
  
  message('pediatric ...', appendLF = F)
  
  # Retrieve pediatric hospitalization associated with influenza
  # from the IMPACT study (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7331140/). 
  # Data provided by Christina Bancej's group. 
  
  # --- Read raw data
  dat = read_sas('data/impactmaster_2020.sas7bdat')
  
  # --- Recode into human readable values
  
  df = dat %>%
    # Derive "province" based on "centre" variable
    # using data dictionary coding as reference
    mutate(province = case_when(centre == 1 ~ "NS",
                                centre %in% c(2, 8, 10) ~ "QC",
                                centre %in% c(3, 6) ~ "ON",
                                centre == 4 ~ "MB",
                                centre == 5 ~ "BC",
                                centre %in% c(7, 11) ~ "AB",
                                centre == 9 ~ "NL",
                                centre == 12 ~ "SK"),
           # Derive variables - label categories
           Subtype2 = case_when(FluType == 2 ~ "B",
                                Subtype %in% 1:3 ~ "H1N1",
                                Subtype == 4 ~ "H3N2",
                                Subtype == 5 ~ "Other",
                                Subtype == 6 ~ "Unknown",
                                Subtype == 7 ~ "Subtyping not done"),
           Subtype = case_when(Subtype == 1 ~ "A/H1N1 2009",
                               Subtype == 2 ~ "A/H1N1 seasonal",
                               Subtype == 3 ~ "H1N1 subtype not known",
                               Subtype == 4 ~ "A/H3N2",
                               Subtype == 5 ~ "Other",
                               Subtype == 6 ~ "Unknown",
                               Subtype == 7 ~ "Subtyping not done"),
           death = case_when(Death == 1 ~ "No",
                             Death == 2 ~ "Yes",
                             TRUE ~ "Unknown/NA"),
           vaccinated = case_when(Vacc_Status == 1 ~ "Yes",
                                  Vacc_Status == 2 ~ "No",
                                  Vacc_Status == NA ~ "NA"),
           gender = case_when(gender == 1 ~ "Male",
                              gender == 2 ~ "Female",
                              gender == NA ~ "NA")) 
    
    # --- Redefine variables
  
    df = df %>% 
           mutate(season  = paste0("20", substr(season, 1, 2), "-", substr(season, 3, 4)),
                  # Note: the "type" variable recoding below could be misclassifying
                  # viruses with FluType = 3 ("Other"), 4 ("Unknown") and 5 ("A and B")
                  # as type B viruses - may want to change the coding
                  # to prevent errors in the future
                  type    = ifelse(FluType==1,'A','B'),
                  subtype = ifelse(Subtype2 %in% c('H1N1','H3N2','B'),Subtype2,NA))
    
    # --- Select & rename relevant variables
    
    df = df %>% 
      select(season, province, adat, age,
             type, subtype, Tot_days_hosp, 
             Days_ICU, death, vaccinated) %>%
      # Rename variable
      rename(date.admission = adat, 
             days.hosp = Tot_days_hosp, days.icu = Days_ICU)
    return(df)
}



get_pop <- function() {
  # Population growth in Canada, by age groups
  
  a = read.csv('data/popcanada.csv')
  
  res = a %>% pivot_longer(-Age.group) %>% 
    mutate(year = as.numeric(substr(name,2,5))) %>% 
    select(-name) %>% rename(agegroup = Age.group)
  
  # age group bounds:
  
  res$agegroup[res$agegroup=='All ages'] = '0 to 999'
  res$agegroup[res$agegroup=='100 years and over'] = '100 to 999'
  
  lo = str_extract(res$agegroup, '^\\d+') %>% as.numeric()
  hi = str_extract(res$agegroup, '\\s\\d+') %>% as.numeric()
  
  res = mutate(res, lo=lo, hi=hi)
  
  return(res)
  df %>% 
    ggplot(aes(x=year, y=value))+ 
    geom_line()+geom_point()+
    facet_wrap(~Age.group, scales = 'free_y')
  
}


hosp_ped_norm_pop <- function(subt = NULL) {
  message('Retrieving hospitalization...', appendLF = FALSE)
  pop  = get_pop()
  hosp = get_hosp_ped()
  
  pop = pop %>% 
    filter(agegroup != '0 to 999') %>% 
    mutate(ag = paste(lo,hi,sep='_'))
  
  # Total population aged 0-19 each year
  pop2 = pop %>% 
    filter(lo < 20) %>%
    group_by(year) %>% 
    summarize(tot = sum(value))
  
  
  bb = unique(c(0, pop$hi))
  ll = unique(paste(pop$lo, pop$hi, sep='_'))
  
  hosp$ag = cut(hosp$age, 
                breaks = bb, 
                labels = ll) 
  hosp$ag[hosp$age <= 0] <- ll[1]
  hosp$year <- as.numeric(substr(hosp$season, 1, 4))
  
  # clean up
  hosp = filter(hosp, !(type=='B' & grepl('^H',subtype)))
  
  # Filter the subtype requested
  if(!is.null(subt)) hosp = filter(hosp, subtype == subt)
  
  # stratified by province, age group and subtype:
  
  df = hosp %>% 
    group_by(year, province, ag, type, subtype) %>% 
    summarize(n.hosp = n(), .groups = 'drop') %>% 
    left_join(pop, by = c('year', 'ag')) %>% 
    # Warning: not quite the right normalization
    # because the _national_ pop is used instead 
    # of the _provincial_ one.
    mutate(hosp.norm = n.hosp / value * 100000)
  
  # National level
  
  df.can =  ungroup(hosp) %>% 
    group_by(year, ag, type, subtype) %>% 
    summarize(n.hosp = n(), .groups = 'drop')  %>% 
    left_join(pop, by = c('year', 'ag')) %>% 
    mutate(hosp.norm = n.hosp / value * 100000)
  
  df.can2 =  ungroup(hosp) %>% 
    group_by(year, type, subtype) %>% 
    summarize(n.hosp = n(), .groups = 'drop')  %>% 
    left_join(pop2, by = c('year')) %>% 
    mutate(hosp.norm = n.hosp / tot * 100000)
  
  
  res = list(
    prov.ag = df,
    can.ag = df.can,
    can = df.can2
  )
  message('done.')
  return(res)
  
  g = df %>% ggplot(aes(x=year, y=hosp.norm, colour=ag))+
    geom_line(size=1) +
    scale_x_continuous(breaks = seq(2010, 2020, by=2))+
    scale_y_log10()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    labs(
      title = 'Pediatric hospitalization for influenza',
      x='', y='hosp. per 100,000')
   
  g.prov = g + facet_grid(type + subtype ~ province, scale='fixed')
  g.prov
  g.can = g %+% df.can + facet_wrap(type ~ subtype, scale='fixed')
  g.can
  
  g.can2 = df.can2 %>% 
    ggplot(aes(x=year, y=hosp.norm))+
    geom_line(size=1)+
    scale_x_continuous(breaks = seq(2010, 2020, by=2))+
    facet_wrap(type ~ subtype)+
    labs(
      title = 'Pediatric hospitalization for influenza',
      x='', y='hosp. per 100,000')
  g.can2  
}



