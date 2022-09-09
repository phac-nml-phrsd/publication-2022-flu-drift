###
###   UTILITIES/HELPER FUNCTIONS
### 



library(lubridate)

#' Return the starting year of the 
#' Northern hemisphere season for 
#' a given date.
#'  
year_season_start <- function(d) {
  w = epiweek(d)
  y = year(d)
  res = y 
  # the dates before week 35 belong to the season 
  # that started the previous calendar year:
  # ZC: only subset values where date is not NA
  # since code wouldn't run due to NA values
  res[!is.na(w)&w<35] <- res[!is.na(w)&w<35] - 1
  return(res)
}

#' Standardize location names across different data sets
#' 
standard_locname <- function(s) {
  x = tolower(s)
  # Remove anything that is not a letter
  # and remove all special characters
  y = stringi::stri_trans_general(stringr::str_remove_all(x, pattern = '\\W+'),
                                  "Latin-ASCII")
  return(y)
}


#' Calculate the mean antigenic distance 
#' for a given season and hemisphere
#' 
calc_mad <- function(df.seqs, y, h) {
  theyears = df.seqs$season.year
  idx = which(theyears == y & df.seqs$hem == h)
  
  # Extract the sub-matrix that is from the given year
  m.i = m[idx,idx]
  # size of the sub-matrix
  n.i = nrow(m.i)
  # We have zeros on the diagonal that 
  # we don't want to include in the mean. 
  # there are n^2 in total in the matrix,
  # and n on the diagonal. So for the mean,
  # we want to divide by (n^2-n), not n^2, 
  # so we adjust the denominator:
  if(!is.null(n.i)){
    if(n.i <= 1) res = NA
    if(n.i >  1) res = mean(m.i, na.rm = TRUE) * n.i^2 / (n.i^2-n.i)
  }
  
  if(is.null(n.i) ) res = NA
  
  return(res)
}


# Calculate mean antigenic distance between seasons
calc_mad_season <- function(m, sobj.d, season.vector, season.dist) {
  
  # loop through the vector of seasons from case dataset (RVDSS)
  for(i in seq_along(season.vector)) {
    
    # get index for each season and the season before it
    idx1 = which(sobj.d$season.year == season.vector[i])
    idx2 = which(sobj.d$season.year == season.vector[i] - 1)
    
    # if there is no comparison between a season and the previous season,
    # note the distance as "NA"
    if (length(idx1) == 0 | length(idx2)== 0) {
      season.dist[i] = NA
    }
    
    #if there is a comparison between a season and the previous one
    else {
      # extract submatrix of the two seasons
      m.i = m[idx1,idx2]
      
      # calculate mean antigenic distance and add to season distance vector
      season.dist[i] = mean(m.i)
    }
  }
  return(season.dist)
}


# Summarize RVDSS cases by season
sum_cases_season <- function(data, subtype) {
  if (subtype == "H1N1") {
    filter.term = "PosFLUAP"
  }
  if (subtype == "H3N2") {
    filter.term = "PosFLUA2"
  }
  cases.by.season <- rvdss %>%
    filter(name == filter.term) %>%
    group_by(season) %>%
    summarise(num.cases = sum(v, na.rm = T))
  
  return(cases.by.season)
}


# Standardize country names to merge datasets
standard_country_name <- function(country_column) {
  # Standardize name format
  country_column_s = standard_locname(country_column)
  
  # Replace country names manually
  country_column = case_when(
    
    # Russia
    country_column_s %in% 
      standard_locname(c("Russian Federation",
                         "The Russian Federation",
                         "Russian Federation, the")) 
    ~ "Russia",
    
    # South Korea
    country_column_s %in%
      standard_locname(c("Korea, Republic of",
                         "Korea, the Republic of",
                         "Republic of Korea",
                         "The Republic of Korea",
                         "Korea, South"))
    ~ "South Korea",
    
    # Macau
    country_column_s == standard_locname("Macao") ~ "Macau",
    
    # Laos
    country_column_s %in% 
      standard_locname(c("Lao People's Democratic Republic",
                         "The Lao People's Democratic Republic",
                         "Lao, People's Democratic Republic",
                         "Lao PDR"))
    ~ "Laos",
    
    # Svalbard and Jan Mayen
    country_column_s %in% 
      standard_locname(c("Svalbard",
                         "Jan Mayen",
                         "Svalbard & Jan Mayen"))
    ~ "Svalbard and Jan Mayen",
    
    # Czech Republic
    country_column_s %in%
      standard_locname(c("Czechia",
                         "The Czech Republic"))
    ~ "Czech Republic",
    
    # Democratic Republic of the Congo
    country_column_s %in%
      standard_locname(c("Congo, the Democatic Republic of",
                         "Congo, Democatic Republic of",
                         "The Democratic Republic of the Congo",
                         "Democratic Republic of Congo",
                         "The Democratic Republic of Congo",
                         "Congo-Kinshasa"))
    ~ "Democratic Republic of the Congo",
    
    # Republic of the Congo
    country_column_s %in%
      standard_locname(c("the Republic of the Congo", 
                         "Republic of the Congo, the",
                         "Congo, Republic of",
                         "Congo, the Republic of",
                         "Congo-Brazzaville"))
    ~ "Republic of the Congo",
    
    # British Virgin Islands
    country_column_s == standard_locname("Virgin Islands, British")
    ~ "British Virgin Islands",
    
    # Gabon
    country_column_s %in%
      standard_locname(c("Gabonese Republic", 
                         "the Gabonese Republic", 
                         "Gabonese Republic"))
    ~ "Gabon",
    
    # Brazil
    country_column_s %in%
      standard_locname(c("the Federative Republic of Brazil", 
                         "Federative Republic of Brazil", 
                         "Federative Republic of Brazil, the",
                         "Brazil, Federative Republic of",
                         "Brazil, the Federative Republic of"))
    ~ "Brazil",
    
    # Colombia
    country_column_s %in%
      standard_locname(c("the Republic of Colombia", 
                         "Republic of Colombia", 
                         "Republic of Colombia, the",
                         "Colombia, Republic of",
                         "Colombia, the Republic of"))
    ~ "Colombia",
    
    # Sao Tome and Principe
    country_column_s %in%
      standard_locname(c("Sao Tome & Principe", 
                         "Democratic Republic of Sao Tome and Principe", 
                         "the Democratic Republic of Sao Tome and Principe", 
                         "Democratic Republic of Sao Tome and Principe, the",
                         "Sao Tome and Principe, Democratic Republic of",
                         "Sao Tome and Principe, the Democratic Republic of",
                         "Republic of Sao Tome and Principe", 
                         "the Republic of Sao Tome and Principe", 
                         "Republic of Sao Tome and Principe, the",
                         "Sao Tome and Principe, Republic of",
                         "Sao Tome and Principe, the Republic of"))
    ~ "Sao Tome and Principe",
    
    # Kenya
    country_column_s %in%
      standard_locname(c("the Republic of Kenya", 
                         "Republic of Kenya", 
                         "Republic of Kenya, the",
                         "Kenya, Republic of",
                         "Kenya, the Republic of"))
    ~ "Kenya",
    
    # Somalia
    country_column_s %in%
      standard_locname(c("Federal Republic of Somalia", 
                         "the Federal Republic of Somalia",
                         "Federal Republic of Somalia, the",
                         "Somalia, Federal Republic of",
                         "Somalia, the Federal Republic of"))
    ~ "Somalia",
    
    # Maldives
    country_column_s %in%
      standard_locname(c("the Maldives",
                         "Maldive Islands",
                         "the Maldive Islands",
                         "Maldive Islands, the",
                         "the Republic of Maldives", 
                         "Republic of Maldives", 
                         "Republic of Maldives, the",
                         "Maldives, Republic of",
                         "Maldives, the Republic of"))
    ~ "Maldives",
    
    # Indonesia
    country_column_s %in%
      standard_locname(c("the Republic of Indonesia", 
                         "Republic of Indonesia", 
                         "Republic of Indonesia, the",
                         "Indonesia, Republic of",
                         "Indonesia, the Republic of"))
    ~ "Indonesia",
    
    # Kiribati
    country_column_s %in%
      standard_locname(c("the Republic of Kiribati", 
                         "Republic of Kiribati", 
                         "Republic of Kiribati, the",
                         "Kiribati, Republic of",
                         "Kiribati, the Republic of"))
    ~ "Kiribati",
    
    # Keep other country names as is
    TRUE ~ country_column
  )
  
  # Remove strings after commas and brackets
  country_column = stringr::str_extract(country_column, "[^,]+")
  country_column = stringr::str_extract(country_column, "[^(]+")

  return(country_column)
}


# Manually assign hemispheres to countries and regions not included
# in the cities database
manual_hem <- function(hem_column, continent_column, country_column) {
  
  hem_column = case_when(
   
    # First assign hemisphere based on continent
    # Countries in Oceania are in southern hemisphere,
    # countries in North America and Europe are in northern hemisphere
    is.na(hem_column) & continent_column == "Oceania" ~ "S",
    is.na(hem_column) & continent_column %in% c("North America", "Europe") ~ "N",
    
    # Assign hemisphere to locations not in the "cities" database
    # Northern hemisphere: Montserrat, Palestine
    standard_locname(country_column) %in%
      standard_locname(c("Montserrat", "Palestine",
                         "Palestinian Territory",
                         "State of Palestine",
                         "Palestine, State of",
                         "The State of Palestine",
                         "Palestine, the State of")) ~ "N",
    
    # Assign hemispheres to countries on the equator
    # based on influenza seasonality (see mad-v2.2 for details)
    standard_locname(country_column) %in%
      standard_locname(c("Gabon", "Republic of the Congo",
                         "Democratic Republic of the Congo")) ~ "N",
    
    standard_locname(country_column) %in%
      standard_locname(c("Ecuador", "Colombia", "Brazil", 
                         "Sao Tome and Principe", "Uganda", "Somalia",
                         "Kenya", "Maldives", "Indonesia", "Kiribati")) ~ "S",
    
    # Keep the rest of the hemisphere assignments as is
    TRUE ~ hem_column
  )
  
  return(hem_column)
}

#' Retrieve the strain name of reference sequences for alignment.
#' @param v String. Virus name (e.g., \code{H1N1})
get_ref_alignment <- function(v) {
  ref = read.csv('data/ref-strains-alignment.csv', strip.white = TRUE)
  if(!v %in% unique(ref$virus)) {
    message(paste('No alignment reference found for virus',v))
    return(NULL)
  }
  res = as.character(ref$reference_strains[ref$virus==v])
  message(paste('Alignment references for ',v,':\n',paste(res,collapse='\n ')))
  return(res)
}

# Code to test (not executed):
if(0){
  d  = c(ymd('2002-02-11'), ymd('2020-08-14'))
  year_season_start(d)
}
