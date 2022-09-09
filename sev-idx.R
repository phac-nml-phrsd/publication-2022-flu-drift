###
###  CALCULATE SEVERITY INDEX
###

logit <- function(x) {
  res = log(x/(1-x))
  res[x <= 0] <- NA
  res[x >= 1] <- NA
  return(res)
}

#' Calculate severity index
#'
#' @param merged Dataframe of merged clinical data 
#' (i.e., reported cases, hospitalizations, etc.) 
#' @param switches Logical vector. Specify the data streams used. 
#' Used for tests only.
#' 
calc_sevidx <- function(merged, switches = c(h=1, d=1, p=1, R=1)) {
  
  a = merged %>% 
    select(subtype, yss, hosp.norm, death.rate, peak.pos, R) %>% 
    # if more than one lag was used,
    # then remove duplicated values
    # not affected by lags:
    distinct()
  
  # log (or logit) transformation
  # except for R (basic reprod. num.)
  # which is already on the log scale.
  # Also apply the switches (used only for sensitivity analysis).
  a2 = a %>% 
    mutate(
      log.h   = log(hosp.norm),
      log.d   = log(death.rate), 
      logit.p = logit(peak.pos)
    ) %>%
    select(subtype, yss, starts_with('log'), R) %>% 
    pivot_longer(cols = -c(subtype, yss))
  
  # mean and std dev needed for subsequent  
  # standardize log-transformed values:
  ss = a2 %>%
    group_by(subtype, name) %>% 
    summarise(m = mean(value, na.rm = TRUE), 
              stdev = sd(value, na.rm = TRUE),
              .groups = 'drop')
  
  # Standardization
  a3 = left_join(a2, ss, by=c('subtype','name')) %>% 
    mutate(stdval = (value - m)/stdev)
  
  # Apply on/off switches (for sensitivity analysis, tests,...)
  a3$stdval[a3$name == 'log.h'] <- a3$stdval[a3$name == 'log.h'] * switches['h']
  a3$stdval[a3$name == 'log.d'] <- a3$stdval[a3$name == 'log.d'] * switches['d']
  a3$stdval[a3$name == 'log.p'] <- a3$stdval[a3$name == 'log.p'] * switches['p']
  a3$stdval[a3$name == 'R']     <- a3$stdval[a3$name == 'R'] * switches['R']
  
  # Number of data type available and temporary   
  # sum for a given subtype and season:
  numdata = a3 %>% 
    group_by(subtype, yss) %>% 
    summarize(n = sum(!is.na(stdval)), 
              s = sum(stdval, na.rm=TRUE),
              .groups='drop')
  
  # Severity index
  sevidx = numdata %>% 
    mutate(si = s/n) %>% 
    rename(si.n = n) %>%
    select(-s)
  
  # Table to save
  y = a2 %>% pivot_wider()
  z = full_join(sevidx, y, by = c('subtype', 'yss'))
  write.csv(z , file = 'sevidx.csv', row.names = FALSE, quote = FALSE)
  
  return(sevidx)
}




