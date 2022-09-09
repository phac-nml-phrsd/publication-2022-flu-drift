###
###   CALCULATE THE BASIC REPRODUCTION NUMBER (R0) FROM REPORTED CASES
###



#' Calculate Ro based on Gamma approximation.
#' 
#' @param Gbar Numeric. Mean of the intrinsic generation interval distribution.
#' @param SD Numeric. Standard deviation of the intrinsic generation interval distribution
#' @param season.year String. The year starting the season.
#' @param cases  Dataframe of cases
#' @param diagnose.plot Display diagnostic plot?
#' 
calc_R <- function(season.year, 
                   cases,
                   Gbar = 3.6, SD=1.6, 
                   diagnose.plot = FALSE) {
  
  # Find week of peak cases for each season
  peak.cases.subtype <- cases$wideformat %>%
    group_by(yss) %>%
    summarise(start_date = date[1],
              peak_date_idx = which.max(pos),
              peak_date = date[which.max(pos)])
  
  
  # Set Gbar (mean generation interval) and kappa (CV^2) values 
  # based on the literature
  # Values from Cowling et al. (2009)
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3057478/
  # Gbar = 3.6
  # SD = 1.6
  kappa = (SD/Gbar)^2
  
  # Determine index of input season in peak cases dataframe
  season.idx = which(peak.cases.subtype$yss == season.year)
  
  # Date of peak cases of given input season based on index
  peak.date = peak.cases.subtype$peak_date[season.idx]
  
  # Set start and end dates based on peak case date
  # End = 3 weeks before peak, start = 7 weeks before end
  date.range.end = seq.Date(peak.date, length.out = 2, by = "-3 weeks")[2]
  date.range.start = seq.Date(date.range.end, length.out = 2, by = "-7 weeks")[2]
  
  # Filter case data between the start and end dates
  df = cases$wideformat %>%
    filter(between(date, date.range.start, date.range.end)) %>%
    mutate(log_cases = log(pos))
  
  # Linear regression model
  model = lm(data = df, formula = log_cases ~ date)
  
  ms = summary(model)
  
  if(diagnose.plot){
    fname = paste(subtype, season.year,sep='-')
    g = ggplot(df, aes(x=date, y=log_cases))+
      geom_point()+
      geom_smooth(method='lm', formula = 'y~x')+
      labs(x='',y='',title=fname, 
           subtitle = round(ms$adj.r.squared,4))
    pdf(paste0(fname,'.pdf'))
    plot(g)
    dev.off()
  }
  
  
  # Slope of regression model (date coefficient)
  # This is the exponential growth rate
  r = as.double(coef(model)[2])
  
  # Calculate R (reproduction number)
  R = (1 + kappa * r * Gbar)^(1/kappa)
  
  res = list(R = R, adj.r.squared = ms$adj.r.squared)
  
  return(res)
}

#' Select season based on subtype 
#' (subtype barely circulates some season, so they are excluded)
season.subtype <- function(subtype) {
  
  if(!subtype %in% c('H1N1', 'H3N2', 'B')) 
    stop(paste('Seasons to calculate R0 not defined for:',subtype))
  
  if (subtype == "H1N1") {
    res <- c(2013, 2015, 2018, 2019)
  }
  if (subtype == "H3N2") {
    res <- c(2010:2012, 2014, 2016:2019)
  }
  if (subtype == "B") {
    res <- c(2007:2008, 2010:2019)
  }
  return(res)
}


#' Calculate R0 for multiple virus subtypes
#'
#' @param subtype String specifying the subtype (e.g., "H1N1")
#' @param cases Dataframe of cases
#' @param all.seasons Is the calculation for all seasons or just selected ones? 
#'
calc_R_multi <- function(subtype, cases, all.seasons=FALSE) {
  
  message(paste('Calculating Ro for',subtype))
  
  # By default, Rt is calculated for seasons
  # where the (sub)type was clearly 
  # circulating. But there is the option to 
  # calculate for all seasons (in that case check diagnostic plots!)
  x = season.subtype(subtype)
  if(all.seasons) x = 2007:2019
  
  rmod = sapply(x, calc_R, 
                cases = cases,
                diagnose.plot = FALSE)
  
  R.estim = as.numeric(rmod[1,])
  adj.r2  = as.numeric(rmod[2,])
  
  res = data.frame(yss = x,
                  R = R.estim, adj.r2 = adj.r2,
                   subtype = subtype)
  
  return(res)
}


