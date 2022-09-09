###
###   CALCULATE GENETIC DISTANCES
###



library(ggplot2)
library(dplyr)
library(angedist)  # devtools::install_github('phac-nml-phrsd/angedist')
library(lubridate)

source('epitope-def.R')


#' Calculate pairwise antigenic distances for all sequences. 
#' 
calc_ad <- function(virus, prms) {
  
  path = paste0('seqs/',virus,'_cleaned_align.fasta')
  sobj = import_seqs(path = path, prms)
  
  sites.list = get_epitope_positions(virus)
  
  m = list()
  for(i in seq_along(sites.list)){
    message(paste('Calculating pairwise genetic distance with `', 
                  names(sites.list)[i],'`'))
    m[[i]] = dist_matrix(sobj, ncores = parallel::detectCores()-1, 
                         sites = sites.list[[i]])
  }
  names(m) = names(sites.list)
  
  fname = paste0('data/dists-', virus,'.rds')
  saveRDS(m, file = fname)
  return(m)
}


#' Calculate antigenic distance between seasons
calc_inter_ad <- function(y, dm, sobj, lag=1, ci=0.95) {
  
  message(paste('calculating interseason distance:', y))
  
  yss = year_season_start(sobj$date.collection)
  idx  = which(yss == y)
  idx1 = which(yss == y-lag)
  
  if(length(idx)==0){
    warning(paste('No sequences found for year',y,': adjust the season range and/or lag. Stopping!'))
  }
  if(length(idx1)==0){
    warning(paste('No sequences found for year',y-lag,': adjust the season range and/or lag. Stopping!'))
  }
  
  if(length(idx) * length(idx1) == 0) 
    return(c(dist.inter.m=NA, 
             dist.inter.sd=NA,
             dist.inter.qlo=NA,
             dist.inter.qhi=NA,
             n.inter=NA))
  
  # Submatrix of the interseason distances:
  a = dm[idx, idx1]
  
  res = c(
    dist.inter.m   = mean(a),
    dist.inter.sd  = sd(a),
    dist.inter.qlo = quantile(a, probs = 0.5 - ci/2, names = FALSE),
    dist.inter.qhi = quantile(a, probs = 0.5 + ci/2, names = FALSE),
    n.inter        = length(idx) * length(idx1)
  )
  return(res)
}


# Modified function to calculate interseason distances, without plots
comp_ad_cases_2 <- function(season.rng, season.lag, dm, sobj) {
  
  # Interseason distances
  df = sapply(X = season.rng, 
              FUN = calc_inter_ad,
              dm=dm, sobj=sobj, lag = season.lag) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(yss = season.rng) %>% 
    # To facilitate comparison with cases
    # display distance at approximately
    # the date where the epidemic peak 
    # usually occurs:
    mutate(date = ymd(paste0(yss+1,'-02-01')))
  return(df)
}

#' Calculate ALL antigenic distances
#' 
calc_interseas_ad <- function(lag.vec, dist.list, 
                              sobj, subtype, 
                              season.rng = 2006:2019) {
  
  if(0){  # DEBUG
    lag.vec = 1
    dist.list = list(hamming = dm.h, blosum=dm.b)
    subtype = 'H1N1'
    season.rng = 2006:2019
  }
  
  message('Calculating interseason antigenic distance...')
  dint = list() ; k = 1
  
  for(i in seq_along(lag.vec)){
    for(j in seq_along(dist.list)){
      message(paste0('Lag = ',lag.vec[i], ' ; dist = ', names(dist.list)[j]))
      
      tmp = comp_ad_cases_2(season.rng, 
                            season.lag = lag.vec[i],
                            dm = dist.list[[j]], 
                            sobj = sobj) %>% 
        mutate(distance.type = names(dist.list)[j],
               subtype = subtype,
               season.lag = lag.vec[i])
      
      dint[[k]] = tmp
      k=k+1
    }
  }
  res = bind_rows(dint)
  return(res)
}

