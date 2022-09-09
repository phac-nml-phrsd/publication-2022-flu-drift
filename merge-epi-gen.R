###
###   MERGE EPIDEMIOLOGICAL AND GENETIC DATA SETS
###

library(ggplot2) ; theme_set(theme_bw())
library(dplyr)
library(angedist)
library(lubridate)
library(stringr)
library(stringi)
library(tidyr)
library(readxl)

source('utils.R')
source('epi-data.R')
source('calc-dist.R') 
source('sev-idx.R')

if(!exists('do.gdist')) do.gdist=0

#' Merge epidemiological and genetic data.
#'
#' @param subtype String. Influenza subtype, i.e., H1N1, H3N2 or B.
#' @param lag.inter.season Integer vector. Lags when calculating
#'  the antigenic distance between seasons.
#' @param epitope.def String. Type of epitope definition (`narrow`, `broad`)
#' @param season.rng Integer vector. Start years of the seasons 
#' (for example 2006:2019)
#' @param prms.gen List. Parameters used to read the genetic sequences.
#'
#' @return A dataframe merging the epidemiological data with the 
#' interseason genetic distances.
#' 
merge_epi_gen <- function(subtype, 
                          lag.inter.season, 
                          epitope.def, 
                          season.rng, 
                          prms.gen,
                          Rt.all.seasons = FALSE) {
  
  message(paste0('\n --- Merging epi & genetic data for ', subtype,'\n'))
  
  # ---- GENETIC DATA
  
  # Sequences object
  sobj = import_seqs(paste0('seqs/',subtype,'_cleaned_align.fasta'), prms)
  
  # Remove reference sequences used in alignment
  refsn = get_ref_alignment(subtype)
  if(!is.null(refsn)){
    idx.rm = sapply(refsn, grep, x=sobj$strain.name)
    if(length(idx.rm[[1]]>0)){
      nn = length(sobj$strain.name)
      idx.keep = c(1:nn)[-idx.rm]
      sobj = filter_index(sobj, idx.keep)
    }
  }
  
  # Antigenic distances
  dist.list = readRDS(paste0("data/dists-",subtype,".rds"))
  dist.list = dist.list[epitope.def]
  # Assuming only Hamming distance.
  # if an other distance, change code accordingly...
  names(dist.list) = paste0('hamming-', names(dist.list))
  
  # Calculate all interseason antigenic distances
  dist.inter = calc_interseas_ad (
    lag.vec    = lag.inter.season, 
    dist.list  = dist.list, 
    sobj       = sobj, 
    subtype    = subtype, 
    season.rng = season.rng)
  
  
  # ---- EPIDEMIOLOGICAL DATA 
  
  epi = get_epi_data(subtype, Rt.all.seasons)
  
  vax.cvg.sivc = epi$vax.cvg$sivc
  vax.cvg.sc   = epi$vax.cvg$statcan
  vax.eff      = epi$vax.eff
  deaths       = epi$deaths
  hosp         = epi$hosp
  R0           = epi$R0
  pkpos        = epi$pkpos
  
  # ---- Merge all 
  
  df.merge = dist.inter %>% 
    left_join(vax.cvg.sivc, by = 'yss') %>% 
    left_join(vax.cvg.sc, by = 'yss') %>% 
    left_join(vax.eff, by = 'yss') %>% 
    left_join(deaths, by ='yss') %>% 
    left_join(hosp, by = 'yss') %>% 
    left_join(R0, by = 'yss') %>% 
    left_join(pkpos, by='yss')
  
  res = df.merge %>%
    select(
      yss,
      starts_with('dist.inter'),
      n.inter,
      starts_with('vax.cvg'),
      starts_with('vax.eff'),
      hosp.norm,
      death.rate,
      distance.type, 
      season.lag,
      R, adj.r2,
      starts_with('peak.pos')) %>% 
    mutate(subtype = subtype)
  
  return(res)
}


# ==== RUN ====

# For the main analysis, 
# select: `broad3`, `full.HA12`, `narrow`

epitope.def = c( 
  # For tests:
                # 'broad',
                # 'full.HA1',
                # 'misc1',
                # 'misc2',
  # For final analysis
                'broad3', 
                'full.HA12',
                'narrow'
                )

prms = list(
  seq.type   = 'AA',
  pathogen   = 'influenza',
  seq.source = 'GISAID'
)

# Perform analysis using two lags (only lag=1 is retained for the final analysis)
lag.inter.season = c(1,2)

# Range of seasons:
season.rng       = 2006:2019

subtype = c('H1N1', 'H3N2', 'B')

# Calculate Rt for all seasons, or just 
# only when virus is clearly circulating: 
Rt.all.seasons   = FALSE   

# Calculate antigenic genetic distance
if(do.gdist){
  distcalc = lapply(subtype, calc_ad, prms=prms)
}

# Merge all the epi data with antigenic distance

merged = lapply(
  X                = subtype,
  FUN              = merge_epi_gen, 
  lag.inter.season = lag.inter.season, 
  epitope.def      = epitope.def, 
  season.rng       = season.rng, 
  prms.gen         = prms.gen,
  Rt.all.seasons   = Rt.all.seasons
) %>% 
  bind_rows()

# Calculate the severity index based 
# on epidemiological data only.

# Switches for sensitivity analysis
# (all equal to 1 for main result)
switches = c(h=1, d=1, p=1, R=1)

sevidx = calc_sevidx(merged, switches)

# Final data frame, formatted
# for downstream stats analysis:

df.final = merged %>% 
  left_join(sevidx, by = c('subtype', 'yss')) %>% 
  select(subtype, season.lag, yss, 
         starts_with('dist'), n.inter, 
         starts_with('si'),
         starts_with('vax.eff'))

save.image('merged-epi.RData')


