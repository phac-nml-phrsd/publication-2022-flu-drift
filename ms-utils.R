###
###   Utilities for manuscript
###


library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(angedist)

# List of accession numbers for a given virus

get_accession <- function(virus) {
  # virus='H1N1'
  prms = list(seq.source='GISAID', pathogen='influenza', seq.type='AA')
  s = angedist::import_seqs(path=paste0('seqs/', virus, '_GISAID_cir.fasta'), prms)
  a = s$accession.num
  return(a)
}

# List all accession number used

get_accessionnumber_gisaid <- function() {
  
  x = sapply(c('H1N1', 'H3N2', 'B'), FUN = get_accession)
  y = unlist(x)
  
  z = y[!duplicated(y)]
  
  write.csv(z, file = 'seqs/accession-number-list.txt', row.names = FALSE, 
            quote = FALSE)
}

