
###
###   IMPORT SEQUENCES AND FILTER OUT LOW QUALITY SEQUENCES
###   (this is typically done before alignment)
###

library(ggplot2)
library(dplyr)
library(angedist)  # devtools::install_github('phac-nml-phrsd/angedist')
library(lubridate)
library(seqinr)

source('utils.R')

# Descriptive parameters needed by `angedist` functions.
prms = list(
  seq.type   = 'AA',
  pathogen   = 'influenza',
  seq.source = 'GISAID'
)


#' Clean the sequences saved in FASTA files.
#'
#' @param v String. The virus considered ("H1N1", "H3N2", "B")
#' @param nsubs Integer. Number of sequences sub-sampled from the FASTA file.
#' @param loc String. Geographical locations to be filtered.
#'
#' @return A FASTA file keeping only the sequences that satisfy quality thresholds.
#' 
clean_seqs <- function(v, nsubs, loc=NULL) {
  
  if(0) { # DEBUG
    v='H3N2'
    v='H1N1'
    v = 'B'
    # nsubs=1e4
  }
  
  path = paste0('seqs/',v,'_GISAID_cir.fasta')
  sobj = import_seqs(path = path, prms)
  
  if(!is.null(sobj)){
    n1   = length(sobj$seq)
    
    message(paste0('Locations selected:\n', paste(loc,collapse = ', ')))
    
    # Strain to exclude from being filtered out
    # to make sure they are present in the alignment
    # (they are references for epitopes positions)
    xsn  = get_ref_alignment(v)
    sobj = filter_location(sobj, loc, except.strain.name = xsn)
    
    n2   = length(sobj$seq)
    message(paste0('Before filtering locations: ',n1,' sequences\n',
                   'After filtering locations:  ',n2,' sequences'))
  }
  
  if(grepl('^H\\dN\\d', v))  
    q = clean_seq_influenza_A(sobj, verbose = T)
  
  if(v == 'B') 
    q = clean_seq_influenza_B(sobj, verbose = T)
  
  if(length(sobj$seq) <= nsubs){
    msg = paste0('Size of subsample requested (',nsubs,
                 ') is larger than total numbers of sequences (',
                 length(sobj$seq),
                 '). Subsampling based on date not performed and returning all sequences.')
    print(msg) ; message(msg)
    qs = q
  }
  
  if(length(sobj$seq) > nsubs){
    qs = subsample_date(q,
                        n.subsample = nsubs,
                        verbose = T)
  }
  
  fname =  paste0('seqs/', v,'_cleaned.fasta')
  
  seqinr::write.fasta(sequences = qs$seq, 
                      names     = qs$header, 
                      file.out  = fname)
  message(paste0('Cleaned FASTA file of ',nsubs,
                 ' subsamples written in: ',fname))
  
  return(list(original=q, subsample=qs))
}

# ---- RUN ----

viruses = c('H1N1', 'H3N2', 'B')

# Only Canadian locations:

loc     = c('britishcolumbia','alberta', 'saskatchewan',
            'manitoba', 'ontario', 'quebec',
            'northwestterritories','nunavut', 'yukon',
            'newfoundland','novascotia','princeedwardisland', 'newbrunswick')

nsubs   = 1e9 # `nsubs` very large means no sub-sampling

# Results are FASTA files
# (variable `a` will not be used)
a = lapply(viruses, clean_seqs, nsubs=nsubs, loc=loc)

  



