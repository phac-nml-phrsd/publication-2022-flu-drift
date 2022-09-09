###
###   MAIN SCRIPT TO RUN ANALYSES PRESENTED 
###   IN THE ASSOCIATED PUBLICATION:
###   
###   Title:   Antigenic drift and epidemiological severity of seasonal influenza in Canada
###   Authors: Zishu Chen, Christina Bancej, Liza Lee, David Champredon (david.champredon@canada.ca)
###   Journal: Scientific Reports (Springer Nature)
###   Year:    2022
### 
###  *******************
###  * * * WARNING * * * 
###  *******************
###  Some scripts in this folder will not run because
###  the data files are missing. We cannot share some of them
###  publicly (data set of pediatric hospitalizations is not public, 
###  and posting the genetic sequences would not comply with 
###  GISAID's sharing agreement). Hence, the scripts are presented 
###  here for methodological transparency. 
###  Please contact the corresponding author (david.champredon@canada.ca)
###  for any enquiries about this code. 
### 


### ==== Switches ====

do.import = 1   # import genetic sequences
do.align  = 1   # align the genetic sequences (Unix-like only)
do.gdist  = 1   # calculate pairwise genetic distances
do.merge  = 1   # merge genetic and epidemiological data
do.stats  = 1   # perform statistical analysis

### ==== Import & clean genetic sequences ====

# The influenza genetic sequences are imported
# and cleaned from FASTA files saved in the `seqs` folder

if(do.import) source('import-clean.R')


### ==== Align genetic sequences ====

# The genetic sequences must be aligned before
# the genetically-based antigenic distance can 
# be calculated. 
# WARNING: this can only be run in a Unix-like
# environment and the software `muscle` must
# be installed.

if(do.align) {
  system('./align.sh H1N1')
  system('./align.sh H3N2')
  system('./align.sh B')
}


### ==== Merge genetic & epi data ====

# The inter-season antigenic distances
# are calculated from the aligned sequences,
# the epidemiological data are retrieved and
# from the `data` folder, and finally both 
# the epi data and genetic distances are merged.

if(do.merge) source('merge-epi-gen.R')


### ==== Statistical analysis ====

# All the steps above led to a dataframe 
# named `df.final` that should be used
# for this statistical analysis.

if(do.stats) {
  source('stats-model.R')
}




