#####################################
#### Creating the mainparams  SEE: https://github.com/labroo2/deStructure
#####################################
#make a temporary output file path
setwd("/blue/mresende/share/Givanildo/Structure_files")
output_file <- "/blue/mresende/share/Givanildo/Structure_files/RUN_STRUCTURE/mainparams.txt"
#run mainparams function; output to temp directory
#install.packages("remotes")
#remotes::install_github("labroo2/deStructure")
library(deStructure)
mainparams(maxpops = NULL, burnin = 1000, numreps = 10000, infile = "structure_output.txt",
           outfile = "my_output_file", numinds = "571", numloci = "85002",
           ploidy = 4, missing = -9, onerowperind = 0, label = 1,
           popdata = 0, popflag = 0, locdata = 0, phenotype = 0, extracols = 0,
           markernames = 1, recessivealleles = 0, mapdistances = 0, phased = 0,
           phaseinfo = 0, markovphase = NULL, notambiguous = NULL, outpath = output_file)
