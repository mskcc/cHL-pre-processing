####################################################################################
args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n")
    cat("   usage: unitTest01.R RDADIR/regex\n")
    cat("\n")
    cat("      - RDADIR/regex          = either a folder with Halo RDA\n")
    cat("                                or a path with a regex\n")
    cat("\n")
    quit()
}

if(len(args)>2) {
    cat("\n  WARNING... More than one FILE/REGEX arg does not do what you expect\n")
    cat("  Only first arg taked as input RDA file\n")
    cat("  You probably want to quote that regEx\n\n")
}

####################################################################################
#
####################################################################################
####################################################################################

####################################################################################

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}

if(PARALLEL) {
    require(furrr)
    require(purrr)
    nCores=floor(availableCores()/4)
    parallel_map=future_map
    plan(multicore, workers = nCores)
    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
    require(purrr)
    parallel_map=map
    cat("Running in SERIAL\n")
}

require(tidyverse, quietly = T, warn.conflicts=F)
require(openxlsx, quietly = T)
require(fs, quietly = T)
require(gtools)
require(readxl)


####################################################################################
####################################################################################
if(interactive()) {
  stop(".INCLUDE")
}
####################################################################################

dataFileURI=args[1]
if(dir.exists(dataFileURI)) {
    dataFileDir=dataFileURI
    fileRegEx=".rda"
} else {
    dataFileDir=dirname(dataFileURI)
    fileRegEx=basename(dataFileURI)
}

rdaFiles=dir_ls(dataFileDir,regex=fileRegEx)
if(len(rdaFiles)==0) {
    cat("\n   No files match pattern\n\n")
    cat("     dataFileDir,fileRegEx=[",dataFileDir,fileRegEx,"]\n\n")
    quit()
}

#pxx=parallel_map(rdaFiles,getPosMarkerTable,cellTypeMarkers) %>% bind_rows

