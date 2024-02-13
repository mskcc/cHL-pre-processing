####################################################################################
args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n")
    cat("   usage: checkHaloTables.R HaloDataURI\n")
    cat("\n")
    cat("      - HaloDataURI           = either a folder with Halo Spreadsheets\n")
    cat("                                or a path with a regex\n")
    cat("\n")
    quit()
}
names(args)=args
####################################################################################
source("HaloX/readHalo.R")
source("HaloX/tools.R")
source("HaloX/VERISON.R")
source("HaloX/VERSION.R")
####################################################################################
####################################################################################

getSampleId<-function(ff) {gsub("\\.csv.*","",basename(ff))}

read_haloFile <- function(ff) {

    sampleId=getSampleId(ff)
    obj=read_halo(ff)

    numAIInputs=len(unique(obj$Analysis_Inputs))

    return(list(sampleId=sampleId,numAIInputs=numAIInputs,numCells=nrow(obj)))

}

####################################################################################

PARALLEL=TRUE
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
####################################################################################
if(interactive()) {
#  stop(".INCLUDE")
}
####################################################################################
####################################################################################

dataFileURI=args[1]

if(dir.exists(dataFileURI)) {
    dataFileDir=dataFileURI
    fileRegEx=".csv.gz"
} else {
    dataFileDir=dirname(dataFileURI)
    fileRegEx=basename(dataFileURI)
}

haloFiles=dir_ls(dataFileDir,regex=fileRegEx)

xx=parallel_map(haloFiles,read_haloFile)

map(xx,as_tibble) %>%
    bind_rows(.id="path") %>%
    write_csv(cc("haloManifest","CheckHaloTable",genTimeUUID(),".csv"))


#parallel_map(haloFiles,readAndSaveHaloObject,sampleManifest,markerNameMap)

