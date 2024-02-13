####################################################################################
args=commandArgs(trailing=T)

#
# Separate out any options arguments
#
optionals=grep("=",args,value=T)

oArgs=list(FILE_TAG=NULL)
if(len(optionals)>0) {
    require(stringr, quietly = T, warn.conflicts=F)
    parseArgs=str_match(optionals,"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){oArgs[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args=grep("=",args,value=T,invert=T)

if(len(args)<1) {
    cat("\n")
    cat("   usage: [FILE_TAG=tag] RDADIR/regex\n")
    cat("\n")
    cat("      - RDADIR/regex          = either a folder with Halo RDA\n")
    cat("                                or a path with a regex\n")
    cat("\n")
    cat("\n")
    quit()
}

if(len(args)>2) {
    cat("\n  WARNING... More than one FILE/REGEX arg does not do what you expect\n")
    cat("  Only first arg taked as input RDA file\n")
    cat("  You probably want to quote that regEx\n\n")
}

####################################################################################
####################################################################################
source("HaloX/VERSION.R")

####################################################################################
####################################################################################

getCellGeom <- function(rdaFile) {

    oo <- readRDS(rdaFile)

    geom <- oo$geom.data %>%
        select(UUID,Exclude,Sample,SPOT,x0,y0,
                matches("Area"),matches("Nucleus_"),
                PixelScaling,matches("[XY](Max|Min)")) %>%
        mutate_at(vars(matches("Area")),~.*PixelScaling**2) %>%
        mutate(Nucleus_Perimeter=Nucleus_Perimeter*PixelScaling) %>%
        mutate_at(vars(matches("[XY](Max|Min)")),~.*PixelScaling)

    geom

}

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

source("HaloX/tools.R")

####################################################################################
####################################################################################

dataFileURI=args[1]

getRDAFiles<-function(dataFileURI) {

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
    } else if(!file_test("-f",rdaFiles[1])) {
        cat("\n   ",rdaFiles[1],"is not a file\n\n")
        cat("    dataFileURI =",paste0("[",dataFileURI,"]"),"\n")
        cat("    dataFileDir =",paste0("[",dataFileDir,"]"),"\n")
        cat("      fileRegEx =",paste0("[",fileRegEx,"]"),"\n")
        quit()
    }

    return(list(dataFileDir=dataFileDir,fileRegEx=fileRegEx,rdaFiles=rdaFiles))

}

rr=getRDAFiles(dataFileURI)
dataFileDir=rr$dataFileDir
fileRegEx=rr$fileRegEx
rdaFiles=rr$rdaFiles

####################################################################################
if(interactive()) {
#  stop(".INCLUDE")
}
####################################################################################
geomTbl=parallel_map(rdaFiles,getCellGeom) %>% bind_rows

####################################################################################
# Save Results
####################################################################################


#
# Dump args, and other parameters for tracking reproducablity
# Also make a signature from them to uniqify the output
#
names(args)=c("arg_1___dataFileURI")

params=c(
    as.list(args),
    list(
        VERSION=VERSION,
        dataFileDir=dataFileDir,
        fileRegEx=fileRegEx,
        rdaFiles=rdaFiles,
        gitTag=names(git2r::tags())[1],
        gitCommit=git2r::commits()[[1]]$sha,
        haloXGitTag=names(git2r::tags("HaloX"))[1],
        haloXGitCommit=git2r::commits("HaloX")[[1]]$sha,
        commandArgs=as.list(commandArgs())
        )
    )

paramSig=substr(digest(params),1,8)
DTS=gsub("[:-]","",Sys.time()) %>% gsub(" ","_",.)

fileBase=cc("cellGeom","Hodgkins",VERSION,"_",paramSig,DTS)

write(as.yaml(params),paste0(fileBase,".params.yaml"))

write_csv(geomTbl,paste0(fileBase,".csv.gz"))
saveRDS(geomTbl,paste0(fileBase,".rda"),compress=T)
