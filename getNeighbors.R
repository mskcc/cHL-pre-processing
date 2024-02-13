####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<1) {
    cat("\n")
    cat("   usage: getNeighbors [THREADS] [RADIUS] RDAFILE\n")
    cat("\n")
    cat("      - RDAFILE                 = HALO OBJ RDAFILE\n")
    cat("\n   Optionals:\n")
    cat("      - THREADS [default=24]    = Number of threads to use\n")
    cat("      - RADIUS  [default=50]    = Neighbor radius (in microns)\n")
    cat("\n")
    quit()
}
names(args)=args

suppressPackageStartupMessages(require(stringr))

opt.args=list(RADIUS=50,THREADS=24)

ii=grep("=",args)
if(len(ii)>0) {
    parseArgs=str_match(args[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){opt.args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

opt.args$RADIUS=as.numeric(opt.args$RADIUS)
opt.args$THREADS=as.numeric(opt.args$THREADS)

if(len(ii)>0) {
    args=args[-ii]
}

####################################################################################
source("HaloX/VERSION.R")
source("HaloX/neighborTools.R")
####################################################################################
####################################################################################

####################################################################################

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}


if(PARALLEL) {
    library(parallel)
    nCores=12
    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
    cat("Running in SERIAL\n")
}

require(tidyverse, quietly = T, warn.conflicts=F)

####################################################################################
if(interactive()) {
#  stop(".INCLUDE")
}
####################################################################################
####################################################################################

rdaFile=args[1]
obj=readRDS(rdaFile)

#
# Hard remove Excluded Cells; do not use in neighbor calcuation
#

obj$geom.data=obj$geom.data %>% filter(!Exclude)
obj$marker.data=obj$marker.data %>% filter(UUID %in% obj$geom.data$UUID)
obj$metadata$exclusion=list(COMMENT="Excluded Cells Deleted")


# ll=obj$geom.data %>%
#     group_split(SPOT) %>%
#     parallel_map(getNeighborhoodsForFOV,opt.args$RADIUS)

os=obj$geom.data %>% group_split(SPOT)

if(PARALLEL) {
    ll=mclapply(os,getNeighborhoodsForFOV,opt.args$RADIUS,mc.cores=nCores)
} else {
    ll=  lapply(os,getNeighborhoodsForFOV,opt.args$RADIUS)
}

if(any(map(ll,"numMaxRltR") %>% unlist != 0)) {
    cat("\n\n  FATAL ERROR::Some FOV's did not have MaxR >= RADIUS\n")
    cat("  Need to fix computation of numNeighborCells\n\n")
    rlang::abort("FATAL ERROR::maxR not satisfied")
}

nnTbl=map(ll,"nnTbl") %>% bind_rows()
numNeighborCells=map(ll,"numNeighborCells") %>% unlist %>% unique

obj$neighbors=nnTbl

obj$metadata$neighborParams=list(
    RADIUS=opt.args$RADIUS,
    numNeighborCells=numNeighborCells
    )



oFile=gsub(".rda",sprintf("_NeighborTbl_%d.rda",opt.args$RADIUS),basename(rdaFile))
saveRDS(obj,oFile,compress=T)
