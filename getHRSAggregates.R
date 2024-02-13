####################################################################################
####################################################################################
if(interactive() && exists("SOURCED") && SOURCED) {
  stop(".INCLUDE")
}
####################################################################################
####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<1) {
    cat("\n")
    cat("   usage: getHRSAggregates [THREADS] [RADIUS] RDAFILE\n")
    cat("\n")
    cat("      - RDAFILE                 = HALO OBJ RDAFILE\n")
    cat("\n   Optionals:\n")
    cat("      - THREADS [default=24]    = Number of threads to use\n")
    cat("\n")
    quit()
}
names(args)=args

suppressPackageStartupMessages(require(stringr))

opt.args=list(THREADS=24)

ii=grep("=",args)
if(len(ii)>0) {
    parseArgs=str_match(args[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){opt.args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

opt.args$THREADS=as.numeric(opt.args$THREADS)

if(len(ii)>0) {
    args=args[-ii]
}

####################################################################################
source("HaloX/VERSION.R")
source("HaloX/delaunay.R")
####################################################################################
####################################################################################

####################################################################################

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}

PARALLEL=TRUE
if(PARALLEL) {
    library(parallel)
    nCores=opt.args$THREADS
    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
    cat("Running in SERIAL\n")
}

require(tidyverse, quietly = T, warn.conflicts=F)
require(RTriangle)

####################################################################################
####################################################################################
####################################################################################

rdaFile=args[1]
obj=readRDS(rdaFile)

obj$neighbors=delaunay_nearest_neighbors(obj)

obj$neighbors=delaunay_second_neighbors(obj$neighbors)

sampleId=gsub("___HaloObj.*","",basename(rdaFile))

obj$metadata$getHRSAggregates=list(VERSION="1")

oFile=gsub(".rda","_WithDN.rda",basename(rdaFile))
saveRDS(obj,oFile,compress=T)

nFile=cc("neighborGraph",sampleId,"HRS",".rda")

ctFile=paste0("rda/v5_CellTypes/",sampleId,"____CellTypes__ExcludeFilt.rda")
ct=readRDS(ctFile)

graphTbl=obj$neighbors %>%
    filter(DelaunayNeighbor | D2Neighbors) %>%
    group_split(SPOT) %>%
    map(left_join,ct %>% select(UUID,Cell_type),by=c(C.UUID="UUID")) %>%
    map(left_join,ct %>% select(UUID,Cell_type),by=c(N.UUID="UUID")) %>%
    map(filter,Cell_type.x==Cell_type.y & Cell_type.x=="HRS") %>%
    map(left_join,select(obj$geom.data,SPOT,C.UUID=UUID,x=x0,y=y0),by = c("SPOT", "C.UUID")) %>%
    map(left_join,select(obj$geom.data,SPOT,N.UUID=UUID,xend=x0,yend=y0),by = c("SPOT", "N.UUID"))

saveRDS(graphTbl,nFile,compress=T)

# graph=igraph::graph_from_edgelist(gi %>% select(C.UUID,N.UUID) %>% as.matrix)
# connected=igraph::components(graph)

# nodes=names(which(connected$membership==11))
# xr=obj$geom.data %>% filter(SPOT==gi$SPOT[1] & UUID %in% nodes) %>% pull(x0) %>% range
# yr=obj$geom.data %>% filter(SPOT==gi$SPOT[1] & UUID %in% nodes) %>% pull(y0) %>% range

# xx %>% filter(SPOT==gi$SPOT[1] & x0>xr[1]-WIN & x0<xr[2]+WIN & y0>yr[1]-WIN & y0<yr[2]+WIN) %>% ggplot(aes(x0,y0,color=Cell_type=="HRS",shape=UUID %in% nodes)) + theme_light() + geom_point(size=3) + scale_color_manual(values=c("grey85","darkred"))

# pg=xx %>% filter(SPOT==gi$SPOT[1] & x0>xr[1]-WIN & x0<xr[2]+WIN & y0>yr[1]-WIN & y0<yr[2]+WIN) %>% ggplot(aes(x0,y0,color=Cell_type=="HRS")) + theme_light() + geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=edgeSegments,color="grey35") + geom_point(size=3) + scale_color_manual(values=c("grey85","darkred"))

SOURCED=T
