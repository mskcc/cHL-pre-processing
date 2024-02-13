####################################################################################
# This version uses the new CellType definitions table
####################################################################################
args=commandArgs(trailing=T)
PROJECT_NAME=basename(getwd())

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
    cat("   usage: [FILE_TAG=tag] PROJParams.yaml RDADIR/regex\n")
    cat("\n")
    cat("      - PROJParams.yaml       = Project paramaters (as paths to metadata files)\n")
    cat("      - RDADIR/regex          = either a folder with Halo RDA\n")
    cat("                                or a path with a regex\n")
    cat("\n")
    cat("        FILE_TAG=tag          = Optional tag to add to file name\n")
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
source("HaloX/tools.R")
source("HaloX/HaloCore/R/assign_cell_types.R")

projParams=read_yaml(args[1])

CELL_TYPE_DEF_VERSION=getCellTypeFileVersion(getMetaDataParameter(projParams,"CELLTYPES"))
cat("\n   CELL_TYPE_DEF_VERSION =",CELL_TYPE_DEF_VERSION,"\n\n")

####################################################################################
####################################################################################

getCellAtlas <- function(rdaFile) {

    cellTypeTbl = read_metadata_as_tibble(projParams,"CELLTYPES")

    cellTypeMarkers = read_metadata_as_tibble(projParams,"MARKERS") %>%
        filter(Identity=="Y") %>%
        pull(Marker_name)

    oo <- readRDS(rdaFile)

    allMarkers <- oo$marker.data %>% distinct(Marker) %>% pull

    oo$cellTypeTableV2=cellTypeTbl
    oo$cellTypeMarkers=cellTypeMarkers

    m_dat=oo$marker.data %>%
        left_join(oo$geom.data %>% select(UUID,Sample,SPOT,Exclude) %>% mutate(FOV=SPOT)) %>%
        select(-Exclude)

    oo$cellTypes=assign_cell_types_by_fov(
                    m_dat,
                    cellTypeTbl,
                    intersect(cellTypeMarkers,allMarkers),
                    flexible=T)

    posMarkers <- getPosMarkerByCell(oo$marker.data,cellTypeMarkers)

    posMarkersOthers <- getPosMarkerByCell(oo$marker.data,setdiff(allMarkers,cellTypeMarkers)) %>%
        rename(PosMarkersOthers=PosMarkers)

    cellTypes <- left_join(oo$cellTypes,posMarkers,by="UUID") %>%
        left_join(posMarkersOthers,by="UUID")

    atlas <- oo$geom.data %>%
        select(UUID,Exclude,Sample,SPOT,x0,y0,XMin,XMax,YMin,YMax) %>%
        left_join(cellTypes,by=c("Sample","UUID")) %>%
        filter(!Exclude)

    list(atlas=atlas,cellTypes=oo$cellTypes %>% filter(UUID %in% atlas$UUID))

}

####################################################################################

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}

if(PARALLEL) {
     library(parallel)
     nCores=ifelse("NCORES" %in% names(Sys.getenv()),as.numeric(Sys.getenv()[["NCORES"]]),12)
     cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
     cat("Running in SERIAL\n")
}

require(tidyverse, quietly = T, warn.conflicts=F)
require(openxlsx, quietly = T)
require(fs, quietly = T)
require(gtools)
require(readxl)

####################################################################################
####################################################################################
####################################################################################

dataFileURI=args[2]
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


####################################################################################
if(interactive()) {
    #stop(".INCLUDE")
}
####################################################################################

if(PARALLEL) {

    aa=mclapply(rdaFiles,getCellAtlas,mc.cores = nCores)

} else {

    aa=lapply(rdaFiles,getCellAtlas)

}

atlas=map(aa,"atlas") %>% bind_rows

if(!("Cell_type" %in% colnames(atlas))) {
    rlang::abort("\n\nCell_type column missing from atlas table\n\n")
}

if(any(is.na(distinct(atlas,Cell_type) %>% pull))) {

    cat("\n\nWARNING: Some Cell_types are == NA in \n\n")

}

####################################################################################
# Save Results
####################################################################################

#
# Dump args, and other parameters for tracking reproducablity
# Also make a signature from them to uniqify the output
#
names(args)=c("arg_1___projParams","arg_2___dataFileURI")

cellTypeMarkers = read_metadata_as_tibble(projParams,"MARKERS") %>%
        filter(Identity=="Y") %>%
        pull(Marker_name)


params=c(
    as.list(args),
    list(
        VERSION=VERSION,
        CELL_TYPE_DEF_VERSION=CELL_TYPE_DEF_VERSION,
        dataFileDir=dataFileDir,
        fileRegEx=fileRegEx,
        rdaFiles=rdaFiles,
        cellTypeMarkers=cellTypeMarkers,
        gitCommit=git2r::commits()[[1]]$sha,
        haloXGitCommit=git2r::commits("HaloX")[[1]]$sha,
        commandArgs=as.list(commandArgs())
        )
    )

paramSig=substr(digest(params),1,8)
DTS=gsub("[:-]","",Sys.time()) %>% gsub(" ","_",.)

fileBase=cc("cellAtlas",PROJECT_NAME,VERSION,"b","CTD",CELL_TYPE_DEF_VERSION,"_",paramSig,DTS)

write(as.yaml(params),paste0(fileBase,".params.yaml"))

write_csv(atlas,paste0(fileBase,".csv.gz"))
saveRDS(atlas,paste0(fileBase,".rda"),compress=T)

posMarkersToTypeMap=atlas %>%
    distinct(PosMarkers,.keep_all=T) %>%
    select(PosMarkers,Category,Cell_type,Subtype)

write_csv(posMarkersToTypeMap,paste0(fileBase,"MarkerTypeMap",".csv"))

#
# Save counts/stats report
#

atlas <- atlas |> filter(!Exclude)

tblG=atlas |>
    group_by(Cell_type,Subtype) |>
    summarize(N=n()) |>
    ungroup() |>
    mutate(PCT=N/sum(N)) |>
    rename(Total=N)

tblS=atlas |>
    group_by(Sample,Cell_type,Subtype) |>
    summarize(N=n()) |>
    group_by(Sample) |>
    mutate(PCT=N/sum(N)) |>
    gather(Metric,Value,N,PCT) |>
    unite(SampleMetric,Sample,Metric,sep=".") |>
    spread(SampleMetric,Value,fill=0)

cellTypesTable=left_join(tblG,tblS) %>% arrange(desc(Total))

tblG=atlas |>
    group_by(Cell_type,Subtype,PosMarkers) |>
    summarize(N=n()) |>
    ungroup() |>
    mutate(PCT=N/sum(N)) |>
    rename(Total=N)

tblS=atlas |>
    group_by(Sample,Cell_type,Subtype,PosMarkers) |>
    summarize(N=n()) |>
    group_by(Sample) |>
    mutate(PCT=N/sum(N)) |>
    gather(Metric,Value,N,PCT) |>
    unite(SampleMetric,Sample,Metric,sep=".") |>
    spread(SampleMetric,Value,fill=0)

markerComboTable=left_join(tblG,tblS) %>% arrange(desc(Total))

paramSig=substr(digest(params),1,8)
DTS=gsub("[:-]","",Sys.time()) %>% gsub(" ","_",.)

fileBase=cc("markerComboCounts",PROJECT_NAME,oArgs$FILE_TAG,VERSION,"c","CTD",CELL_TYPE_DEF_VERSION,"_",paramSig,DTS)

write(as.yaml(params),paste0(fileBase,".params.yaml"))

param.df=data.frame(Value=unlist(params)) %>% rownames_to_column("Param")
write.xlsx(
    list(
        CellTypes=cellTypesTable,
        MarkerCombo=markerComboTable,
        Params=param.df),
    cc(fileBase,".xlsx")
    )

cat("\nSaving CellType RDA files\n")
for(ii in seq(len(aa))) {
    ctFile=gsub("___Halo.*","____CellTypes__ExcludeFilt.rda",basename(rdaFiles[ii]))
    cat("  ",ctFile,"\n")
    saveRDS(aa[[ii]]$cellTypes,ctFile,compress=T)
}


