####################################################################################
args=commandArgs(trailing=T)

#
# Separate out any options arguments
#
optionals=grep("=",args,value=T)

oArgs=list(MARKER_RENAME=NULL,SAMPLE_METADATA=NULL)
if(len(optionals)>0) {
    require(stringr, quietly = T, warn.conflicts=F)
    parseArgs=str_match(optionals,"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){oArgs[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args=grep("=",args,value=T,invert=T)

if(len(args)<1) {
    cat("\n")
    cat("   usage: markerScan.R [MARKER_RENAME=MARKER_RENAME.yaml] [SAMPLE_METADATA=SampleMetaData.xlsx] haloFile1.csv [haloFile2.csv ...]\n\n")
    quit()
}

names(args)=args

####################################################################################
####################################################################################
source("HaloX/readHalo.R")
#source("HaloX/tools.R")
####################################################################################
####################################################################################

getHaloMetaData <- function(ff){

    df=read_halo(ff)
    markers=df %>%
        select(matches("_Cell_Intensity$")) %>%
        gather(Marker,Value) %>%
        distinct(Marker) %>%
        mutate(Marker=gsub("_Cell_Intensity$","",Marker))

    list(markers=markers)

}

path2sampleName <- function(pp) {
    basename(pp) %>% gsub(".csv.*$","",.)
}

getMarkerMetaData <- function(ai1) {
    map(ai1$thresholds,as_tibble) %>% bind_rows()
}

PARALLEL=TRUE
require(purrr)
if(PARALLEL) {
    library(parallel)
    nCores=floor(availableCores()/4)
    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
    cat("Running in SERIAL\n")
}

require(tidyverse)
require(openxlsx)
require(gtools)
####################################################################################
if(interactive()) {
#   stop(".INCLUDE")
}
####################################################################################
####################################################################################

cat("# Invalid Halo Files\n", file="markerScanErrorLog.txt")

if(PARALLEL) {
    x.meta=mclapply(args,getHaloMetaData,mc.cores=nCores)
} else {
    x.meta=  lapply(args,getHaloMetaData)
}


samplesWithErrors=map(x.meta,~is.null(.)) %>% unlist

x.meta=x.meta[!samplesWithErrors]

x.markers=map(x.meta,1) %>%
    bind_rows(.id="PATH") %>%
    mutate(sample=path2sampleName(PATH)) %>%
    mutate(sample=factor(sample,levels=mixedsort(unique(sample)))) %>%
    select(-PATH)

if(!is.null(oArgs$MARKER_RENAME)) {

    renameMarkers=read_yaml(oArgs$MARKER_RENAME)
    renameMarkers=unlist(renameMarkers$markerRenames)

    x.markers=mutate(x.markers,Marker=recode(Marker,!!!renameMarkers))

}

uniqueMarkers=x.markers %>% distinct(Marker) %>% arrange(Marker)

markerTbl=x.markers %>% mutate(INC="x") %>% spread(Marker,INC) %>% arrange(sample)

missingMarkerTbl=markerTbl %>%
    gather(Marker,Found,-sample) %>%
    mutate(Found=ifelse(is.na(Found),F,T)) %>%
    rename(Sample=sample) %>%
    filter(!Found) %>%
    group_by(Sample) %>%
    summarize(Observed_missing=paste(sort(Marker),collapse=","))

sampleErrorTable=tibble(Samples=names(samplesWithErrors),Valid=!samplesWithErrors) %>% filter(!Valid)

if(!is.null(oArgs$SAMPLE_METADATA)) {

    expectedMissing=read_xlsx(oArgs$SAMPLE_METADATA) %>%
        select(Sample=CellDive_ID,matches("Final_missing")) %>%
        gather(Metric,Markers,-Sample) %>%
        mutate(Metric="ExpectedMissing") %>%
        separate_rows(Markers,sep=",") %>%
        filter(!is.na(Markers))

    mm=missingMarkerTbl %>% gather(Metric,Markers,-Sample) %>% separate_rows(Markers,sep=",")


}

tbl=list(
    missingMarkers=missingMarkerTbl,
    uniqueMarkers=uniqueMarkers,
    sampleMarkers=markerTbl,
    invalid.halo.format=sampleErrorTable
    )

uuid=genTimeUUID()

write.xlsx(tbl,cc("markerInfo",uuid,".xlsx"))
write_csv(mutate(uniqueMarkers,NewName=""),cc("markerNames",uuid,".csv"))

