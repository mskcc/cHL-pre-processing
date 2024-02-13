####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<4) {
    cat("\n")
    cat("   usage: [MARKER_INFO=markerMetaData.yaml] haloParameters.yaml haloMarkerRenames.yaml sampleManifest.yaml listOfHaloFiles\n")
    cat("\n")
    cat("      - haloParamters.yaml      = analysis parameters\n")
    cat("      - haloMarkerRenames.yaml  = marker name fixes\n")
    cat("      - sampleManifest.yaml     = meta data for samples\n")
    cat("      - listOfHaloFiles         = is a file of halo spreadsheet files\n")
    cat("\n")
    cat("      - MARKER_INFO             = markerMetaData.yaml\n")
    cat("\n")
    quit()
}
names(args)=args

suppressPackageStartupMessages(require(stringr))

opt.args=list(MARKER_INFO="")

ii=grep("=",args)
if(len(ii)>0) {
    parseArgs=str_match(args[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){opt.args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

if(len(ii)>0) {
    args=args[-ii]
}

####################################################################################
source("HaloX/readHalo.R")
source("HaloX/tools.R")
source("HaloX/VERSION.R")
####################################################################################
####################################################################################

getSampleId<-function(ff) {gsub("\\.csv.*","",basename(ff))}

#:DEBUG markerMetaDataFile=opt.args$MARKER_INFO
#:DEBUG markerNameMap=unlist(read_yaml("data/meta/haloMarkerRenames.yaml")$markerRenames)

read_haloObject <- function(ff,markerNameMap,markerMetaDataFile) {

    sampleId=getSampleId(ff)

    #
    # If in Debug mode (envar DEBUG=="Y")
    # only load 10,000 cells per sample
    #
    if("DEBUG" %in% names(Sys.getenv()) && Sys.getenv()[["DEBUG"]]=="Y") {
        obj=parse_halo(ff,n_max=10000)
    } else {
        obj=parse_halo(ff)
    }

    obj$geom.data=obj$geom.data %>% mutate(Sample=sampleId)

    if(!is.null(markerNameMap)) {
        obj$marker.data=obj$marker.data %>% mutate(Marker=recode(Marker,!!!markerNameMap))
    }

    #
    # Join in any marker metadata
    #

    if(file.exists(markerMetaDataFile)) {

        markerMetaData=read_yaml(markerMetaDataFile) %>%
            map(as_tibble) %>%
            bind_rows %>%
            rename(Marker=Marker_name)

        validMarkers=markerMetaData$Marker
        markersRead=obj$marker.data %>% distinct(Marker) %>% pull
        missingMarkers=setdiff(validMarkers,markersRead)

        if(len(missingMarkers)>0) {
            #
            # Backfill missing markers
            #
            uu=obj$marker.data %>% distinct(UUID)
            obj$marker.data=bind_rows(
                obj$marker.data,
                map(missingMarkers,~bind_cols(uu,Marker=.)) %>% bind_rows
            ) %>% arrange(UUID,Marker)
        }

        #
        # Then filter out any markers not in marker manifest
        #
        obj$marker.data=filter(obj$marker.data,Marker %in% validMarkers)

        mm=obj$marker.data %>%
            select(UUID,Marker,matches("_Intensity$")) %>%
            left_join(markerMetaData %>% select(Marker,Threshold_compartment),by="Marker") %>%
            gather(IntensityCompartment,TIntensity,matches("_Intensity")) %>%
            mutate(IntensityCompartment=gsub("_Intensity","",IntensityCompartment)) %>%
            filter(
                Threshold_compartment==IntensityCompartment
                | (is.na(Threshold_compartment)&IntensityCompartment=="Cell")
                )

        obj$marker.data=left_join(obj$marker.data,mm,by = c("UUID", "Marker"))

    }

    obj

}

readAndSaveHaloObject <- function(ff,sampleManifest,markerNameMap,markerMetaDataFile,outDir=".") {

    dir.create(outDir,recursive=T,showWarnings=F)

    obj=read_haloObject(ff,markerNameMap,markerMetaDataFile)

    obj$sample.data=sampleManifest %>%
        filter(File==ff) %>%
        select(
            CellDive_ID,Patient_ID,PixelScaling,Microscope,
            Pixel_to_micron,Image_type,Staining_batch,
            DAPI_first,DAPI_last,File
            ) %>%
        as.list

    obj$geom.data = obj$geom.data %>%
        mutate(PixelScaling=obj$sample.data$PixelScaling) %>%
        mutate(X=(XMax+XMin)/2,Y=(YMax+YMin)/2) %>%
        mutate(x0=PixelScaling*X,y0=-PixelScaling*Y)

    obj$metadata=list(
        markerMetaDataFile=markerMetaDataFile,
        markerNameMap=markerNameMap
        )

    saveRDS(obj,file.path(outDir,cc(obj$sample.data$CellDive_ID,"__HaloObj",VERSION,".rda")),compress=T)

    #
    # Return some basic stats
    # - count of marker positive states
    # - count of cells by FOV/SPOT
    #

    mpTbl=obj$marker.data %>%
        filter(Positive_Classification==1) %>%
        group_by(UUID) %>%
        summarize(MarkerPos=paste0(sort(Marker),collapse=",")) %>%
        ungroup %>%
        count(MarkerPos) %>%
        arrange(desc(n)) %>%
        mutate(Patient_ID=obj$sample.data$Patient_ID) %>%
        select(Patient_ID,everything())

    cellCounts=obj$geom.data %>%
        group_by(SPOT) %>%
        count() %>%
        ungroup %>%
        mutate(Patient_ID=obj$sample.data$Patient_ID) %>%
        select(Patient_ID,everything())

    list(cellCounts=cellCounts,mpTbl=mpTbl)

}

####################################################################################

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}

#
# This code needs R-3.x will not work with R-4.x
#
# if(PARALLEL) {
#     require(furrr)
#     require(purrr)
#     nCores=floor(availableCores()/4)
#     parallel_map=future_map
#     plan(multicore, workers = nCores)
#     cat("Running in PARALLEL nCores=",nCores,"\n")
# } else {
#     require(purrr)
#     parallel_map=map
#     cat("Running in SERIAL\n")
# }

if(PARALLEL) {
    library(parallel)

    nCores=ifelse("NCORES" %in% names(Sys.getenv()),as.numeric(Sys.getenv()[["NCORES"]]),8)

    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
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


haloParameters=read_yaml(args[1])
scalingMap=unlist(haloParameters$scalingPx2Microns)

markerRenames=read_yaml(args[2])
if(!is.null(markerRenames$markerRenames)) {
    markerNameMap=unlist(markerRenames$markerRenames)
} else {
    markerNameMap=NULL
}

sampleManifest.o=read_yaml_as_tibble(args[3])

haloFiles=scan(args[4],"")

fileManifest=tibble(File=haloFiles) %>% mutate(CellDive_ID=getSampleId(File))

sampleManifest=left_join(fileManifest,sampleManifest.o,by="CellDive_ID") %>%
    mutate(PixelScaling=as.numeric(scalingMap[Microscope])) %>%
    select(CellDive_ID,Patient_ID,PixelScaling,File,everything())

missingManifestSamples=sampleManifest %>% filter(is.na(PixelScaling))
if(nrow(missingManifestSamples)>0){
    cat("\n   Samples missing metadata Nsamples =",nrow(missingManifestSamples),"\n")
    cat("   ",paste0(missingManifestSamples$CellDive_ID,collapse=", "),"\n")
    cat("    Will not be processed\n\n")
    sampleManifest=sampleManifest %>% filter(!is.na(PixelScaling))
    haloFiles=sampleManifest$File
}

#
# This code needs R-3.x will not work with R-4.x
#
#ret=parallel_map(haloFiles,readAndSaveHaloObject,sampleManifest,markerNameMap,opt.args$MARKER_INFO)

if(PARALLEL) {
    ret=mclapply(haloFiles,readAndSaveHaloObject,sampleManifest,markerNameMap,opt.args$MARKER_INFO,mc.cores=nCores)
} else {
    ret=  lapply(haloFiles,readAndSaveHaloObject,sampleManifest,markerNameMap,opt.args$MARKER_INFO)
}



tbl=ret %>%
    map("mpTbl") %>%
    bind_rows %>%
    group_by(MarkerPos) %>%
    mutate(Total=sum(n)) %>%
    spread(Patient_ID,n,fill=0) %>%
    ungroup %>%
    arrange(desc(Total)) %>%
    mutate(CPCT=cumsum(Total)/sum(Total)) %>%
    select(MarkerPos,Total,CPCT,everything()) %>%
    head(1000)

cellCounts=ret %>%
    map("cellCounts") %>%
    bind_rows %>%
    group_by(Patient_ID) %>%
    mutate(Total=sum(n)) %>%
    spread(SPOT,n)

uuid=genTimeUUID()

ll=list(cellCounts=cellCounts,markerPos=tbl)
write.xlsx(ll,cc("loadStates",uuid,".xlsx"))

