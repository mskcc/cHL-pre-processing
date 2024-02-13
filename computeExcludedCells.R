drawRegion <- function(rr) {
    lines(rr$x0,rr$y0,col="#BA000088",lwd=3)
}

source("HaloX/parseHaloCoordFiles.R")
source("HaloX/spatialTools.R")

boundaryLength=20

suppressPackageStartupMessages({
    require(fs)
    require(yaml)
})

args=commandArgs(trailingOnly = T)

if(len(args)<2) {
    cat("\n")
    cat("   usage: computeExcludedCells.R PROJParams.yaml RDAFILE\n")
    cat("\n")
    quit()
}

projParams=read_yaml(args[1])

boundaryCoorDirectory=projParams$datapaths$EXCLUSION_DIR
driftCoorDirectory=projParams$datapaths$DRIFTDATA_DIR

if(!dir.exists(boundaryCoorDirectory)) {
    rlang::abort(paste(
        "\n\nFATAL ERROR: boundaryCoorDirectory does not exist [",boundaryCoorDirectory,"]"
        ))
}

rdaFile=args[2]
obj=readRDS(rdaFile)

sample=obj$geom.data %>% distinct(Sample) %>% pull(Sample)
driftFile=file.path(driftCoorDirectory,paste0(sample,"_drift_summary.txt.gz"))

if(file.exists(driftFile)){
    cat("Found driftFile for sample [",sample,"]",driftFile,"\n")
    obj$geom.data = mark_drift_cells(obj$geom.data,driftFile)
}

fovMetaDataFile=dir_ls("data/meta",regex="_FOVs") %>% dir_ls
if(len(fovMetaDataFile)==0) {
    cat("\n\nERROR: Can not find FOV meta data file\n\n")
    stop("ERROR::Missing meta data")
}

spots=obj$geom.data %>% distinct(SPOT) %>% pull(SPOT) %>% as.character(.)

fovs.exclude=read_yaml_as_tibble(fovMetaDataFile) %>%
    filter(CellDive_ID==sample) %>%
    select(CellDive_ID,SPOT=FOV_number,FOV_exclusion) %>%
    mutate(Exclude=!is.na(FOV_exclusion))

fovsToKeep=fovs.exclude %>% filter(!Exclude) %>% pull(SPOT)

fovsToExcludeAndMissing=as.numeric(setdiff(spots,fovsToKeep))

obj$geom.data$Exclude.FOV=F
if(len(fovsToExcludeAndMissing)!=0) {
    cat("\n\nFOVs to exclude=",paste(sort(fovsToExcludeAndMissing)),"\n\n")
    obj$geom.data=obj$geom.data %>% mutate(Exclude.FOV=SPOT %in% fovsToExcludeAndMissing)
}

geom.data.bySPOTs=group_split(obj$geom.data,SPOT)

exclusionRegions=list()
fovAreas=list()

excludedGeomData=list()

for(ii in seq(len(geom.data.bySPOTs))) {


    print(ii)

    geom.data=geom.data.bySPOTs[[ii]]

    spot=geom.data %>% distinct(SPOT) %>% pull(SPOT)
    if(len(spot)>1) {
        cat("\n    Too many spots in geom.data data.frame [Num Spots =",len(spot),"\n\n")
        stop("FATAL_ERROR")
    }
    sample=geom.data %>% distinct(Sample) %>% pull(Sample)
    if(len(sample)>1) {
        cat("\n    Too many samples in geom.data data.frame [Num Samples =",len(sample),"\n\n")
        stop("FATAL_ERROR")
    }

    geom.data=mark_boundary_cells(geom.data,boundaryLength)

    ps=geom.data$PixelScaling[1]
    boundaryFile=file.path(boundaryCoorDirectory,paste0(sample,"_Spot",spot,".annotations"))
    if(file.exists(boundaryFile)) {

        cat(boundaryFile,"\n")

        geom.data=mark_excluded_cells_by_spot(geom.data,boundaryCoorDirectory)
        regions=read_halo_coordinate_file(boundaryFile) %>% mutate(x0=X*ps,y0=-Y*ps)

        exclusionRegions[[cc("spot",spot)]]=group_split(regions,Region)

    } else {
        geom.data$Exclude.Manual=FALSE
    }

    pDir=file.path("tmp/plots",sample)
    dir.create(pDir,recursive = T, showWarnings = F)
    pFile=file.path(pDir,paste0("exclusionPlot__",sample,"_Spot",spot,".png"))

    png(filename=pFile,
    type="cairo",
    units="in",
    width=11,
    height=7.7,
    pointsize=12,
    res=150)

    ptCls=c("#BEBEBE88","#fc8d62")

    plot(geom.data$x0,geom.data$y0,cex=.5,pch=19,col=ptCls[as.factor(geom.data$Exclude.Boundary)],
            main=paste("Sample",sample,"Spot",spot,paste0("\nBoundary ",boundaryLength,"um")),
        );
    ex1=geom.data$Exclude.Manual
    points(geom.data$x0[ex1],geom.data$y0[ex1],cex=.5,pch=19,col="#8b0000AA")

    ex2=geom.data$Exclude.Drift
    points(geom.data$x0[ex2],geom.data$y0[ex2],cex=.7,pch=19,col="#008b00AA")


    if(file.exists(boundaryFile)) {
        dum=group_split(regions,Region) %>% map(drawRegion)
    }

    dev.off()

    geom.data$Exclude = geom.data %>% select(matches("^Exclude")) %>% mutate(Exclude=rowSums(.)>0) %>% pull(Exclude)

    # if(file.exists(boundaryFile)) {
    #     area=computeFOVArea(geom.data,boundaryLength,group_split(regions,Region))
    # } else {
    #     area=computeFOVArea(geom.data,boundaryLength,NULL)
    # }

    # area$Spot=spot

    # fovAreas[[len(fovAreas)+1]]=as_tibble(area)

    excludedGeomData[[len(excludedGeomData)+1]] = geom.data

}

obj$geom.data=excludedGeomData %>% bind_rows

obj$sample.data$Exclusions=c("Boundary","Drift","Manual","FOVs")
obj$sample.data$Exclusions.Boundary.Length=boundaryLength
obj$sample.data$Exclusions.Drift=driftFile

obj$exclusionRegions=exclusionRegions

#
# Now do marker level exclusions from SAMPLE MetaData file
#

sampleMetaData=transpose(read_yaml_as_tibble(projParams$datapaths$SAMPLES))
names(sampleMetaData)=map(sampleMetaData,"CellDive_ID") %>% unlist

markersToExclude=c(
    strsplit(sampleMetaData[[obj$sample.data$CellDive_ID]]$Final_missing_markers_functional,",")[[1]],
    strsplit(sampleMetaData[[obj$sample.data$CellDive_ID]]$Final_missing_markers_identity,",")[[1]]
    )

markersToExclude=markersToExclude[!is.na(markersToExclude)]

mii=obj$marker.data$Marker %in% markersToExclude

colsToExclude=c("Cell_Intensity","Cytoplasm_Intensity","Nucleus_Intensity","Positive_Classification","TIntensity")

for(ci in colsToExclude) {
    obj$marker.data[[ci]][mii]=NA
}

saveRDS(obj,gsub(".rda","Exclusions_.rda",basename(rdaFile)),compress=T)

system2(
    "convert",
    c(
        file.path(pDir,paste0("exclusionPlot__",sample,"_Spot","*",".png")),
        paste0("exclusionPlot__",sample,".pdf")
        )
    )

#dum=map(dir(pattern="exclusionPlot.*.png"),file.remove)

