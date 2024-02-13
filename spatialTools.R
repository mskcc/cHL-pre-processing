suppressPackageStartupMessages({
    library(sp)
    library(readr)
    library(randtoolbox)
})

source("HaloX/tools.R")

intersect_region_with_geom.data<-function(regionI,gd) {
    point.in.polygon(gd$X,gd$Y,regionI$X,regionI$Y)
}

mark_excluded_cells_by_spot <- function(geom.data.spoti,boundaryCoorDirectory) {

    spot=geom.data.spoti %>% distinct(SPOT) %>% pull(SPOT)
    sample=geom.data.spoti %>% distinct(Sample) %>% pull(Sample)

    if(len(spot)>1 | len(sample)>1) {
        cat("\n    Too many spots or samples in geom.data data.frame","Spots",spot,"Samples",samples,"\n\n")
        stop("FATAL_ERROR")
    }

    geom.data.spoti = geom.data.spoti %>% mutate(Exclude.Manual=FALSE)
    boundaryFile=file.path(boundaryCoorDirectory,paste0(sample,"_Spot",spot,".annotations"))
    if(file.exists(boundaryFile)) {

        regions=read_halo_coordinate_file(boundaryFile) %>% group_split(Region)
        geom.data.spoti$Exclude.Manual = map(regions,intersect_region_with_geom.data,geom.data.spoti) %>%
                                            map(~.>0) %>% reduce(`|`)

    }

    geom.data.spoti

}

getFOVRect<-function(gd){
    c(min(gd$XMin),min(gd$YMin),max(gd$XMax),max(gd$YMax))
}

getRectArea<-function(rect) {
    (rect[3]-rect[1])*(rect[4]-rect[2])
}

drawRect<-function(b0,...) {
    #rect(xleft, ybottom, xright, ytop)
    #rect(b0[1],-b0[2],b0[3],-b0[4],...)
    lines(c(b0[1],b0[3],b0[3],b0[1],b0[1]),-c(b0[2],b0[2],b0[4],b0[4],b0[2]),...)
}

getBoundaryRect <- function(gd,boundaryLength) {
    ps=gd$PixelScaling[1]
    (
        c(min(gd$XMin),min(gd$YMin),max(gd$XMax),max(gd$YMax)) +
        c(1,1,-1,-1)*boundaryLength/ps
    )
}

rectToPolygon<-function(rect) {
    pp=matrix(
        c(
            rect[1],rect[2],
            rect[1],rect[4],
            rect[3],rect[4],
            rect[3],rect[2],
            rect[1],rect[2]
            ),ncol=2,byrow=T
        )
    colnames(pp)=c("X","Y")
    data.frame(pp)
}

getBoundaryPolygon <- function(gd,boundaryLength) {
    bb=rectToPolygon(getBoundaryRect(gd,boundaryLength))
}

mark_boundary_cells <-function(gd,boundaryLength) {
    boundaryPolygon <- getBoundaryPolygon(gd,boundaryLength)

    inside=intersect_region_with_geom.data(boundaryPolygon,gd)
    gd$Exclude.Boundary=inside==0
    gd
}

mark_drift_cells <- function(gd,driftFile,PCT.DRIFT.CUTOFF=0.10) {

    driftCellUUID=read_tsv(driftFile) %>%
        filter(drift_loss_pixel_pct>PCT.DRIFT.CUTOFF) %>%
        mutate(UUID=generateCellUUID(image_location,x_min,x_max,y_min,y_max)) %>%
        pull(UUID)
    gd$Exclude.Drift=FALSE
    gd$Exclude.Drift[gd$UUID %in% driftCellUUID]=TRUE

    gd

}

randSobolInRectangle<-function(n,rect,...) {
    width=rect[3]-rect[1]
    height=rect[4]-rect[2]
    ss=sobol(n,2,...)
    data.frame(X=ss[,1]*width+rect[1],Y=ss[,2]*height+rect[2])
}

computeFOVArea<-function(gdS,boundaryLength,exclusionRegions,Nss=10000000) {

#gdS=obj$geom.data %>% filter(SPOT==spotI)
#exclusionRegions=obj$exclusionRegions[[cc("spot",spotI)]]

    #
    # gdS is the geom.data for one SPOT in one Sample. Check that
    #

    SPOTS=unique(gdS$SPOT)
    samps=unique(gdS$Sample)

    if(len(SPOTS)>1 | len(samps)>1) {
        cat("\n    FATAL ERROR: more than one spot or sample in geom.data\n\n")
        stop("spatialTools::computeFOVArea")
    }

    ps=gdS$PixelScaling[1]
    areas=list()

    #
    # Areas are in sq-millimeters
    #

    areas$Full=getRectArea(getFOVRect(gdS))*ps^2/1e6

    paddedFov=getBoundaryRect(gdS,boundaryLength)

    areas$withPadding=getRectArea(paddedFov)*ps^2/1e6

    excludedCellPolygons=gdS %>%
        filter(Exclude.Drift) %>%
        transpose %>%
        map(~rectToPolygon(c(.$XMin,.$YMin,.$XMax,.$YMax)))

    pts=randSobolInRectangle(Nss,paddedFov)
    system.time({integral=map(c(exclusionRegions,excludedCellPolygons),intersect_region_with_geom.data,pts) %>% map(~.>0) %>% reduce(`|`)})

    areas$withExclusions=areas$withPadding*(1-mean(integral))

    areas

}


