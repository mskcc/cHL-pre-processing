library(dplyr)

source("HaloX/parseHaloCoordFiles.R")

# boundaryCoorDirectory="data/raw/Hodgkin_Lymphoma_ALLcases/HaloCoordinates/ExclusionCoordinates"


exclude_points_within_regions <- function(obj,boundaryCoorDirectory) {

    geom.data.bySPOTs=group_split(obj$geom.data,SPOT)

}

exclude_points_within_regions_bySpot <- function(geom.data,boundaryCoorDirectory) {

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

    boundaryFile=file.path(boundaryCoorDirectory,paste0(sample,"_Spot",spot,".annotations"))

    if(!file.exists(boundaryFile)) {
        return(geom.data %>% mutate(ExcludeRegion=FALSE))
    }

    regions=read_halo_coordinate_file(boundaryFile)



}


## ggplot(regions,aes(X,-Y)) + geom_point(color="darkred") + geom_point(data=(geom.data),aes(XMin,-YMin),size=1,color="#BEBEBE",alpha=.75)
