suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(yaml)
})

source("HaloX/tools.R")

#
# This version is for reading Halo files version
#    HALO v3.2.1851.354
#

read_halo <- function(ff,...) {

    colTypes=cols(
      .default = col_double(),
      `Image File Name` = col_character()
    )

    read_csv(ff,col_types=colTypes,progress=F,...) %>% rename_all(~fixColNames(.))

}

fixColNames <- function(ss) {
    gsub(" ","_",ss) %>% gsub("_\\(.*\\)$","",.)
}

extractSpot <- function(Image_File_Name) {
    as.numeric(str_match(Image_File_Name,"_Spot(\\d+)")[,2])
}

# markerFields=c(
#     "Positive_Classification",
#     "Positive_Nucleus_Classification",
#     "Nucleus_Intensity",
#     "Cell_Intensity",
#     "Positive_Cytoplasm_Classification",
#     "Cytoplasm_Intensity",
#     "%_Cytoplasm_Completeness",
#     "%_Nucleus_Completeness"
# )

markerFields=c(
    "Positive_Classification",
    "Nucleus_Intensity",
    "Cell_Intensity",
    "Cytoplasm_Intensity",
    "%_Cytoplasm_Completeness",
    "%_Nucleus_Completeness"
)

markerFields.re=paste0("_",markerFields,"$",collapse="|")

parse_halo <- function(ff,...) {

    dd=read_halo(ff,...) %>%
        mutate(SPOT=extractSpot(Image_File_Name)) %>%
        mutate(UUID=generateCellUUID(Image_File_Name,XMin,XMax,YMin,YMax))

    #
    # Cell Intensity always seems to be in HaloFile even
    # if the thresholding failed. So use it to find all
    # markers
    #

    markerNames=tibble(Colnames=colnames(dd)) %>%
        filter(grepl("_Cell_Intensity$",Colnames)) %>%
        mutate(Marker=gsub("_Cell_Intensity$","",Colnames)) %>%
        pull(Marker)

    #
    # Need to be able to work with markers that have (_) in their name
    # explicity find all columns that start with a marker name
    #

    markerCols=map(markerNames,function(mm){grep(paste0("^",mm,"_"),colnames(dd))}) %>% unlist

    geom.data=dd %>% select(-all_of(markerCols))

    markerNameFields=grep(markerFields.re,colnames(dd)[markerCols],value=T)

    marker.data=dd %>%
        select(UUID,all_of(markerNameFields)) %>%
        gather(MarkerMetric,Value,-UUID) %>%
        mutate(Marker=gsub(markerFields.re,"",MarkerMetric)) %>%
        mutate(Metric=str_match(MarkerMetric,markerFields.re)[,1] %>% gsub("^_","",.)) %>%
        select(-MarkerMetric) %>%
        spread(Metric,Value)

    return(list(geom.data=geom.data,marker.data=marker.data))

}

cvtToFullHaloObject <- function(obj) {
    right_join(obj$geom.data,gather(obj$marker.data,ValueType,Value,-UUID,-Marker))
}
