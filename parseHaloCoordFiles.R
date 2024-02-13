suppressPackageStartupMessages({
    library(xml2)
    library(tibble)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(stringr)
})

read_halo_coordinate_file <- function(xmlFile) {

    root=read_xml(xmlFile)

    vertexTable = xml_find_all(root,".//Region") %>%
        map(extract_vertices_from_region) %>%
        bind_rows(.id="Region")

    parseFileData=str_match(basename(xmlFile),"(.*)_Spot(\\d+)\\.")
    sample=parseFileData[,2]
    spot=as.numeric(parseFileData[,3])

    vertexTable %>%
        mutate(Sample=sample,SPOT=spot)

}

extract_vertices_from_region <- function(region) {
    xml_find_all(region,".//V") %>%
        map(xml_attrs) %>%
        map(enframe) %>%
        map(spread,name,value) %>%
        bind_rows(.id="V") %>%
        mutate(X=as.numeric(X),Y=as.numeric(Y))
}

