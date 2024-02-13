####################################################################################
args=commandArgs(trailing=T)

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
    cat("   usage: [FILE_TAG=tag] CellTypeTable.yaml RDADIR/regex\n")
    cat("\n")
    cat("      - CellTypeTable.xlsx    = Cell Type Definition Table\n")
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

CELL_TYPE_DEF_VERSION=getCellTypeFileVersion(basename(args[1]))
cat("\n   CELL_TYPE_DEF_VERSION =",CELL_TYPE_DEF_VERSION,"\n\n")

####################################################################################
####################################################################################

getCellAtlas <- function(rdaFile,cellTypeTbl) {

    cellTypeMarkers <- getCellTypeMarkers(cellTypeTbl)
    oo <- readRDS(rdaFile)

    allMarkers <- oo$marker.data %>% distinct(Marker) %>% pull

    oo <- getCellTypeByCell(oo,cellTypeTbl)

    cellTypes = oo$cellTypes %>% separate(CellType,c("Category","CellType","SubType"),sep=";;")

    posMarkersOthers <- getPosMarkerByCell(oo$marker.data,setdiff(allMarkers,cellTypeMarkers)) %>%
        rename(PosMarkersOthers=PosMarkers)

    cellTypes <-  left_join(cellTypes,posMarkersOthers,by="UUID")

    atlas <- oo$geom.data %>%
        select(UUID,Exclude,Sample,SPOT,x0,y0) %>%
        left_join(cellTypes,by="UUID")

}

####################################################################################

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}

if(PARALLEL) {
    require(furrr)
    require(purrr)
    nCores=floor(availableCores()/4)
    parallel_map=future_map
    plan(multicore, workers = nCores)
    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
    require(purrr)
    parallel_map=map
    cat("Running in SERIAL\n")
}

require(tidyverse, quietly = T, warn.conflicts=F)
require(openxlsx, quietly = T)
require(fs, quietly = T)
require(gtools)
require(readxl)


####################################################################################

source("HaloX/tools.R")

####################################################################################
####################################################################################

cellTypeTable=args[1]
cellTypes=loadCellTypes(cellTypeTable)

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
#  stop(".INCLUDE")
}
####################################################################################
atlas=parallel_map(rdaFiles,getCellAtlas,cellTypes) %>% bind_rows

# cat("\n\n  Get rid of this once you fix Sample in the RDA files\n")
# sampleTable=read_yaml("meta/HodgkinLymphoma_Samples/HodgkinLymphoma_Samples___Sheet1.yaml") %>%
#     map(as_tibble) %>%
#     bind_rows
# renameIDs=sampleTable %>% select(Sample=CellDive_ID,Patient_ID)
# atlas=atlas %>% left_join(renameIDs) %>% rename(CellDive_ID=Sample,Sample=Patient_ID)

####################################################################################
# Save Results
####################################################################################


#
# Dump args, and other parameters for tracking reproducablity
# Also make a signature from them to uniqify the output
#
names(args)=c("arg_1___cellTypeTable","arg_2___dataFileURI")

params=c(
    as.list(args),
    list(
        VERSION=VERSION,
        CELL_TYPE_DEF_VERSION=CELL_TYPE_DEF_VERSION,
        dataFileDir=dataFileDir,
        fileRegEx=fileRegEx,
        rdaFiles=rdaFiles,
        cellTypeMarkers=getCellTypeMarkers(cellTypes),
        gitTag=names(git2r::tags())[1],
        gitCommit=git2r::commits()[[1]]$sha,
        haloXGitTag=names(git2r::tags("HaloX"))[1],
        haloXGitCommit=git2r::commits("HaloX")[[1]]$sha,
        commandArgs=as.list(commandArgs())
        )
    )

paramSig=substr(digest(params),1,8)
DTS=gsub("[:-]","",Sys.time()) %>% gsub(" ","_",.)

fileBase=cc("cellAtlas","Hodgkins",VERSION,"b","CTD",CELL_TYPE_DEF_VERSION,"_",paramSig,DTS)

write(as.yaml(params),paste0(fileBase,".params.yaml"))

write_csv(atlas,paste0(fileBase,".csv.gz"))
saveRDS(atlas,paste0(fileBase,".rda"),compress=T)
