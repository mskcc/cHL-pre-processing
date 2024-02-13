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
source("HaloX/VERSION.R")
source("HaloX/tools.R")

CELL_TYPE_DEF_VERSION=getCellTypeFileVersion(basename(args[1]))
cat("\n   CELL_TYPE_DEF_VERSION =",CELL_TYPE_DEF_VERSION,"\n\n")

####################################################################################
####################################################################################

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
if(interactive()) {
#  stop(".INCLUDE")
}
####################################################################################

cellTypeTable=args[1]
cellTypes=loadCellTypes(cellTypeTable)
cellTypeMarkers=getCellTypeMarkers(cellTypes)

cat("CellTypeMarkers =",paste0(cellTypeMarkers,collapse=","),"\n")

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

pxx=parallel_map(rdaFiles,getPosMarkerTable,cellTypeMarkers) %>% bind_rows

# cat("\n\n  Get rid of this once you fix Sample in the RDA files\n")
# sampleTable=read_yaml("data/meta/HodgkinLymphoma_Samples/HodgkinLymphoma_Samples___Sheet1.yaml") %>%
#     map(as_tibble) %>%
#     bind_rows
# renameIDs=sampleTable %>% select(Sample=CellDive_ID,Patient_ID)
# pxx=pxx %>% left_join(renameIDs) %>% select(-Sample) %>% rename(Sample=Patient_ID)

computeTotalCounts<-function(pxx) {

    totalCounts=pxx %>%
        group_by(Sample,PosMarkers) %>%
        summarize(N=sum(n)) %>%
        spread(Sample,N,fill=0) %>%
        mutate(Total=rowSums(.[-1])) %>%
        arrange(desc(Total))

    tbl=pxx %>%
        group_by(Sample,PosMarkers) %>%
        summarize(N=sum(n)) %>%
        mutate(Total=sum(N),PCT=N/Total) %>%
        select(-Total) %>%
        gather(Stat,Value,N,PCT) %>%
        unite(SampleStat,Sample,Stat,sep=".") %>%
        group_by(PosMarkers) %>%
        mutate(Total=sum(ifelse(grepl("\\.N$",SampleStat),Value,0))) %>%
        ungroup %>%
        mutate(SampleStat=factor(SampleStat,levels=mixedsort(unique(SampleStat)))) %>%
        arrange(SampleStat) %>%
        spread(SampleStat,Value,fill=0) %>%
        mutate(PCT.Total=Total/sum(Total)) %>%
        arrange(desc(PCT.Total))

    tbl$FullType=map(tbl$PosMarkers,getCellType,cellTypes) %>% unlist
    tbl=tbl %>%
        separate(FullType,c("Category","CellType","SubType"),sep=";;") %>%
        select(-Category) %>%
        select(CellType,SubType,PosMarkers,PCT.Total,Total,everything())

    tblS=tbl %>%
        select(-matches("PCT")) %>%
        group_by(CellType,SubType) %>%
        summarize_at(vars(matches("Total|\\.N$")),sum) %>%
        ungroup %>%
        gather(Sample,N,-matches("Type")) %>%
        group_by(Sample) %>%
        mutate(PCT=N/sum(N)) %>%
        gather(Stat,Value,N,PCT) %>%
        unite(SampleStat,Sample,Stat,sep=".") %>%
        mutate(SampleStat=gsub("\\.N$","",SampleStat)) %>%
        spread(SampleStat,Value) %>%
        select(CellType,SubType,Total,Total.PCT,everything()) %>%
        arrange(desc(Total.PCT))

    ret=list(tbl=tbl,tblS=tblS,totalCounts=totalCounts)

    ret

}

#halt("BREAK")

ret=computeTotalCounts(pxx)

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
        cellTypeMarkers=cellTypeMarkers,
        gitTag=names(git2r::tags())[1],
        gitCommit=git2r::commits()[[1]]$sha,
        haloXGitTag=names(git2r::tags("HaloX"))[1],
        haloXGitCommit=git2r::commits("HaloX")[[1]]$sha,
        commandArgs=as.list(commandArgs())
        )
    )

paramSig=substr(digest(params),1,8)
DTS=gsub("[:-]","",Sys.time()) %>% gsub(" ","_",.)

fileBase=cc("markerComboCounts","Hodgkins",oArgs$FILE_TAG,VERSION,"b","CTD",CELL_TYPE_DEF_VERSION,"_",paramSig,DTS)

write(as.yaml(params),paste0(fileBase,".params.yaml"))

param.df=data.frame(Value=unlist(params)) %>% rownames_to_column("Param")
write.xlsx(
    list(CellTypes=ret$tblS,MarkerCombo=ret$tbl,Totals=ret$totalCounts,Params=param.df),
    cc(fileBase,".xlsx")
    )

# #
# # Restrict to thresholded FOV's
# #

# sampleMetaData=read_yaml("data/meta/HodgkinLymphoma_Samples/HodgkinLymphoma_Samples___Sheet1.yaml") %>% map(as_tibble) %>% bind_rows

# thresholdedFOVs=sampleMetaData %>%
#     select(Sample=CellDive_ID,Thresholded_FOV) %>%
#     mutate(Thresholded_FOV=gsub(" ","",Thresholded_FOV)) %>%
#     separate_rows(Thresholded_FOV) %>%
#     arrange(Sample,Thresholded_FOV) %>%
#     filter(!is.na(Thresholded_FOV)) %>%
#     mutate(SampleFOVs=paste0(Sample,".f_",Thresholded_FOV)) %>%
#     pull(SampleFOVs)

# pxx.flt=pxx %>%
#     mutate(SampleFOVs=paste0(Sample,".f_",SPOT)) %>%
#     filter(SampleFOVs %in% thresholdedFOVs) %>%
#     select(-SampleFOVs)

# rm(ret)
# ret.flt=computeTotalCounts(pxx.flt)

# write.xlsx(
#     list(CellTypes=ret.flt$tblS,MarkerCombo=ret.flt$tbl,Totals=ret.flt$totalCounts,Params=param.df),
#     cc(fileBase,"ThresholdedFOVs",".xlsx")
#     )
