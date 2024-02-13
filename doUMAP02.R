stop("NOT READY TO RUN")
#
# Try to copy the UMAP done in
#
#   The immunoregulatory landscape of human tuberculosis granulomas
#   Erin F. McCaffrey, Michele Donato, …Michael Angelo Show authors
#   Nature Immunology (2022)Cite this article
#   [https://www.nature.com/articles/s41590-021-01121-x]
#
#   GitHub Repo: [https://github.com/angelolab/publications/tree/master/2022-McCaffrey_etal_HumanTB]
#   File: [https://github.com/angelolab/publications/blob/master/2022-McCaffrey_etal_HumanTB/R/TBMIBIrunUMAP.R]
#
#umap configuration parameters
#           n_neighbors: 15
#          n_components: 2
#                metric: euclidean
#              n_epochs: 200
#                 input: data
#                  init: spectral
#              min_dist: 0.5
#      set_op_mix_ratio: 1
#    local_connectivity: 1
#             bandwidth: 1
#                 alpha: 1
#                 gamma: 1
#  negative_sample_rate: 5
#                     a: 0.583029805174352
#                     b: 1.33416738413742
#                spread: 1
#          random_state: 628998288
#       transform_state: NA
#                   knn: NA
#           knn_repeats: 1
#               verbose: FALSE
#       umap_learn_args: NA
#                method: naive
#       metric.function: [function]

suppressPackageStartupMessages(require(stringr))

###############################################################################
###############################################################################

cArgs=commandArgs(trailing=T)

#
# This code will parse command line args in the form of
#    KEY=VAL
# and sets
#    args[[KEY]]=VAL
#

# Set defaults first

args=list(
            n_neighbors=15,metric="euclidean",
            min_dist=0.01,spread=1,paramTag="default",DEBUG=F,a=NULL,b=NULL,
            DOWNSAMPLE=NULL,N_THREADS=12,
            MARKERS="TOTAL"
            )

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args$n_neighbors=as.numeric(args$n_neighbors)
args$min_dist=as.numeric(args$min_dist)
args$spread=as.numeric(args$spread)
args$DEBUG=as.logical(args$DEBUG)
if(!is.null(args$a)) {
    args$a=as.numeric(args$a)
    args$b=as.numeric(args$b)
}
pArgs=grep("=",cArgs,invert=T,value=T)
normIntenFile=pArgs[1]
atlasFile=pArgs[2]

suppressPackageStartupMessages({
    require(tidyverse)
    require(umap)
    require(yaml)
    require(gplots)
})

#atlas=readRDS(atlasFile)

MARKERMETADATAFILE="data/meta/HodgkinLymphoma_Markers__V0/HodgkinLymphoma_Markers__V0___Sheet1.yaml"
markerMetaData=read_yaml(MARKERMETADATAFILE) %>%
    map(as_tibble) %>%
    bind_rows

if(args$MARKERS=="TOTAL") {

    umapMarkers=markerMetaData %>%
        filter(Total_UMAP=="Y") %>%
        pull(Marker_name)

} else if(args$MARKERS=="IDENTITY") {

        umapMarkers=markerMetaData %>%
            filter(Total_UMAP=="Y") %>%
            filter(Identity=="Y") %>%
            pull(Marker_name)


} else {

        cat("\n\n    Not implemented, args$MARKERS ==", args$MARKERS,"\n\n")
        quit()

}

CELLSET="ALL"
dd=read_csv(normIntenFile,col_types=cols(.default = "n",UUID=col_character()))
if(!is.null(args$DOWNSAMPLE)) {
    args$DOWNSAMPLE=as.numeric(args$DOWNSAMPLE)
    cat("Downsample by",1/args$DOWNSAMPLE,"\n")
    dd=dd %>% sample_frac(1/args$DOWNSAMPLE)
    CELLSET=cc("Dn",args$DOWNSAMPLE)
}

check=dd %>% select(UUID,all_of(umapMarkers)) %>% filter_all(~is.na(.))
if(nrow(check)>0) {
    rlang::abort("FATAL: NA's in intensity matrix")
}

mm=dd %>% 
    select(UUID,all_of(umapMarkers)) %>% 
    distinct(across(umapMarkers),.keep_all=T) %>% 
    filter_all(~!is.na(.)) %>% 
    as.data.frame %>% 
    column_to_rownames("UUID") %>% 
    as.matrix

#if(interactive()) halt("\n\n\tSTOP\n\n")

ncores=as.numeric(args$N_THREADS)

halt("DDDDD")

if(is.null(args$a)) {
    uu=umap(mm,
            verbose=T,n_threads=ncores,
            min_dist=args$min_dist,
            spread=args$spread,
            metric=args$metric,
            n_neighbors=args$n_neighbors
        )
} else {
    args$spread=NULL
    args$min_dist=NULL
    uu=umap(mm,
            verbose=T,n_threads=ncores,
            metric=args$metric,
            n_neighbors=args$n_neighbors,
            a=args$a,
            b=args$b
        )
}
uu1=list(embedding=uu)

params=list(
    DTS=as.character(Sys.time()),
    commandArgs=commandArgs(),
    cellSet=CELLSET,
    markers=args$MARKERS,
    normIntenFile=normIntenFile,
    args=args
    )


obj=list(
    uu=uu1,
    mm=mm,
    params=params
    )

uBase=cc("umap",
            "normFile",substr(digest::digest(normIntenFile),1,8),
            CELLSET,paste0("Markers-",args$MARKERS),
            "args",substr(digest::digest(params),1,8)
            )

saveRDS(obj,paste0(uBase,"___UUObj.rda"),compress=T)
write_yaml(params,paste0(uBase,"___PARAMS.yaml"))
