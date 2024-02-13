require(tidyverse, quietly = T, warn.conflicts=F)
require(yaml, quietly = T, warn.conflicts=F)

#DROOT="/Users/socci/Desktop/Halo/dev/juno"

#RFILE="H16_5562___HaloObj_v10.7__210714_Exclusions_.rda"
#obj=readRDS(RFILE)

parse_AnalysisInput<-function(obj,DROOT="") {

    summaryFile=file.path(DROOT,
        dirname(obj$sample.data$File),
        "Summary_data",
        gsub(".csv.*","_Summary.csv.gz",basename(obj$sample.data$File))
        )

    xx=read_csv(summaryFile) %>% filter(!is.na(`Analysis Inputs`))

    ai=xx[["Analysis Inputs"]]

    #
    # Check that "Analysis Inputs" same for all spots
    #

    if(len(unique(ai))>1) {
        cat("\n\tFATAL ERROR::parseAI.R 'Analysis Inputs' different for different SPOTS\n")
        cat("\tSummaryFile =",summaryFile,"\n")
        cat("\n\n")
        stop("ERROR")
    }

    aid=tibble(aif=strsplit(ai,";")[[1]]) %>%
        separate(aif,c("Key","Value"),sep=":")

    markerRenames=obj$metadata$markerNameMap
    markers=aid %>%
        filter(grepl("Dye",Key)) %>%
        mutate(Marker=recode(Value,!!!markerRenames)) %>%
        mutate(MarkerNo=as.numeric(gsub("Dye ","",Key)))

    extractFieldByRegex<-function(aid,regex) {
        aid %>%
            filter(grepl(regex,Key)) %>%
            bind_cols(markers %>% select(Marker,MarkerNo)) %>%
            mutate(Value=as.numeric(Value))
    }

    nucTheta=extractFieldByRegex(aid,"nuclei.*thres.*weak")
    cytoTheta=extractFieldByRegex(aid,"cyto.*thres.*weak")
    nucCompTheta=extractFieldByRegex(aid,"dye_nuclei_completeness_threshold")
    cytoCompTheta=extractFieldByRegex(aid,"dye_cyto_completeness_threshold")

    thetas=bind_rows(list(nucTheta,cytoTheta,nucCompTheta,cytoCompTheta)) %>%
        spread(Key,Value) %>%
        arrange(MarkerNo)
        # %>%
        # mutate(Theta=ifelse(is.na(dye_nuclei_positive_threshold_weak),dye_cyto_positive_threshold_weak,dye_nuclei_positive_threshold_weak)) %>%
        # mutate(Theta.Comp=ifelse(is.na(dye_nuclei_positive_threshold_weak),"cyto","nuclei"))

    obj$thetas=thetas

    obj

}
