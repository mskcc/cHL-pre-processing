readAndExtractNormIntensities<-function(rr) {

    cat(">",rr,"\n")
    oo=readRDS(rr)

    oo$marker.data %>%
        select(UUID,Marker,NIntensity) %>%
        spread(Marker,NIntensity)

}

require(tidyverse)

if(interactive()) {
    PARALLEL=FALSE
} else {
    PARALLEL=TRUE
}

if(PARALLEL) {
    require(furrr)
    require(purrr)
    if(Sys.getenv()[["HOSTNAME"]]=="terra") {
        nCores=8
    } else {
        nCores=availableCores()
    }
    parallel_map=future_map
    plan(multicore, workers = nCores)
    cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
    require(purrr)
    parallel_map=map
    cat("Running in SERIAL\n")
}

rdaFiles=commandArgs(trailing=T)

tbl=parallel_map(rdaFiles,readAndExtractNormIntensities) %>% bind_rows

write_csv(tbl,"normalizedIntensities.csv.gz")
