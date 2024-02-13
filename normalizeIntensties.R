####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<1) {
    cat("\n")
    cat("   usage: [DATAROOT] RDAFILE\n")
    cat("\n")
    cat("      - RDAFILE                 = HALO OBJ RDAFILE\n")
    cat("\n   Optionals:\n")
    cat("      - DATAROOT [default=\"\"]   = Root of data file folder\n")
    cat("\n")
    quit()
}
names(args)=args

suppressPackageStartupMessages(require(stringr))

opt.args=list(DATAROOT="")

ii=grep("=",args)
if(len(ii)>0) {
    parseArgs=str_match(args[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){opt.args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

if(len(ii)>0) {
    args=args[-ii]
}

####################################################################################
source("HaloX/parseAI.R")
source("HaloX/tools.R")
source("HaloX/VERSION.R")
source("HaloX/HaloCore/utils.R")
####################################################################################
####################################################################################

####################################################################################

require(tidyverse, quietly = T, warn.conflicts=F)

####################################################################################
if(exists("SOURCED") && SOURCED) halt(".INCLUDE")
####################################################################################
####################################################################################

rdaFile=args[1]
obj=readRDS(rdaFile)
oo=parse_AnalysisInput(obj,opt.args$DATAROOT)
mmd=read_yaml(obj$metadata$markerMetaDataFile) %>% map(as_tibble) %>% bind_rows
comp=mmd %>% select(Marker=Marker_name,Threshold_compartment)
oo$thetas=left_join(oo$thetas,comp) %>%
    mutate(Theta = case_when(
                        Threshold_compartment=="Cytoplasm" ~ dye_cyto_positive_threshold_weak,
                        Threshold_compartment=="Nucleus" ~ dye_nuclei_positive_threshold_weak,
                        T ~ as.numeric(NA)
                        )
            )

thetas=oo$thetas %>% select(Marker,Theta)

nonExcluded.UUID=oo$geom.data %>% filter(!Exclude) %>% pull(UUID)


#    mutate(NIntensity=TIntensity) %>%

normI=oo$marker.data %>%
    filter(UUID %in% nonExcluded.UUID) %>%
    left_join(thetas,by="Marker") %>%
    mutate(NIntensity=TIntensity/Theta) %>%
    mutate(NIntensity=ifelse(
        NIntensity<1 | Positive_Classification==0,
        1,
        NIntensity)
    ) %>%
    mutate(NIntensity=log(NIntensity)) %>%
    group_by(Marker) %>%
    mutate(Scale=quantile(NIntensity[NIntensity>0],.95,na.rm=T)) %>%
    mutate(Scale=ifelse(is.na(Scale),1,Scale)) %>%
    mutate(NIntensity=NIntensity/Scale) %>%
    ungroup

#ranks=getRanksPerMarkerByFOV(oo)
#normI=left_join(normI,ranks,by=c("UUID","Marker"))

oo$marker.data=normI
oo$geom.data=oo$geom.data %>% filter(UUID %in% nonExcluded.UUID)

oo$metadata$normalizeParams=list(
    ExcludedFilter=TRUE,
    Type="NoiseFlatten,ThetaScale",
    Transform="log",
    Scale="95%-tile"
    )

oFile=gsub(".rda","_Normalized.rda",basename(rdaFile))
saveRDS(oo,oFile,compress=T)

SOURCED=T

