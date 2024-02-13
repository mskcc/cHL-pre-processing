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
####################################################################################
####################################################################################

####################################################################################

library(tidyverse)

####################################################################################
if(interactive() && exists(".INCLUDE") && .INCLUDE) {
    halt("INCLUDE")
}
.INCLUDE=TRUE
####################################################################################
####################################################################################

rdaFile=args[1]
oo=readRDS(rdaFile)

nonExcluded.UUID=oo$geom.data %>% filter(!Exclude) %>% pull(UUID)

geom.data=oo$geom.data %>% filter(UUID %in% nonExcluded.UUID) %>% select(UUID,Sample,SPOT)
marker.pseudo=oo$marker.data %>% filter(UUID %in% nonExcluded.UUID) %>% group_by(Marker) %>% summarize(MinIP=min(TIntensity[TIntensity>0]))

#
# Get thresholds
#
ai=parse_AnalysisInput(oo,opt.args$DATAROOT)
mmd=read_yaml(ai$metadata$markerMetaDataFile) %>% map(as_tibble) %>% bind_rows
comp=mmd %>% select(Marker=Marker_name,Threshold_compartment)
ai$thetas=left_join(ai$thetas,comp) %>%
    mutate(Theta = case_when(
                        Threshold_compartment=="Cytoplasm" ~ dye_cyto_positive_threshold_weak,
                        Threshold_compartment=="Nucleus" ~ dye_nuclei_positive_threshold_weak,
                        T ~ as.numeric(NA)
                        )
            )
thetas=ai$thetas %>% select(Marker,Theta)
thetas=thetas %>% left_join(marker.pseudo) %>% mutate(ThetaS=Theta/MinIP)

thetasS=thetas$ThetaS
names(thetasS)=thetas$Marker
thetasS=as.list(thetasS)

#
# Get intensity table
#

#nonExcluded.UUID=sample(nonExcluded.UUID,10000)

dd=oo$marker.data %>%
    filter(UUID %in% nonExcluded.UUID) %>%
    left_join(thetas,by="Marker") %>%
    left_join(geom.data,by="UUID") %>%
    mutate(Pos=factor(Positive_Classification),SPOT=factor(SPOT)) %>%
    mutate(SIntensity=TIntensity/MinIP)

FOVMetaDataFile="data/meta/HodgkinLymphoma_FOVs__V0/HodgkinLymphoma_FOVs__V0___Sheet1.yaml"
thresholdFOV=read_yaml_as_tibble(FOVMetaDataFile) %>% filter(FOV_used_for_thresholding=='y' & CellDive_ID==dd$Sample[1]) %>% pull(FOV_number)

fSpotLevels=levels(dd$SPOT)
recodeSpot=ifelse(fSpotLevels==thresholdFOV,paste0(fSpotLevels,"*"),fSpotLevels)
names(recodeSpot)=fSpotLevels
dd$SPOT=recode_factor(dd$SPOT,!!!recodeSpot)

#markersToUse=c("CD4","CD8","CD30","MUM1")
markersToUse=c("CD4","CD8")

thetaD=filter(thetas,Marker %in% markersToUse)

# logBreaks=c(0,1,10,100,1000,10000,100000,1000000,10000000)
# logBreaksLabels=c("0","1","10","100","1,000","10,000","100,000","1,000,000","10,000,000")

logBreaksS=c("0","1","10",paste0("10^",(2:7)))
logBreaks=map(logBreaksS,str2expression)%>%map(eval)%>%unlist
logBreaksLabels=str2expression(logBreaksS)

qTbl=dd %>%
    group_by(Marker,SPOT) %>%
    mutate(Thres=ifelse(TIntensity>Theta,1,0)) %>%
    select(Pos,Thres) %>%
    mutate(N=n()) %>%
    group_by(Marker,SPOT,N) %>%
    count(Pos,Thres) %>%
    unite(PosThres,Pos,Thres) %>%
    spread(PosThres,n) %>%
    mutate(PCT.Cond=(`0_0`+`1_1`)/N) %>%
    mutate(Q=-10*log10(1-PCT.Cond)) %>%
    mutate(Sample=dd$Sample[1]) %>%
    select(Sample,SPOT,Marker,everything())

qTblD=qTbl %>% filter(Marker %in% markersToUse) %>% mutate(Ql=substr(sprintf("%.2f",round(Q,2)),1,3)) %>% mutate(yl=(as.numeric(SPOT)-1)%%2+.5)

# pg=dd %>% filter(Marker %in% markersToUse) %>%
#     ggplot(aes(SPOT,SIntensity,color=Pos)) + theme_light() + geom_violin() +
#     geom_hline(data=thetaD,aes(yintercept=ThetaS)) + facet_wrap(~Marker,ncol=1,scales="free") +
#     scale_y_continuous(trans="log1p",breaks=logBreaks,labels=logBreaksLabels) +
#     ggtitle(dd$Sample[1]) + geom_text(data=qTblD,aes(x=SPOT,y=yl,label=Ql),color="black",hjust=0)

# pm=list()
# for(mi in markersToUse) {
#     cat(mi)
#     pm[[mi]]=dd %>% filter(Marker==mi) %>%
#         ggplot(aes(SIntensity,fill=Pos)) + theme_light() + geom_density(alpha=.5) + geom_vline(xintercept=thetasS[[mi]]) + facet_wrap(~SPOT,scale="free") + scale_x_continuous(trans="log1p",breaks=logBreaks,labels=logBreaksLabels) + ggtitle(paste(dd$Sample[1],mi))
#     cat("\n")
# }

# pfile=cc("qcPlt","Threshold","01",dd$Sample[1],".pdf")
# pdf(file=pfile,width=11,height=8.5)
# print(pg)
# print(pm)
# dev.off()

# tfile=cc("qcTbl","Threshold","Q",dd$Sample[1],".csv")
# write_csv(qTbl,tfile)
