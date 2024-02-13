plotFitTest <- function(pDat,lORCut) {

  cls=RColorBrewer::brewer.pal(5,"PRGn")
  pgf1=pDat %>% ggplot(aes(CD20__TIntensity,CD3__TIntensity,color=ifelse(abs(lOR)<2,NA,lOR))) +
    coord_fixed() +
    geom_point(alpha=0.5) +
    geom_abline(slope=1,intercept=0,alpha=.1,size=3) +
    scale_color_gradientn(colors=cls,name=expression(log[10](OR)))

  pgf1

}

suppressPackageStartupMessages(require(tidyverse))

theme_set(theme_light(base_size=16))

args=commandArgs(trailing=T)
sampleObjFile=args[1]

cellAtlasFile=grep("cellAtlas.*.rda$",dir(),value=T)
atlas=readRDS(cellAtlasFile)

#sampleObjFile="rda/v5_Exclusion/H16_4923___HaloObj_v10.7__210714_Exclusions_.rda"
#cellAtlasFile="cellAtlas_Hodgkins_v10.7__210714_b_CTDv____ba5cf549_20210718_195242.rda"

#sampleObjFile="../HodgkinsV10c/rda/v5_Exclusion/H16_4923___HaloObj_v10.7__210714_Exclusions_.rda"
#cellAtlasFile="../HodgkinsV10c/cellAtlas_Hodgkins_v10.7__210714_b_CTDv____c19fcbab_20210718_182340.rda"

mA="CD3"
mB="CD20"

#
# Load object data & get markerTbl
# - Filter out excluded data
# - join cell data

obj=readRDS(sampleObjFile)

sampleName=paste0(obj$sample.data$CellDive_ID," / ",obj$sample.data$Patient_ID)
sampleId=obj$sample.data$CellDive_ID

if(!"Exclude" %in% colnames(obj$geom.data)) {

  cat("\n\n  ERROR: Not implemented; Exclude column missing\n")
  cat("  You should be using the HaloObjects Files with Exclude column set\n\n")
  stop("Not implemented")

}


uuids=obj$geom.data %>% filter(!Exclude) %>% pull(UUID)

geom.data = obj$geom.data %>% filter(UUID %in% uuids)
marker.data = obj$marker.data %>% filter(UUID %in% uuids) %>%
  left_join(select(geom.data,UUID,Sample,SPOT))

source("HaloX/plotTools.R")
pdf(file=cc("dblPosMarkerDist",sampleId,"_",mA,mB,".pdf"),width=11,height=8.5)
ppg=plotMarkerDist(marker.data,c(mA,mB),sampleName)
print(ppg)
dev.off()

mTbl=marker.data %>%
  filter(Marker %in% c(mA,mB)) %>%
  mutate(Marker=factor(Marker,levels=c(mA,mB))) %>%
  mutate(SPOT=factor(SPOT)) %>%
  select(UUID,Marker,Positive_Classification,TIntensity,Sample,SPOT) %>%
  left_join(atlas %>% select(UUID,CellType,PosMarkers)) %>%
  gather(Metric,Value,Positive_Classification,TIntensity) %>%
  unite(MarkerMetric,Marker,Metric,sep="__") %>%
  spread(MarkerMetric,Value)
mTbl$DPos=cc(mTbl[[paste0(mA,"__Positive_Classification")]],mTbl[[paste0(mB,"__Positive_Classification")]])

if(1) {
  source("HUtils/loadColorPalettes.R")
  pg=mTbl %>% ggplot(aes(CD3__TIntensity,CD20__TIntensity,color=CellType)) + geom_point(alpha=.75) + facet_wrap(~DPos) + scale_color_manual(values=hColors$cellTypes) + ggtitle(paste("Sample",sampleName))
  pfile=cc("dblPosMarker2DDensity",sampleId,"_",mA,mB,".png")
  png(pfile,width=14,height=11)
  print(pg)
  dev.off()
  convertPNGtoPDF(pfile)
}

mTbl=split(mTbl,mTbl$SPOT)

source("HaloX/modelTools.R")

#
# Need to split by spots
#

#
# abs(lOR) needs to be great then this to reassign
#
lORCut=2

replacePos=list()
pngFiles=list()
statsL=list()

for(si in seq(mTbl)) {

  cat("FOV",names(mTbl)[si],"si =",si,"of",len(mTbl),"\n")

  mDat=mTbl[[si]]
  posN=count(mDat,DPos) %>% spread(DPos,n)

  training.dat=mDat %>% filter(CellType!="UNKNOWN" & DPos != "0_0")
  posN.train=count(training.dat,DPos) %>% spread(DPos,n)


  stats=tibble(
    SampleId=sampleId,
    FOV=mDat$SPOT[1],
    NumCells=nrow(mDat),
    N.DblPos=sum(mDat$DPos=="1_1"),
    PCT.DblPos=mean(mDat$DPos=="1_1"),
    N.TraininSet=nrow(training.dat),
    lORCut=lORCut)

  stats[[paste0("PCT.",mA)]]=posN[["1_0"]]/(posN[["1_0"]]+posN[["0_1"]])

  if( nrow(training.dat)<500 | !all(dim(posN.train)==c(1,2)) ) {
    cat("WARNING:: Not enough training data for FOV",mDat$SPOT[1],"Sample",mDat$Sample[1],"\n")
    statsL[[len(statsL)+1]]=stats

    next
  }



  stats[[paste0("PCT.",mA,".Train")]]=posN.train[["1_0"]]/(posN.train[["1_0"]]+posN.train[["0_1"]])

  modelA=logisticFitter(paste0(c(mA,mB),"__TIntensity"),paste0(mA,"__Positive_Classification"),training.dat)
  modelB=logisticFitter(paste0(c(mA,mB),"__TIntensity"),paste0(mB,"__Positive_Classification"),training.dat)
  converged=modelA$model$converged & modelB$model$converged

  stats$Converged.A=modelA$model$converged
  stats$Converged.B=modelB$model$converged

  stats$Converged=converged

  if(!converged) {
    cat("WARNING:: Models did not converge, can not reassign markers Pos for FOV",
      mDat$SPOT[1],"Sample",mDat$Sample[1],"\n")
    statsL[[len(statsL)+1]]=stats

    next
  }

  pDat=mDat %>%
    mutate(probA=modelA$predict(mDat),probB=modelB$predict(mDat)) %>%
    mutate(lOA=log10(probA/(1-probA)),lOB=log10(probB/(1-probB))) %>%
    mutate(lOR=lOA-lOB) %>%
    rename_at(
      vars(matches("__Positive_Classification")),
      ~gsub("__Positive_Classification","__Positive_Classification.orig",.))

  pgf1=plotFitTest(pDat,lORCut)
  fovNo=pDat$SPOT[1]
  pfile=cc("dblPosMarkerModelTest",sampleId,sprintf("%02d",fovNo),"_",mA,mB,".png")
  png(pfile,width=14,height=11)
  print(pgf1)
  dev.off()
  pngFiles=c(pngFiles,pfile)
  #convertPNGtoPDF(pfile)


  posA=paste0(mA,"__Positive_Classification")
  posB=paste0(mB,"__Positive_Classification")

  pDat[[posB]]=pDat[[paste0(posB,".orig")]]
  pDat[[posA]]=pDat[[paste0(posA,".orig")]]

  pDat[[posA]][pDat$DPos=="1_1" & pDat$lOR < -lORCut]=0
  pDat[[posB]][pDat$DPos=="1_1" & pDat$lOR > lORCut]=0

  pDat$DPos.new=cc(
    pDat[[paste0(mA,"__Positive_Classification")]],
    pDat[[paste0(mB,"__Positive_Classification")]]
    )

  replacePos[[si]]=pDat %>%
    select(UUID,matches("Positive_Classification")) %>%
    gather(MarkerMetric,Value,-UUID) %>%
    separate(MarkerMetric,c("Marker","Metric"),sep="__") %>%
    arrange(UUID) %>%
    spread(Metric,Value)

  count.New=pDat %>% filter(DPos=="1_1") %>% count(DPos.new) %>% spread(DPos.new,n)

  try({stats[[paste0("n.New.",mA)]]=count.New[["1_0"]]})
  try({stats[[paste0("n.New.",mB)]]=count.New[["0_1"]]})
  try({stats[["n.Unresolved"]]=count.New[["1_1"]]})
  try({stats[["PCT.Unresolved"]]=count.New[["1_1"]]/nrow(pDat)})

  posN.new=pDat %>% count(DPos.new) %>% spread(DPos.new,n)
  stats[[paste0("PCT.",mA,".Fix")]]=posN.new[["1_0"]]/(posN.new[["1_0"]]+posN.new[["0_1"]])

  statsL[[len(statsL)+1]]=stats

}

replacePos=bind_rows(replacePos)

if(nrow(replacePos)>0) {
  obj$marker.data=obj$marker.data %>%
    rename(Positive_Classification.orig=Positive_Classification) %>%
    left_join(replacePos) %>%
    mutate(Positive_Classification=ifelse(
                  is.na(Positive_Classification),
                  Positive_Classification.orig,
                  Positive_Classification)
    )
    system2("convert",c(unlist(pngFiles),cc("dblPosMarkerModelTest",sampleId,"","_",mA,mB,".pdf")))
}

gsub(".rda",cc("fixM",mA,mB,".rda"),basename(sampleObjFile))
saveRDS(obj,file=gsub(".rda",cc("fixM",mA,mB,".rda"),basename(sampleObjFile)),compress=T)

write_csv( bind_rows(statsL),gsub(".rda",cc("fixM",mA,mB,".csv"),basename(sampleObjFile)))



