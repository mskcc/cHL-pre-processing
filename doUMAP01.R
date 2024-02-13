cat("TimeStamp:Start",Sys.time(),"\n\n")
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
            SEED=42,
            DOWNSAMPLE=NULL,N_THREADS=12,
            MARKERS="TOTAL",
            REMOVE_UNKNOWNS=TRUE,
            CELLTYPES="ALL",
            CACHE=NULL
        )

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args$n_neighbors=as.numeric(args$n_neighbors)
args$min_dist=as.numeric(args$min_dist)
args$spread=as.numeric(args$spread)
args$SEED=as.numeric(args$SEED)
args$DEBUG=as.logical(args$DEBUG)
args$REMOVE_UNKNOWNS=as.logical(args$REMOVE_UNKNOWNS)
if(!is.null(args$DOWNSAMPLE)) args$DOWNSAMPLE=as.numeric(args$DOWNSAMPLE)

if(!is.null(args$a)) {
    args$a=as.numeric(args$a)
    args$b=as.numeric(args$b)
}
pArgs=grep("=",cArgs,invert=T,value=T)
normIntenFile=pArgs[1]
atlasFile=pArgs[2]

################################################################################
#
# Set random seed

set.seed(args$SEED)


# ###############################################################################
# defaultUMAP=list(
#   n_neighbors = 15,
#   n_components = 2,
#   metric = "euclidean",
#   n_epochs = NULL,
#   learning_rate = 1,
#   scale = FALSE,
#   init = "spectral",
#   init_sdev = NULL,
#   spread = 1,
#   min_dist = 0.01,
#   set_op_mix_ratio = 1,
#   local_connectivity = 1,
#   bandwidth = 1,
#   repulsion_strength = 1,
#   negative_sample_rate = 5,
#   a = NULL,
#   b = NULL,
#   nn_method = NULL,
#   n_trees = 50,
#   search_k = 2 * n_neighbors * n_trees,
#   approx_pow = FALSE,
#   y = NULL,
#   target_n_neighbors = n_neighbors,
#   target_metric = "euclidean",
#   target_weight = 0.5,
#   pca = NULL,
#   pca_center = TRUE,
#   pcg_rand = TRUE,
#   fast_sgd = FALSE,
#   ret_model = FALSE,
#   ret_nn = FALSE,
#   ret_extra = c(),
#   n_threads = NULL,
#   n_sgd_threads = 0,
#   grain_size = 1,
#   tmpdir = tempdir(),
#   verbose = getOption("verbose", TRUE)
# )
# ###############################################################################
# # S3 method for Seurat
# seurateUMAP=list(
#   n_neighbors = 30L,
#   metric = "cosine",
#   min_dist = 0.3,
#   paramTag="seurate"
# )
# ###############################################################################

###############################################################################

suppressPackageStartupMessages({
    require(tidyverse)
    require(uwot)
    require(yaml)
    require(gplots)
    require(parallel)
})

source("HaloX/tools.R")

# MARKERMETADATAFILE="data/meta/HodgkinLymphoma_Markers__V2/HodgkinLymphoma_Markers__V2___Sheet1.yaml"
# markerMetaData=read_yaml(MARKERMETADATAFILE) %>%
#     map(as_tibble) %>%
#     bind_rows

projParams=read_yaml("config/study.yaml")

markerMetaData=read_metadata_as_tibble(projParams,"MARKERS")

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


##############################################################################
#
# Load data
#

if(is.null(args$CACHE)) {
    CELLSET="ALL"
    dd=read_csv(normIntenFile,col_types=cols(.default = "n",UUID=col_character()))

    if(args$REMOVE_UNKNOWNS) {

        if(!exists("atlas")) {
            cat("\nReading atlas file ... ")
            atlas=readRDS(atlasFile)
            cat("done\n")
        }

        known.uuid=atlas %>%
            filter(
                !Exclude &
                Category!="UNKNOWN" &
                Category!="superNeg"
                ) %>% pull(UUID)

        dd=dd %>% filter(UUID %in% known.uuid)
        CELLSET=cc("RmUnk")


    if(args$CELLTYPES!="ALL") {

        if(!exists("atlas")) {
            cat("\nReading atlas file ... ")
            atlas=readRDS(atlasFile)
            cat("done\n")
        }

        aa=atlas %>% filter(UUID %in% dd$UUID)

        if(toupper(args$CELLTYPES)=="IMMUNE") {

            ct.uuid=aa %>% filter(Category=="Immune_All") %>% pull(UUID)

        } else if(toupper(args$CELLTYPES) %in% c("IMMUNEB","IMMUNE.B","IMMUNE-B")) {

            ct.uuid=aa %>% filter(Category=="Immune_All" & Cell_type!="Leukocyte_Other") %>% pull(UUID)

        } else if(toupper(args$CELLTYPES) %in% c("NOLEUK","NO-LEUK")) {

            ct.uuid=aa %>% filter(Cell_type!="Leukocyte_Other") %>% pull(UUID)

        } else if(toupper(args$CELLTYPES)=="HRS") {

            ct.uuid=aa %>% filter(Cell_type=="HRS") %>% pull(UUID)

        }

        dd=dd %>% filter(UUID %in% ct.uuid)
        CELLSET=cc(CELLSET,args$CELLTYPES)

    }



    }

    if(!is.null(args$DOWNSAMPLE)) {

        #
        # For debugging
        #
        #dd=dd %>% sample_frac(1/100)

        if(!exists("atlas")) {
            cat("\nReading atlas file ... ")
            atlas=readRDS(atlasFile)
            cat("done\n")
        }

        aa=atlas %>% filter(UUID %in% dd$UUID)

        Wsamp=aa %>% count(Subtype) %>% arrange(desc(n)) %>% mutate(Ws=1/sqrt(n))

        cat("Downsample by",1/args$DOWNSAMPLE,"\n")

        xx=left_join(aa,Wsamp) %>% mutate(BLOCK=substr(UUID,1,2)) %>% group_split(BLOCK)
        yy=mclapply(xx,function(x){sample_frac(x,size=1/args$DOWNSAMPLE,weight=Ws)},mc.cores=24)
        aw=bind_rows(yy)

        typeCountTbl=count(aw,Subtype) %>% arrange(desc(n)) %>% full_join(Wsamp,by="Subtype")

        dd=dd %>% filter(UUID %in% aw$UUID)

        CELLSET=cc(CELLSET,"Dn.SQRT",args$DOWNSAMPLE)

        write_csv(typeCountTbl,cc("stats_UMAP01_cellTypeCounts",CELLSET,".csv"))

    }

    check=dd %>% select(UUID,all_of(umapMarkers)) %>% filter_all(~is.na(.))
    if(nrow(check)>0) {
        rlang::abort("FATAL: NA's in intensity matrix")
    }

    cacheFile=cc("umapCacheFile",basename(normIntenFile),CELLSET,".rda")
    saveRDS(
        list(
                dd=dd,
                normIntenFile=normIntenFile,
                atlasFile=filter(atlas,UUID %in% dd$UUID),
                CELLSET=CELLSET
            ),
            cacheFile,
            compress=T
            )

} else {

    lc=readRDS(args$CACHE)
    dd=lc$dd
    normIntenFile=lc$normIntenFile
    atlasFile=lc$atlasFile
    CELLSET=lc$CELLSET

}

##############################################################################
##############################################################################

CELLSET=cc("Markers",args$MARKERS,"_",CELLSET)

##############################################################################
if(exists("aa")) rm(aa)
if(exists("atlas")) rm(atlas)
if(exists("xx")) rm(xx)
if(exists("yy")) rm(yy)
if(exists("aw")) rm(aw)
##############################################################################

mm=dd %>% 
    select(UUID,all_of(umapMarkers)) %>% 
    distinct(across(umapMarkers),.keep_all=T) %>% 
    filter_all(~!is.na(.)) %>% 
    as.data.frame %>% 
    column_to_rownames("UUID") %>% 
    as.matrix

#if(interactive()) halt("\n\n\tSTOP\n\n")

ncores=as.numeric(args$N_THREADS)

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
            paste0("Tag-",args$paramTag),
            "args",substr(digest::digest(params),1,8)
            )

saveRDS(obj,paste0(uBase,"___UUObj.rda"),compress=T)
write_yaml(params,paste0(uBase,"___PARAMS.yaml"))

cat("\n\n\tUMAP DONE\n\n")
cat("TimeStamp:End",Sys.time(),"\n\n")
quit()
#rlang::abort("\n\n\tUMAP DONE\n\n")

source("HaloX/tools.R")
source("toolsUMAP.R")

cellTypeTbl=loadCellTypes("data/meta/HodgkinLymphoma_CellTypes__V1/HodgkinLymphoma_CellTypes__V1___Sheet1.yaml")

rownames(uu)=rownames(mm)

dg=as.data.frame(uu) %>%
    rownames_to_column("UUID") %>%
    as_tibble %>%
    left_join(atlas) %>%
    left_join(cellTypeTbl,by=c(CellSubType="Subtype"))

umap.var="Subtype_full_name"
colors=read_yaml("hodgkinCellTypeColors_NDS.yaml")
cls=unlist(colors$cell_class_colors)
ncls=names(cls)
ii=grep("^#",cls,invert=T)
cls[ii]=col2hex(cls[ii])
cls=paste0(cls,"0a")
names(cls)=ncls
dg[[umap.var]]=factor(dg[[umap.var]],levels=names(cls))

#pt.size=sqrt(1e4/nrow(dg))
pt.size=sqrt(1e5/nrow(dg))

paramStr=gsub("\n",";",as.yaml(args)) %>% gsub(" ","",.) %>% gsub(";"," ",.)
nCells=nrow(dg)

require(aricode)

nSubTypes=len(unique(dg$Subtype_figures))
kk=kmeans(uu,centers=nSubTypes)
uu.AMI=AMI(kk$cluster,dg$Subtype_figures)

params$metrics=list(AMI=uu.AMI)

X0=min(uu[,1])
YN=max(uu[,2])

pg0=dg %>%
    ggplot(aes_string("V1","V2",color=umap.var)) +
    geom_point(size=pt.size) +
    scale_color_manual(values=cls,na.value="#FEFEFE01",labels=cellTypeTbl$Subtype_figures) +
    theme_bw() +
    guides(color = guide_legend(col=3,override.aes = list(alpha = 1, size=4))) +
    ggtitle(paste(paramStr,"nC =",nCells)) +
    annotate("text",x=X0,y=YN,label=paste("AMI =",sprintf("%.4f",uu.AMI)),size=8, hjust=0)

png(filename=paste0(uBase,"___Plt_",args$paramTag,".png"),type="cairo",units="in",width=17,height=11,pointsize=12,res=150)
print(pg0)
dev.off()

if(nrow(dg)>1e6) {
    dg=dg %>% sample_n(1e6)

    pg0=dg %>%
        ggplot(aes_string("V1","V2",color=umap.var)) +
        geom_point(size=3*pt.size) +
        scale_color_manual(values=cls,na.value="#FEFEFE01",labels=cellTypeTbl$Subtype_figures) +
        theme_bw() +
        guides(color = guide_legend(col=3,override.aes = list(alpha = 1, size=4))) +
        ggtitle(paste(paramStr,"nC =",nCells)) +
        annotate("text",x=X0,y=YN,label=paste("AMI =",sprintf("%.4f",uu.AMI)),size=8, hjust=0)


    png(filename=paste0(uBase,"___Plt1Mb_",args$paramTag,".png"),type="cairo",units="in",width=17,height=11,pointsize=12,res=150)
    print(pg0)
    dev.off()
}

write_yaml(params,paste0(uBase,"___PARAMS.yaml"))
