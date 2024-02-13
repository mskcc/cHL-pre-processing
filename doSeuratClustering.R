if(R.Version()$major<4) {
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}

require(tidyverse)
require(Seurat)
require(fs)
require(yaml)

#NCells=1000000
args=commandArgs(trailing=T)
NCells=as.numeric(args[1])
cat("\n\nNCells =",NCells,"\n")

atlas=readRDS(dir_ls(regex="cellAtlas_Hodgkins.*rda")) 

markers=read_yaml("data/meta/HodgkinLymphoma_Markers__V0/HodgkinLymphoma_Markers__V0___Sheet1.yaml") %>% bind_rows
cMarkers=markers %>% filter(Identity=="Y") %>% pull(Marker_name)

#uu=readRDS(dir_ls(regex="umap_normFile_.*Dn.*rda"))
uu=readRDS(dir_ls(regex="umap_normFile_.*ALL.*rda"))

umap=uu$uu$embedding
colnames(umap)=c("V1","V2")

dd=data.frame(uu$mm) %>% 
    rownames_to_column("UUID") %>% 
    bind_cols(as_tibble(umap)) %>% 
    left_join(atlas) %>% 
    tibble

dups=duplicated(uu$mm)
dd=dd[!dups,]

dd=dd %>% filter(CellType!="superNeg")

if(nrow(dd)>NCells) {
    cat("\n\n")
    cat("   Downsampling from",nrow(dd),"to",NCells,"\n\n")
    ds=dd %>% sample_n(NCells)
} else {
    cat("\n\n")
    cat("   Running all cells ",nrow(dd),"\n\n")
    ds=dd
}

xx=ds %>% 
    select(UUID,all_of(cMarkers)) %>% 
    column_to_rownames("UUID") %>% 
    as.matrix

# need to increase future.globals.maxSize. Was getting error
#     Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) :
#       The total size of the 7 globals exported for future expression (‘FUN()’) is 1.94 GiB..
#       This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
#       The three largest globals are ‘query’ (1.94 GiB of class ‘numeric’), ‘index’
#
# Setting to +Inf disables checking
#    https://www.rdocumentation.org/packages/future/versions/1.22.1/topics/future.options

options(future.globals.maxSize=Inf)

fn=FindNeighbors(xx)

md=ds %>% select(-all_of(colnames(uu$mm)))

saveRDS(list(fn=fn,xx=xx,md=md,NCells=NCells),cc("seuratClusterV4",NCells,"FindNeigh.rda"),compress=T)

resolutions=c(0.1,0.2,0.4,0.8,1.2)

# PARALLEL=FALSE
# if(PARALLEL) {

#     require(furrr)
#     require(purrr)
#     nCores=floor(availableCores()/4)
#     parallel_map=future_map
#     plan(multicore, workers = nCores)
#     cat("Running in PARALLEL nCores=",nCores,"\n")

#     fc=parallel_map(resolutions,~FindClusters(fn$snn,resolution=.),.options = furrr_options(seed = TRUE))

# } else {

fc=list()
for(ri in resolutions) {
    cat("running ri =",ri,"...")
    fc[[as.character(ri)]]=FindClusters(fn$snn,resolution=ri)
    cat("\n")
}

#}

fc=bind_cols(fc)

md=ds %>% select(-all_of(colnames(uu$mm))) 

clust=fc %>% rownames_to_column("UUID") %>% left_join(md) %>% tibble

write_csv(clust,file=cc("seuratClusterV4",NCells,".csv.gz"))
