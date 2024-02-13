if(R.Version()$major<4) {
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}

args=commandArgs(trailing=T)
if(len(args)<2) {
    cat("\n")
    cat("    usage: doSeuratClusteringPartII.R neighborObj.rda res1 [res2 ... resN]\n\n")
    quit()
}

neighObjFile=args[1]
resolutions=as.numeric(args[-1])

require(tidyverse)
require(Seurat)

nn=readRDS(neighObjFile)
fn=nn$fn

fc=list()
for(ri in resolutions) {
    cat("running ri =",ri,"...")
    fc[[as.character(ri)]]=FindClusters(fn$snn,resolution=ri)
    cat("\n")
}

fc=bind_cols(fc)

clust=fc %>% rownames_to_column("UUID") %>% left_join(nn$md) %>% tibble

rTag=substr(digest::digest(resolutions),1,8)
write_csv(clust,file=cc("seuratClusterV4ii",nn$NCells,rTag,".csv.gz"))
