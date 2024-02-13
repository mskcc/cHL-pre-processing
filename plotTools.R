suppressPackageStartupMessages(require(ggforce))

plotMarkerDist <- function(mm,markers,sampleName="") {

    if(sampleName=="") {
        sampleName = mm$Sample[1]
    }

    pg = filter(mm,Marker %in% markers) %>%
        mutate(Pos=factor(Positive_Classification)) %>%
        ggplot(aes(Marker,TIntensity,color=Pos)) +
        geom_violin() +
        ggtitle(paste("Sample",sampleName)) +
        facet_wrap_paginate(~SPOT,nrow=2,ncol=3)

    pgl=list()
    nPages=n_pages(pg)
    for(ni in seq(nPages)) {
        pgl[[ni]]=pg+facet_wrap_paginate(~SPOT,nrow=2,ncol=3,page=ni)
    }

    pgl

}

convertPNGtoPDF <- function(pfile,keep=F) {
    system2("convert",c(pfile,gsub(".png",".pdf",pfile)))
    if(!keep)
        unlink(pfile)
}

png <- function(pfile,width,height,type="cairo",units="in",pointsize=12,res=150,...) {
    grDevices::png(filename=pfile,width=width,height=height,type="cairo",units="in",pointsize=12,res=150,...)
}

require(patchwork)
paginatePlots<-function(plts,pRows,pCols) {

    nPlots=pRows*pCols

    pp=list()
    page=1
    currPlot=1
    pa=NULL

    for(ii in seq(len(plts))) {
        if(is.null(pa)) {
            pa=plts[[ii]]
        } else {
            pa=pa+plts[[ii]]
            if(len(pa$patches$plots)==(nPlots-1)) {
                pp[[len(pp)+1]]=pa + plot_layout(ncol=pCols,nrow=pRows)
                pa=NULL
            }
        }
    }

    if(!is.null(pa))
        pp[[len(pp)+1]]=pa + plot_layout(ncol=pCols,nrow=pRows)

    pp

}

