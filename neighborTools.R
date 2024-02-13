require(SearchTrees)

distGeomType <- function(g1,g2) {
    as.numeric(sqrt((g1$x0-g2$x0)^2+(g1$y0-g2$y0)^2))
}

getNeighborhoodsForFOV<-function(gdi,radius) {

    geom=gdi %>% filter(!Exclude)
    cellGeom=cbind(X=geom$x0,Y=geom$y0)
    nnTree=createTree(cellGeom)

    dx=abs(max(gdi$x0)-min(gdi$x0))
    dy=abs(max(gdi$y0)-min(gdi$y0))
    numNeighborCells=radius**2/(dx*dy/nrow(gdi))
    numNeighborCells=5*ceiling(numNeighborCells/100)*100

    cat("numNeighborCells=",numNeighborCells,"\n")

    #
    # kNN crashes hard if you try to find more neighbors then available cells
    #
    if(nrow(geom)<numNeighborCells) {
        cat("\n\n  WARNING estimated number of neighbors =",numNeighborCells,"too large\n")
        cat("  Only",nrow(geom),"cell in FOV",gdi$FOV[1],gdi$Sample[1],"\n\n")
        numNeighborCells=nrow(geom)-1
    }

    neighborTable=knnLookup(nnTree,geom$x0,geom$y0,k=numNeighborCells)
    yy=map(
            seq(nrow(neighborTable)),
            getNeighborTableForCell,
            neighborTable,
            geom
        ) %>%
        bind_rows %>%
        arrange(C.UUID,Dij.micron) %>%
        mutate(SPOT=gdi$SPOT[1]) %>%
        select(SPOT,C.UUID,N.UUID,Dij.micron)

    maxR=yy %>% group_by(C.UUID) %>% summarize(maxR=max(Dij.micron))

    cat("\n\n  ",mean(maxR$maxR<radius)*100,"PCT of cells with maxR<radius,\n")

    list(
        nnTbl=yy %>% filter(Dij.micron<=radius),
        maxR=maxR,
        FOVArea=(dx*dx),
        cellDensity=(dx*dy/nrow(gdi)),
        numNeighborCells=numNeighborCells,
        numMaxRltR=sum(maxR$maxR<radius)
        )

}

getNeighborTableForCell <- function(iCenter,knnTbl,cg) {

    nCells=purrr::transpose(cg[knnTbl[iCenter,],] %>% select(UUID,x0,y0))
    cCell=purrr::transpose(cg[iCenter,] %>% select(UUID,x0,y0))[[1]]

    dij=lapply(nCells,function(x){distGeomType(x,cCell)}) %>% unlist
    nij=knnTbl[iCenter,]

    tibble(
            Center=iCenter,
            C.UUID=cCell$UUID,
            Neighbor=nij,
            N.UUID=cg$UUID[nij],
            Dij.micron=round(dij,6)
            ) %>% filter(C.UUID!=N.UUID)

}
