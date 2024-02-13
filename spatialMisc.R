findCleanRects<-function(gd,bb) {
    cat(bb,"\n")
    drawRect(bb*ps,col=1)
    area=getRectArea(bb)
    cat("Area=",getRectArea(bb),"\n")
    if(area<.1) return(bb)

    gd1=gd %>% filter(X>bb[1] & X<bb[3] & Y>bb[2] & Y<bb[4])
    if(nrow(gd1)<1) {
        return(bb)
    }

    exclusions=gd1 %>% filter(Exclude)

    if(nrow(gd1)==nrow(exclusions)) {
        drawRect(bb*ps,col=2,lwd=2)
        return(NULL)
    }

    if(nrow(exclusions)>0) {
        if(bb[4]-bb[2]>bb[3]-bb[1]) {
            mid=(bb[4]+bb[2])/2
            b1=bb;b2=bb
            b1[4]=mid
            b2[2]=mid
        } else {
            mid=(bb[3]+bb[1])/2
            b1=bb;b2=bb
            b1[3]=mid
            b2[1]=mid
        }
        return(rbind(findCleanRects(gd,b1),findCleanRects(gd,b2)))
    } else {
        return(bb)
    }
}

getBoundaryPolygon.0 <- function(gd,boundaryLength) {
    ps=gd$PixelScaling[1]
    bb=(
        c(min(gd$XMin),min(gd$YMin),max(gd$XMax),max(gd$YMax)) +
        c(1,1,-1,-1)*boundaryLength/ps
    )

    bb=matrix(c(bb[1],bb[2],bb[1],bb[4],bb[3],bb[4],bb[3],bb[2],bb[1],bb[2]),ncol=2,byrow=T)
    colnames(bb)=c("X","Y")

    data.frame(bb)

}

