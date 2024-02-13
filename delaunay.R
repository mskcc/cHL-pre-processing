delaunay_nearest_neighbors<-function(oo) {

    edges_by_spot=oo$geom.data %>%
        select(UUID,SPOT,x0,y0) %>%
        group_split(SPOT) %>%
        map(delaunay_edges_by_spot)

    nn_by_spot=oo$neighbors %>% group_split(SPOT)

    DAT=transpose(list(edges_by_spot,nn_by_spot))

    if(PARALLEL) {
        ll=mclapply(DAT,assign_nearest_neighbors,mc.cores=nCores)
    } else {
        ll=  lapply(DAT,assign_nearest_neighbors)
    }

    bind_rows(ll)

}

delaunay_edges_by_spot<-function(geom_data_by_spot) {

    #
    # D=T ==> Delaunay Triangulation
    # S=F ==> No Steiner Points on Boundry
    #
    rt=triangulate(pslg(select(geom_data_by_spot,x0,y0)),D=T,S=F)

    edges=bind_rows(
        tibble(
            C1=rt$E[,1],C1.UUID=geom_data_by_spot$UUID[rt$E[,1]],
            C2=rt$E[,2],C2.UUID=geom_data_by_spot$UUID[rt$E[,2]]
        ),
        tibble(
            C1=rt$E[,2],C1.UUID=geom_data_by_spot$UUID[rt$E[,2]],
            C2=rt$E[,1],C2.UUID=geom_data_by_spot$UUID[rt$E[,1]]
        )
    )

    edges

}

assign_nearest_neighbors<-function(args) {
    #
    # edges=args[[1]]
    # neighbors=args[[2]]

    cat("assign_nearest_neighbors",args[[2]]$SPOT[1],"\n")

    edge_UUIDs=args[[1]] %>% mutate(EDGE.UUID=cc(C1.UUID,C2.UUID)) %>% pull(EDGE.UUID)

    args[[2]] %>%
        mutate(E.UUID=cc(C.UUID,N.UUID)) %>%
        mutate(DelaunayNeighbor=E.UUID %in% edge_UUIDs) %>%
        select(-E.UUID)

}


delaunay_second_neighbors<-function(nn) {

    ns=nn %>% group_split(SPOT)

    if(PARALLEL) {
        ll=mclapply(ns,delaunay_second_neighbors_by_spot,mc.cores=nCores)
    } else {
        ll=  lapply(ns,delaunay_second_neighbors_by_spot)
    }

    bind_rows(ll)

}

delaunay_second_neighbors_by_spot<-function(ns) {

    cat("delaunay_second_neighbors_by_spot",ns$SPOT[1],"\n")

    d2Edges=full_join(
                    ns %>% filter(DelaunayNeighbor),
                    ns %>% filter(DelaunayNeighbor),
                    by="N.UUID"
                ) %>%
            filter(C.UUID.x!=C.UUID.y) %>%
            distinct(C.UUID.x,C.UUID.y) %>%
            mutate(E.UUID=cc(C.UUID.x,C.UUID.y)) %>%
            select(E.UUID) %>%
            pull

    ns %>% mutate(D2Neighbors=!DelaunayNeighbor & cc(C.UUID,N.UUID) %in% d2Edges)

}

