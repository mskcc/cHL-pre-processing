getRanksPerMarkerByFOV<-function(oo) {
    oo$marker.data %>%
        left_join(oo$geom.data,by="UUID") %>%
        filter(!Exclude) %>%
        filter(!is.na(TIntensity)) %>%
        group_by(SPOT,Marker) %>%
        mutate(R=rank(TIntensity),Pr=(R+1)/(n()+2)) %>%
        ungroup %>%
        select(UUID,Marker,R,Pr)
}
