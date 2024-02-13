dropPositiveMarker<-function(mStr,m1) {
    paste0(setdiff(strsplit(mStr,",")[[1]],m1),collapse=",")
}
