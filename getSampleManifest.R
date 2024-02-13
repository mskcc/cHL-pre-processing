require(tidyverse)
require(openxlsx)
args=commandArgs(trailing=T)

getLineCounts<-function(ff) {
    read_csv(ff) %>%
        rename_all(~gsub(" ","_",.)) %>%
        select(Image_File_Name) %>%
        mutate(STag=gsub("_Spot.*","",Image_File_Name)) %>%
        mutate(Spot=as.numeric(gsub(".*Spot","",Image_File_Name) %>% gsub(".afi","",.))) %>%
        count(Image_File_Name, STag, Spot)
}

xx=map(args,getLineCounts)
tbl=xx %>% bind_rows %>% select(-Image_File_Name) %>% spread(Spot,n) %>% arrange(STag)
write.xlsx(tbl,"sampleManifest.xlsx")
