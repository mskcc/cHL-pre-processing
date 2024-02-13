require(tidyverse)
require(openxlsx)

atlas=readRDS("cellAtlas_Hodgkins_v10.9__220718_b_CTD_V1___xx_20220718_000000.rda")
tbl1=atlas %>% filter(CellType=="T_All") %>% group_by(Sample,SPOT) %>% mutate(N=n()) %>% group_by(Sample,SPOT,N) %>% count(SubType) %>% mutate(PCT=n/N) %>% select(-n) %>% spread(SubType,PCT,fill=0)

write.xlsx(list(T_All=tbl1),"qcCellType.xlsx")
