require(tidyverse)

args=commandArgs(trailingOnly = T)
P=as.numeric(args[1])
hFile=args[2]

oDir=file.path("SubSample",P)
dir.create(oDir,recursive=T,showWarnings=F)

read_csv(hFile) %>% sample_frac(P) %>% write_csv(file.path(oDir,basename(hFile)))
