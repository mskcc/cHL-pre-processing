args=commandArgs(trailing=T)

suppressPackageStartupMessages({
    require(yaml)
    require(purrr)
    require(readxl)
})

xls2YAML <- function(xfile) {
    #cat("processing",basename(xfile))

    ydir=basename(xfile) %>% gsub(".xls.*","",.)
    dir.create(ydir,showWarnings=F)

    sheets=excel_sheets(xfile)
    for(si in sheets) {
        #cat(" sheet",si)
        xx=read_xlsx(xfile,sheet=si)
        ofile=file.path(ydir,gsub(".xlsx",paste0("___",make.names(si),".yaml"),basename(xfile)))
        write_yaml(xx,ofile,column.major=F)
    }
    #cat("\n")
}


dum=map(args,xls2YAML)


