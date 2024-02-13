cat("\n   Starting tests ...\n")
suppressPackageStartupMessages(require(testthat))
source("HaloX/readHalo.R")
obj=parse_halo("data/TS_0001.csv.gz")

expect_true(all(c("geom.data","marker.data") %in% names(obj)))

req.cols.geom.data=c("Image_File_Name", "Object_Id", "XMin", "XMax", "YMin", "YMax",
"Cell_Area", "Cytoplasm_Area", "Nucleus_Area", "Nucleus_Perimeter",
"Nucleus_Roundness", "SPOT", "UUID")

req.cols.marker.data=c("UUID", "Marker", "Positive_Classification")

expect_true(all(req.cols.geom.data %in% colnames(obj$geom.data)))
expect_true(all(req.cols.marker.data %in% colnames(obj$marker.data)))

sig.test.geom.data="17a9ed5fee618bfc1180755e3699602c"
sig.test.marker.data="0a409658927c6e0361dfe9f02a99ac71"

sig.geom.data=obj$geom.data %>%
    arrange(Object_Id) %>%
    select(Image_File_Name,Object_Id,XMin,XMax,YMin,YMax,SPOT,UUID) %>%
    digest::digest(.)

sig.marker.data=obj$marker.data %>%
    arrange(UUID,Marker) %>%
    select(UUID,Marker,Positive_Classification) %>%
    digest::digest(.)

expect_equal(sig.test.geom.data,sig.geom.data)
expect_equal(sig.test.marker.data,sig.marker.data)

cat("   DONE\n\n")

