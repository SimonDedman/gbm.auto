mygrids <- gbm.auto::grids # load grids
mysamples <- gbm.auto::samples # load samples
setwd("/home/simon/Desktop/gbm temp/sepParamTestbench/")
source('~/Dropbox/Galway/Analysis/R/gbm.auto/R/gbm.utils.R')
library("shapefiles")
Crop_Map <- read.shapefile("/home/simon/Desktop/gbm temp/CroppedMap/Crop_Map")
source('~/Dropbox/Galway/Analysis/R/gbm.auto/Gbm.auto_extras/gbm.auto.binGausSepParams.R')

testlist <- list(c(1,2),c(3,4))
testlist
testlist[[1]]
testlist[[2]]
testved <- c(5,6)
testved
testlist2 <- list(c(7,8),9)
testlist2
testlist2[[1]]
testlist2[[2]]
testlist2 <- list(10,c(11,12))
testlist2
testlist2[[1]]
testlist2[[2]]
testlist3 <- list(c(1,2),c(3,4),c(5,6))
testlist3
testlist3[[1]]
testlist3[[2]]
testlist3[[3]]
if (length(testlist3) > 2) {stop("Only 2 list items allowed: 1 bin 1 Gaus")}

gbm.auto2(samples = mysamples,
          grids = mygrids,
          expvar = c(4:10),
          resvar = 11,
          tc = 2,
          lr = list(c(0.05,0.04),0.03),
          bf = 0.5,
          mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE, BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
