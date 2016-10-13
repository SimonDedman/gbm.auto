# gbm.auto variance idea
# variance: just need to re-run it to different folders, so feasibly just have a loop with
# all unecessary stuff switched off and a new folder each time?

library("gbm.auto")
mygrids <- gbm.auto::grids # load grids
mysamples <- gbm.auto::samples # load samples
setwd("/home/simon/Desktop/gbm temp/variance/")
source('~/Dropbox/Galway/Analysis/R/gbm.auto/R/gbm.utils.R')
# library("shapefiles")
# Crop_Map <- read.shapefile("/home/simon/Desktop/gbm temp/CroppedMap/Crop_Map")
#source('~/Dropbox/Galway/Analysis/R/gbm.auto/Gbm.auto_extras/gbm.auto.binGausSepParams.R')

gbmlooptest <- gbm.loop(loops = 2,
         savecsv = T, #unnecessary
         samples = mysamples,
         grids = mygrids,
         expvar = c(4:10),
         resvar = 11,
         simp = FALSE,
         RSB = F)
