# Script to Run/test gbm/basemap
source('~/Dropbox/Galway/Analysis/R/gbm.auto/gbm.basemap.R')
setwd("~/Desktop")
setwd("GSHHS_shp")
getwd()
gbm.basemap(getzip = FALSE,
            bounds = c(-122.586557, -122.196168, 37.393566, 37.659590),
            res = 1, savename = "SF_Bay")


####
bounds = c(range(grids[,1]),range(grids[,2])) # defaults for gbm.map

# could instead have option in gbm.map to set shape as "CALC", or even default to CALC
# whereby it runs gbm.basemap using bounds derived from:
bounds = c(range(x),range(y))