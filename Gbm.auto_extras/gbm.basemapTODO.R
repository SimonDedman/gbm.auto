# gbm.basemap todo
grids <- read.csv("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", header = TRUE)
gridslat = 2
gridslon = 1
xbounds <- range(grids[,gridslon], na.rm = TRUE)
ybounds <- range(grids[,gridslat], na.rm = TRUE)
bounds <- c(xbounds, ybounds)

####TODO####
# Make it work
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds) # Error in res == 1 :
# comparison (1) is possible only for atomic and list types
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = 1)
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = 2)
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = 3)
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = 4)
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = 5)
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = "c")
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = "l")
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = "i")
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = "h")
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = "f")
mymap <- gbm.basemap(getzip = "GSHHS_shp", bounds = bounds, res = "CALC")

# Can I unpack only the relevant folder? Or at least only GSHGG?
setwd("/home/simon/Desktop/")
ifelse(getzip == TRUE, { # download & unzip GSHGG if getzip = TRUE
  download.file(paste("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-", zipvers, ".zip", sep = ""), "GSHHG.zip")
  unzip("GSHHG.zip")
  unzip("GSHHG.zip", files = "GSHHS_shp") #fails
  unzip("GSHHG.zip", files = "LICENSE.TXT") #works
  setwd("./GSHHS_shp")}
  , setwd(getzip)) # else just setwd to there
####end todo####