####Shiela Heymans' chat @ WKIRISH 14&15/09/2015####
# Added "linesfiles" to gbm.auto to export CSVs of line plots
# WIth the hope that they can be used in EcoPath w/ EcoSim
# If it works:
#1. Need to write code for processing DATRAS files, potentially set lat/long extent
#   then filter by species, year?; then vlookup and pivot table to end up with:
#   site/lat/long/explanatory variables/response variable(s)
#   ideally 1 sheet, lat/long down all sites, species as columns across.
#   assumes summed over X years, but used can set that.
#   see /home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/DATRAS/DATRAS All Combined/AllSurveysCombined.xlsx
# haul: datras "tab" 1, concatenate Survey/StNo/HaulNo/Year, need HaulDur, ShootLat, ShootLon, HaulLat, HaulLon, StatRec?, Depth?
# sp len wgt: datras "tab" 2, concatenate Survey/StNo/HaulNo/Year, need SpecCode, Latin Name (I looked up?), Sex, TotalNo, NoMeas?, SubWgt?, CatCatchWgt, LngtClass, HLNoAtLngt
# sp len sex age mat wgt: datras "tab" 3, concatenate Survey/StNo/HaulNo/Year, need SpecCode, Latin Name (I looked up?), LngtClass, Sex, Maturity, Age, NoAtLngt, IndWgt
#
#
#2. Also need to combine bin & gaus linesfiles CSVs somehow. Multiply? Mean?

#prep
# set & test working directory
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles")
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles")
getwd()

# load gmb.utils, which contains various of Elith's packages not contained in dismo
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")

# load gbm functions
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R")

source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R")

# load samples
mysamples<-read.csv("Hauls&J&Preds&Enviros_Trimmed_ISonly.csv", header = TRUE, row.names=NULL)

setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/testbench")

gbm.auto(expvar=c(4:11,15,17,21,25,29,33,37),
         resvar=c(44),
         grids=mygrids,
         samples=mysamples,
         #tc=c(2,5,15),
         lr=c(0.005, 0.001),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = FALSE,
         RSB= FALSE,
         BnW = FALSE,
         savegbm = FALSE,
         VARINT = FALSE,
         linesfiles = TRUE,
         legendtitle = "CPUE")