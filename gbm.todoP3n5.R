####Load scripts & data####
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")
source('C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source('C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')

####run gbm.auto without E####
# Load data
mysamples<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Data/Samples_allRays_Env_F_E.csv", header = TRUE, row.names=NULL)
mygrids<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", header = TRUE)

# set directory, without fishing E as an expvar
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E")
#cuckoo
gbm.auto(expvar=c(4:9),resvar=c(12),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#thornback
gbm.auto(expvar=c(4:9),resvar=c(13),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#blonde
gbm.auto(expvar=c(4:9),resvar=c(14),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#spotted
gbm.auto(expvar=c(4:9),resvar=c(15),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#all at once
gbm.auto(expvar=c(4:9),resvar=c(12:15),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

####run gbm.auto with E####
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E")
#cuckoo
gbm.auto(expvar=c(4:9,11),resvar=c(12),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#thornback
gbm.auto(expvar=c(4:9,11),resvar=c(13),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#blonde
gbm.auto(expvar=c(4:9,11),resvar=c(14),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#spotted
gbm.auto(expvar=c(4:9,11),resvar=c(15),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

#all at once
gbm.auto(expvar=c(4:9,11),resvar=c(12:15),grids=mygrids,tc=c(2,6),lr=c(0.005, 0.001),ZI = TRUE)


####run valuemaps####
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/ValueMaps")
preddata <- read.csv(file="AllPredsWithE_F_E.csv")

setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E/ValueMaps")
preddata <- read.csv(file="AllPredsWithE_F_E.csv")

# Cuckoo
gbm.valuemap(data=preddata,
             scalerange=c(3,8),
             goodcols=c(3),
             badcols=c(8),
             goodname="Cuckoo",
             #goodweight = c(1),
             #badweight = c(1),
             steps=2,
             latcolno = 1,
             loncolno = 2)


gbm.valuemap(data=preddata,
               #loncolno = 1,
               #latcolno = 2,
               scalerange = c(3,4),  # need to scale for the moment since I don't have bpa and therefore the values are all over the place
               goodcols = c(4), # only one: resvar
               badcols = c(3), # only one: E. For this loop could make badcols=badcols since it's in the environment
               #species = "TESTSPECIES", #or set manually. For gbm.auto: names(samples[i])  ? 
# Esteps now created in function
               goodname = "Sept02run03", #not required, has default but test w input
#               goodname = paste("GoodSpecies E_wgt ",rev(which(rev(Esteps) == q))-1,"pct",sep=""),
#               badname = paste("BadSpecies E_wgt ",rev(which(rev(Esteps) == q))-1,"pct",sep=""),
#               bothname = paste("BothSpecies E_wgt ",rev(which(rev(Esteps) == q))-1,"pct",sep=""),
               #legendtitleV = "Conservation Value",
               steps = 2)
               #goodweight = 1, #not weighting resvar in proper run?
               #badweight = q,  #.
               #mapmain = paste("FishingE weighting: ",rev(which(rev(Esteps) == q))-1,"%, Total Bpa: ",bpasum=sum(bothdata),sep=""), # Print bothdata sum on map
  
  
  
  ####running this (in gbm.auto to begin with)####
  source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")  #SD specific, remove @ end
  ####MatF####
  # set wd
  setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/ConsValMaps")
  # load data. Expecting Longitude, Latitude, data columns: predicted abundances, fishing effort etc.. E.g.: Abundance_Preds_All.csv from gbm.auto
  data<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Both_Abundance_Preds_plusE.csv", header = TRUE)
  
  gbm.valuemap(data,
               loncolno = 1,
               latcolno = 2,
               scalerange = c(3:11),
               goodcols = c(4:11),
               badcols = c(3),
               goodweight = c(1.5,1,3,1,1.5,1,3,1), #CTBS mf & j
               badweight = 4,  #fishing 4 times as important, for example.
               species = names(samples[i]),
               ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
  
  
  # preddata[,5] <- preddata[,badcols] - Esteps[q]
  