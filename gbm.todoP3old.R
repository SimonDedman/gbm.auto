####Load scripts & data####
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")
source('C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source('C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')

####Load scripts & linux####
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')

####run gbm.auto without E####
####FULL DATASET GBM.AUTO RUNS FOR FUTURE REFERENCE####
# Load data
mysamples<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Data/Samples_allRays_Env_F_E.csv", header = TRUE, row.names=NULL)
mygrids<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", header = TRUE)

# Load linux
mysamples<-read.csv("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Data/Samples_allRays_Env_F_E.csv", header = TRUE, row.names=NULL)
mygrids<-read.csv("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", header = TRUE)

# set directory, without fishing E as an expvar
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E")
setwd("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E")
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

####SINGLE LINE FULL DATASET GBM.AUTO RUN####
#all at once
gbm.auto(expvar=c(4:9),resvar=c(12:15),grids=mygrids,samples=mysamples,tc=c(2,6),lr=c(0.005, 0.001),bf=c(0.5),
         gridslat = 2,gridslon = 1,ZI = TRUE,map = TRUE,RSB= TRUE)

####everything below now defunct####





####run gbm.auto with E####
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E")
setwd("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E")
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
# with E
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/ValueMaps")
preddata <- read.csv(file="C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/AllPredsWithE_F_E.csv")

# without E
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E/ValueMaps")
preddata <- read.csv(file="C:/Users/Simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E/AllPredsWithE_F_E.csv")

#linux with E
setwd("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/ValueMaps")
preddata <- read.csv(file="/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/AllPredsWithE_F_E.csv")

#linux without E
setwd("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E/ValueMaps")
preddata <- read.csv(file="/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/Without E/AllPredsWithE_F_E.csv")

# Cuckoo
gbm.valuemap(data=preddata,
             scalerange=c(3,8),
             goodcols=c(3),
             badcols=c(8),
             goodname="Cuckoo66scaled",
             #goodweight = c(3),
             #badweight = c(1),
             steps=10,
             limitline = 2000,
             #savethis = c("data","report"),
             savethis = NULL,
             #plotthis = c("good","bad","both","line","fish","fishermen"),
             plotthis = c("both"),
             #plotthis = NULL,
             latcolno = 1,
             loncolno = 2,
             m = 1,
             zero = FALSE)
             
# fishermen map takes ages, but now works. Weird that it takes so long...
# 34 effort map (bad) inverted. Put inversion lower down. Only good/bad/both, no line, no "fish"?
# 35 effort map same. All same. Good fine. Bad wrong, scaling bad, somehow inverted.
# Bothdata using same name as baddata, overwriting it. m=0.001 + scaling a problem?
# combo map: black overwriting 1st layer? 1st layer not working in any case.
# Was using a bad mix of data[,] and bothdatadf[,], fixd, black bit commented out
# line currently commented out until I work it out!
# 36 goodmap works. badmap works. bothmap seems to be subtractive. 
# Error in .subset(x, j) : only 0's may be mixed with negative subscripts
# related to baddata: turned scaling off, inversion process breaks? Built for 0:1, is now 0:100000 or whatever.
# change 1s to max(badddata)? testbench it. Fixed.
# 37 good bad both works.
# Error in `[.data.frame`(data, , bothdata) : undefined columns selected 
# bothdata's a vector unatached to data. Needs to be linked else how to know where to plot the points?
# bothdata = gooddata+baddata. Both good & bad are in original data[,] order so both should be too?
# Order doesn't happen til line 244 
# 38 what's the difference between both & fish? Fish is both + closed area, just closed area inactive at moment.
# good bad both combo work. line inactive. closed area inactive. Activate it now.
# 39 Error in xj[i] : only 0's may be mixed with negative subscripts
# closed area z problem. Test with dud z to make sure it IS that
# 40 still fails... could it be x? global assign bothdatadf to view
# 41 ok, problrem is limitline value is 0.5 you fucking goon. So it stops at 1 point. What SHOULD it be?
# Should it be sorted by good+inverted_bad but adding good, i.e. sum the Bs as you go down the list of best combo site?
# need to create df with lon lat goodcol combocol DONE, sort by combocol DONE, add by goodcol DONE. Globally assigned bothdatadf twice
# changed liimtline to 15000. Guess. Update after checking global assigns.
# 42 good bad both combo work. bothdatapreloop looks fine. Goodcol order v similar to bothdata as expected
# bothdatadf looks good, 4221 rows (out of 378570). Now map it. Removed preloop global
# 43 combo: whole thing is black. also white lines issue. Lower limitline to 1500 to test
# 44 still bad. Scale goes to 1 only... check. You're still using the z=1 temp you doughnut. Changed.
# 45 black legend scales to 46 but whole thing still black. Limitline 1500 -> 500
# 46 still all black. 500 -> 100
# 47 same. 100 -> 50
# 48 same. Still out of 47 though... change to 20
# 49 still. 5
# 50 same. Try something else. bothdatafr[,4] is now only 1 row. Try mapping on its own
# 51 mapped together still... Comment out map2
# 52 STILL BLACK! HOW?!
# 53 same
# 54 same. manually set heatcol to colourramp
# 55 STILL black!! Moved commented section to another sheet. Is it running commented out sections?!
# 56 no change. And still globally assigning bothdatadf! Restarting.
# 57 works. fucking finally. So how/why did this occur?! Conclusion has to be that it was running commented code. Remove heatcol
# 58 fine. Uncommenting next block, with heatcol commented out.
# 59 closed area effect not noticeable. Earlier dev.off kills it? Removed. Black commented out still.
# 60 no noticeable change. Move closed area to separate image to test. YOU HAD 'SOURCE ON SAVE' OFF YOU FUCKING MARSHMALLOW!!
# 61 closed areas generates points but it's area is constrained to the limited points i.e. is smaller. Also no coast: middle of IS, coast out of range.
# can use points() ? Before I do, try adding them together again. Changed limitline to 2000, heatcol to black.
# 62 black works; underlying map is deleted. Added min & max to stretch closed areas map to full extent
# 63 extents did nothing. Changed from NA to 0
# 64 works. try par(new=TRUE) to overlay
# 65 nope.
# 66 dunno what I changed. Potentially manually drew in the blacks w/ GIMP colour replace. Need to run gbm.map elements separately, see Hans' email.

goodweight = NULL
if(!is.null(goodweight)) print("shouldn't evaluate")

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
  