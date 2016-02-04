# Work out why new (Dec15) maps look weird. Compared to Aug-Oct maps:
# Oct maps: CPUE numbers wrong (too low) - samples dbase wrong, rays caught & shown in DATAS defo absent
 # distributions reasonably as expected from prior knowledge
 # patterning fractured: fewer samples (dec:oct): c:159:146, t:286:193, b:200:47, s:285:151; less trees: 700&1250:250, lower TD correlation: 78:61
# Dec maps: #s "right", patterning normal (diffuse) but
 # blonde line anomaly,
 # distributions omit expected areas & previous hotspots
 # extreme points obscure nuances (relates to 2)

setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/map discrepancy")

# re-run oct data with dec code to see if code changes are involved or if it's just the data.
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R")
library(beepr)
mygrids<-read.csv("JUN_grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA.csv", header = TRUE)
mysamples<-read.csv("AUG_Hauls&J&Preds&Enviros_Trimmed_ISonly.csv", header = TRUE, row.names=NULL)
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/map discrepancy/modeloutput")

gbm.auto(expvar=c(4:11,15,19,23,27,31,35),
         resvar=c(46),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,14),  #15 added 6.12.15 no improvement in output line issue
         lr=c(0.005),  #0.01, 
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         varint = FALSE) #zero removed 6.12.15 no improvement in output line issue

# oct data, dec code -> oct maps. Code is fine. Grids are similar.. but not identical!
# Run oct samples with dec grids
mygrids<-read.csv("DEC_grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA.csv", header = TRUE)
mysamples<-read.csv("AUG_Hauls&J&Preds&Enviros_Trimmed_ISonly.csv", header = TRUE, row.names=NULL)
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/map discrepancy/modeloutput")
# oct samples, dec grids, dec code -> oct maps

# Run dec samples with oct grids
mygrids<-read.csv("JUN_grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA.csv", header = TRUE)
mysamples<-read.csv("DEC_Hauls&J&Preds&Enviros_Trimmed_ISonly.csv", header = TRUE, row.names=NULL)
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/map discrepancy/modeloutput")
# dec samples, jun grids, dec code -> Dec maps

# Dec samples to blame, but why the line anomaly? Data should be MORE right? Look at samples again.
# row 1389 IE-IGFS/MB38/143/2004: ctbs dec: 2,0,0,2, aug: 2,2,16,34 = aug has samples in places dec doesn't.
# could dec data be in wrong order?? check above row against datras. 13 blondes all >36cm i.e. not juves. Should be 0. Dec's correct

# then plot both samples on map DONE "nonzero blondes"

# then compare histograms of nonzero samples? why does NIGFS/86/24/2011 & 2 others have <1 ray scores? >60 min hauls, 1 ray caught.
 # done, both look fine.


####Manually maps outputs####
grids<-read.csv("Abundance_Preds_only_cheat.csv", header = TRUE)

png(filename = "MapTmp.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
# run gbm.map function with generated parameters
gbm.map(x = grids[,2],
        y = grids[,1],
        z = grids[,3],
        species = "Blonde",
        legendtitle = "CPUE")
# byx, byy, mapmain, heatcol, shape, mapback, landcol, lejback, legendloc, grdfun, zero, quantile, heatcolours, colournumber
dev.off()

png(filename = "MapTmpBnW.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
gbm.map(x = grids[,2],
        y = grids[,1],
        z = grids[,3],
        species = "Blonde",
        legendtitle = "CPUE",
        landcol = grey.colors(1, start=0.8, end=0.8), #light grey. 0=black 1=white
        mapback = "white",
        heatcolours = grey.colors(8, start=1, end=0))
dev.off()