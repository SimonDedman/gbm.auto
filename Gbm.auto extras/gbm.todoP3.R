####Script to Run P3 analysis:
# valuemaps, scaled data, weighted scaled data, scales (&weighted?) + inverted effort combo data, maps & closed areas
# test

####todo####
# Label perspeciesareamap legend items w/ species names not 1:4 numbers (L387 & 402)
# legend <- c("Open", goodname), in legend.grid(), L61 of gbm.map. Hmm...
# gbm.map allows ... optionals to be passed to it
# but legend.grid is currently sealed
# ok end of gbm.map L61 add "..." to see if I'm allowed to add that to legend.grid DONE OK
# then run something which will call map & see if it breaks, e.g. baddata DONE OK
# if not, then try running baddata again with legend specified & see if it changes. valuemap L69
# fails: if I pass 'legend' to gbm.map then in legend.grid it doesn't know wether to interpret the legend as legend the function
# or legend the parameter of legend the function. UNBELIEVABLY stupid parameter naming by legend function builders.
# Asked hans c("Open", goodname, recursive = TRUE) for later, if it's possible
# decompile hans code for legend.grid and build it back up manually thius passing legend(parameter) to legend (function) directly
# in legend.grid: 'type' populates 'legend' argument to legend function.
# change legend.grid so there's a type option to allow character strings
# change gbm.map L61 to call legend.grid with type = type
# change gbm.map arguments: add type = 2,
# in gbm.valuemap L387 & 402 add "type = goodname"
goodname <- c("Blonde", "Cuckoo", "Spotted", "Thornback")
class(goodname) # "character"
if (class(goodname) == "character") print("success")
goodname <- 1
# insert into legend.grid @ L17:
if (class(type) == "character") legend <- type
# asked hans


# L177 do I use k for anything? Is that a problem?

# Improve closed area calculation algorithm. Use some kind of sorting algorithm rather than one by one counting. Clever loop system or something?
# Closing the INVERSE/worst bottom-up counted-to HRMSY% is closing the BEST 1-HRMSY%.
# I'm counting to up to 8% instead of down to 92%. Fine. How to make this quicker though?
# asked on stackoverflow, no reply
dbase <- data.frame(CPUE = rep(2,10), row.names = NULL)
dbase <- data.frame(CPUE = runif(378570, min = 0, max = 1), row.names = NULL)

HRMSY = 0.4
HRMSY = 0.08

n <- nrow(dbase)

CPUEMSY <- (sum(dbase[,"CPUE"]) * HRMSY)

for (k in 1:nrow(dbase)) {
  print(paste(n,", ",round((sum(dbase[n:nrow(dbase),"CPUE"])/CPUEMSY) * 100, 3),"%",sep = ""))  # progress printer
  if (sum(dbase[n:nrow(dbase),"CPUE"]) < CPUEMSY) {n = n - 1} else { #if sum of rows from end to this point < CPUEMSY aka HRMSY% add another row and sum again, ELSE:
    print(paste("close rows 1:", n - 1, sep = ""))
    break
    }
  }
# Asked Coilin. See email, potentially investigate halving with rounding.
# 18 splits for my 378570 rows. 18 loops would phenomenally reduce processing time.
# Coilin answer:
which.min((cumsum(dbase[nrow(dbase):1,"CPUE"]) - CPUEMSY) ^ 2)
# This asks which index (-1) has the cumulative CPUE closest to the CPUEMSY.
# Maybe not what you are after but cumsum could assist.
coilins <- n - which.min((cumsum(dbase[n:1,"CPUE"]) - CPUEMSY) ^ 2)
# cumsum from end (n) upwards towards 1 of CPUE until the row X where end:X = CPUEMSY i.e. 'open' cells
# close the rest i.e. inverse i.e. n-X
#edited Ls 186 & 313. try it out.
# edited L193&5 changing n+1 to n. Not sure why it's n+1. Doubt 1 row has much effect in any case but check it out later.


# search this: ####remove xlab & ylab above for general code####
# actually a change in gbm.map L57. Put xlab & ylab as gbm.map parameters?

# Only 1 badcols allowed.
# e.g. L77 invert baddata meh.

# PerSpeciesClosedAreaMpas: need to make the scale from 1:4 (0:4?) not auto. There are 1+4 items in the legend, correctly.
# L309 set breaks. gbm.map allows them to be parsed in "..."? Didn't change anything.
# tweaked colours, added breaks1 & 2 globals to gbm.map, manually set valuemap to only run the effort loop
# breaks1: 0:4 correct. breaks2: didn't appear, so unhelpful run! But at least means the breaks I want are passed into the system...
# breaks is being changed somewhere near bottom of gbm.map, but by what?
# instead of 0-1,1-2,2-3,3-4, have it be C T B S! Still need to sort the numbering problem first.
# numberwise should be breaks = c(0,0:length(goodcols)), L381 & 396
# ALREADY DONE???

gbm.valuemap(dbase = mydata, loncolno = 2, latcolno = 1, goodcols = c(5,3,6,4), badcols = 7,
             conservecol = 8, HRMSY = c(0.14,0.08,0.08,0.15), plotthis = c("bad"),
             maploops = "Effort")

plotthis = c("both","close")

# PerSpeciesClosedArea colours kinda suck. How to change? It goes through an awkard route.

# In cumulative closures, sum the HRMSY% for each species, e.g. 4th one will be ~100%, how much is the 1st? 150%?
# where to put that, in legend?

# Metrics: add %area covered
 # dispersion factor / standard deviation / clustering /contagious-uniform distribution. Number of closed areas created. Don’t need all that right now.


####Load scripts & data####
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.basemap.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')

####run gbm.auto with E####
# Load linux
mysamples <- read.csv("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Data/Samples_allRays_Env_F_E.csv", header = TRUE, row.names = NULL)
mygrids <- read.csv("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", header = TRUE)
# Load linux no grain
mysamples <- read.csv("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Data/Samples_allRays_Env_F_E_NoGrain.csv", header = TRUE, row.names = NULL)
mygrids <- read.csv("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE_NoGrain.csv", header = TRUE)
# set directory, with fishing E as an expvar
setwd("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E")

# load mapplots coast and set as shape
library(mapplots)
data(coast)
shape = coast
# run gbm.autos; cuckoo
gbm.auto(expvar = c(4:9,11), resvar = 12, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)
#thornback
gbm.auto(expvar = c(4:9,11), resvar = 13, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)
#blonde
gbm.auto(expvar = c(4:9,11), resvar = 14, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)
#spotted
gbm.auto(expvar = c(4:9,11), resvar = 15, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)
#all at once
gbm.auto(expvar = c(4:9,11), resvar = c(12:15), grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)

# run gbm.autos, NoGrain; cuckoo
gbm.auto(expvar = c(4:7, 9, 11), resvar = 12, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE, mapshape = coast)
#thornback
gbm.auto(expvar = c(4:7, 9, 11), resvar = 13, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE, mapshape = coast)
#blonde
gbm.auto(expvar = c(4:7, 9, 11), resvar = 14, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE, mapshape = coast)
#spotted
gbm.auto(expvar = c(4:7, 9, 11), resvar = 15, grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE, mapshape = coast)
#all at once
gbm.auto(expvar = c(4:7, 9, 11), resvar = c(12:15), grids = mygrids, lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE, mapshape = coast)

####From here to run all automatically####
#(this is an overnight job!)
# Load data & set WD
conserve <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/ConservationMaps/Combo/AllScaledData.csv", header = TRUE, row.names = NULL)
mydata <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/AllPreds_E.csv", header = TRUE, row.names = NULL)
mydata <- cbind(mydata, conserve = conserve[,3]) #add conservation data as a column to mydata
setwd("/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/ValueMaps")
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')

# Create target folders
dir.create("Blonde 0.08")
dir.create("Goodweight 10s")
dir.create("Badweight 10")
dir.create("Goodweight 4.17 3.5 2.33 1")

####Run gbm.valuemap####
data(coast, package = "mapplots")

# normal 1:1 values
setwd("Blonde 0.08")
gbm.valuemap(dbase = mydata,  # data.frame to load. Expects Lon, Lat & data columns: predicted abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
             loncolno = 2, # column number in data which has longitudes
             latcolno = 1, # column number in data which has latitudes
             goodcols = c(5,3,6,4),  # which column numbers are abundances (where higher = better)? C B S T
             badcols = 7,  # which column numbers are 'negative' elements e.g. fishing (where higher = worse)?
             conservecol = 8, #conservation column
             HRMSY = c(0.14,0.08,0.08,0.15),
             mapshape = coast)

# run with effort weight as 10 "badweight"
rm(list = ls()) # remove everything to clear workspace, free memory & reload
conserve <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/ConservationMaps/Combo/AllScaledData.csv", header = TRUE, row.names = NULL)
mydata <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/AllPreds_E.csv", header = TRUE, row.names = NULL)
mydata <- cbind(mydata, conserve = conserve[,3])
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')
setwd("../Badweight 10")
gbm.valuemap(dbase = mydata,
             loncolno = 2,
             latcolno = 1,
             goodcols = c(5,3,6,4),
             badcols = 7,
             conservecol = 8,
             HRMSY = c(0.14,0.08,0.08,0.15),
             badweight = 10) # effort weight as 10

# run with species weights all as 10 "goodweight"
rm(list = ls())
conserve <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/ConservationMaps/Combo/AllScaledData.csv", header = TRUE, row.names = NULL)
mydata <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/AllPreds_E.csv", header = TRUE, row.names = NULL)
mydata <- cbind(mydata, conserve = conserve[,3])
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')
setwd("../Goodweight 10s")
gbm.valuemap(dbase = mydata,
             loncolno = 2,
             latcolno = 1,
             goodcols = c(5,3,6,4),
             badcols = 7,
             conservecol = 8,
             HRMSY = c(0.14,0.08,0.08,0.15),
             goodweight = c(10,10,10,10)) # species weights all 10s

# run with species weights set individually "Goodweight 4 3.5 1.5 1"
rm(list = ls())
conserve <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/ConservationMaps/Combo/AllScaledData.csv", header = TRUE, row.names = NULL)
mydata <- read.csv(file = "/home/simon/Dropbox/Galway/Project Sections/3b. BRT plus Bpa Sam & Dave/Analysis/Model Outputs/With E/AllPreds_E.csv", header = TRUE, row.names = NULL)
mydata <- cbind(mydata, conserve = conserve[,3])
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R')
source('/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.valuemap.R')
setwd("../Goodweight 4.17 3.5 2.33 1")
gbm.valuemap(dbase = mydata,
             loncolno = 2,
             latcolno = 1,
             goodcols = c(5,3,6,4),
             badcols = 7,
             conservecol = 8,
             HRMSY = c(0.14,0.08,0.08,0.15),
             goodweight = c(4.17,3.5,2.33,1)) # species weights individualised

####Run Notes####
#1. plotthis=NULL, savethis=data: ...) used in an incorrect context
#2. beepr didn't load
#3.  Error in `[<-.data.frame`(`*tmp*`, , datacoln + i, value = c(0.00480244254949487,  : new columns would leave holes after existing columns
  #L34. L33 was i in c(goodcols,badcols), would have meant col# data[,datacoln+i] was 7+3 not 7+1
#4.  Error in `[.data.frame`(data, , datacoln + goodcols, drop = FALSE) : undefined columns selected
  #L39,40,42,43: colnumbers logical error, maths was wrong
#5. Error in `[.data.frame`(bioboth, , SortCol0 + 1) : undefined columns selected
  #L183: "close" switched off but bioboth processing after j loop finished was done in top environment
  #i.e. irrespective of "close" switch status, would expect to be processing bioboth
#6. COMPLETED: plotthis=NULL & savethis=data. Scaling wrong: L34 not fully fixed, missed bits
#7. COMPLETED: scaling columns worked, bothdata max is 3.6077 : wrong. Should be max 2, right? And should be more than 1 column. L62 or earlier.
  # baddata can only be 1 col. goodata should be 4. Global save good/bad/both to inspect them. Gooddata is 1 col vector i.e. culprit. L40.
  # L39: if goodweight given, gooddata is goodcols matrix multiplied by goodweight. If NOT given, it's DELIBERATELY the row sums i.e. 1 col, so so L40&43 are wrong.
  # if goodweight not given, gooddata is simply data[,datacoln+seq(1,length(goodcols)). Edited.
#8. Works but need to create new colnames for gooddata before they get added to data via bothdata (L62 & 63). Done, L64 & 65
#9. works. Minor tweak from _wc to _swc. RUN 1 COMPLETE: "data".
#10. plotthis=null savethis=null test. Works
#11. plotthis="bad".  Error in seq.default(xlim[1], xlim[2], by = byx) : wrong sign in 'by' argument
  #L50, gbm.map L50, make.grid // lat & lon numbers were the wrong way round, might not change anything but re-running
#12. COMPLETED bad. try good
#13. COMPLETED good. try both (needs good & bad first? nope)
#14. COMPLETED both. try close
#15. Error: object 'byx' not found. Where? L136 was !is.null but that assumes byx is in the function parameters. True for gbm.map but not here. Changed to !exists("byx")
#16. Error in `[.data.frame`(z, i) : undefined columns selected. L159, z = bioboth[,get(paste("sort",j,sep=""))] probably the problem.
  # From L156, L110. Is just last col, use that
#17. j loop overlay maps work! Yay! Error in `[.data.frame`(bioboth, order(-bioboth[, SortCol0 + 1], -bioboth[,  : undefined columns selected
  # L186. Not sure if bioboth is available; probably a less messy way to set the column numbers also. Global assign sortcols & bioboth to inspect.
  # Just needed a comma at the end to indicate it was a row sort. Apparently defaults to column sort.
#18. Error in rowSums(bioboth[, SortCol0 + 1:SortCol0 + 4], drop = FALSE) : unused argument (drop = FALSE). L194. Removed.
#19.  Error in `[.data.frame`(bioboth, , SortCol0 + 1:SortCol0 + 4) : undefined columns selected L194. Dunno why; copied pmax L190 info instead of range
#20. Error in rowSums(bioboth[, SortCol0 + 1], bioboth[, SortCol0 + 2], bioboth[,  : unused argument (bioboth[, SortCol0 + 4])
  # can't use rowSums as a list, has to be one x for the array. Went back to range, needed to put (SortCol0+1):(SortCol0+4) in brackets otherwise it evaluates wrongly
#21. Error in MPAgrow[, 3] : subscript out of bounds. Was cbinding to zeroes each time instead of growing MPAgrow.
#22 COMPLETED! Black cells are an odd size... Do I need blue background? Not for those at least. Changed. Ideally turn off legend.
#23. Goodcols reordering attempt, 5,3,4,6. COMPLETED. White background is worse, changed back. Might be better if cells were black. L232
#24. plotthis=close savethis=close. Fresh run new day, Error: object 'coast' not found. mapplots was loaded... for some reason L136 data(coast) now  needs ",package="mapplots""
#25. coast STILL not found! Saved in gbm.map as it should have been, trying again
#26. STILL!!! not found. Moved to near-top of gbm.valuemap
#27. COMPLETED. ClosuresData.csv 1st 2 cols named data..loncolno, why? Assumedly Closures object? Which was bioboth. Test @ L70, add colnames?
#28. Fixed. Assigned by lat/loncolno L70
#29. COMPLETED. ClosuresData.csv looks good. Now switch everything on.
#30. COMPLETED. NOICE. Fix black cell size from #22. byx size is different in ClosedValueComboMap & CumulativeClosedAreaMap?
  # added lines 171 & 241
#31. COMPLETED. Nope, byxs & byys are the same size. Hmm. is byx/byy the only control of cell size? Maybe check the plotting data?
  # make.grid & draw.grid are the grid points. grd is the actual object. Swapped 171 & 241 for grd export. Not 100% sure it's not
  # just using the same grd/byx values both times though...This will tell me. byxs MIGHT be the same but grds MUST NOT be. They're identical
  # this isn't getting the two different grds. How to do? Do in gbm.map AND valuemap (as tmps)?
  # L171 in valuemap ^ last line of gbm.map, is from closed area grd2 and the final loop of j so BCTS closed areas.
#32. Now different. ClosedValueComboMapGRD (works) 4.9mb, CumulativeClosedAreaMapGRD (doesn't) 27.3mb
  # different second lat & lon: = different breaks? = different byx/byy? Check that again, editde L172 & end of gbm.map
#33. ClosedValueComboMapBYX works = 0-00417. CumulativeClosedAreaMapBYX (doesn't) = 0.0017
  # ClosedValueComboMapBYX / CumulativeClosedAreaMapBYX = 2.369137 # !! strange
  # find out why this is going wrong or just force it?
  # whats wrong. What's the difference between the lat/long of
  # ClosedValueComboMap (data[,latcolno & loncolno]: works): 38750
  # & CumulativeClosedAreaMap (Closures:doesn't)
  # L217 Closurestmp
#34. length: 378570, same as data. Odd! How does byx work? Completely based on x i.e. lon
  # does Closurestmp = data[,1]? Maybe reordered?
  # Closures first isn't reordered since bioboth in L189
  # maybe just reorder it back to "normal" @ L218??
  # data is ordered by latitude (descending) then longitude (descending). So maybe gbm.map NEEDS it to be like that.
  # Yep: byx calc works on a next-cell basis
  # so I COULD do the ordering OR I could force CumulativeClosedAreaMap to use the same byx?
  # Maybe reordering will work NOW but not generally?
  # fuck it, do it first then see. Done, L218 area, using cols 2 & 1 since: bioboth<-data.frame(data[,loncolno],data[,latcolno]
#35. Didn't put a motherfucking comma at the end of the motherfucking ORDER line, AGAIN. Argh!!
  # Remove notes in valuemap afterwards. byx failure. Why!? Lat & lon numbers the wrong way around before
#36. Data's lat/lon, closure's is lon/lat. Data ISN'T sorted by -lat -lon! How IS is sorted then? IDK.
  # Make a sort index for data, append it at the start, carry it all the way through to Closures...then sort by that?
  # First try not sorting closures but swapping columns? Tried
#37. Nope same. I guess try the sort index? No wait, first try forcing it to use data's byx. what was the last working byx?
  # non-global objects set @ L171, called as gbm.map parameters @ L244&5
#38. COMPLETED! Blacks working but lat/longs flipped! Blanked out column flip @ L220. Delete these notes when it works.
#39. COMPLETED. NoooooooICE! Swap blue for white background, blue too close to the grey
#40. COMPLETED. UUUULLLLLLTTTTRRRAAA WWWWWWIIIIIINNNNNN

####Dave Chat####
#A# Add maps: sort by -biomass, closedcombomap
#B# Add maps: sort by +effort, closedcombomap
#C# What % of effort reduced by each closure? Print on map??

#A#
#41. export pre-sort bioboth L109&110 to view it
 # bioboth doesn't contain Hans_E so can't sort by it
 #L70 add Hans_E. Need to name it? see if that breaks anything?
#42. COMPLETED but closed areas are huge - much bigger than before, basically everything but the effort peaks, and cum.closed takes up everything.
 #put 0.15 for all, closedvaluecombomaps are normal again but CumClosed still broken. Restart R in case they're using saved objcets
#43. nope still broken.
  # ading baddata has thrown everyithng off by a row. Add "+1" @ L182
#44. add biomass2 & 3 sections
#45. Fail undefined columns. I changed final section to bioboth1, now changed back to bioboth
#46. Nothing changed, same error. Export SortCol0 & bioboth for inspection
  #L320 onwards changed bioboht to bioboth1: bioboth hasn't got any sorters. Should I do this section thrice, for bioboth 1 2 3? See if it works first.
#47. same problem! biobothtmp has sort4 ONLY not 1 2 3. L110 the problem: @start of each new J it creates a fresh bioboth1 from bioboth.
  # Fixed L75:77 aded bioboth to bioboth1,2,3. Also changed HRMSY to real values
#48. Still not working! Why?! Change to global assign
#49. Still! bioboth1/2/3 all have 11. Sorts aren't added. Because it's globally assigned?
#50. bioboth1 to bioboth, comment out new sections, get it working again.
#51. works again; HRMSY values causing allblacks? Must be an error though, spotted ray should be the same! Rerunning with 0.15s
#52. works, looks good, rerunning with real values, spotted should be the same.
#53. changes to 0.15s, stil the same. Baddata screwed up somehow in biobothtmp, why?
#54. L115 added +1 after j to put sort rows after baddata
#55. baddata still wrong in biobothtmp.
#56. clear environment, re-run with 0.15s
#57.   switched off some plotting. ClosuresData.csv, species_s cols scaled to 2, should be 1?
  # data: lat lon goodcols badcols
  # scaled to 1 l36 &37, called _s, added to data. (+length goolcols)
  # L65 bothdata added to data. (+length(goodcols). L67 _swc named onto bothdata
  # =lat lon goodcols badcols goodcols_s goodcols_swc
  # bioboth should just be data but replace bothdata with bothdata*m
  # SHITLOADS of changes. Data is now the only object running through the code & being processed. Column references updated for that.
  # bioboth gone. Closures no longer in savethis, only data
#58. Error in if (sum(data[1:n, 2 + j]) < CPUEMSY[j]) { : missing value where TRUE/FALSE needed
  # data completely FUCKED, what happened? Create intermediary saveouts
  # L120 culprit. Was using bothdatarange not bothdatarange[1]
#59. Closures replaced with data
#60. Error in if (sum(data[1:n, 2 + j]) < CPUEMSY[j]) { : missing value where TRUE/FALSE needed
  # used n as an object in a separate sheet, clean env reload
#61. STILL fails, why?!?!
  # nightmarish. data is an R core term which means I was trying to sum a function. Doh!!! Changed to x
#62. CPUEMSY too high? Removed CPUEMSY, calculated in place
#63. closedvaluecombo cell sizes weird. cumclosedarea maps lat/lons are reversed and cell sizes are also weird
#64. better but cells are still weird. Reorder dbase before mapping? L125, then use dbase colrefs in L133:135
#65. Added byxout to gbm.map to export byx then set in environment as byx & byy for other maps to use. CumClosedArea still broken, why?
#66. SortCol0 must be wrong, referencing wrong place? AllClosed & SumClosed are wrong. Closure 1:4 wrong too: all 6 are real numbers
#67. SORTED! SortCol0 was out by 1. Now, do the loop thingy.
#68. Working. Closed areas all look suspiciously like thornback abundance though... Closed value effort maps: is it closing high effort FIRST? Looks suspicious
  # Fixed, changed maploopnames to include -data[,etc]. Still need to fix thornback suspiciousness.
#69. Not working again, fails a@ N=1 L125 again
#70. both CPUEMSYtmp & sumtmp = NA_real_
  # maploopcode doesn't caputre enough, needed to add -data etc. But now it doesn't work, 1st loop evaluates to -0.8569
  # mlc is working, it's just the wrong code, when I type it out (mlcdbase) it's wrong. But it's the same as working one from last night!!
  # could it be that last night eval(parse(text=())) evaluated to a simple integer cos it was inside a section.
  # now it includes data. Expand ther eval(parse) to include the whole line
#71. Error in dbase[, bothdatarange[1] - 1 + j] : incorrect number of dimensions
#72 rebuilt the eval c(). WORKS. NICE! Re-run with correct HRMSYs
#73. using correct HSMSY, closed value effort map doesn't really work, is the same for each? Sorted the same but biomass should be different?
  # it's different species, differnet HRMSYs, it must SURELY be different!
  # l197 & 201 not working, colnames not evaluating. Added get().
  # maploopnames added to sort L124:126, ditto L195 to 215.
  # the problem is the SORTS not the closurs (closures are from sorts)
  # n<-1 moved inside J loop so it resets to 1 for each species.
#74. cLOSEDvALUEcOMBO maps are closing max biomass first, not max combo
  # closed value biomass only works for spotted. It's using spotted for all of them. Hasn't reset the J on the fresh O loop??
  # is it not resetting J at the start of the loop (doesn't the for loop do that automatically?) or is it overwriting the results?
  # added reset j @ L111.
#75. Same. So it's not sorting properly? l119?
  # combo is sorted by & closing on biomass, j correct, so basically combo is what biomass should be
  # biomass is closing on spotted only, o loop #2
  # effort is sorted by spotted biomass & closing species biomass number
  # after o loop 1, j is locked to 4 still
  # fixed colnames
  # global o & j counter L120 & dbase 121
#76. print percentage counter not working properly. Round to 0dp, sum or calculation is wrong
  # NA%s for something, what is that?! If one species is aligning with another then that might explain it
  # 1: In sum(1:n) : integer overflow - use sum(as.numeric(.))
  # Os & Js are CORRECT in the loop, i think?
#77. The plus hans e biomass combo maps are slightly wrong. The basemaps are abundance NOT combo
  # are the goodcols_s being affected by hansE?!
#78. Problem: _swc, bothdatarange*m. Has baddata been scaled properly? baddata not inverted?
  # Max rays = 1. No E = 1. Both = 2. I think it's actually working.
  # badcut/all/pct L126:8. Extra legend L184
#79. #legend(x="topright",title=paste("Effort reduced: ",badpct,"%",sep="")) #Error in as.graphicsAnnot(legend):argument "legend" is missing, with no default
#80. legend.grid(legendloc="topright", breaks=NA, bg=lejback, title=paste("Effort reduced: ",badpct,"%",sep="")) error unused argument legendloc
  # L183 title changed. 128 round & digits. WORKS
  # bothdata colour scale. L101, 141, heatcolours
#81. Looks good. L249 changed 1 = Closed to closed E pct. How, less simple. L240
#82. Added effort percentage calculator for cum.closed.areas L235:238






####DONE####
# allow different image device for Mac OSX
# add alert pings for maps and loops and such
# perspeciesclosedarea maps colours get fucked
# Option to set which loops you run (Combo (aka both), Biomass (aka good), Effort (aka bad), Conservation)
# does stuff break if you set plotthis to null or switch certain stuff off? Documented
# L30 data(coast) & search "shape = coast": change to universal
# L349: do algorithmically so I don't have to know how many species there are
# Need B&W outputs as well for per species closed area DONE
# note down where I've used fixed values instead of algorithmic: just the last section? yes. DONE
# improve notes DONE
# L22 object alerts not found DONE
# Add processs printers DONE checkem DONE
# edit savethis since it can only be one thing. DONE

# Thornback closed areas based on P2 conservation
# cbind conservation map predabund column, initial latlongs are the same? YES DONE here L57
# ProcessedData.csv from valuemap run as normal, compare to
# see if valuemap still works now!
##NOT QUITE: ClosedValueMapXXXX_Species are wrong for 1st species i.e. blonde: basemap is effort,
# I.e. it's using a relative column which isn't expecting the presence of conserve.
# L133. L140 is the z which is the problem. CHANGED
# Fixed for 1st species but now differently weird for 3rd & 4th. CTBS problem?
# ClosedValueMapCombo_Blonde has blonde closed area w/ cuckoo basemap, i.e.:
# it's observing the original mydata column order not goodcols
# z = dbase[,ncol(dbase)-length(goodcols)+j]  L140
# colnames are Blonde_swc	Cuckoo_swc	Thornback_swc	Spotted_swc
# colnames are correct... data are wrong?
# L33 scaling using correct columns, so species_s col data are CORRECT
# L64 species_swc data look to be CORRECT also. Simply calling them wrongly?
# maploopnames for combo (L111) not doing the same as hans_e_map, which is correct? (bothdata map L97)
# L97: z = dbase[,bothdatarange[j]]  #12:16 original, 13:17 new
# L70: dbase[,bothdatarange] <- bothdata*m
# L69: bothdatarange <- (ncol(dbase)-length(goodcols)+1):ncol(dbase)
# L111 "dbase[order(-dbase[,bothdatarange[1]-1+j]),]"
# L141 z = dbase[,ncol(dbase)-length(goodcols)+j]. FIXED
#
# ClosedValueMapBiomass & ClosedValueMapEffort_species is CTBS not BCTS order for closed areas (blacks!)
# L176 z = dbase[,ncol(dbase)] # it's ONLY spotted i.e. the last one, doesn't move with the loop.
# Or is the last column necessarily the loop?
# [L125:7 added species name instead of number for sort.]
# if the last column is wrong, called @ L179, it's made @ L129
# so is J pointing to the wrong place in maploopcodes?
# After o loop (maploopcode) 1, the blacks go from BCTS to CTBS. It looks like they're summing the correct amounts,
# but are distributed in the wrong place, due to the wrong sort?
# L114 o2 & 3?
# L115: reverts back to the top of the o loop, j is still 4, L116 j set back to 1.
#
# biomass is CTBS, effort is just S. What the actual FUCK?!
# biomass: L114 2+j is defo wrong, this results in the CTBS. Should be what? Simply goodcols[j]? CHANGED
# biomass FIXED. Effort I think is CTBS rather than BCTS? Defo. So what, it's ordered by least effort.
# A-HA! Least effort value is 0. LOADS of values are 0. Thus for MOST of the dataset, for all of the zeroes,
# it'll just be using the previous sort? But that would be spotted, NOT CTBS...?
# is there something in the effort sort which would pseudo-default biomass sorts to ctbs?
# Needed to sort by effort THEN adundance. Which makes sense if you think about it. DONE
#
# So: closed areas based on P2 abundance i.e. "conserve" column, colno=8
# added to L112 & 118. 119 changed from 1:3 to 1:length(maploopnames)
# DONE, WORKED!


# 4 species cumulative closed area based on combo (4 pic panel; different colour per growth)
# Changed 4s to length(goodcols) in a few places for generalisability later
# L236 removed zeroes in the cbind since they're a space waste
# L275 down: coded SpeciesGrow but commented out for 1st run to check eveything working after minor tweaks.
# L293&294 need to work out colours.
# tweaks work. Uncomment out that bit, sort out colours and see what happens. DONE RUN
# done!

# Change text in L175: main=paste("Predicted Abundance + Fishing Effort:   ?
# DONE

# weighting simulations: 10:1 rays:E, 1:10, (1-10):1 mixed ray weightings
# already coded, just try it out
# works for multiple goodweights
# fails for single badweight: trying to array multiple a length1 item:
# "Error in as.matrix(dbase[, dbasecoln + length(goodcols) + seq(1, length(badcols))]) %*%: non-conformable arguments" L44

# Sam email/chat HRMSY blonde (?) changed from 8% to 5% to make more conservative? LESS.
# nope. HRMSY is max that can be sustainably removed, NOT the minimum you need to conserve.
# Added invHRMSY. Need to re-run combos.

# Re-running now super slow because it's counting up to e.g. 95% instead of up to 5%
# also the column referencing for the counts was wrong: was using CTBS_original not BCTS_userset. CHANGED
# Reverse the ordering & the counting? removed invHRMSY
# same order, counting up from bottom to get to the 5% not down from top to get to the 95% DONE

# L134: add "1 of 16" counter for each line so you can see where you are in the grand scheme DONE






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