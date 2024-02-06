# gbm.auto variance idea
# bootstrapping & variance estimates: how to run gbm.auto or gbm.step in a bootstrap?

# variance: just need to re-run it to different folders, so feasibly just have a loop with
# all unecessary stuff switched off and a new folder each time?
#
# A looping function which repeats the same param combo a ton of times and
# takes the variance of the predicted abundance at each site, for each loop. This
# will hopefully end up with a variance surface and/or variance score, since
# variance is kinda a missing metric in these analyses.


# "It does seem that the Gaussian models stop working reliably (I got individual
# runs to work for bull and sandbar sharks, but could never get the same
# parameters to work more than once) somewhere between 44 and 33 “positive”
# sets.  I wonder if it might be worth a separate paper doing some kind of
# sensitivity analysis to figure out where that line actually is?" Chuck
# Need to work out how to do code for
# if (fails) record info, change params, re-run.

# A bootstrapping function. Essentially looping the same params, but removing
# random single/multiple rows/columns of data to test for e.g. time series effect
# even if single year splits aren't powerful enough to run a BRT on their own
# because of insufficient data.
#
# So maybe this kind of analysis could fit into the coding for one of these? Or
# all 3 together. They're all clearly related. Repeating, sometimes taking stuff
# out, and collating answers at the end.

# source('~/Dropbox/Galway/Analysis/R/gbm.auto/R/gbm.utils.R')
library("gbm.auto")
mygrids <- gbm.auto::grids # load grids
mysamples <- gbm.auto::samples # load samples
setwd("/home/simon/Dropbox/Galway/Project Sections/5. Intro & Conclusion/Extra graphics/CofV/")

# Crop_Map <- st_read(dsn = paste0("Crop_Map", ".shp"), layer = savename, quiet = TRUE)
# library(mapplots)
data(coast)
shape <- coast

gbmlooptest <- gbm.loop(loops = 2,
         savecsv = T, #unnecessary
         samples = mysamples,
         grids = mygrids,
         expvar = c(4:10),
         resvar = 11,
         simp = FALSE,
         RSB = F)

gbmlooptest <- gbm.loop(loops = 10, # the number of loops required, integer
                        savecsv = T, # save the variances in simple & extended format
                        calcpreds = T,
                        varmap = T, # create a map of the variance outputs?
                        measure = "CPUE", # map legend, variance of what? Default CPUE
                        cleanup = F, # remove gbm.auto-generated directory each loop?
                        grids = mygrids,         # explantory data to predict to. Import with (e.g.)
                        # read.csv and specify object name.
                        samples = mysamples,  # explanatory and response variables to predict from.
                        # Keep col names short, no odd characters, starting numerals or terminal periods
                        # Spaces may be converted to periods in directory names, underscores won't.
                        # Can be a subset
                        expvar = c(4:10),               # list of column numbers of explanatory variables in
                        # 'samples', expected e.g. c(1,35,67,etc.). No default
                        resvar = 11,               # column number of response variable (e.g. CPUE) in
                        # samples. Expected, e.g. 12. No default. Column name should be species name
                        tc = 2,            # permutations of tree complexity allowed, can be a
                        # vector with the largest sized number no larger than the number of
                        # explanatory variables e.g. c(2,7), or a list of 2 single numbers or vectors,
                        # the first to be passed to the binary BRT, the second to the Gaussian, e.g.
                        # tc = list(c(2,6), 2) or list(6, c(2,6))
                        lr = 0.01,   # permutations of learning rate allowed. Can be a
                        # vector or a list of 2 single numbers or vectors, the first to be passed to
                        # the binary BRT, the second to the Gaussian, e.g.
                        # lr = list(c(0.01,0.02),0.0001) or list(0.01,c(0.001, 0.0005))
                        bf = 0.5,             # permutations of bag fraction allowed, can be single
                        # number, vector or list, per tc and lr
                        ZI = "CHECK",         # are data zero-inflated? TRUE/FALSE/"CHECK".
                        # TRUE: delta BRT, log-normalised Gaus, reverse log-norm and bias corrected.
                        # FALSE: do Gaussian only, no log-normalisation.
                        # CHECK: Tests data for you. Default is TRUE.
                        simp = F,          # try simplfying best BRTs?
                        gridslat = 2,         # column number for latitude in 'grids'
                        gridslon = 1,         # column number for longitude in 'grids'
                        cols = grey.colors(1,1,1), # barplot colour vector. Assignment in order of
                        # explanatory variables. Default 1*white: white bars black borders. '1*' repeats
                        linesfiles = T,   # save individual line plots' data as csv's?
                        savegbm = F,       # save gbm objects and make available in environment after running? Open with load("Bin_Best_Model")
                        varint = F,        # calculate variable interactions? Default:TRUE, FALSE
                        # for error "contrasts can be applied only to factors with 2 or more levels"
                        map = T,           # save abundance map png files?
                        shape = shape,      # set coast shapefile, else downloaded and autogenerated
                        RSB = F,           # run Unrepresentativeness surface builder?
                        BnW = F,           # repeat maps in black and white e.g. for print journals
                        alerts = T,        # play sounds to mark progress steps
                        pngtype = "cairo-png"
                        ) # filetype for png files, alternatively try "quartz"

gbmlooptest <- gbm.loop(grids = mygrids,
                        samples = mysamples,
                        expvar = c(4:10),
                        resvar = 11,
                        simp = F,
                        shape = shape)
