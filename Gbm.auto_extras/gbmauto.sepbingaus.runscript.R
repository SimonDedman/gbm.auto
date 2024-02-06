library("gbm.auto")
mygrids <- gbm.auto::grids # load grids
mysamples <- gbm.auto::samples # load samples
# setwd("/home/simon/Desktop/gbm temp/sepParamTestbench/")
# source('~/Dropbox/Galway/Analysis/R/gbm.auto/R/gbm.utils.R')
Crop_Map <- st_read(dsn = paste0("Crop_Map", ".shp"), layer = savename, quiet = TRUE)
#source('~/Dropbox/Galway/Analysis/R/gbm.auto/Gbm.auto_extras/gbm.auto.binGausSepParams.R')

# test1 basic clean run
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:10), resvar = 11,
         tc = 2,
         lr = 0.01,
         bf = 0.5,
         mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE,
         BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
# worked fine, note most stuff switched off. Try full run.

# test run 2 basic w/ all optionals on. Note expvars was 4:10 before, wrongly.
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = 0.01,
         bf = 0.5,
         mapshape = Crop_Map)
# works fine, try basic length2 vectors for params
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = c(2,3),
         lr = c(0.01,0.009),
         bf = c(0.5,0.6),
         mapshape = Crop_Map)
# fine. single numbers with one list, turn off optionals
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = 0.01,
         bf = 0.5,
         mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE,
         BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
# fine. tc lr bf list. try with no grids
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(c(0.5,0.6), 0.7),
         mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE,
         BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
# Failed, report error, edited bin not gaus to test it, retrying. Failed again.
# reduce bf to list of 2
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE,
         BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
# Failed again. L640 & 656 should be (simp)! Turn simp on, retry.
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map, simp = TRUE, varint = FALSE, map = FALSE,
         BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
# finally worked. Try with simp off
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map, simp = F, varint = FALSE, map = FALSE,
         BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
# fine tho decide whether I wanna remove the simp columns if not simplifying
# try with everything on
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map)
# says it worked but no maps produced...
# singles & vectors work. Basic list:
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = 0.01,
         bf = 0.5,
         mapshape = Crop_Map)
# tc list worked. trying lr
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = list(0.01, c(0.01, 0.005)),
         bf = 0.5,
         mapshape = Crop_Map)
# worked. trying bf
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = 0.01,
         bf = list(0.5, 0.7),
         mapshape = Crop_Map)
# worked. tc & lr
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = 0.5,
         mapshape = Crop_Map)
# worked. lr & bf
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map)
# worked. tc & bf
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = 0.01,
         bf = list(0.5, 0.7),
         mapshape = Crop_Map)
# worked. all on.
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map)
# worked. BFs WEREN'T MULTI LISTS!
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = 0.01,
         bf = list(c(0.5, 0.7),0.5),
         mapshape = Crop_Map)
# worked. Try opposite side lists.
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = list(c(0.01, 0.02), 0.005),
         bf = list(0.5, c(0.5, 0.7)),
         mapshape = Crop_Map)
# FAILED: ran through but no maps produced.
# So what combo of list bf with vector, and list tc/lr causes failure? tc?
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = 0.01,
         bf = list(0.5, c(0.5, 0.7)),
         mapshape = Crop_Map)
# FAILED. BF & TC fail. BF & LR?
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = list(c(0.01, 0.02), 0.005),
         bf = list(0.5, c(0.5, 0.7)),
         mapshape = Crop_Map)
# FAILED. So is it just BF? Non-BF-multilist worked...
gbm.auto(samples = mysamples, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = list(c(0.01, 0.02), 0.005),
         bf = 0.5,
         mapshape = Crop_Map)
# FAILED. Weird. This worked when bf was a list of 2... Try that again
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(c(2,3),2),
         lr = list(0.01, c(0.01, 0.005)),
         bf = list(0.5, 0.7),
         mapshape = Crop_Map)
# Worked.
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Variable interactions done XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Binomial predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Gaussian predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Final abundance calculated  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Output CSVs written     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Report CSV written      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       RSB CSV written       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Reticulating splines     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour map generated     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Black & white map generated XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Colour RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     B&W RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Grids/maps/everything done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# So... when
# (bf is a single AND lr OR tc are lists)
# OR
# (bf is a list with a vector AND lr OR tc are lists)
# if fails. Purple lines not printed in fail runs.
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Variable interactions done XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Binomial predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Gaussian predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Final abundance calculated  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Output CSVs written     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Report CSV written      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       RSB CSV written       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Reticulating splines     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour map generated     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Black & white map generated XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Colour RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     B&W RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# [1] "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Grids/maps/everything done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# > gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = list(c(0.01, 0.02), 0.005),
         bf = list(0.5, c(0.5, 0.7)),
         mapshape = Crop_Map)
# worked
# So what combo of list bf with vector, and list tc/lr causes failure? tc?
dir.create("./2")
setwd("./2")
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = 0.01,
         bf = list(0.5, c(0.5, 0.7)),
         mapshape = Crop_Map)
# FINE
dir.create("./3")
setwd("./3")
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = 2,
         lr = list(c(0.01, 0.02), 0.005),
         bf = list(0.5, c(0.5, 0.7)),
         mapshape = Crop_Map)
# FINE
dir.create("./4")
setwd("./4")
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = list(c(0.01, 0.02), 0.005),
         bf = 0.5,
         mapshape = Crop_Map)
# FINE
dir.create("./4")
setwd("./4")
gbm.auto(samples = mysamples, grids = mygrids, expvar = c(4:9), resvar = 11,
         tc = list(2,c(2,3)),
         lr = list(c(0.01, 0.02), 0.005),
         bf = 0.5,
         mapshape = Crop_Map, simp = F)
