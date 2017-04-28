## Run gbm.cons
library("devtools")
install_github("SimonDedman/gbm.auto") # update gbm.auto to latest
library("gbm.auto")
mygrids <- gbm.auto::grids
Juveniles <- gbm.auto::Juveniles #NoGrain
Adult_Females <- gbm.auto::Adult_Females #NoGrain

setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Testbench")

# load mapplots coast and set as shape
library(mapplots)
data(coast)
shape <- coast

# NoGrain Runs
gbm.cons(mygrids = mygrids, # csv file (/ & location) of gridded lat long data to predict to
         subsets = c("Juveniles","Adult_Females"), # Subset names
         alerts = TRUE,                # play sounds to mark progress steps
         map = TRUE,                   # produce maps
         BnW = TRUE,
         resvars = c(43:46,10:13), #NoGrain dataset values
         gbmautos = TRUE,
         expvars = list(c(4:10,14,16,20,24,28,36), #NoGrain dataset values
                        c(4:10,14,17,21,25,29,37),
                        c(4:10,14,18,22,26,30),
                        c(4:10,14,19,23,27,31,38),
                        4:9,
                        4:9,
                        4:9,
                        4:9),
         tcs = list(c(2,13),
                    c(2,13),
                    12,
                    c(2,13),
                    c(2,5),
                    c(2,5),
                    5,
                    c(2,5)),
         lrs = list(c(0.01,0.005),
                    c(0.01,0.005),
                    0.005,
                    c(0.01,0.005),
                    0.005,
                    0.005,
                    0.001,
                    0.005),
         ZIs = rep(TRUE,8),
         savegbms = rep(FALSE, 8),
         varints = rep(FALSE, 8),
         RSBs = rep(FALSE, 8),
         BnWs = rep(FALSE, 8),
         shape = shape,
         zeroes = rep(FALSE,8))

# With Grain
mygrids <- read.csv("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto data csvs/All Sheets/Grain/grids.csv")
Juveniles <- read.csv("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto data csvs/All Sheets/Grain/Juveniles.csv")
Adult_Females <- read.csv("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto data csvs/All Sheets/Grain/Adult_Females.csv") # single or vector of samples data csv files corresponding to subsets
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Testbench")

gbm.cons(mygrids = mygrids, # csv file (/ & location) of gridded lat long data to predict to
         subsets = c("Juveniles","Adult_Females"), # Subset names
         alerts = TRUE,                # play sounds to mark progress steps
         map = TRUE,                   # produce maps
         BnW = TRUE,
         resvars = c(44:47,11:14), #Grain dataset values
         gbmautos = TRUE,
         expvars = list(c(4:11,15,17,21,25,29,37), #Grain dataset values
                        c(4:11,15,18,22,26,30,38),
                        c(4:11,15,19,23,27,31),
                        c(4:11,15,20,24,28,32,39),
                        4:10,
                        4:10,
                        4:10,
                        4:10),
         tcs = list(c(2,14),
                    c(2,14),
                    13,
                    c(2,14),
                    6,
                    c(2,6),
                    6,
                    6),
         lrs = list(c(0.01,0.005),
                    c(0.01,0.005),
                    0.005,
                    c(0.01,0.005),
                    0.001,
                    0.001,
                    0.00001, #still fails
                    0.00001), # refuses to resolve irrespective of LR: list(0.001,0.000000001),
         ZIs = rep(TRUE,8),
         savegbms = rep(FALSE, 8),
         varints = rep(FALSE, 8),
         RSBs = rep(FALSE, 8),
         BnWs = rep(FALSE, 8),
         shape = shape,
         zeroes = rep(FALSE,8))
# failed @ matF C, gaus. Fixed. Failing on blonde.
# Failing on spotted: bin simpified model:L rror in round(gbm.object$cv.statistics$deviance.mean, 4) :
#   non-numeric argument to mathematical function
# Spotted binonly (gaus=F, lr=0.0005): Error in x[[jj]]: attempt to select less than one element in get1index <real>
# > traceback()
# 3: `[<-.data.frame`(`*tmp*`, 1:(length(Bin_Best_Simp_Check$final.drops$preds) -
#                                   dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1]),
#                     (reportcolno - 11), value = c("Depth", "Current_Speed"))
# 2: `[<-`(`*tmp*`, 1:(length(Bin_Best_Simp_Check$final.drops$preds) -
#                        dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1]),
#          (reportcolno - 11), value = c("Depth", "Current_Speed"))
# 1: gbm.auto(grids = mygrids, samples = Adult_Females, resvar = 14,
#             expvar = 4:10, tc = c(6), lr = 5e-04, savegbm = FALSE, varint = FALSE,
#             RSB = FALSE, BnW = FALSE, simp = T, shape = shape, gaus = F)
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Testbench/Adult_Females/")
gbm.auto(grids = mygrids,
         samples = Adult_Females,
         resvar = 14, #spotted
         expvar = 4:10,
         tc = c(6),
         lr = 0.0005,
         savegbm = FALSE,
         varint = FALSE,
         RSB = FALSE,
         BnW = FALSE,
         simp = F,
         shape = shape,
         gaus = F) # worked, binonly, no simp

gbm.auto(grids = mygrids,
         samples = Adult_Females,
         resvar = 13, # blonde
         expvar = 4:10,
         tc = c(6),
         lr = 0.0001,
         savegbm = FALSE,
         varint = FALSE,
         RSB = FALSE,
         BnW = FALSE,
         simp = F,
         shape = shape,
         gaus = F)

# without doing gmb.auto runs
# make it clear how to do this in P4
gbm.cons(mygrids = mygrids,
         subsets = c("Juveniles","Adult_Females"), # Subset names
         gbmautos = FALSE, resvars = c(44:47,11:14), shape = shape)
