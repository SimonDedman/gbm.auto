### BRT code todo ####
setwd("C:/Users/Simon Dedman/Dropbox/Galway/Analysis/R/P1 Map analysis")
setwd("/home/simon/Dropbox/Galway/Analysis/R/P1 Map analysis")
setwd("/media/Windows7_OS/Users/Simon Dedman/Dropbox/Galway/Analysis/R/P1 Map analysis")
#source('C:/Users/Simon Dedman/Dropbox/Galway/Analysis/R/Elith 2008 BRT/gbm.predict.grids.R')
# mygrids<-read.csv("ALL_ENV_MGT.csv", header = TRUE)
mygrids<-read.csv("grids_essential.csv", header = TRUE)
# mygrids$BOTTEMP <- mygrids$SBT_SEP #do I need this still?
mysamples<-read.csv("DATRAS_ENV_MGT2.csv", header = TRUE, row.names=NULL)
load("Bin_Best_Model")
load("Gaus_Best_Model")
source('../gbm.auto/gbm.utils.R')


###remove these###
expvar <- c(115,124,136,147,153,154)   # list of column numbers of explanatory variables in SAMPLES.
resvar <- c(111)  # column number response variable (e.g. CPUE) in SAMPLES. resvar=c(111,60,87,75,102), all,c,t,b,s
tc <- c(2,5)   # list permutations of tree complexity allowed c(2,5)
lr <- c(0.01,0.005)   # list permutations of learning rate allowed c(0.01,0.005)
bf <- c(0.5)   # list permutations of bag fraction allowed
gridslat = 2
gridslon = 1
ZI = TRUE
map = TRUE
RSB = TRUE
cols = grey.colors(6,0.3,1)
###
i <- resvar
j <- 2
k <- 0.01
l <- 0.5
predcpue=3
samples<-mysamples
grids<-mygrids
###
rm(list = ls())
###remove these###

####run gbm.auto####
gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(111),  #,60,87,75,102
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = FALSE,   #TRUE
         RSB=FALSE,     #TRUE
         legendtitle = "CPUE")


#### run gbm.rsb ####
gbm.rsb(samples,
        grids,
        expvarnames,
        gridslat=2,
        gridslon=1,
        rsbres=resvar)


####run gbm.map####
# setup output image parameters
png(filename = paste("./",names(samples[i]),"/PredAbundanceMap.png",sep=""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
# run function
gbm.map(x = grids[,gridslon],
        y = grids[,gridslat],
        z = grids[,predabund],
        byx = byx,
        byy = byy,
        mapmain = "Predicted abundance: ",
        species = names(samples[i]),
        shape = coast,
        landcol = "darkgreen",
        legendloc = "bottomright",
        legendtitle = "CPUE")
dev.off()


#####BRT Subsets:####
####Juveniles####
# cuckoo
c_sample <- subset(mysamples, CRAY_AV.L > 0 & CRAY_AV.L < 400)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(60),
         grids=mygrids,
         samples=c_sample,
         tc=c(5),
         lr=c(0.00000000001),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# thornback
t_sample <- subset(mysamples, TRAY_AV.L > 0 & TRAY_AV.L < 450)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(87),
         grids=mygrids,
         samples=t_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# blonde
b_sample <- subset(mysamples, BRAY_AV.L > 0 & BRAY_AV.L < 550)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(75),
         grids=mygrids,
         samples=b_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# spotted
s_sample <- subset(mysamples, SRAY_AV.L > 0 & SRAY_AV.L < 400)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(102),
         grids=mygrids,
         samples=s_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)


####Mature Females####
# cuckoo
c_sample <- subset(mysamples, CRAY_F_AV.L > 400)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(60),
         grids=mygrids,
         samples=c_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# thornback
t_sample <- subset(mysamples, TRAY_F_AV.L > 450)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(87),
         grids=mygrids,
         samples=t_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# blonde
b_sample <- subset(mysamples, BRAY_F_AV.L > 550)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(75),
         grids=mygrids,
         samples=b_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# spotted
s_sample <- subset(mysamples, SRAY_F_AV.L > 400)

gbm.auto(expvar=c(115,124,136,147,153,154),
         resvar=c(102),
         grids=mygrids,
         samples=s_sample,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB=TRUE)

# Allow user to call extra parameters with '...' which will be parsed to gbm.step in places?
# function (data,                             # the input dataframe
#           gbm.x,                                    # the predictors
#           gbm.y,                                    # and response
#           offset = NULL,                            # allows an offset to be specified
#           fold.vector = NULL,                       # allows a fold vector to be read in for CV with offsets,
#           tree.complexity = 1,                      # sets the complexity of individual trees
#           learning.rate = 0.01,                     # sets the weight applied to inidivudal trees
#           bag.fraction = 0.75,                      # sets the proportion of observations used in selecting variables
#           site.weights = rep(1, nrow(data)),        # allows varying weighting for sites
#           var.monotone = rep(0, length(gbm.x)),     # restricts responses to individual predictors to monotone
#           n.folds = 10,                             # number of folds
#           prev.stratify = TRUE,                     # prevalence stratify the folds - only for p/a data
#           family = "bernoulli",                     # family - bernoulli (=binomial), poisson, laplace or gaussian
#           n.trees = 50,                             # number of initial trees to fit
#           step.size = n.trees,                      # numbers of trees to add at each cycle
#           max.trees = 10000,                        # max number of trees to fit before stopping
#           tolerance.method = "auto",                # method to use in deciding to stop - "fixed" or "auto"
#           tolerance = 0.001,                        # tolerance value to use - if method == fixed is absolute,
#           # if auto is multiplier * total mean deviance
#           keep.data = FALSE,                        # keep raw data in final model
#           plot.main = TRUE,                         # plot hold-out deviance curve
#           plot.folds = FALSE,                       # plot the individual folds as well
#           verbose = TRUE,                           # control amount of screen reporting
#           silent = FALSE,                           # to allow running with no output for simplifying model)
#           keep.fold.models = FALSE,                 # keep the fold models from cross valiation
#           keep.fold.vector = FALSE,                 # allows the vector defining fold membership to be kept
#           keep.fold.fit = FALSE,                    # allows the predicted values for observations from CV to be kept
#           ...)                                      # allows for any additional plotting parameters
#
# include eliths BRT_ALL within this code?

Bag fraction optimiser
# Trial & error iterative approach to determine the optimal bag fraction for a data set? How? Stop @ whole percentages?
# Try OPTIM function & see http://r.789695.n4.nabble.com/Optimization-in-R-similar-to-MS-Excel-Solver-td4646759.html
# Possibly have an option to start with this in the R function.

ZI TEST

OPTIMISE PARAMETERS
# Trial & error iterative approach to determine the optimal bag fraction (/ all parameters?)
# for a data set? How? Stop @ whole percentages. Try OPTIM function & see
# http://r.789695.n4.nabble.com/Optimization-in-R-similar-to-MS-Excel-Solver-td4646759.html
# Possibly have an option to start with this in the R function? Separate function?
# Maybe do as separate function then feed into this so the outputs are jkl?
# or make one uber function but can use all 3 separately. Uberfunction doesn't need the loops?
# optim: use Method "L-BFGS-B"
# require(optimx)
# see: https://stats.stackexchange.com/questions/103495/how-to-find-optimal-values-for-the-tuning-parameters-in-boosting-trees/105653#105653
## The caret package in R is tailor made for this.
# Its train function takes a grid of parameter values and evaluates the performance using various flavors of cross-validation or the bootstrap. The package author has written a book, Applied predictive modeling, which is highly recommended. 5 repeats of 10-fold cross-validation is used throughout the book.
# For choosing the tree depth, I would first go for subject matter knowledge about the problem, i.e. if you do not expect any interactions - restrict the depth to 1 or go for a flexible parametric model (which is much easier to understand and interpret). That being said, I often find myself tuning the tree depth as subject matter knowledge is often very limited.
# I think the gbm package tunes the number of trees for fixed values of the tree depth and shrinkage.
# https://www.youtube.com/watch?v=7Jbb2ItbTC4
# see gbm.fixed in BRT_ALL.R - having processed the optimal BRT, might as well just use those details going forward rather than re-running the best one again.

Processing time estimate
#	option to have R poll the computer and, when you say go, popup a box saying
# "Based on the parameters you've selected to try, this will run X models on a Y-item-sized dataset (Y = variables * count)
# which may take about Z minutes based on your processor. You have a multicore processor so it will take Z/#processors (time)
# if you have multicore processing enabled - see here LINK" OK/Cancel
# Edit the BRT progrss counter to account for multiple resvars: currently does e.g. n/8 then loops back to 1. Not useful.
# Also could print the current resvar name. Maybe do a running time thing? "This is BRT N of X, Y% complete, took time Z, total
# expected time AA, time remaining AB
PROCESSING TIMER - nice to have, non essential
# option to have R poll the computer and, when you say go, popup a box saying
# "Based on the parameters you've selected to try, this will run X models on a
# Y-item-sized dataset (Y = variables * count) which may take about Z minutes based
# on your processor. You have a multicore processor so it will take Z/#processors (time)
# if you have multicore processing enabled - see here LINK" OK

3D PLOT
# what to do? Bother with it?

Audio notifications
# beepr for completed run and error (different noises)? Already in the code, works?

Clean workspace - dump fails: TRY rm() rather than dump() - dont need this for function
# rm(list = ls()) removes everything
# running function the dumping everyithng means R is still using 1.7gb of RAM!

MULTICORE PROCESSING - already supposedly incorporated; doesnt work
# Re-investigate multicore R
# https://stackoverflow.com/questions/4775098/r-with-a-multi-core-processor
# http://www.google.com/url?q=http%3A%2F%2Fcran.r-project.org%2Fweb%2Fviews%2FHighPerformanceComputing.html&sa=D&sntz=1&usg=AFQjCNEshsLTzWAF9g4cUGzvn5zIfIywsA
# https://groups.google.com/forum/#!topic/davis-rug/VLIXz5i3vZI
# n.cores in gbm
# The number of CPU cores to use. The cross-validation loop will attempt to send different CV folds off to different cores.
# If n.cores is not specified by the user, it is guessed using the detectCores function in the parallel package.
# Note that the documentation for detectCores makes clear that it is not failsave and could return a spurious number of available cores.
# for gbm.step gbm.simplify gbm.plot gbm.plot.fit gbm.predict.grids gbm.percpec (3d plot). Try all of them?
# 'parallel' package loaded by gmb?
# n.cores = detectCores(all.tests = FALSE, logical = FALSE)
# Warning messages:
#   1: In plot.window(...) : "n.cores" is not a graphical parameter
# 2: In plot.xy(xy, type, ...) : "n.cores" is not a graphical parameter
# 3: In axis(side = side, at = at, labels = labels, ...) :
#   "n.cores" is not a graphical parameter
# 4: In axis(side = side, at = at, labels = labels, ...) :
#   "n.cores" is not a graphical parameter
# 5: In box(...) : "n.cores" is not a graphical parameter
# 6: In title(...) : "n.cores" is not a graphical parameter
# 7: "n.cores" is not a graphical parameter
# So: fails and is supposedly done by default anyway.