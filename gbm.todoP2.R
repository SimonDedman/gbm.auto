## gbm.todo.P2
# working script for paper 2
# 1/6/2015
# could I incorporate gbm.utils, gbm.rsb and gbm.map into gbm.auto?
# gbm.auto: add in second map & rsb generator for greyscale versions with appropriate names
##

####prep####
# set & test working directory
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles")
getwd()

# load gmb.utils, which contains various of Elith's packages not contained in dismo
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")

# load gbm functions
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R")

# load grids
mygrids<-read.csv("grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA.csv", header = TRUE)

# load saved models if re-running aspects from a previous run
# load("Bin_Best_Model")
# load("Gaus_Best_Model")

####model:juve individual preds####
# load samples
mysamples<-read.csv("Hauls&J&Preds&Enviros_Trimmed.csv", header = TRUE, row.names=NULL)

# Samples Response variables (CTBS) columns 44-47 inc.
# colnames(mysamples), I need 4:9 (enviros), 10 (hans lpue), 11 (combined scallop icesrects), 15 (whelk CPUE icesrects)
# Samples Explanatory variables columns 4-11,15 (all), plus
# Cuckoo: 17, 21, 25, 29, 33, 37 // 40 (all together),
# Thornback: 18, 22, 26, 30, 34, 38 // 41 (all together),
# Blonde: 19, 23, 27, 31, 35 // 42 (all together),
# Spotted: 20, 24, 28, 32, 36, 39 // 43 (all together),

# Grids Explanatory variables columns 3-15 (all), plus
# Cuckoo: 17, 21, 25, 29, 33, 37, 40 (all together),
# Thornback: 18, 22, 26, 30, 34, 38, 41 (all together),
# Blonde: 19, 23, 27, 31, 35, 42 (all together),
# Spotted: 20, 24, 28, 32, 36, 39, 43 (all together),
# (Model doesn't require you to put these in since it just looks up the
# column names from #samples in #grids)

# Set WD so I can run all these at once
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators")
#cuckoo
gbm.auto(expvar=c(4:11,15,17,21,25,29,33,37),
         resvar=c(44),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,15),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#thornback
gbm.auto(expvar=c(4:11,15,18,22,26,30,34,38),
         resvar=c(45),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,15),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#blonde without 15 interactions
gbm.auto(expvar=c(4:11,15,19,23,27,31,35),
         resvar=c(46),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#re run blonde w/ smaller LR for 15 interactions (1 combo), then test against scores for 2 & 5 interactions
gbm.auto(expvar=c(4:11,15,19,23,27,31,35),
         resvar=c(46),
         grids=mygrids,
         samples=mysamples,
         tc=c(15),
         lr=c(0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#spotted
gbm.auto(expvar=c(4:11,15,20,24,28,32,36,39),
         resvar=c(47),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,15),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

####model:juve combined preds####
# Set WD so I can run all these at once
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Combined Predators")
#cuckoo
gbm.auto(expvar=c(4:11,15,40),
         resvar=c(44),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,10),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#thornback
gbm.auto(expvar=c(4:11,15,41),
         resvar=c(45),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,10),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#blonde without 10 interactions
gbm.auto(expvar=c(4:11,15,42),
         resvar=c(46),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#re run blonde w/ smaller LR for 10 interactions (1 combo), then test against scores for 2 & 5 interactions
gbm.auto(expvar=c(4:11,15,42),
         resvar=c(46),
         grids=mygrids,
         samples=mysamples,
         tc=c(10),
         lr=c(0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#spotted
gbm.auto(expvar=c(4:11,15,43),
         resvar=c(47),
         grids=mygrids,
         samples=mysamples,
         tc=c(2,5,10),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

# individual - blonde: no prediced abundance & rsb figures, why? Not enough data.

####model: mat F####
# set wd for mature female samples sheets
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F")
# load samples
mysamples<-read.csv("F_Mat_plus_LPUE_plus_Enviro_IS_AllSp.csv", header = TRUE, row.names=NULL)

# run models: cuckoo
gbm.auto(expvar=4:10,
         resvar=11,
         grids=mygrids,
         samples=mysamples,
         tc=c(2,6),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#thornback
gbm.auto(expvar=4:10,
         resvar=12,
         grids=mygrids,
         samples=mysamples,
         tc=c(2,6),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

#blonde
gbm.auto(expvar=4:10,
         resvar=13,
         grids=mygrids,
         samples=mysamples,
         tc=c(6), #2, #report shows 6 is best
         lr=c(0.001), #0.01, 0.005
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

# restart model with a smaller learning rate or smaller step size...
# Error in if (get(paste("Bin_BRT", ".tc", j, ".lr", k * 100, ".bf", l, :argument is of length zero
# removed lr=0.01
#
# Warning messages:
# 1: In cor(y_i, u_i) : the standard deviation is zero
# 2: In cor(y_i, u_i) : the standard deviation is zero
#
# Samples Explanatory variables columns 4:10 inclusive
# Samples Response variables (CTBS) columns 11:14 inclusive
#
# blonde has only 20 positives in 2137 samples, 0.9%, compared to 5.6 - 7.8 for the other species.
# try higher learning rate?
# tried 0.001, lasted longer, failed @ simplified gaus gbm.interactions:
#
# 1: In cor(y_i, u_i) : the standard deviation is zero
# 2: In cor(y_i, u_i) : the standard deviation is zero
# 3: glm.fit: algorithm did not converge 
# 4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
#
# so y.data == u_i which is (from gbm.step)
##73 u_i <- sum(y.data * site.weights)/sum(site.weights)
##3 site.weights = rep(1, nrow(data))
##74 u_i <- rep(u_i, length(y_i))
#
# if line73 u_i calcaulation is one number repeated across the length of y.data,
# and y_i = u_i, then y_i - u_i for each row must be ~2200 repeats of 0-0, so
# variance (&SD) is BASICALLY 0 (but not actually, surely?)
# BUT!
# All outputs appear to be in the folder. How? I guess it finished then showed warnings.
# So is this a problem? Maybe discuss in results. Need to understand the problem though.
# test the data
y_i2 <- mysamples[,13]
site.weights2 <- rep(1, length(y_i2))
u_i2 <- sum(y_i2 * site.weights2)/sum(site.weights2)
u_i2 <- rep(u_i2, length(y_i2))
u_i2 # show data
y_i2 # show data
cor(y_i2, u_i2)
# In cor(y_i2, u_i2) : the standard deviation is zero
var(y_i2, u_i2) # = 0
# do variance/sd manually
z_i2 <- u_i2 - y_i2
x_i2 <- z_i2^2
var_i2 <- mean(x_i2) #0.211
sqrt(var_i2) #0.459 so why isn't this working?!
# alternative tack: try working data to see the difference
y_i3 <- mysamples[,12] #thornback
site.weights3 <- rep(1, length(y_i3))
u_i3 <- sum(y_i3 * site.weights3)/sum(site.weights3)
u_i3 <- rep(u_i3, length(y_i3))
u_i3 # show data: values ~10X larger than u_i2
y_i3 # show data
cor(y_i3, u_i3) # In cor(y_i3, u_i3) : the standard deviation is zero
# Huh. So it fails for that too...!
# but not when run as part of the full model.

#spotted
gbm.auto(expvar=4:10,
         resvar=14,
         grids=mygrids,
         samples=mysamples,
         tc=c(2,6),
         lr=c(0.01, 0.005),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")

####test section####
#cuckoo
mysamples<-read.csv("F_Mat_plus_LPUE_plus_Enviro_Cuckoo_Filtered_AllLPUE.csv", header = TRUE, row.names=NULL)
gbm.auto(expvar=5:11,
         resvar=4,
         grids=mygrids,
         samples=mysamples,
         tc=c(5),
         lr=c(0.00001),
         bf=c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         legendtitle = "CPUE")


# create binary (0/1) response variable, for bernoulli BRTs
mysamples$brv <- ifelse(mysamples[4] > 0, 1, 0)
brvcol <- which(colnames(mysamples)=="brv") # brv column number for BRT

gbm.step(data = mysamples, gbm.x = 5:11, gbm.y = brvcol, family = "bernoulli", 
         tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5)

# Error in while (delta.deviance > tolerance.test & n.fitted < max.trees) { : 
# missing value where TRUE/FALSE needed

gbm.step

#see gbm.step.debug
#159: code after { doesn't evaluate to TRUE/FALSE ?
# no,  The condition must have either a TRUE or FALSE result. The { afterwards is what to do.
# So "delta.deviance > tolerance.test & n.fitted < max.trees" doesn't result in TRUE or FALSE
# what are the terms?

# delta.deviance: 1

# tolerance.test:
##78 tolerance.test <- tolerance = 0.001 (default, I set this by learning rate, I think)
##79 tolerance.test <- mean.total.deviance * tolerance
##77 mean.total.deviance <- total.deviance/n.cases
##75 total.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)
##72 y_i <- y.data
##18 y.data <- data[, gbm.y]  #resvar?
##74 u_i <- rep(u_i, length(y_i)
##73 u_i <- sum(y.data * site.weights)/sum(site.weights)
##3 site.weights = rep(1, nrow(data))
##23 n.cases <- nrow(data)  == 119

##71 n.fitted <- n.trees
##5 n.trees = 50 (set by user? Think this is default

##5 max.trees = 10000 (set by user? think this is default)

# "delta.deviance > tolerance.test & n.fitted < max.trees"
# "n.fitted < max.trees" = 50<10000 = true

# "delta.deviance > tolerance.test" = 1 > tolerance.test
# tolerance.test: 0.001 * total.deviance/119

## calc.deviance(obs, pred, weights = rep(1,length(obs)), family="binomial", calc.mean = TRUE)
# calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)
calc.deviance(mysamples[,4], rep(sum(mysamples[,4] * rep(1,119)/119), length(mysamples[,4])), weights = rep(1,119), family = "bernoulli", calc.mean = FALSE)
# [1] NaN    Warning message: In log(1 - u_i) : NaNs produced
# why is u_i being called? it's native to gbm.step not calc.deviance & I don't see how i've called it...
calc.deviance
##12 log(1 - u_i)
##8 u_i <- pred
log(1-rep(sum(mysamples[,4] * rep(1,119)/119), length(mysamples[,4])))
##error confirmed
1-rep(sum(mysamples[,4] * rep(1,119)/119), length(mysamples[,4]))
rep(sum(mysamples[,4] * rep(1,119)/119), length(mysamples[,4]))

# log(1-u_i) in calc.deviance: u_i, which is preds, can't be positive, because when inverted (1-u_i)
# it's then negative, & can't be logged. So calc.deviance::pred is the culprit. Probably.
# Which is u_i in gbm.step. To reiterate:
##73 u_i <- sum(y.data * site.weights)/sum(site.weights)
##74 u_i <- rep(u_i, length(y_i)

##73 u_i <- sum(y.data * site.weights)/sum(rep(1,119))
##73 u_i <- sum(y.data * rep(1,119))/119
##73 u_i <- sum(mysamples[,brvcol] * rep(1,119))/119

###All binary results are 1: mysamples[,brvcol]
### So there are no zeroes
### see also original CPUE: mysamples[,4]
### So these are pre-filtered for positives
### you twat###
# Go back & get the original data inc zeroes.
# Done: mothersheet already done, I was using child sheets I never shoulda made.
####end test section####