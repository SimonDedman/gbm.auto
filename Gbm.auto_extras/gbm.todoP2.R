## gbm.todo.P2
# working script for paper 2
# 1/6/2015
##
####TODO####
# have per-species cons maps underpin the Conservation Sort in p3, rather than
# all 4 species underpinned by the 4-sp gbm.cons map.

####prep####
# set & test working directory
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles")
getwd()
install_github("SimonDedman/gbm.auto")

# load grids
mygrids <- gbm.auto::grids

# load saved models if re-running aspects from a previous run
# load("Bin_Best_Model")
# load("Gaus_Best_Model")

####model:juve individual preds####
# load samples
samples <- gbm.auto::Juveniles

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
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators")

#cuckoo original without ComSkt_C since it's empty
names(mysamples)[32]
gbm.auto(expvar = c(4:10,14,16,20,24,28,36),
         resvar = c(43),
         grids = mygrids,
         samples = mysamples,
         tc = c(2,13),
         lr = c(0.005, 0.001),
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB = TRUE,
         varint = FALSE,
         zero = FALSE)

#thornback
max((mysamples)[33]) # ComSkt_T max = 0, remove it
gbm.auto(expvar = c(4:10,14,17,21,25,29,37),
         resvar = c(44),
         grids = mygrids,
         samples = mysamples,
         tc = c(2,13),
         lr = c(0.01, 0.005),
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB = TRUE,
         varint = FALSE,
         zero = FALSE)
# 27/2/16 worked

# Blonde
gbm.auto(expvar = c(4:10,14,18,22,26,30),
         resvar = c(45),
         grids = mygrids,
         samples = mysamples,
         tc = c(12),
         lr = c(0.005),
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB = TRUE,
         varint = FALSE,
         zero = FALSE)
# 27/2/16 worked

#17/08/15 succeeded
####12/08/2015 FAIL####
##Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) :
##contrasts can be applied only to factors with 2 or more levels In addition: There were 50 or more warnings (use warnings() to see the first 50)
# commented out lines 361 & 362, find interactions in gbm.auto, causing the problem.
# also 413-416, which includes the results in the report

#spotted
max((mysamples)[35]) # ComSkt_B max = 0, remove it
gbm.auto(expvar = c(4:10,14,19,23,27,31,38),
         resvar = c(46),
         grids = mygrids,
         samples = mysamples,
         tc = c(2,13),
         lr = c(0.01, 0.005),
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB = TRUE,
         varint = FALSE,
         zero = FALSE)
# FAIL: contrasts can be applied only to factors with 2 or more levels

####model: mat F####
# set wd for mature female samples sheets
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F")
# load samples
mysamples <- gbm.auto::Adult_Females

# run models: cuckoo
gbm.auto(expvar = 4:9,
         resvar = 10,
         grids = mygrids,
         samples = mysamples,
         tc = c(2,6),
         lr = c(0.005),
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB = TRUE,
         varint = FALSE,
         zero = FALSE)

#thornback: fails @lr=0.01
gbm.auto(expvar = 4:9,
         resvar = 11,
         grids = mygrids,
         samples = mysamples,
         tc = c(2,6),
         lr = c(0.005),
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB= TRUE,
         varint = FALSE,
         zero = FALSE)

#blonde
gbm.auto(expvar = 4:9,
         resvar = 12,
         grids = mygrids,
         samples = mysamples,
         tc = c(6), #2, #report shows 6 is best
         lr = c(0.001), #0.01, 0.005
         bf = c(0.5),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         map = TRUE,
         RSB = TRUE,
         varint = FALSE,
         zero = FALSE)

# 12/08/2015:
1: In cor(y_i, u_i) : the standard deviation is zero
2: glm.fit: algorithm did not converge
3: glm.fit: fitted probabilities numerically 0 or 1 occurred
4: glm.fit: fitted probabilities numerically 0 or 1 occurred
# seems to have worked though?

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

#spotted #failed @lr=0.01
gbm.auto(expvar = 4:9,
         samples = Adult_Females,
         resvar = 13,
         grids = mygrids,
         tc = c(2,6),
         lr = c(0.005),
         gridslat = 2,
         gridslon = 1,
         ZI = TRUE,
         varint = FALSE,
         zero = FALSE)

####Conservation maps####
# simply add Abundance_Preds_only.csv[,Predabund] for juve individual preds & mat Fs
# then use as z in gbm.map
# where it is currently: setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F")
setwd("../") # go up to /Maps
dir.create("ConservationMaps") # create conservation maps directory

for (i in c("Cuckoo","Blonde","Thornback","Spotted")) {
  juves <- read.csv(paste("./Juveniles/Individual Predators/",i,"/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  matfs <- read.csv(paste("./Mature Females plus Hans' F/", i,"/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  dir.create(paste("./ConservationMaps/", i, sep = ""))
  setwd(paste("./ConservationMaps/", i, "/", sep = ""))

# map juve abundance + matf abundance
png(filename = paste("./Conservation_Map_",i,".png", sep = ""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
gbm.map(x = juves[,2],
        y = juves[,1],
        z = juves[,3] + matfs[,3],
        mapmain = "Predicted CPUE (numbers per hour): ",
        species = i,
        zero = FALSE,
        legendtitle = "CPUE")
dev.off()

# again in B&W
png(filename = paste("./Conservation_Map_BnW",i,".png", sep = ""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
gbm.map(x = juves[,2],
        y = juves[,1],
        z = juves[,3] + matfs[,3],
        mapmain = "Predicted CPUE (numbers per hour): ",
        species = i,
        zero = FALSE,
        colournumber = 5,
        heatcolours = grey.colors(5, start = 1, end = 0),
        mapback = "white",
        legendtitle = "CPUE")
dev.off()

### Scale both to 1
  juvescale <- juves[,3] / max(juves[,3], na.rm = TRUE)
  matfscale <- matfs[,3] / max(matfs[,3], na.rm = TRUE)
  png(filename = paste("./Scale1-1_Conservation_Map_",i,".png",sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = juves[,2],
          y = juves[,1],
          z = (juvescale + matfscale) * 50, #why *50? 1+1=2 max so this makes it 100 max
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = i,
          zero = FALSE,
          breaks = c(0,20,40,60,80,100),
          colournumber = 5,
          legendtitle = "CPUE (Scaled %)")
  dev.off()

  # Do again in B&W
  png(filename = paste("./Scale1-1_Conservation_Map_BnW_",i,".png", sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = juves[,2],
          y = juves[,1],
          z = (juvescale + matfscale)*50,
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = i,
          zero = FALSE,
          breaks = c(0,20,40,60,80,100),
          colournumber = 5,
          heatcolours = grey.colors(5, start = 1, end = 0),
          mapback = "white",
          legendtitle = "CPUE (Scaled %)")
  dev.off()
  setwd("../") # go back up to ConservationMaps
  setwd("../")} # go back up to Maps & close loop
beep(8)


#### Glue scaled outputs together ####
#Can do this all algorithmically?
  juve_C <- read.csv(paste("./Juveniles/Individual Predators/Cuckoo/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  juve_T <- read.csv(paste("./Juveniles/Individual Predators/Thornback/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  juve_B <- read.csv(paste("./Juveniles/Individual Predators/Blonde/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  juve_S <- read.csv(paste("./Juveniles/Individual Predators/Spotted/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  matf_C <- read.csv(paste("./Mature Females plus Hans' F/Cuckoo/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  matf_T <- read.csv(paste("./Mature Females plus Hans' F/Thornback/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  matf_B <- read.csv(paste("./Mature Females plus Hans' F/Blonde/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  matf_S <- read.csv(paste("./Mature Females plus Hans' F/Spotted/Abundance_Preds_only.csv", sep = ""), header = TRUE)
  dir.create(paste("./ConservationMaps/Combo/", sep = ""))
  setwd(paste("./ConservationMaps/Combo/", sep = ""))

  juve_C_scale <- juve_C[,3] / max(juve_C[,3], na.rm = TRUE)
  juve_T_scale <- juve_T[,3] / max(juve_T[,3], na.rm = TRUE)
  juve_B_scale <- juve_B[,3] / max(juve_B[,3], na.rm = TRUE)
  juve_S_scale <- juve_S[,3] / max(juve_S[,3], na.rm = TRUE)
  matf_C_scale <- matf_C[,3] / max(matf_C[,3], na.rm = TRUE)
  matf_T_scale <- matf_T[,3] / max(matf_T[,3], na.rm = TRUE)
  matf_B_scale <- matf_B[,3] / max(matf_B[,3], na.rm = TRUE)
  matf_S_scale <- matf_S[,3] / max(matf_S[,3], na.rm = TRUE)

  allscaled <-  {juve_C_scale +
                juve_T_scale +
                juve_B_scale +
                juve_S_scale +
                matf_C_scale +
                matf_T_scale +
                matf_B_scale +
                matf_S_scale}

  #save as csv for later
  allscaleddf <- data.frame(LATITUDE = juve_C[,1], LONGITUDE = juve_C[,2], allscaled)
  write.csv(allscaleddf, row.names = FALSE, file = paste("./AllScaledData.csv", sep = ""))

  png(filename = "Scaled_Conservation_Map.png",
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = juves[,2],
          y = juves[,1],
          z = allscaled * 12.5, #8 max, this raises to 100
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = "All Species",
          zero = FALSE,
          legendtitle = "CPUE (Scaled %)")
  dev.off()

  # again in B&W
  png(filename = "Scaled_Conservation_Map_BnW.png",
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = juves[,2],
          y = juves[,1],
          z = allscaled * 12.5,
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = "All Species",
          zero = FALSE,
          heatcolours = grey.colors(5, start = 1, end = 0),
          mapback = "white",
          legendtitle = "CPUE (Scaled %)")
  dev.off()
  beep(8)










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
# 1 > 0.001 * total.deviance/119
# 119 > 0.001 * total.deviance
# 119000 > total.deviance
# where's total.deviance?
# gbm.step:
# total.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

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
