"gbm.auto" <-
  function (grids,                        # explantory data to predict to. Import with (e.g.) read.csv and specify object name. Default is mygrids
            samples,                      # explanatory & response variables to predict from. Keep col names short, no odd characters, starting numerals or terminal periods. Spaces may be converted to periods in directory names, underscores won't. Can be a subset. Default is mysamples
            expvar,                       # list of column numbers of explanatory variables in 'samples', expected e.g. c(1,35,67,etc.). No default
            resvar,                       # column number of response variable (e.g. CPUE) in samples. Expected, e.g. 94. No default
            tc = c(2,5),                  # list of permutations of tree complexity allowed, expected e.g. and default: c(2,5)
            lr = c(0.01,0.005),           # list of permutations of learning rate allowed, expected e.g. and default: c(0.01,0.005)
            bf = 0.5,                     # list permutations of bag fraction allowed, expected e.g. & default: 0.5
            ZI = TRUE,                    # are data zero-inflated? TRUE/FALSE/"CHECK". If TRUE do delta BRT, log-normalised Gaussian, later reverse log-normalised & bias corrected. If FALSE do Gaussian only, no log-normalisation. CHECK: Tests data for you. Default is TRUE.
            gridslat = 2,                 # column number for latitude in 'grids'
            gridslon = 1,                 # column number for longitude in 'grids'
            cols = grey.colors(6,1,1),    # vector of colours for plots. Assignment in order of explanatory variables.
            savegbm = TRUE,               # save the gbm objects externally? Can reopen later with (e.g.) load("Bin_Best_Model")
            RSB = FALSE,                  # run representativeness surface builder?
            map = FALSE,                  # save abundance map png files?
            legendtitle = "CPUE",         # gbm.map metric of abundance for legend title, from legend.grid in mapplots
            ...)                          # additional parameters, e.g for gbm.map (mainmap, heatcol, shape, landcol, legendcol, bg). VALID?
{
# Generalised Boosting Model / Boosted Regression Tree process chain automater.
# Simon Dedman, 2014, simondedman@gmail.com, https://github.com/SimonDedman/gbm.auto

# Function to automate the many steps required to use boosted regression trees to predict abundances in a delta process,
# i.e. binary (0/1) proportion prediction coupled with presence-only abundance prediction to give total prediction.
# Loops through all permutations of parameters provided (learning rate, tree complexity, bag fraction), chooses the best,
# then tries to simplify that. Generates line, dot & bar plots, and outputs these and the predictions and a report of all
# variables used, statistics for tests, variable interactions, predictors used & dropped, etc.. If selected, generates
# predicted abundance maps, and representativeness surfaces.
#
# Underlying functions are from packages gbm and dismo, functions from Elith et al. 2008 (bundled as gbm.utils.R), mapplots,
# and my own functions gbm.map and gbm.rsb

####1. Check packages, start loop####
if (!require(gbm)) {stop("you need to install the gbm package to run this function")}
if (!require(dismo)) {stop("you need to install the dismo package to run this function")}
if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
if (!require(labeling)) {stop("you need to install the labeling package to run this function")}
if(map==TRUE) if (!require(mapplots)) {stop("you need to install the mapplots package to run this function")}
if(RSB==TRUE) if (!exists("gbm.rsb")) {stop("you need to install the gbm.rsb function to run this function")}
if(RSB==TRUE) if (!exists("gbm.map")) {stop("you need to install the gbm.map function to run this function")}
if(!exists("gbm.predict.grids")) {stop("you need to install the gbm.predict.grids function from gbm.utils.R to run this function")}
if(!exists("roc")) {stop("you need to install the roc function from gbm.utils.R to run this function")}
if(!exists("calibration")) {stop("you need to install the calibration function from gbm.utils.R to run this function")}
require(gbm)
require(dismo)
require(beepr)
require(labeling)
#options(error = function() {beep(9)})
expvarnames<-names(samples[expvar]) # list of explanatory variable names
expvarcols<-cbind(cols[1:length(expvarnames)],expvarnames) # assign explanatory variables to colours

for(i in resvar){
m=1 # jkl combo loop counter to allow best bin/gaus BRT choice
n=1   # Print counter for all loops of BRT combos (i.e. counter for l)
if(!all(expvarnames %in% names(grids))) {stop("Not all expvar column names found as column names in grids")}

####2. ZI check & log TODO####
# logs of resvar for ZI data & reverse lognormalise & bias correct later (line 117 gaus BRT gbm.y & sections 22 & 23)
# asked here: https://stats.stackexchange.com/questions/112292/r-are-data-zero-inflated
# bother with this or just have an element in the function function(x,...,zeroinflated) where zeroinflated can be Y or N?
# ifelse(samples[resvar]=zeroinflated, #test
#  samples[paste("grv_",names(samples[resvar]),"4model",sep="")] <- log1p(samples[resvar]), #yes
#  samples[paste("grv_",names(samples[resvar]),"4model",sep="")] <- samples[resvar]) #no
# pscl package

# Zero-inflation is about the shape of the distribution. Therefore, you will have to specify the distribution for the non-zero part (Poisson, Negative Binomial, etc), if you want a formal test. Then you can use a likelihood ratio test to see if the zero-inflated parameters can be dropped from the model. This can be done in R.
# In cruder terms, zero inflation is defined not only by proportion of zeros but also by the total number of observations. Say, if you assume a zero-inflated Poisson model and your data contain 50% of zeros, you still won't be able to say with certainty that it's zero inflated if the total number of points is only 4. On the other hand, 10% of zeros in 1000 observations can result in a positive test for zero-inflation.
# Zero-inflated property is associated with count-based data, so I haven't heard of "zero-inflated normal". E.g. in this package:

# fit Poisson model where the zero and non-zero components contain only the intercept
# then
# check if the intercept from the zero component has a significant p-value.

# remember: if NOT ZI then it's ONLY the gaussian BRTs and NOT log-normalised. Need to work out how to do that
# maybe just, for BIN BRT: if(ZI) {run as is} (no else, don't run it, move on)
# for Gaus: samples have already been logged or not depending on ZI-ness so no change
# at end: if(ZI) {unlog procedure} else {maybe do nothing? check}
# ZI entered in function above, this is auto entered by the test code only if check enabled
# if(ZI=="CHECK") {run ZI test resulting in ZI=TRUE/FALSE} (else nothing: continue with ZI as is)

# test to make sure response variable has zeroes (prevent user wasting hours debugging nested functions
# only to find it's failing because they made a stupid mistake. Like I did!)
if(min(samples[i])>0) print("No zeroes in response variable. Method expects unsuccessful, as well as successful, samples")
# create binary (0/1) response variable, for bernoulli BRTs
samples$brv <- ifelse(samples[i] > 0, 1, 0)
brvcol <- which(colnames(samples)=="brv") # brv column number for BRT

# create logged response variable, for gaussian BRTs when data is zero-inflated (otherwise just use resvar directly)
logem <- log1p(samples[,i])
dont  <- samples[,i]
if (ZI==TRUE){samples$grv <- logem} else {samples$grv <- dont}
grvcol <- which(colnames(samples)=="grv") # grv column number for BRT

# presence-only subsets for gaussian element of delta BRTs
grv_yes <- subset(samples,grv > 0)

####3. Begin Report####
reportcolno = (3+(5*(length(tc)*length(lr)*length(bf)))+14)  # calculate number of columns for report: standard elements (3&12) + 5 elements for each loop
Report <- data.frame(matrix(NA, nrow = (max(6,length(expvar))), ncol = (reportcolno)))  # build blank df, rows=biggest of 6 (max static row number of stats) or n of exp. vars
colnames(Report) <- c("Explanatory Variables","Response Variables","Zero Inflated?") # populate static colnames
colnames(Report)[(reportcolno-13):reportcolno] <- c("Best Binary BRT", "Best Gaussian BRT", "Bin_BRT_simp predictors kept (ordered)",
                                                     "Bin_BRT_simp predictors dropped", "Gaus_BRT_simp predictors kept (ordered)",
                                                     "Gaus_BRT_simp predictors dropped", "Simplified Binary BRT stats", "Simplified Gaussian BRT stats",
                                                     "Best Binary BRT variables", "Relative Influence (Bin)", "Best Gaussian BRT variables",
                                                     "Relative Influence (Gaus)", "Biggest Interactions (Bin)", "Biggest Interactions (Gaus)")
Report[1:length(expvar),1] <- names(samples[expvar])
Report[1,2] <- names(samples[i])
Report[1,3] <- ZI

####4. Binomial BRT####
for(j in tc){   # list permutations of tree complexity allowed
 for(k in lr){   # list permutations of learning rate allowed
  for(l in bf){   # list permutations of bag fraction allowed

(assign(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""),gbm.step(data=samples,
    gbm.x = expvar, gbm.y = brvcol, family = "bernoulli", tree.complexity = j, learning.rate = k, bag.fraction = l)))

####5. Select best Bin model####
# Makes an object with BRT training data correlation score of the 1st variable combo in the loop, & the name of that BRT combo,
# per species (m is reset at the end of each species loop so m=1 is the first loop per species)
# Then if subsequent combos have higher values (more correlation), it replaces the value, and name of combo it came from

# create blanks for best results
Bin_Best_Score <- 0
Bin_Best_Model <- 0
Gaus_Best_Score <- 0
Gaus_Best_Model <- 0

if(m==1)
  {Bin_Best_Score <- get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
  Bin_Best_Model <- paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
}  else if(get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]>Bin_Best_Score)
      {Bin_Best_Score <- get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
      Bin_Best_Model <- paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep="")}

# progress printer, right aligned
print(paste("
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   Completed BRT ",n," of ",2*length(i)*length(tc)*length(lr)*length(bf),"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
n <- n+1   # Add to print counter

####6. Gaussian BRT####
assign(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""),gbm.step(data=samples,
    gbm.x = expvar, gbm.y = grvcol, family = "gaussian", tree.complexity = j, learning.rate = k, bag.fraction = l))

####7. Select best Gaus model####
if(m==1)
  {Gaus_Best_Score <- get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
  Gaus_Best_Model <- paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
} else if(get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]>Gaus_Best_Score)
      {Gaus_Best_Score <- get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
      Gaus_Best_Model <- paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep="")}

####8. Add BRT stats to report####
Report[1:3,((m*5)-1)] <- c(paste("tree complexity: ",j,sep=""), paste("learning rate: ",k,sep=""), paste("bag fraction: ",l,sep=""))
Report[1:6,(m*5)] <- c(paste("trees: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
                       paste("Training Data Correlation: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]],sep=""),
                       paste("CV Mean Deviance: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.mean,sep=""),
                       paste("CV Deviance SE: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.se,sep=""),
                       paste("CV Mean Correlation: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.mean,sep=""),
                       paste("CV Correlation SE: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.se,sep=""))
Report[1,(m*5)+1] <- paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
Report[1:6,((m*5)+2)] <- c(paste("trees: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
                            paste("Training Data Correlation: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]],sep=""),
                            paste("CV Mean Deviance: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.mean,sep=""),
                            paste("CV Deviance SE: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.se,sep=""),
                            paste("CV Mean Correlation: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.mean,sep=""),
                            paste("CV Correlation SE: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.se,sep=""))
Report[1,(m*5)+3] <- paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
colnames(Report)[((m*5)-1):((m*5)+3)] <- c(paste("Parameter Combo ",m,sep=""),
                                           paste("Binary BRT ",m," stats",sep=""),
                                           paste("Binary BRT ",m," name",sep=""),
                                           paste("Gaussian BRT ",m," stats",sep=""),
                                           paste("Gaussian BRT ",m," name",sep=""))

# progress printer, right aligned for visibility
print(paste("
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ",2*length(i)*length(tc)*length(lr)*length(bf),"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
n <- n+1   # Add to print counter: 2 per loop, 1 bin 1 gaus BRT
m <- m+1   # Add to loop counter: 1 per loop, used for bin/gaus_best model selection
}}}        # close loops, producing all BRT/GBM objects then continuing through model selection

####9. Test simplification benefit, do so if better####
samples<<-samples # globally assign samples: bad practice but fixes an odd problem whereby the code runs if done manually
# but crashes when in a huge looping function, saying it's unable to access samples.
Bin_Best_Simp_Check <- gbm.simplify(get(Bin_Best_Model)) # run the simplification tester on the best model

# if best number of variables to remove isn't 0 (i.e. it's worth simplifying), re-run the best model (Bin_Best_Model, using gbm.call to get its values)
# with just-calculated best number of variables to remove, removed. gbm.x asks which number of drops has the minimum mean (lowest point on the line)
# & that calls up the list of predictor variables with those removed, from $pred.list

if(min(Bin_Best_Simp_Check$deviance.summary$mean) < 0)
  assign("Bin_Best_Simp", gbm.step(data = samples,
                                 gbm.x = Bin_Best_Simp_Check$pred.list[[which.min(Bin_Best_Simp_Check$deviance.summary$mean)]],
                                 gbm.y = get(Bin_Best_Model)$gbm.call$gbm.y,
                                 tree.complexity = get(Bin_Best_Model)$gbm.call$tree.complexity,
                                 learning.rate = get(Bin_Best_Model)$gbm.call$learning.rate,
                                 family = get(Bin_Best_Model)$gbm.call$family,
                                 bag.fraction = get(Bin_Best_Model)$gbm.call$bag.fraction))

# progress printer, right aligned for visibility
print(paste("
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Bin model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))

# Same for Gaus
Gaus_Best_Simp_Check <- gbm.simplify(get(Gaus_Best_Model))
if(min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0)
  assign("Gaus_Best_Simp", gbm.step(data = samples,
                                 gbm.x = Gaus_Best_Simp_Check$pred.list[[which.min(Gaus_Best_Simp_Check$deviance.summary$mean)]],
                                 gbm.y = get(Gaus_Best_Model)$gbm.call$gbm.y,
                                 tree.complexity = get(Gaus_Best_Model)$gbm.call$tree.complexity,
                                 learning.rate = get(Gaus_Best_Model)$gbm.call$learning.rate,
                                 family = get(Gaus_Best_Model)$gbm.call$family,
                                 bag.fraction = get(Gaus_Best_Model)$gbm.call$bag.fraction))

# progress printer, right aligned for visibility
print(paste("
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Gaus model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))

####10. Select final best models####
# if Bin_Best has a simplified model:
if(min(Bin_Best_Simp_Check$deviance.summary$mean) < 0)
# & if the simplified model has better correlation than Bin_Best itself
  if(Bin_Best_Simp$self.statistics$correlation > Bin_Best_Score[1])
# then replace Bin_Best score/model values with those from the simplified model
  {Bin_Best_Score <- Bin_Best_Simp$self.statistics$correlation
   Bin_Best_Model <- "Bin_Best_Simp"}

Bin_Best_Model<<-get(Bin_Best_Model) # globally assign final model for external testing later

# Same for Gaus:
if(min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0)
  if(Gaus_Best_Simp$self.statistics$correlation > Gaus_Best_Score[1])
    {Gaus_Best_Score <- Gaus_Best_Simp$self.statistics$correlation
     Gaus_Best_Model <- "Gaus_Best_Simp"}

Gaus_Best_Model<<-get(Gaus_Best_Model) # globally assign final model for external testing later

####11. Line plots####
dir.create(names(samples[i])) # create resvar-named directory for outputs
opar<-par() # save par defaults
# par(mar=c(10,2,2,2), fig=c(0,1,0.5,1), cex.lab=0.8) # grow margins so labels fit. This bit redundant: needs to be between png() & gbm.plot()

# All plots on one image for Bin & Gaus
png(filename = paste("./",names(samples[i]),"/Bin_Best_line.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
gbm.plot(get(Bin_Best_Model),
         n.plots=length(get(Bin_Best_Model)$contributions$var),
         write.title = F, y.label = "Marginal Effect",
         plot.layout = c(ceiling(sqrt(length(get(Bin_Best_Model)$contributions$var))),
             ifelse(sqrt(length(get(Bin_Best_Model)$contributions$var))
              - floor(sqrt(length(get(Bin_Best_Model)$contributions$var))) <0.5,
              floor(sqrt(length(get(Bin_Best_Model)$contributions$var))),
              floor(sqrt(length(get(Bin_Best_Model)$contributions$var)))+1)))
dev.off()

png(filename = paste("./",names(samples[i]),"/Gaus_Best_line.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
gbm.plot(get(Gaus_Best_Model),
         n.plots=length(get(Gaus_Best_Model)$contributions$var),
         write.title = F, y.label = "Marginal Effect",
         plot.layout = c(ceiling(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
             ifelse(sqrt(length(get(Gaus_Best_Model)$contributions$var))
              - floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))) <0.5,
              floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
              floor(sqrt(length(get(Gaus_Best_Model)$contributions$var)))+1)))
dev.off()

# All plots individually, named by explanatory variable, bin & gaus
for (o in 1:length(get(Bin_Best_Model)$contributions$var)){
  png(filename = paste("./",names(samples[i]),"/Bin_Best_line_",as.character(get(Bin_Best_Model)$contributions$var[o]),".png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(1.35,3.4,0.4,0.5), fig=c(0,1,0,1), las=1, lwd=8, bty="n", mgp=c(2,0.5,0), xpd=NA)  #mar=c(1.35,2.4,0.4,0.5)
    # bg=expvarcols[match(get(Bin_Best_Model)$contributions$var[o],expvarcols[,2]),1]) #changed margin to hide label #XPD YPD ALLOWS AXES TO EXTEND FURTHER TO ENCOMPASS ALL DATA? #colour removed
plotgrid<-plot.gbm(get(Bin_Best_Model),match(get(Bin_Best_Model)$contributions$var[o], get(Bin_Best_Model)$gbm.call$predictor.names),lwd=8,return.grid=TRUE)
xx <- labeling::extended(min(plotgrid[1]), max(plotgrid[1]),7, only.loose=TRUE) # sets range & ticks
yy <- labeling::extended(min(plotgrid[2]), max(plotgrid[2]),7, only.loose=TRUE) # sets range & ticks
plot(range(xx),range(yy), t="n", xaxt="n", yaxt="n", bty="n", ylab=NA)
lines(plotgrid,type="l")
axis(1, lwd.ticks=8, lwd=8, at=xx) # is providing only the thick line & downticks
axis(2, lwd.ticks=8, lwd=8, at=yy)
#rug(quantile(samples[as.character(get(Bin_Best_Model)$contributions$var[o])], probs=seq(0,1,0.01), na.rm=TRUE), side=1, lwd=5, ticksize=0.03) #n of ticks probs seq arg 3: 0.1, 0.05, 0.01
rug(samples[as.character(get(Bin_Best_Model)$contributions$var[o])][,1], side=1, lwd=5, ticksize=0.03) # all points rug
dev.off() }

for (p in 1:length(get(Gaus_Best_Model)$contributions$var)){
  png(filename = paste("./",names(samples[i]),"/Gaus_Best_line_",as.character(get(Gaus_Best_Model)$contributions$var[p]),".png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = "cairo-png")
  #par(mar=c(2.9,2,0.4,0.5), fig=c(0,1,0,1), las=1, lwd=8, bty="n", mgp=c(2,0.5,0))
par(mar=c(1.35,3.4,0.4,0.5), fig=c(0,1,0,1), las=1, lwd=8, bty="n", mgp=c(2,0.5,0), xpd=NA)
    # bg=expvarcols[match(get(Gaus_Best_Model)$contributions$var[p],expvarcols[,2]),1]) #changed margin to hide label #XPD YPD ALLOWS AXES TO EXTEND FURTHER TO ENCOMPASS ALL DATA? #colour removed
plotgrid<-plot.gbm(get(Gaus_Best_Model),match(get(Gaus_Best_Model)$contributions$var[p], get(Gaus_Best_Model)$gbm.call$predictor.names),lwd=8,return.grid=TRUE)
xx <- labeling::extended(min(plotgrid[1]), max(plotgrid[1]),7, only.loose=TRUE) # sets range & ticks
yy <- labeling::extended(min(plotgrid[2]), max(plotgrid[2]),7, only.loose=TRUE) # sets range & ticks
plot(range(xx),range(yy), t="n", xaxt="n", yaxt="n", bty="n", ylab=NA)
lines(plotgrid,type="l")
axis(1, lwd.ticks=8, lwd=8, at=xx) # is providing only the thick line & downticks
axis(2, lwd.ticks=8, lwd=8, at=yy)
#rug(quantile(samples[as.character(get(Gaus_Best_Model)$contributions$var[p])], probs=seq(0,1,0.01), na.rm=TRUE), side=1, lwd=5, ticksize=0.03) #n of ticks probs seq arg 3: 0.1, 0.05, 0.01
rug(samples[as.character(get(Gaus_Best_Model)$contributions$var[p])][,1], side=1, lwd=5, ticksize=0.03) # all points rug
dev.off() }

####12. Dot plots####
png(filename = paste("./",names(samples[i]),"/Bin_Best_dot.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
gbm.plot.fits(get(Bin_Best_Model))
#   unused? , plot.layout = c(ceiling(sqrt(length(get(Bin_Best_Model)$contributions$var))), #plot.layout unused argument?
#              ifelse(sqrt(length(get(Bin_Best_Model)$contributions$var))
#               - floor(sqrt(length(get(Bin_Best_Model)$contributions$var))) <0.5,
#               floor(sqrt(length(get(Bin_Best_Model)$contributions$var))),
#               floor(sqrt(length(get(Bin_Best_Model)$contributions$var)))+1))
dev.off()

png(filename = paste("./",names(samples[i]),"/Gaus_Best_dot.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
gbm.plot.fits(get(Gaus_Best_Model))
#   unused? , plot.layout = c(ceiling(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
#              ifelse(sqrt(length(get(Gaus_Best_Model)$contributions$var))
#               - floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))) <0.5,
#               floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
#               floor(sqrt(length(get(Gaus_Best_Model)$contributions$var)))+1))

dev.off()

####13. 3D plot TODO####
# gbm.perspec(Bin_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)
# gbm.perspec(Gaus_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)

####14. Bar plots of variable influence####
# create tables
Bin_Bars <- summary(get(Bin_Best_Model),
        cBars=length(get(Bin_Best_Model)$var.names),
        n.trees=get(Bin_Best_Model)$n.trees,
        plotit=FALSE, order=TRUE, normalize=TRUE, las=1, main=NULL)
write.csv(Bin_Bars, file=paste("./",names(samples[i]),"/Binary BRT Variable contributions.csv",sep=""), row.names=FALSE)

Gaus_Bars <- summary(get(Gaus_Best_Model),
        cBars=length(get(Gaus_Best_Model)$var.names),
        n.trees=get(Gaus_Best_Model)$n.trees,
        plotit=FALSE, order=TRUE, normalize=TRUE, las=1, main=NULL)
write.csv(Gaus_Bars, file=paste("./",names(samples[i]),"/Gaussian BRT Variable contributions.csv",sep=""), row.names=FALSE)

# produce graphics
png(filename = paste("./",names(samples[i]),"/Bin_Bars.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "",
    type = "cairo-png")
par(mar=c(2.5,0.3,0,0.5), fig=c(0,1,0,1), cex.lab=0.8,mgp=c(1.5,0.5,0), cex=1.3)
midpoints<-barplot(rev(Bin_Bars[,2]), cex.lab=1.2, las=1, horiz=TRUE, cex.names=0.8, xlab="Influence %", col=rev(expvarcols[match(Bin_Bars[,1],expvarcols[,2]),1]), xlim=c(0,2.5+ceiling(max(Bin_Bars[,2]))))
text(0.1, midpoints,labels=rev(Bin_Bars[,1]),adj=0,cex=1.5)
axis(side = 1, lwd = 6, outer=TRUE, xpd=NA)
dev.off()

png(filename = paste("./",names(samples[i]),"/Gaus_Bars.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "",
    type = "cairo-png")
par(mar=c(2.5,0.3,0,0.5), fig=c(0,1,0,1), cex.lab=0.8,mgp=c(1.5,0.5,0), cex=1.3)
midpoints<-barplot(rev(Gaus_Bars[,2]), cex.lab=1.2, las=1, horiz=TRUE, cex.names=0.8, xlab="Influence %", col=rev(expvarcols[match(Gaus_Bars[,1],expvarcols[,2]),1]), xlim=c(0,2.5+ceiling(max(Gaus_Bars[,2]))))
text(0.1, midpoints,labels=rev(Gaus_Bars[,1]),adj=0,cex=1.5)
axis(side = 1, lwd = 6, outer=TRUE, xpd=NA)
dev.off()

####15. Variable interactions####
find.int_Bin <- gbm.interactions(get(Bin_Best_Model))
find.int_Gaus <- gbm.interactions(get(Gaus_Best_Model))

####16. Binomial predictions####
gbm.predict.grids(get(Bin_Best_Model), grids, want.grids = F, sp.name = "Bin_Preds")
grids$Bin_Preds <- Bin_Preds

####17. Gaussian predictions####
gbm.predict.grids(get(Gaus_Best_Model), grids, want.grids = F, sp.name = "Gaus_Preds")
grids$Gaus_Preds <- Gaus_Preds

####18. Backtransform logged Gaus to unlogged####
grids$Gaus_Preds_Unlog <- exp(Gaus_Preds + 1/2 * sd(get(Gaus_Best_Model)$residuals,na.rm=FALSE)^2)

####19. BIN*positive abundance = final abundance####
grids$PredAbund <- grids$Gaus_Preds_Unlog * grids$Bin_Preds # add column in grids for predicted abundance
predabund <- which(colnames(grids)=="PredAbund") # predicted abundance column number for writecsv

####20. Final saves####
# should names(samples[i]) be something else? This is currently "CPUE" but should be e.g. "Blonde Ray"...
# CSV of Predicted values at each site inc predictor variables' values.
write.csv(grids,row.names=FALSE, file = paste("./",names(samples[i]),"/Abundance_Preds_All.csv",sep=""))
# CSV of Predicted values at each site without predictor variables' values.
write.csv(grids[c(gridslat,gridslon,predabund)], row.names=FALSE, file = paste("./",names(samples[i]),"/Abundance_Preds_only.csv",sep=""))
if (savegbm==TRUE) {save(Bin_Best_Model,file = paste("./",names(samples[i]),"/Bin_Best_Model",sep=""))
                    save(Gaus_Best_Model,file = paste("./",names(samples[i]),"/Gaus_Best_Model",sep=""))}

####21. Finalise & Write Report####
Report[1:2,(reportcolno-13)] <- c(paste("Model combo: ",Bin_Best_Model,sep=""),paste("Model CV score: ",Bin_Best_Score,sep=""))
Report[1:2,(reportcolno-12)] <- c(paste("Model combo: ",Gaus_Best_Model,sep=""),paste("Model CV score: ",Gaus_Best_Score,sep=""))
Report[1:dim(subset(Bin_Best_Simp_Check$final.drops,order>0))[1],(reportcolno-11)] <- as.character(subset(Bin_Best_Simp_Check$final.drops,order>0)$preds)
Report[1:(length(Bin_Best_Simp_Check$final.drops$preds)-dim(subset(Bin_Best_Simp_Check$final.drops,order>0))[1]),(reportcolno-10)] <-
  as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order>0))[1])+1):length(Bin_Best_Simp_Check$final.drops$preds)])
Report[1:dim(subset(Gaus_Best_Simp_Check$final.drops,order>0))[1],(reportcolno-9)] <- as.character(subset(Gaus_Best_Simp_Check$final.drops,order>0)$preds)
Report[1:(length(Gaus_Best_Simp_Check$final.drops$preds)-dim(subset(Gaus_Best_Simp_Check$final.drops,order>0))[1]),(reportcolno-8)] <-
  as.character(Gaus_Best_Simp_Check$final.drops$preds[((dim(subset(Gaus_Best_Simp_Check$final.drops,order>0))[1])+1):length(Gaus_Best_Simp_Check$final.drops$preds)])
Report[1:6,(reportcolno-7)] <-   c(paste("trees: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
                                   paste("Training Data Correlation: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]],sep=""),
                                   paste("CV Mean Deviance: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.mean,sep=""),
                                   paste("CV Deviance SE: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.se,sep=""),
                                   paste("CV Mean Correlation: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.mean,sep=""),
                                   paste("CV Correlation SE: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.se,sep=""))
Report[1:6,(reportcolno-6)] <- c(paste("trees: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
                                 paste("Training Data Correlation: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]],sep=""),
                                 paste("CV Mean Deviance: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.mean,sep=""),
                                 paste("CV Deviance SE: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.se,sep=""),
                                 paste("CV Mean Correlation: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.mean,sep=""),
                                 paste("CV Correlation SE: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.se,sep=""))
Report[1:(length(Bin_Bars[,1])),(reportcolno-5)] <- as.character(Bin_Bars$var)
Report[1:(length(Bin_Bars[,2])),(reportcolno-4)] <- as.character(Bin_Bars$rel.inf)
Report[1:(length(Gaus_Bars[,1])),(reportcolno-3)] <- as.character(Gaus_Bars$var)
Report[1:(length(Gaus_Bars[,2])),(reportcolno-2)] <- as.character(Gaus_Bars$rel.inf)
Report[1:2,(reportcolno-1)] <- c(paste(find.int_Bin$rank.list$var1.names[1]," and ",find.int_Bin$rank.list$var2.names[1],". Size: ",find.int_Bin$rank.list$int.size[1],sep=""),
                                 paste(find.int_Bin$rank.list$var1.names[2]," and ",find.int_Bin$rank.list$var2.names[2],". Size: ",find.int_Bin$rank.list$int.size[2],sep=""))
Report[1:2,(reportcolno)] <- c(paste(find.int_Gaus$rank.list$var1.names[1]," and ",find.int_Gaus$rank.list$var2.names[1],". Size: ",find.int_Gaus$rank.list$int.size[1],sep=""),
                               paste(find.int_Gaus$rank.list$var1.names[2]," and ",find.int_Gaus$rank.list$var2.names[2],". Size: ",find.int_Gaus$rank.list$int.size[2],sep=""))
write.csv(Report, row.names=FALSE, na="", file = paste("./",names(samples[i]),"/Report.csv",sep=""))

####22. Representativeness surface builder####
# does NOT plot the surface, just builds it. If built, it will be plotted by map maker.
if(RSB==TRUE){rsbdf <- gbm.rsb(samples,grids,expvarnames,gridslat,gridslon,rsbres=i)}

####23. Map maker####
if(map==TRUE) {
  # generate output image & set parameters
  png(filename = paste("./",names(samples[i]),"/PredAbundMap_",names(samples[i]),".png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
  # run gbm.map function with generated parameters
  gbm.map(x = grids[,gridslon],
               y = grids[,gridslat],
               z = grids[,predabund],
#               byx = byx, #default now null so shouldn't need to provide anything here
#               byy = byy,
               mapmain = "Predicted abundance: ",
               species = names(samples[i]),
               shape = coast,
               landcol = "darkgreen",
               legendloc = "bottomright",
               legendtitle = legendtitle,
               grids=grids,
               gridslon=gridslon,
               gridslat=gridslat,
               predabund=predabund)  # hopefully parses grids dataset to gbm.map to use
  dev.off()

  # if RSB called, plot that surface separately
  if(RSB==TRUE){
    png(filename = paste("./",names(samples[i]),"/RSB_Map_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf[,"Unrepresentativeness"],
#            byx = byx,  # both will have been created by gbm.map if they weren't provided by the user
#            byy = byy,
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            #heatcol= rgb(seq(from=1,to=0,length.out=12),seq(from=1,to=0,length.out=12),seq(from=1,to=0,length.out=12)), # rgb(0,0,0,seq(from=0,to=0.5,length.out=12)) black, from opaque to transparent
            shape = coast,
            landcol = "darkgreen",
            legendloc = "bottomright",
            legendtitle = "UnRep 0-1",
            grids=grids,
            gridslon=gridslon,
            gridslat=gridslat,
            predabund=predabund)
    dev.off()
   } # close RSB mapper
  } # close Map Maker
 } # close response variable (resvar) loop
####END####
beep(8) # notify the user with a noise, since this process can take a long time.
} # close the function