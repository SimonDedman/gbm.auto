"gbm.auto" <-
  function (grids = NULL,                 # explantory data to predict to. Import with (e.g.) read.csv and specify object name. Defaults to NULL (won't predict to grids)
            samples = mysamples,          # explanatory & response variables to predict from. Keep col names short, no odd characters, starting numerals or terminal periods. Spaces may be converted to periods in directory names, underscores won't. Can be a subset. Default is mysamples
            expvar,                       # list of column numbers of explanatory variables in 'samples', expected e.g. c(1,35,67,etc.). No default
            resvar,                       # column number of response variable (e.g. CPUE) in samples. Expected, e.g. 94. No default. Column name should be species name.
            #tc = c(2,5),                  # list of permutations of tree complexity allowed, expected e.g. and default: c(2,5) SHOULD DEFAULT TO 2,length(expvar)
            lr = c(0.01,0.005),           # list of permutations of learning rate allowed, expected e.g. and default: c(0.01,0.005)
            bf = 0.5,                     # list permutations of bag fraction allowed, expected e.g. & default: 0.5
            ZI = "CHECK",                 # are data zero-inflated? TRUE/FALSE/"CHECK". TRUE? do delta BRT, log-normalised Gaussian, later reverse log-normalised & bias corrected. FALSE: do Gaussian only, no log-normalisation. CHECK: Tests data for you. Default is TRUE.
            gridslat = 2,                 # column number for latitude in 'grids'
            gridslon = 1,                 # column number for longitude in 'grids'
            cols = grey.colors(6,1,1),    # barplot colour vector. Assignment in order of explanatory variables. Default: 6*white i.e. blank bars (has border)
            linesfiles = FALSE,           # save individual line plots' data as csv's?
            savegbm = TRUE,               # save the gbm objects externally? Can reopen later with (e.g.) load("Bin_Best_Model")
            varint = TRUE,                # calculate variable interactions? Default:TRUE, set FALSE if code fails with "contrasts can be applied only to factors with 2 or more levels"
            map = TRUE,                   # save abundance map png files?
            RSB = TRUE,                   # run representativeness surface builder?
            BnW = TRUE,                   # repeat map (& RSB if TRUE) in black & white for print journals
            alerts = TRUE,                # play sounds to mark progress steps
            ...)                          # additional parameters, esp for gbm.map (byx, byy, mapmain, heatcolours, colournumber, shape, mapback, landcol, legendtitle, lejback, legendloc, grdfun, zero, quantile)
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
# and my own functions gbm.map, gbm.rsb, gbm.valuemap

####1. Check packages, start loop####
if (!require(gbm)) {stop("you need to install the gbm package to run this function")}
if (!require(dismo)) {stop("you need to install the dismo package to run this function")}
if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
if (!require(labeling)) {stop("you need to install the labeling package to run this function")}
if(map==TRUE) if(!require(mapplots)) {stop("you need to install the mapplots package to run this function")}
if(RSB==TRUE) if(!exists("gbm.rsb")) {stop("you need to install the gbm.rsb function to run this function")}
if(RSB==TRUE) if(!exists("gbm.map")) {stop("you need to install the gbm.map function to run this function")}
if(!is.null(grids)) if(!exists("gbm.predict.grids")) {stop("you need to install the gbm.predict.grids function from gbm.utils.R to run this function")}
if(!exists("roc")) {stop("you need to install the roc function from gbm.utils.R to run this function")}
if(!exists("calibration")) {stop("you need to install the calibration function from gbm.utils.R to run this function")}
require(gbm)
require(dismo)
require(beepr)
require(labeling)

if(alerts) options(error = function() {beep(9)})  # give warning noise if it fails
    
expvarnames<-names(samples[expvar]) # list of explanatory variable names
expvarcols<-cbind(cols[1:length(expvarnames)],expvarnames) # assign explanatory variables to colours
if(!exists("tc")) tc <- c(2,length(expvar)) # if tc not set by user, default to 2,length(expvar)

for(i in resvar){
m=1 # jkl combo loop counter to allow best bin/gaus BRT choice
n=1   # Print counter for all loops of BRT combos (i.e. counter for l)
if(!is.null(grids)) if(!all(expvarnames %in% names(grids))) {stop("Not all expvar column names found as column names in grids")}

####2. ZI check & log####
# if user has asked code to check for ZI, check it & set new ZI status
if(ZI=="CHECK") if(sum(samples[,resvar]==0,na.rm=TRUE)/length(samples[,resvar])>=0.5) ZI=TRUE else ZI=FALSE

# ensure resvar has zeroes (expects mix of successful & unsuccessful samples)
if(min(samples[i])>0) print("No zeroes in response variable. Method expects unsuccessful, as well as successful, samples")

# create binary (0/1) response variable, for bernoulli BRTs
samples$brv <- ifelse(samples[i] > 0, 1, 0)
brvcol <- which(colnames(samples)=="brv") # brv column number for BRT

# create logged response variable, for gaussian BRTs when data is zero-inflated (otherwise just use resvar directly)
logem <- log1p(samples[,i])
dont  <- samples[,i]
if(ZI) {samples$grv <- logem} else {samples$grv <- dont}
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

# Begin loops
for(j in tc){   # list permutations of tree complexity allowed
 for(k in lr){   # list permutations of learning rate allowed
  for(l in bf){   # list permutations of bag fraction allowed

####4. Binomial BRT####
if(ZI) {  # don't do if ZI=FALSE
(assign(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""),gbm.step(data=samples,
    gbm.x = expvar, gbm.y = brvcol, family = "bernoulli", tree.complexity = j, learning.rate = k, bag.fraction = l)))

####5. Select best Bin model####
# Makes an object with BRT training data correlation score of the 1st variable combo in the loop, & the name of that BRT combo,
# per species (m is reset at the end of each species loop so m=1 is the first loop per species)
# Then if subsequent combos have higher values (more correlation), it replaces the value, and name of combo it came from

# create blanks for best results
Bin_Best_Score <- 0
Bin_Best_Model <- 0

if(m==1)
  {Bin_Best_Score <- get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
  Bin_Best_Model <- paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
}  else if(get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]>Bin_Best_Score)
      {Bin_Best_Score <- get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
      Bin_Best_Model <- paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep="")}

# progress printer, right aligned
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXX  Completed BRT ",n," of ",2*length(i)*length(tc)*length(lr)*length(bf),"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)
n <- n+1   # Add to print counter
} # close ZI option

####6. Gaussian BRT####
assign(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""),gbm.step(data=samples,
    gbm.x = expvar, gbm.y = grvcol, family = "gaussian", tree.complexity = j, learning.rate = k, bag.fraction = l))

####7. Select best Gaus model####
# create blanks for best results
Gaus_Best_Score <- 0
Gaus_Best_Model <- 0

if(m==1)
  {Gaus_Best_Score <- get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
  Gaus_Best_Model <- paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
} else if(get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]>Gaus_Best_Score)
      {Gaus_Best_Score <- get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]]
      Gaus_Best_Model <- paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep="")}

####8. Add BRT stats to report####
Report[1:3,((m*5)-1)] <- c(paste("tree complexity: ",j,sep=""), paste("learning rate: ",k,sep=""), paste("bag fraction: ",l,sep=""))
# don't do if ZI=FALSE
if(ZI) {Report[1:6,(m*5)] <- c(paste("trees: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
                       paste("Training Data Correlation: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]],sep=""),
                       paste("CV Mean Deviance: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.mean,sep=""),
                       paste("CV Deviance SE: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.se,sep=""),
                       paste("CV Mean Correlation: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.mean,sep=""),
                       paste("CV Correlation SE: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.se,sep=""))
Report[1,(m*5)+1] <- paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep="")}
Report[1:6,((m*5)+2)] <- c(paste("trees: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
                            paste("Training Data Correlation: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$self.statistics$correlation[[1]],sep=""),
                            paste("CV Mean Deviance: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.mean,sep=""),
                            paste("CV Deviance SE: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$deviance.se,sep=""),
                            paste("CV Mean Correlation: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.mean,sep=""),
                            paste("CV Correlation SE: ",get(paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$cv.statistics$correlation.se,sep=""))
Report[1,(m*5)+3] <- paste("Gaus_BRT",".tc",j,".lr",k*100,".bf",l,sep="")
# different colnames depending on ZI state
if(ZI) {colnames(Report)[((m*5)-1):((m*5)+3)] <- c(paste("Parameter Combo ",m,sep=""),
                                           paste("Binary BRT ",m," stats",sep=""),
                                           paste("Binary BRT ",m," name",sep=""),
                                           paste("Gaussian BRT ",m," stats",sep=""),
                                           paste("Gaussian BRT ",m," name",sep=""))
} else {
        colnames(Report)[((m*5)-1):((m*5)+3)] <- c(paste("Parameter Combo ",m,sep=""),
                                             paste("No Binary BRT"),
                                             paste("No Binary BRT"),
                                             paste("Gaussian BRT ",m," stats",sep=""),
                                             paste("Gaussian BRT ",m," name",sep=""))
}

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Completed BRT ",n," of ",2*length(i)*length(tc)*length(lr)*length(bf),"    XXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)
n <- n+1   # Add to print counter: 2 per loop, 1 bin 1 gaus BRT
m <- m+1   # Add to loop counter: 1 per loop, used for bin/gaus_best model selection
}}}        # close loops, producing all BRT/GBM objects then continuing through model selection

####9. Test simplification benefit, do so if better####
samples<<-samples # globally assign samples: bad practice but fixes an odd problem whereby the code runs if done manually
# but crashes when in a huge looping function, saying it's unable to access samples.
# only do if ZI=TRUE
if(ZI) {Bin_Best_Simp_Check <- gbm.simplify(get(Bin_Best_Model)) # run the simplification tester on the best model

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
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Bin model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)
} # close ZI

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
if(alerts) beep(2)
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Gaus model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))

####10. Select final best models####
if(ZI) {  # don't do if ZI=FALSE
# if Bin_Best has a simplified model:
if(min(Bin_Best_Simp_Check$deviance.summary$mean) < 0)
# & if the simplified model has better correlation than Bin_Best itself
  if(Bin_Best_Simp$self.statistics$correlation > Bin_Best_Score[1])
# then replace Bin_Best score/model values with those from the simplified model
  {Bin_Best_Score <- Bin_Best_Simp$self.statistics$correlation
   Bin_Best_Model <- "Bin_Best_Simp"}
# globally assign final model for external testing later
Bin_Best_Model<<-get(Bin_Best_Model)
} # close ZI

# Same for Gaus:
if(min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0)
  if(Gaus_Best_Simp$self.statistics$correlation > Gaus_Best_Score[1])
    {Gaus_Best_Score <- Gaus_Best_Simp$self.statistics$correlation
     Gaus_Best_Model <- "Gaus_Best_Simp"}
# globally assign final model for external testing later
Gaus_Best_Model<<-get(Gaus_Best_Model)

# progress printer, right aligned for visibility
if(alerts) beep(2)
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Best models selected    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))

####11. Line plots####
dir.create(names(samples[i])) # create resvar-named directory for outputs

# All plots on one image for Bin & Gaus
if(ZI) {  # don't do if ZI=FALSE
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
} # close ZI

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
if(ZI) {  # don't do if ZI=FALSE
for (o in 1:length(get(Bin_Best_Model)$contributions$var)){
  png(filename = paste("./",names(samples[i]),"/Bin_Best_line_",as.character(get(Bin_Best_Model)$contributions$var[o]),".png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(1.35,3.4,0.4,0.5), fig=c(0,1,0,1), las=1, lwd=8, bty="n", mgp=c(2,0.5,0), xpd=NA)  #mar=c(1.35,2.4,0.4,0.5)
    # bg=expvarcols[match(get(Bin_Best_Model)$contributions$var[o],expvarcols[,2]),1]) #changed margin to hide label #XPD YPD ALLOWS AXES TO EXTEND FURTHER TO ENCOMPASS ALL DATA? #colour removed
plotgrid<-plot.gbm(get(Bin_Best_Model),match(get(Bin_Best_Model)$contributions$var[o], get(Bin_Best_Model)$gbm.call$predictor.names),lwd=8,return.grid=TRUE)
if(linesfiles) write.csv(plotgrid, row.names=FALSE, na="", file = paste("./",names(samples[i]),"/Bin_Best_line_",as.character(get(Bin_Best_Model)$contributions$var[o]),".csv",sep=""))
xx <- labeling::extended(min(plotgrid[1]), max(plotgrid[1]),7, only.loose=TRUE) # sets range & ticks
yy <- labeling::extended(min(plotgrid[2]), max(plotgrid[2]),7, only.loose=TRUE) # sets range & ticks
plot(range(xx),range(yy), t="n", xaxt="n", yaxt="n", bty="n", ylab=NA)
lines(plotgrid,type="l")
axis(1, lwd.ticks=8, lwd=8, at=xx) # is providing only the thick line & downticks
axis(2, lwd.ticks=8, lwd=8, at=yy)
#rug(quantile(samples[as.character(get(Bin_Best_Model)$contributions$var[o])], probs=seq(0,1,0.01), na.rm=TRUE), side=1, lwd=5, ticksize=0.03) #n of ticks probs seq arg 3: 0.1, 0.05, 0.01
rug(samples[as.character(get(Bin_Best_Model)$contributions$var[o])][,1], side=1, lwd=5, ticksize=0.03) # all points rug
dev.off() }
} # close ZI option
  
for (p in 1:length(get(Gaus_Best_Model)$contributions$var)){
  png(filename = paste("./",names(samples[i]),"/Gaus_Best_line_",as.character(get(Gaus_Best_Model)$contributions$var[p]),".png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(1.35,3.4,0.4,0.5), fig=c(0,1,0,1), las=1, lwd=8, bty="n", mgp=c(2,0.5,0), xpd=NA)
    # bg=expvarcols[match(get(Gaus_Best_Model)$contributions$var[p],expvarcols[,2]),1]) #changed margin to hide label #XPD YPD ALLOWS AXES TO EXTEND FURTHER TO ENCOMPASS ALL DATA? #colour removed
plotgrid<-plot.gbm(get(Gaus_Best_Model),match(get(Gaus_Best_Model)$contributions$var[p], get(Gaus_Best_Model)$gbm.call$predictor.names),lwd=8,return.grid=TRUE)
if(linesfiles) write.csv(plotgrid, row.names=FALSE, na="", file = paste("./",names(samples[i]),"/Gaus_Best_line_",as.character(get(Gaus_Best_Model)$contributions$var[p]),".csv",sep=""))
xx <- labeling::extended(min(plotgrid[1]), max(plotgrid[1]),7, only.loose=TRUE) # sets range & ticks
yy <- labeling::extended(min(plotgrid[2]), max(plotgrid[2]),7, only.loose=TRUE) # sets range & ticks
plot(range(xx),range(yy), t="n", xaxt="n", yaxt="n", bty="n", ylab=NA)
lines(plotgrid,type="l")
axis(1, lwd.ticks=8, lwd=8, at=xx) # is providing only the thick line & downticks
axis(2, lwd.ticks=8, lwd=8, at=yy)
#rug(quantile(samples[as.character(get(Gaus_Best_Model)$contributions$var[p])], probs=seq(0,1,0.01), na.rm=TRUE), side=1, lwd=5, ticksize=0.03) #n of ticks probs seq arg 3: 0.1, 0.05, 0.01
rug(samples[as.character(get(Gaus_Best_Model)$contributions$var[p])][,1], side=1, lwd=5, ticksize=0.03) # all points rug
dev.off() }

# progress printer, right aligned for visibility
if(alerts) beep(2)
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Line plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))

####12. Dot plots####
if(ZI) {  # don't do if ZI=FALSE
png(filename = paste("./",names(samples[i]),"/Bin_Best_dot.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
gbm.plot.fits(get(Bin_Best_Model))
dev.off()
} # close ZI

png(filename = paste("./",names(samples[i]),"/Gaus_Best_dot.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
gbm.plot.fits(get(Gaus_Best_Model))
dev.off()

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Dot plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

####13. 3D plot TODO####
# gbm.perspec(Bin_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)
# gbm.perspec(Gaus_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)

####14. Bar plots of variable influence####
# create tables
if(ZI) {  # don't do if ZI=FALSE
Bin_Bars <- summary(get(Bin_Best_Model),
        cBars=length(get(Bin_Best_Model)$var.names),
        n.trees=get(Bin_Best_Model)$n.trees,
        plotit=FALSE, order=TRUE, normalize=TRUE, las=1, main=NULL)
write.csv(Bin_Bars, file=paste("./",names(samples[i]),"/Binary BRT Variable contributions.csv",sep=""), row.names=FALSE)
} # close ZI

Gaus_Bars <- summary(get(Gaus_Best_Model),
        cBars=length(get(Gaus_Best_Model)$var.names),
        n.trees=get(Gaus_Best_Model)$n.trees,
        plotit=FALSE, order=TRUE, normalize=TRUE, las=1, main=NULL)
write.csv(Gaus_Bars, file=paste("./",names(samples[i]),"/Gaussian BRT Variable contributions.csv",sep=""), row.names=FALSE)

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Bar plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

# produce graphics
if(ZI) {  # don't do if ZI=FALSE
png(filename = paste("./",names(samples[i]),"/Bin_Bars.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "",
    type = "cairo-png")
par(mar=c(2.5,0.3,0,0.5), fig=c(0,1,0,1), cex.lab=0.8,mgp=c(1.5,0.5,0), cex=1.3)
midpoints<-barplot(rev(Bin_Bars[,2]), cex.lab=1.2, las=1, horiz=TRUE, cex.names=0.8, xlab="Influence %", col=rev(expvarcols[match(Bin_Bars[,1],expvarcols[,2]),1]), xlim=c(0,2.5+ceiling(max(Bin_Bars[,2]))))
text(0.1, midpoints,labels=rev(Bin_Bars[,1]),adj=0,cex=1.5)
axis(side = 1, lwd = 6, outer=TRUE, xpd=NA)
dev.off()
} # close ZI

png(filename = paste("./",names(samples[i]),"/Gaus_Bars.png",sep=""),
    width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "",
    type = "cairo-png")
par(mar=c(2.5,0.3,0,0.5), fig=c(0,1,0,1), cex.lab=0.8,mgp=c(1.5,0.5,0), cex=1.3)
midpoints<-barplot(rev(Gaus_Bars[,2]), cex.lab=1.2, las=1, horiz=TRUE, cex.names=0.8, xlab="Influence %", col=rev(expvarcols[match(Gaus_Bars[,1],expvarcols[,2]),1]), xlim=c(0,2.5+ceiling(max(Gaus_Bars[,2]))))
text(0.1, midpoints,labels=rev(Gaus_Bars[,1]),adj=0,cex=1.5)
axis(side = 1, lwd = 6, outer=TRUE, xpd=NA)
dev.off()

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Bar plots plotted      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

####15. Variable interactions####
# only do them if varint=TRUE, the default. Only do bin if ZI=TRUE
if(ZI) if(varint) find.int_Bin <- gbm.interactions(get(Bin_Best_Model))
if(varint) find.int_Gaus <- gbm.interactions(get(Gaus_Best_Model))
# progress printer, right aligned for visibility
if(varint) {print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXX    Variable interactions calculated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))}
if(alerts) beep(2)

#avoid sections 16-20 if not predicting to grids
if(!is.null(grids)) {
####16. Binomial predictions####
if(ZI) {  # don't do if ZI=FALSE
gbm.predict.grids(get(Bin_Best_Model), grids, want.grids = F, sp.name = "Bin_Preds")
grids$Bin_Preds <- Bin_Preds
} # close ZI
  
# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Binomial predictions calculated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

####17. Gaussian predictions####
gbm.predict.grids(get(Gaus_Best_Model), grids, want.grids = F, sp.name = "Gaus_Preds")
if(ZI) {grids$Gaus_Preds <- Gaus_Preds

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Gaussian predictions calculated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

####18. Backtransform logged Gaus to unlogged####
grids$Gaus_Preds_Unlog <- exp(Gaus_Preds + 1/2 * sd(get(Gaus_Best_Model)$residuals,na.rm=FALSE)^2)

####19. BIN*positive abundance = final abundance####
grids$PredAbund <- grids$Gaus_Preds_Unlog * grids$Bin_Preds} else {grids$PredAbund <- Gaus_Preds} #if ZI=TRUE, unlog gaus & multiply by bin. Else just use gaus preds.
predabund <- which(colnames(grids)=="PredAbund") # predicted abundance column number for writecsv

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Final abundance calculated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

####20. Final saves####
# CSV of Predicted values at each site inc predictor variables' values.
write.csv(grids,row.names=FALSE, file = paste("./",names(samples[i]),"/Abundance_Preds_All.csv",sep=""))
# CSV of Predicted values at each site without predictor variables' values.
write.csv(grids[c(gridslat,gridslon,predabund)], row.names=FALSE, file = paste("./",names(samples[i]),"/Abundance_Preds_only.csv",sep=""))
if(savegbm) {save(Gaus_Best_Model,file = paste("./",names(samples[i]),"/Gaus_Best_Model",sep=""))
  if(ZI) {save(Bin_Best_Model,file = paste("./",names(samples[i]),"/Bin_Best_Model",sep=""))}} #only save bin if ZI=TRUE
} #close grids option from above section 16

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Output CSVs written   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

####21. Finalise & Write Report####
# only do final variable interaction lines if varint=TRUE
if(ZI) Report[1:2,(reportcolno-13)] <- c(paste("Model combo: ",Bin_Best_Model,sep=""),paste("Model CV score: ",Bin_Best_Score,sep=""))
Report[1:2,(reportcolno-12)] <- c(paste("Model combo: ",Gaus_Best_Model,sep=""),paste("Model CV score: ",Gaus_Best_Score,sep=""))
if(ZI) Report[1:dim(subset(Bin_Best_Simp_Check$final.drops,order>0))[1],(reportcolno-11)] <- as.character(subset(Bin_Best_Simp_Check$final.drops,order>0)$preds)
if(ZI) Report[1:(length(Bin_Best_Simp_Check$final.drops$preds)-dim(subset(Bin_Best_Simp_Check$final.drops,order>0))[1]),(reportcolno-10)] <-
  as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order>0))[1])+1):length(Bin_Best_Simp_Check$final.drops$preds)])
Report[1:dim(subset(Gaus_Best_Simp_Check$final.drops,order>0))[1],(reportcolno-9)] <- as.character(subset(Gaus_Best_Simp_Check$final.drops,order>0)$preds)
Report[1:(length(Gaus_Best_Simp_Check$final.drops$preds)-dim(subset(Gaus_Best_Simp_Check$final.drops,order>0))[1]),(reportcolno-8)] <-
  as.character(Gaus_Best_Simp_Check$final.drops$preds[((dim(subset(Gaus_Best_Simp_Check$final.drops,order>0))[1])+1):length(Gaus_Best_Simp_Check$final.drops$preds)])
if(ZI) Report[1:6,(reportcolno-7)] <-   c(paste("trees: ",get(paste("Bin_BRT",".tc",j,".lr",k*100,".bf",l,sep=""))$n.trees,sep=""),
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
if(ZI) Report[1:(length(Bin_Bars[,1])),(reportcolno-5)] <- as.character(Bin_Bars$var)
if(ZI) Report[1:(length(Bin_Bars[,2])),(reportcolno-4)] <- as.character(Bin_Bars$rel.inf)
Report[1:(length(Gaus_Bars[,1])),(reportcolno-3)] <- as.character(Gaus_Bars$var)
Report[1:(length(Gaus_Bars[,2])),(reportcolno-2)] <- as.character(Gaus_Bars$rel.inf)
if(ZI) if(varint) Report[1:2,(reportcolno-1)] <- c(paste(find.int_Bin$rank.list$var1.names[1]," and ",find.int_Bin$rank.list$var2.names[1],". Size: ",find.int_Bin$rank.list$int.size[1],sep=""),
                                paste(find.int_Bin$rank.list$var1.names[2]," and ",find.int_Bin$rank.list$var2.names[2],". Size: ",find.int_Bin$rank.list$int.size[2],sep=""))
if(varint) Report[1:2,(reportcolno)] <- c(paste(find.int_Gaus$rank.list$var1.names[1]," and ",find.int_Gaus$rank.list$var2.names[1],". Size: ",find.int_Gaus$rank.list$int.size[1],sep=""),
                              paste(find.int_Gaus$rank.list$var1.names[2]," and ",find.int_Gaus$rank.list$var2.names[2],". Size: ",find.int_Gaus$rank.list$int.size[2],sep=""))
write.csv(Report, row.names=FALSE, na="", file = paste("./",names(samples[i]),"/Report.csv",sep=""))

# progress printer, right aligned for visibility
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Report CSV written    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
if(alerts) beep(2)

#avoid sections 22&23 if not predicting to grids
if(!is.null(grids)) {
  
####22. Representativeness surface builder####
# does NOT plot the surface, just builds it. If built, it will be plotted by map maker.
  if(RSB==TRUE){rsbdf_bin <- gbm.rsb(samples,grids,expvarnames,gridslat,gridslon,rsbres=i)
  pos_samples <- subset(samples, brv >0)
  rsbdf_gaus <- gbm.rsb(pos_samples,grids,expvarnames,gridslat,gridslon,rsbres=i)}
  
####23. Map maker####
if(!exists("mainlegendtitle")) mainlegendtitle = "CPUE" # create mainlegendtitle default if absent. Else will cause error if unset.
if(map==TRUE) {
  # generate output image & set parameters
  png(filename = paste("./",names(samples[i]),"/PredAbundMap_",names(samples[i]),".png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
  # run gbm.map function with generated parameters
  gbm.map(x = grids[,gridslon],
          y = grids[,gridslat],
          z = grids[,predabund],
          species = names(samples[i]),
          legendtitle = mainlegendtitle,
          ...)  # allows gbm.auto's optional terms to be passed to subfunctions:
  # byx, byy, mapmain, heatcol, shape, mapback, landcol, lejback, legendloc, grdfun, zero, quantile, heatcolours, colournumber
  dev.off()
  
  # progress printer, right aligned for visibility
  print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour map generated     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
  if(alerts) beep(2)
  
  # if BnW=TRUE, run again in black & white for journal submission
  if(BnW){
  png(filename = paste("./",names(samples[i]),"/PredAbundMap_BnW_",names(samples[i]),".png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
  gbm.map(x = grids[,gridslon],
          y = grids[,gridslat],
          z = grids[,predabund],
          species = names(samples[i]),
          legendtitle = mainlegendtitle,
          landcol = grey.colors(1, start=0.8, end=0.8), #light grey. 0=black 1=white
          mapback = "white",
          heatcolours = grey.colors(8, start=1, end=0), #default 8 greys; setting heatcolours & colournumber overrides this
          ...) # allows gbm.auto's optional terms to be passed to subfunctions
  dev.off()} # close & save plotting device & close BnW optional
  
  # progress printer, right aligned for visibility
  print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Black & white map generated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
  if(alerts) beep(2)
  
  # if RSB called, plot that surface separately
  if(RSB==TRUE){
    png(filename = paste("./",names(samples[i]),"/RSB_Map_Bin_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf_bin[,"Unrepresentativeness"],
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            legendtitle = "UnRep 0-1",
            ...)
    dev.off()
    
    # progress printer, right aligned for visibility
    print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour RSB binary map generated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
    if(alerts) beep(2)
    
    png(filename = paste("./",names(samples[i]),"/RSB_Map_Gaus_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf_gaus[,"Unrepresentativeness"],
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            legendtitle = "UnRep 0-1",
            ...)
    dev.off()
    
    # progress printer, right aligned for visibility
    print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour RSB Gaussian map generated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
    if(alerts) beep(2)
    
    png(filename = paste("./",names(samples[i]),"/RSB_Map_Both_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf_bin[,"Unrepresentativeness"]+rsbdf_gaus[,"Unrepresentativeness"],
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            legendtitle = "UnRep 0-2",
            ...)
    dev.off()
    
    # progress printer, right aligned for visibility
    print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXX    Colour RSB combination map generated    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
    if(alerts) beep(2)
    
    # if BnW=TRUE, do again for b&w
    if(BnW){
    png(filename = paste("./",names(samples[i]),"/RSB_Map_BnW_Bin_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf_bin[,"Unrepresentativeness"],
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            heatcolours = grey.colors(8, start=1, end=0), #default 8 greys; setting heatcolours & colournumber overrides this
            landcol = grey.colors(1, start=0.8, end=0.8), #light grey. 0=black 1=white
            legendtitle = "UnRep 0-1",
            ...)
    dev.off()
    
    # progress printer, right aligned for visibility
    print(paste("XXXXXXXXXXXXXXXXXXXXXXXXX    Black & white RSB binary map generated    XXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
    if(alerts) beep(2)
    
    png(filename = paste("./",names(samples[i]),"/RSB_Map_BnW_Gaus_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf_gaus[,"Unrepresentativeness"],
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            heatcolours = grey.colors(8, start=1, end=0), #default 8 greys; setting heatcolours & colournumber overrides this
            landcol = grey.colors(1, start=0.8, end=0.8), #light grey. 0=black 1=white
            legendtitle = "UnRep 0-1",
            ...)
    dev.off()
    
    # progress printer, right aligned for visibility
    print(paste("XXXXXXXXXXXXXXXXXXXXXXX    Black & white RSB Gaussian map generated    XXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
    if(alerts) beep(2)
    
    png(filename = paste("./",names(samples[i]),"/RSB_Map_BnW_Both_",names(samples[i]),".png",sep=""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
    gbm.map(x = grids[,gridslon], # add representativeness alpha surface
            y = grids[,gridslat],
            z = rsbdf_bin[,"Unrepresentativeness"]+rsbdf_gaus[,"Unrepresentativeness"],
            mapmain = "Unrepresentativeness: ",
            species = names(samples[i]),
            heatcolours = grey.colors(8, start=1, end=0), #default 8 greys; setting heatcolours & colournumber overrides this
            landcol = grey.colors(1, start=0.8, end=0.8), #light grey. 0=black 1=white
            legendtitle = "UnRep 0-2",
            ...)
    dev.off()}

    # progress printer, right aligned for visibility
    print(paste("XXXXXXXXXXXXXXXXXXXXXX    Black & white RSB Combination map generated    XXXXXXXXXXXXXXXXXXXXXXXXXX",sep=""))
    if(alerts) beep(2)
    
    } # close RSB mapper
   } # close Map Maker
  } #close grids option from above section 22
 } # close response variable (resvar) loop
if(alerts) beep(8) # notify the user with a noise, since this process can take a long time.
} # close the function
####END####