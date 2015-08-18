####Line plots begin @ 0####
# they shouldn't: salinity & temp (& current speed?) don't lend themselves to this.
# check whether data ranges include zeroes in place of NAs?
# see dotplots: they do have zeroes. Aren't these (x) just the raw values though? Fitted (y) are generated
summary(mysamples[,"Salinity"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00   11.51   33.99   35.08 
hist(mysamples[,"Salinity"])
sum(samples[,"Salinity"]==0,na.rm=TRUE)
#4452

####mapmaker defaults @ toplevel####
# see gbm.auto section23, landcol mapback legendloc legendtitle??
# mainmap, heatcol, shape, landcol, legendcol, bg, mapback, grdfun (line17)
# can the user edit these, ACTUALLY? Try changing all of them and see which work

####ZI TEST####
# See paper by Tu in Qiqqa, ZI data.
if(sum(samples[,resvar]==0,na.rm=TRUE)/length(samples[,resvar])>=0.5) ZI=TRUE else ZI=FALSE
if(ZI=TRUE) #run bin analysis
  # "else" do nothing i.e. only do gaus. Also have to change final multiplication and folder creation etc]

  ZI="CHECK"
if(ZI=="CHECK") if(sum(samples[,resvar]==0,na.rm=TRUE)/length(samples[,resvar])>=0.5) ZI=TRUE else ZI=FALSE
  

# logs of resvar for ZI data & reverse lognormalise & bias correct later (line 117 gaus BRT gbm.y & sections 22 & 23)
# asked here: https://stats.stackexchange.com/questions/112292/r-are-data-zero-inflated

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

# Lines to edit:
# 91-119
# 115-118 move/edit progress printer?
# 137-156 sort report
# 165-186
# 203-212
# 223-224 delete?
# 227-237
# 252-266
# 286-294
# 313-317
# 326-333
# 345
# 349-350
# 357: unlog
# 360-361: multiply
# 369
# 373-401 sort report


####Bag fraction optimiser####
# Trial & error iterative approach to determine the optimal bag fraction for a data set? How? Stop @ whole percentages?
# Try OPTIM function & see http://r.789695.n4.nabble.com/Optimization-in-R-similar-to-MS-Excel-Solver-td4646759.html
# Possibly have an option to start with this in the R function.

####OPTIMISE PARAMETERS####
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

####MULTICORE PROCESSING####
#already supposedly incorporated; doesnt work

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
#
# See:
# C:\Users\Simon\Dropbox\Galway\Analysis\R\Coilin R code\rmpi_example.R

####Processing time estimate####
#  option to have R poll the computer and, when you say go, popup a box saying
# "Based on the parameters you've selected to try, this will run X models on a Y-item-sized dataset (Y = variables * count)
# which may take about Z minutes based on your processor. You have a multicore processor so it will take Z/#processors (time)
# if you have multicore processing enabled - see here LINK" OK/Cancel
# Edit the BRT progrss counter to account for multiple resvars: currently does e.g. n/8 then loops back to 1. Not useful.
# Also could print the current resvar name. Maybe do a running time thing? "This is BRT N of X, Y% complete, took time Z, total
# expected time AA, time remaining AB

####Clean workspace####
#dump fails: TRY rm() rather than dump() - dont need this for function
# rm(list = ls()) removes everything
# running function the dumping everyithng means R is still using 1.7gb of RAM!
# rm all model-generated files after each run? List what they are here.

####3D PLOT####
# what to do? Bother with it?

####Allow deep subfunction calls####
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