####TO DO####
# Diya:
# Error in FUN(X[[i]], ...) : only defined on a data frame with all numeric variables
# 3 env vars categorical. Changed to numeric.
# If data need to be numeric I need to reflect that in the documentation AOR
# create a converter.

# separate bin & gaus parameters?
# bin/gaus jkl loops 114:195
# currently:
# bin brt, best bin model, progress printer
# gaus brt, best gaus model,
# add stats to report,
# progres printer
# counters
#
# could be:
# LOOP1: bin brt, best bin model, add bin stats to report, counters, progress printer
# LOOP2: gaus brt, best gaus model, add gaus stats to report, counters, progress printer
# fix report L434
#
# but then what? have lrbin & lrgaus, tcbin & tcgaus, bfbin & bfgaus?
# lr is always lrbin but lrgaus is optional, lr used unless lrgaus supplied?
if (isnull(lrgaus)) lr <- lrgaus
for (k in lrgaus){}

lr = c(0.01,0.005)
lr = list(bin = c(0.01, 0.05), gaus = c(0.005, 0.0001))
class(lr) # numeric / list
lr
lr$bin
length(lr$bin)
as.vector(lr)
for (k in lr) {
  }
gbm.step(data = samples, gbm.x = expvar, gbm.y = brvcol, family = "bernoulli", tree.complexity = j, learning.rate = k, bag.fraction = l, ...)


# Process & map bin only? nightmarish. requires reordering of all bin/gaus/bin/gaus to bin/bin/bin
 # gaus = TRUE
# if (gaus)
#Ls: 94:99,

# improve how functions are checked/loaded.
# see this from gbm:
if (!requireNamespace("gbm")) {
  stop("you need to install the gbm package to run this function")

####gbm.auto heatcols for BnW gbm.map####
# manually setting heatcolours in the gbm.auto(call) passes heatcolours to gbm.auto for gbm.map
# but it'll override the grey colours for the B&W maps
# See gbm.auto L578

####Multicore Processing####
# already supposedly incorporated; doesnt work
# Re-investigate multicore R
# https://stackoverflow.com/questions/4775098/r-with-a-multi-core-processor
# http://www.google.com/url?q=http%3A%2F%2Fcran.r-project.org%2Fweb%2Fviews%2FHighPerformanceComputing.html&sa=D&sntz=1&usg=AFQjCNEshsLTzWAF9g4cUGzvn5zIfIywsA
# https://groups.google.com/forum/#!topic/davis-rug/VLIXz5i3vZI
# n.cores in gbm
# The number of CPU cores to use. The cross-validation loop will attempt to send different CV folds off to different cores.
# If n.cores is not specified by the user, it is guessed using the detectCores function in the parallel package.
# Note that the documentation for detectCores makes clear that it is not failsafe and could return a spurious number of available cores.
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
# See:
# C:\Users\Simon\Dropbox\Galway\Analysis\R\Coilin R code\rmpi_example.R
# mpi doesn't work on my laptop? due to daily rstudio build, or broken laptop, or neither?
gbm.auto(expvar = c(4:9,11), resvar = 12, grids = mygrids, lr = c(0.02), ZI = TRUE, map = FALSE, RSB = FALSE, tc = 2, varint = FALSE, savegbm = FALSE)
# posted to stack exchange, no answer

####OPTIMISE PARAMETERS####
# Bag fraction optimiser
# Trial & error iterative approach to determine the optimal bag fraction for a data set? How? Stop @ whole percentages?
# Try OPTIM function & see http://r.789695.n4.nabble.com/Optimization-in-R-similar-to-MS-Excel-Solver-td4646759.html
# Possibly have an option to start with this in the R function.
#
# Trial & error iterative approach to determine the optimal bag fraction & lr concurrently
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


####Processing time estimate####
#  option to have R poll the computer and, when you say go, popup a box saying
# "Based on the parameters you've selected to try, this will run X models on a Y-item-sized dataset (Y = variables * count)
# which may take about Z minutes based on your processor. You have a multicore processor so it will take Z/#processors (time)
# if you have multicore processing enabled - see here LINK" OK/Cancel
# Edit the BRT progrss counter to account for multiple resvars: currently does e.g. n/8 then loops back to 1. Not useful.
# Also could print the current resvar name. Maybe do a running time thing? "This is BRT N of X, Y% complete, took time Z, total
# expected time AA, time remaining AB
# see proc.time()
# see http://www.ats.ucla.edu/stat/r/faq/timing_code.htm
ptm <- proc.time() # Start the clock!
proc.time() - ptm # Stop the clock
ptm$elapsed # is the time taken in seconds
# maybe add to gbm.auto's report. Plus:
   #  size (total cells, dim) of 'active' database,
   #  whether maps generated,
   #  RSB,
   #  BnW,
   #  savegbm,
   #  varint,
   #  sizes of tc,lr,bf, (smaller numbers take longer)
   #  length of tc,lr,bf,
   #  CPU speed,
   #  number of CPU cores (assuming multithreading sorted!),
   #  RAM size
   #
   #  Assumedly something like time=
   #  dbasesize*sizes(tc,lr,bf)*length(tc*lr*bf) #sizes would need a reciprocal, smaller takes longer
   #  +savegbm
   #  +varint
   #   all*(2 if ZI=true?)
   #  +map*prod(lengths(tc,lr,bf))
   #  +RSB*prod(lengths(tc,lr,bf))
   #  +BnW*prod(lengths(tc,lr,bf))
   # All divided by CPU*cores*RAM?
# Run gbm.auto loads of times, changing one parameter through its range each time, and calculate relationship of
# model parameters to processing time on MY laptop. Then can logically work out how long it should take on another
# machine? Try it on my home one. RAM doesn't matter unless it's limiting? Limitingness is a function of dbasesize
# and parameters, especially if I'm not cleaning the workspace within the function run?

####3D PLOT####
# what to do? Bother with it?

# include eliths BRT_ALL within this code?


####DONE####
# remove log1p since grv_yes subset now has no zeroes L104 & 8

# incorporate gbm.basemap as a gbm.auto default

# gbm.basemap: calculate default res based on size of bounds?

# if nrows (either original for bin or grv_yes for gaus) are <= 42 it'll crash @ 0.5 bf.
# make a function which solves for your row number?? gbm.bfcheck

# Allow user to call extra parameters with '...' which will be parsed to gbm.step
# Ls 124 & 149.
#           offset = NULL,                            # allows an offset to be specified
#           fold.vector = NULL,                       # allows a fold vector to be read in for CV with offsets,
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

# allow different image device for Mac OSX

# user setting tc overrides auto-calc default?

# gbm.rsb plots: log+1 & have scale go to 2 (for 'both' else 1)

# clean workspace

# gbm.basemap

# unrep plots: Change title back to unrep

####grids to be optional####
# if grids=FALSE, ignore:
# lines 10,11, (don't need to change, can be set but unused objects) 40, 54,
# sections 16-20, 22-23
#WORKED!

####variable interactions switch####
# option to switch off variable interactions if they fail
# topline switch defaults to varint=TRUE
# false wraps if(varint) {existing lines}
# "contrasts can be applied only to factors with 2 or more levels"
# worked!

####mapmaker defaults @ toplevel####
# see gbm.auto section23, landcol mapback legendloc legendtitle??
# mapmain, heatcol, shape, landcol, legendcol, lejback, mapback, grdfun (line17)
# can the user edit these, ACTUALLY? Try changing all of them and see which work
# works!

####make b&w an option?####
# default to true
# WORKS!

####win linux mac####
# does directory creation work the same? aren't the slashes different ways around?
# dir.create, l237
# probably works fine!

####gbm.map terms####
# x/y/z are the same as grids[,gridslon]/grids[,gridslat]/grids[,predabund]  ??
# if so: replace grids/gridslon/gridslat/predabund with x/y/z etc in byx/byy generator
# done
# and remove those those terms in gbm.map function : done
# and references to them in gbm.auto: done

####gbm.auto: BnW####
# heatcol defaulting to colours not b&w, why? Multiple arguments to same parameter. Fixed now I've set a default?
# in gbm.map, if people define heatcols, and zero=TRUE, L51 will overwrite heatcol with (alpha + actual heatcolours)
# zero defaults to true. it will always override heatcols.
# DONE!

####gbm.auto tc default####
# should be 2,length(expvar) not 2,5. How to do though? Make not default, if not set by user, sets to that?
# DONE!