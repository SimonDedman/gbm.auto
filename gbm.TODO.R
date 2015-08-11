### BRT code todo ####

ZI TEST
# so far it's proving tricky to automate this but humans do it easily.
# Heuristic IF statements? If >x% of data == 0 then isZI=TRUE? See ILL journal article.

Bag fraction optimiser
# Trial & error iterative approach to determine the optimal bag fraction for a data set? How? Stop @ whole percentages?
# Try OPTIM function & see http://r.789695.n4.nabble.com/Optimization-in-R-similar-to-MS-Excel-Solver-td4646759.html
# Possibly have an option to start with this in the R function.

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
#
# See:
# C:\Users\Simon\Dropbox\Galway\Analysis\R\Coilin R code\rmpi_example.R

Processing time estimate
#  option to have R poll the computer and, when you say go, popup a box saying
# "Based on the parameters you've selected to try, this will run X models on a Y-item-sized dataset (Y = variables * count)
# which may take about Z minutes based on your processor. You have a multicore processor so it will take Z/#processors (time)
# if you have multicore processing enabled - see here LINK" OK/Cancel
# Edit the BRT progrss counter to account for multiple resvars: currently does e.g. n/8 then loops back to 1. Not useful.
# Also could print the current resvar name. Maybe do a running time thing? "This is BRT N of X, Y% complete, took time Z, total
# expected time AA, time remaining AB

Clean workspace - dump fails: TRY rm() rather than dump() - dont need this for function
# rm(list = ls()) removes everything
# running function the dumping everyithng means R is still using 1.7gb of RAM!
# rm all model-generated files after each run? List what they are here.

3D PLOT
# what to do? Bother with it?

Allow deep subfunction calls. So deep.
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