#' Automated Boosted Regression Tree modelling and mapping suite
#'
#' Automates delta log normal boosted regression trees abundance prediction.
#' Loops through all permutations of parameters provided (learning
#' rate, tree complexity, bag fraction), chooses the best, then simplifies it.
#' Generates line, dot and bar plots, and outputs these and the predictions
#' and a report of all variables used, statistics for tests, variable
#' interactions, predictors used and dropped, etc. If selected, generates
#' predicted abundance maps, and Unrepresentativeness surfaces.
#' See www.github.com/SimonDedman/gbm.auto for issues, feedback, and development
#' suggestions. See SimonDedman.com for links to walkthrough paper, and papers
#' and thesis published using this package.
#'
#' @param grids Explantory data to predict to. Import with (e.g.) read.csv and
#' specify object name. Defaults to NULL (won't predict to grids)
#' @param samples Explanatory and response variables to predict from. Keep col
#' names short (~17 characters max), no odd characters, spaces, starting
#' numerals or terminal periods. Spaces may be converted to periods in directory
#' names, underscores won't. Can be a subset of a large dataset.
#' @param expvar List of names or column numbers of explanatory variables in
#' 'samples': c(1,3,6) or c("Temp","Sal"). No default
#' @param resvar Name or column number(s) of response variable in samples: 12,
#' c(1,4), "Rockfish". No default. Column name is ideally species name
#' @param tc Permutations of tree complexity allowed, can be vector with
#' the largest sized number no larger than the number of explanatory variables
#' e.g. c(2,7), or a list of 2 single numbers or vectors, the first to be passed
#' to the binary BRT, the second to the Gaussian, e.g. tc = list(c(2,6), 2) or
#' list(6, c(2,6))
#' @param lr Permutations of learning rate allowed. Can be a vector or a list of
#' 2 single numbers or vectors, the first to be passed to the binary BRT, the
#' second to the Gaussian, e.g. lr = list(c(0.01,0.02),0.0001) or
#' list(0.01,c(0.001, 0.0005))
#' @param bf Permutations of bag fraction allowed, can be single number, vector
#' or list, per tc and lr. Defaults to 0.5
#' @param n.trees from gbm.step, number of initial trees to fit. Can be
#' single or list but not vector i.e. list(fam1,fam2)
#' @param ZI are data zero-inflated? TRUE FALSE "CHECK". TRUE: delta BRT,
#' log-normalised Gaus, reverse log-norm and bias corrected. FALSE: do Gaussian
#' only, no log-normalisation. CHECK: Tests data for you. Default is CHECK.
#' @param fam1 probability distribution family for 1st part of delta process,
#' defaults to "bernoulli"
#' @param fam2 probability distribution family for 2nd part of delta process,
#' defaults to "gaussian"
#' @param simp Try simplfying best BRTs?
#' @param gridslat Column number for latitude in 'grids'
#' @param gridslon Column number for longitude in 'grids'
#' @param multiplot create matrix plot of all line files? Default true
#' turn off if big n of exp vars causes an error due to margin size problems.
#' @param cols Barplot colour vector. Assignment in order of explanatory
#' variables. Default 1*white: white bars black borders. '1*' repeats
#' @param linesfiles Save individual line plots' data as csv's? Default TRUE
#' @param smooth Apply a smoother to the line plots? Default FALSE
#' @param savegbm Save gbm objects and make available in environment after
#' running? Open with load("Bin_Best_Model") Default TRUE
#' @param loadgbm Relative or absolute location of folder containing
#' Bin_Best_Model and Gaus_Best_Model. If set will skip BRT calculations and do
#' predicted maps and csvs. Avoids re-running BRT models again (the slow bit),
#' can run normally once with savegbm=T then multiple times with new grids &
#' loadgbm to predict to multiple grids e.g. different seasons, areas, etc.
#' Default NULL, character vector, "./" for working directory
#' @param varint Calculate variable interactions? Default:TRUE, FALSE for error:
#' "contrasts can be applied only to factors with 2 or more levels"
#' @param map Save abundance map png files?
#' @param shape Set coast shapefile, else bounds calculated by gbm.map which
#' then calls gbm.basemap to download and autogenerate the base map. Read in
#' existing files by installing the shapefiles package then
#' DesiredMapName <- read.shapefile("ShapeFileName")
#' omitting the .shp extension
#' @param RSB Run Unrepresentativeness surface builder? Default TRUE
#' @param BnW Repeat maps in black and white e.g. for print journals. Default TRUE
#' @param alerts Play sounds to mark progress steps. Default TRUE but running
#' multiple small BRTs in a row (e.g. gbm.loop) can cause RStudio to crash
#' @param pngtype Filetype for png files, alternatively try "quartz"
#' @param gaus Do Gaussian runs as well as Bin? Default TRUE
#' @param MLEvaluate do machine learning evaluation metrics & plots? Default TRUE
#' @param ... Optional arguments for zero in breaks.grid in gbm.map, legend in
#' legend.grid in gbm.map, and gbm.step (dismo package) arguments n.trees and
#' max.trees, both of which can be added in list(1,2) format to pass to fam1 and 2
#'
#' @return Line, dot and bar plots, a report of all variables used, statistics
#' for tests, variable interactions, predictors used and dropped, etc. If
#' selected generates predicted abundance maps, and Unrepresentativeness surface
#'
#' @details Errors and their origins:
#'
#' 0. install ERROR: dependencies ‘rgdal’, ‘rgeos’ are not available for package ‘gbm.auto’
#' for linux/*buntu systems, in terminal, type
#' sudo apt install libgeos-dev
#' sudo apt install libproj-dev
#' sudo apt install libgdal-dev
#'
#' 1. Error in FUN(X[[i]], ...) : only defined on a data frame with all numeric variables
#' > Check your variable types are correct, e.g. numerics haven't been imported
#' as factors because there's an errant first row of text information before the
#' data. Remove NA rows from the response variable if present: convert blank
#' cells to NA on import with read.csv(x, na.strings = "") then
#' samples2 <- samples1[-which(is.na(samples[,resvar_column_number])),]
#'
#' 2. At bf=0.5, if nrows <= 42 gbm.step will crash
#' > Use gbm.bfcheck to determine optimal viable bf size
#'
#' 3. Maps/plots dont work/output
#' > If on a Mac, try changing pngtype to "quartz"
#'
#' 4. Error in while (delta.deviance > tolerance.test AMPERSAND n.fitted < max.trees)  :
#'  missing value where TRUE/FALSE needed
#' > If running a zero-inflated delta model (bernoilli/bin & gaussian/gaus),
#' Data are expected to contain zeroes (lots of them in zero-inflated cases),
#' have you already filtered them out, i.e. are only testing the positive cases?
#' Or do you only have positive cases? If so only run (e.g.) Gaussian: set ZI to FALSE
#'
#' 5. Error in round(gbm.object$cv.statistics$deviance.mean, 4) : non-numeric
#' argument to mathematical function
#' > LR or BF probably too low in earlier BRT (normally Gaus run with highest TC)
#'
#' 6. Error in if (n.trees > x$n.trees) { : argument is of length zero}
#' > LR or BF probably too low in earlier BRT (normally Gaus run with highest TC)
#'
#' 7. Error in gbm.fit(x, y, offset = offset, distribution = distribution, w = w)
#' The dataset size is too small or subsampling rate is too large:
#' nTrain*bag.fraction <= n.minobsinnode
#' > LR or BF probably too low in earlier BRT (normally Gaus run with highest TC)
#' It may be that you don't have enough positive samples to run BRT modelling
#' Run gbm.bfcheck to check recommended minimum BF size
#'
#' 8. Warning message: In cor(y_i, u_i) : the standard deviation is zero
#' > LR or BF probably too low in earlier BRT (normally Gaus run with highest TC)
#' It may be that you don't have enough positive samples to run BRT modelling
#' Run gbm.bfcheck to check recommended minimum BF size
#'
#' 9. Anomalous values can obfuscate clarity in line plots e.g. salinity range
#' 32:35ppm  but dataset has errant 0 value: plot axis will be 0:35 and 99.99%
#' of the data will be in the tiny bit at the right. Clean your data beforehand
#'
#' 10. Error in plot.new() : figure margins too large:
#' > In RStudio, adjust plot frame (usually bottom right) to increase its size
#' Still fails? Set multiplot=FALSE
#'
#' 11. Error in dev.print(file = paste0("./", names(samples[i]), "/pred_dev_bin.jpeg"),
#' : can only print from a screen device
#' > An earlier failed run (e.g. LR/BF too low) left a plotting device open.
#' Close it with: dev.off()
#'
#' 12. RStudio crashed: set alerts=F and pause cloud sync programs if outputting
#' to a synced folder
#'
#' 13. Error in grDevices::dev.copy(device = function (filename = "Rplot%03d.jpeg",
#' could not open file './P_PECTINATA..../pred_dev_bin.jpeg' (or similar)
#' > Your resvar column name contains an illegal character e.g. / & ' _
#' Fix with colnames(samples)[n] <- "BetterName"
#'
#' @examples gbm.auto(expvar = c(4:8, 10), resvar = 11, grids = mygrids,
#' tc = c(2,7), lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
#' @export
#' @import dismo
#' @importFrom beepr beep
#' @importFrom gbm plot.gbm
#'
gbm.auto <- function(
  grids = NULL,         # explantory data to predict to. Import with (e.g.)
  # read.csv and specify object name. Defaults to NULL (won't predict to grids)
  samples,  # explanatory and response variables to predict from.
  # Keep col names short, no odd characters, starting numerals or terminal periods
  # Spaces may be converted to periods in directory names, underscores won't.
  # Can be a subset
  expvar,               # list of column numbers of explanatory variables in
  # 'samples', expected e.g. c(1,35,67,etc.). No default
  resvar,               # column number(s) of response variable (e.g. CPUE) in
  # samples, e.g. 12 or c(4,5,6). No default. Column name should be species name
  tc = c(2),            # permutations of tree complexity allowed, can be a
  # vector with the largest sized number no larger than the number of
  # explanatory variables e.g. c(2,7), or a list of 2 single numbers or vectors,
  # the first to be passed to the binary BRT, the second to the Gaussian, e.g.
  # tc = list(c(2,6), 2) or list(6, c(2,6))
  lr = c(0.01,0.005),   # permutations of learning rate allowed. Can be a
  # vector or a list of 2 single numbers or vectors, the first to be passed to
  # the binary BRT, the second to the Gaussian, e.g.
  # lr = list(c(0.01,0.02),0.0001) or list(0.01,c(0.001, 0.0005))
  bf = 0.5,             # permutations of bag fraction allowed, can be single
  # number, vector or list, per tc and lr
  n.trees = 50,         # from gbm.step, number of initial trees to fit. Can be
  # single or list but not vector i.e. list(fam1,fam2)
  ZI = "CHECK",         # are data zero-inflated? TRUE/FALSE/"CHECK".
  # TRUE: delta BRT, log-normalised Gaus, reverse log-norm and bias corrected.
  # FALSE: do Gaussian only, no log-normalisation.
  # CHECK: Tests data for you. Default is TRUE.
  fam1 = "bernoulli",   # probability distribution family for 1st part of delta
  # process, defaults to "bernoulli",
  fam2 = "gaussian",   # probability distribution family for 2nd part of delta
  # process, defaults to "gaussian",
  simp = TRUE,          # try simplfying best BRTs?
  gridslat = 2,         # column number for latitude in 'grids'
  gridslon = 1,         # column number for longitude in 'grids'
  multiplot = TRUE,     # create matrix plot of all line files? Default true
  # turn off if big n of exp vars causes an error due to margin size problems.
  cols = grey.colors(1,1,1), # barplot colour vector. Assignment in order of
  # explanatory variables. Default 1*white: white bars black borders. '1*' repeats
  linesfiles = TRUE,    # save individual line plots' data as csv's?
  smooth = FALSE,       # apply a smoother to the line plots? Default FALSE
  savegbm = TRUE,       # save gbm objects and make available in environment after running? Open with load("Bin_Best_Model")
  loadgbm = NULL,       # relative or absolute location of folder containing
  # Bin_Best_Model and Gaus_Best_Model. If set will skip BRT calculations and do
  # predicted maps and csvs. Default NULL, character vector, "./" for working directory
  varint = TRUE,        # calculate variable interactions? Default:TRUE, FALSE
  # for error "contrasts can be applied only to factors with 2 or more levels"
  map = TRUE,           # save abundance map png files?
  shape = NULL,         # set coast shapefile, else bounds calculated by gbm.map
  # which then calls gbm.basemap to download and autogenerate the base map.
  RSB = TRUE,           # run Unrepresentativeness surface builder?
  BnW = TRUE,           # repeat maps in black and white e.g. for print journals
  alerts = TRUE,        # play sounds to mark progress steps. Running many small
  # BRTs e.g. gbm.loop can cause RStudio to crash, if so set this to FALSE
  pngtype = "cairo-png",# filetype for png files, alternatively try "quartz"
  gaus = TRUE,          # do Gaussian runs as well as Bin? Default TRUE.
  MLEvaluate = TRUE,    # do machine learning evaluation metrics & plots? Default TRUE
  ...)                  # Optional arguments for zero in breaks.grid in gbm.map,
# legend in legend.grid in gbm.map, and gbm.step (dismo
# package) argument max.trees and others.

{
  # Generalised Boosting Model / Boosted Regression Tree process chain automater.
  # Simon Dedman, 2012-6 simondedman@gmail.com github.com/SimonDedman/gbm.auto

  # Function to automate the many steps required to use boosted regression trees
  # to predict abundances in a delta process, i.e. binary (0/1) proportion
  # prediction coupled with presence-only abundance prediction to give total
  # prediction. Loops through all permutations of parameters provided (learning
  # rate, tree complexity, bag fraction), chooses the best, then tries to simplify
  # that. Generates line, dot and bar plots, and outputs these and the predictions
  # and a report of all variables used, statistics for tests, variable
  # interactions, predictors used and dropped, etc.. If selected, generates
  # predicted abundance maps, and Unrepresentativeness surfaces.
  #
  # Underlying functions are from packages gbm and dismo, functions from Elith
  # et al. 2008 (bundled as gbm.utils.R), mapplots, and my own functions gbm.map,
  # gbm.rsb, gbm.valuemap, gbm.cons, gbm.basemap

  ####1. Check packages, start loop####
  if (!require(gbm)) {stop("you need to install the gbm package to run this function")}
  if (!require(dismo)) {stop("you need to install the dismo package to run this function")}
  if (alerts) if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
  if (map) if (!require(mapplots)) {stop("you need to install the mapplots package to run this function")}
  if (map) if (!exists("gbm.map")) {stop("you need to install the gbm.map function to run this function")}
  if (RSB) if (!exists("gbm.rsb")) {stop("you need to install the gbm.rsb function to run this function")}
  if (RSB) if (!exists("gbm.map")) {stop("you need to install the gbm.map function to run this function")}
  if (!is.null(grids)) if (!exists("gbm.predict.grids")) {stop("you need to install the gbm.predict.grids function from gbm.utils.R to run this function")}
  if (!exists("roc")) {stop("you need to install the roc function from gbm.utils.R to run this function")}
  if (!exists("calibration")) {stop("you need to install the calibration function from gbm.utils.R to run this function")}
  require(gbm)
  require(dismo)
  if (alerts) require(beepr)

  # create basemap using gbm.basemap & these bounds, else basemap will be called for every map
  if (!is.null(grids)) if (map) { # create basemap grids not null, map requested, basemap not provided
    if (is.null(shape)) {
      if (!exists("gbm.basemap")) {stop("you need to install gbm.basemap to run this function")}
      bounds = c(range(grids[,gridslon]),range(grids[,gridslat]))
      #create standard bounds from data, and extra bounds for map aesthetic
      xmid <- mean(bounds[1:2])
      ymid <- mean(bounds[3:4])
      xextramax <- ((bounds[2] - xmid) * 1.6) + xmid
      xextramin <- xmid - ((xmid - bounds[1]) * 1.6)
      yextramax <- ((bounds[4] - ymid) * 1.6) + ymid
      yextramin <- ymid - ((ymid - bounds[3]) * 1.6)
      extrabounds <- c(xextramin, xextramax, yextramin, yextramax)
      shape <- gbm.basemap(bounds = extrabounds)
    }}

  wd <- getwd() #store original working directory
  if (alerts) options(error = function() {beep(9)  # give warning noise if it fails
    graphics.off() # kill all graphics devices
    setwd(wd)}) # reinstate original working directory

  expvarnames <- names(samples[expvar]) # list of explanatory variable names
  expvarcols <- cbind(cols[1:length(expvarnames)],expvarnames) # assign explanatory variables to colours

  if (is.list(tc)) { # if lists entered for tc lr or bf, split them to bin and gaus
    if (length(tc) > 2) {stop("Only 2 tc list items allowed: 1 per family")}
    tcgaus <- tc[[2]]
    tc <- tc[[1]]
  } else {tcgaus <- tc} # else make the gaus object the same as the bin

  if (is.list(lr)) {
    if (length(lr) > 2) {stop("Only 2 lr list items allowed: 1 per family")}
    lrgaus <- lr[[2]]
    lr <- lr[[1]]
  } else {lrgaus <- lr}

  if (is.list(bf)) {
    if (length(bf) > 2) {stop("Only 2 bf list items allowed: 1 per family")}
    bfgaus <- bf[[2]]
    bf <- bf[[1]]
  } else {bfgaus <- bf}

  if (is.list(n.trees)) { # if list entered n.trees, split to fam1 and fam2
    if (length(n.trees) > 2) {stop("Only 2 n.trees list items allowed: 1 per family")}
    ntf1 <- n.trees[[1]]
    ntf2 <- n.trees[[2]]
  } else {
    ntf1 <- n.trees
    ntf2 <- n.trees} # else make fam1 and fam2 the same

  for (i in resvar) { # loop everything for each response variable (e.g. species)
    dir.create(names(samples[i])) # create resvar-named directory for outputs
    m = 1 # Gaus only loop counter to allow best gaus BRT choice
    n = 1   # Print counter for all loops of BRT combos & best bin BRT choice
    if (!is.null(grids)) if (!all(expvarnames %in% names(grids))) stop(print("Expvar column names in samples but missing from grids:"), print(expvarnames[which(!expvarnames %in% names(grids))]))
    if (anyNA(samples[i])) stop("Response variable range contains NA values, please filter out these rows with: mysamples <- mysamples[-which(is.na(mysamples[resvar])),]")

    ####2. ZI check & log####
    # if user has asked code to check for ZI, check it & set new ZI status
    if (ZI == "CHECK") if (sum(samples[,i] == 0,na.rm = TRUE)/length(samples[,i]) >= 0.5) ZI = TRUE else ZI = FALSE
    # ensure resvar has zeroes (expects mix of successful & unsuccessful samples for bernoulli/binary runs)
    if (ZI == F) if (min(samples[i]) > 0) print("No zeroes in response variable. If using a zero inflated model, Method expects unsuccessful, as well as successful, samples")

    # create binary (0/1) response variable, for bernoulli BRTs
    samples$brv <- ifelse(samples[i] > 0, 1, 0)
    brvcol <- which(colnames(samples) == "brv") # brv column number for BRT

    # create logged response variable, for gaussian BRTs when data are zero-inflated (otherwise just use resvar directly)
    logem <- log(samples[,i]) # logs resvar i.e. containing zeroes
    dont  <- samples[,i]
    if (ZI) {samples$grv <- logem} else {samples$grv <- dont}
    grvcol <- which(colnames(samples) == "grv") # grv column number for BRT
    grv_yes <- subset(samples, grv >= 0) # nonzero subset for gaussian BRTs
    # actually not nonzero but 'not -Inf' since zeroes logged to "-Inf"
    # Change this to grv_yes <- samples if using a hurdle model including zeroes
    # need to keep it as log1p in that case?
    # ">=" same as just ">" since any nonzero in samples will be logged to a
    # nonzero & any 0 will be logged to -Inf so there'll be no zeroes

    if (is.null(loadgbm)) { #if loadgbm is NULL i.e. you're running BRTs not
      # predicting from existing models. Skip to L863

      ####3. Begin Report####
      #reportcolno = (3 + (5*(length(tc)*length(lr)*length(bf))) + 14)
      if (ZI) {
        reportcolno = 3 + (length(tc)*length(lr)*length(bf)) + (length(tcgaus)*length(lrgaus)*length(bfgaus)) + 14
        # if only 1 permutation, = 19
      } else {
        reportcolno = 3 + (length(tcgaus)*length(lrgaus)*length(bfgaus)) + 7
        # if only 1 permutation = 11
      }
      if (!gaus) reportcolno = 3 + (length(tc)*length(lr)*length(bf)) + 7
      # if only 1 permutation = 11

      # calculate number of columns for report:
      # 3: expvar names, resvar name, ZI state
      # 14: best bin brt, best gaus brt,
      # Bin_BRT_simp predictors kept (ordered), Bin_BRT_simp predictors dropped,
      # Gaus_BRT_simp predictors kept (ordered),Gaus_BRT_simp predictors dropped,
      # Simplified Binary BRT stats, Simplified Gaussian BRT stats,
      # Best Binary BRT variables, Relative Influence (Bin),
      # Best Gaussian BRT variables, Relative Influence (Gaus),
      # Biggest Interactions (Bin), Biggest Interactions (Gaus)
      # + 5 elements for each loop: parameter combo n (tc lr & bf values),
      # Bin BRT n stats, Bin BRT n name
      # Gaus BRT n stats, Gaus BRT n name
      Report <- data.frame(matrix(NA, nrow = (max(6,length(expvar))), ncol = (reportcolno)))
      # build blank df, rows=biggest of 6 (max static row number of stats) or n of exp. vars
      colnames(Report) <- c("Explanatory Variables","Response Variables","Zero Inflated?") # populate static colnames 1:3
      # name bin columns if ZI
      if (!gaus) {colnames(Report)[(reportcolno - 6):reportcolno] <- c("Best Binary BRT",
                                                                       "Bin_BRT_simp predictors dropped",
                                                                       "Bin_BRT_simp predictors kept",
                                                                       "Simplified Binary BRT stats",
                                                                       "Best Binary BRT variables",
                                                                       "Relative Influence (Bin)",
                                                                       "Biggest Interactions (Bin)")
      } else {if (ZI) {colnames(Report)[(reportcolno - 13):(reportcolno - 7)] <- c("Best Binary BRT",
                                                                                   "Bin_BRT_simp predictors dropped",
                                                                                   "Bin_BRT_simp predictors kept",
                                                                                   "Simplified Binary BRT stats",
                                                                                   "Best Binary BRT variables",
                                                                                   "Relative Influence (Bin)",
                                                                                   "Biggest Interactions (Bin)")}
        colnames(Report)[(reportcolno - 6):reportcolno] <- c("Best Gaussian BRT",
                                                             "Gaus_BRT_simp predictors dropped",
                                                             "Gaus_BRT_simp predictors kept",
                                                             "Simplified Gaussian BRT stats",
                                                             "Best Gaussian BRT variables",
                                                             "Relative Influence (Gaus)",
                                                             "Biggest Interactions (Gaus)")}
      # populate the final 14 column names
      Report[1:length(expvar),1] <- names(samples[expvar]) # put expvar names in first column
      Report[1,2] <- names(samples[i]) # put resvar in col 2
      Report[1,3] <- ZI # ZI in col 3

      Bin_Best_Score <- 0 # create blanks for best results to use in loops
      Bin_Best_Model <- 0
      Gaus_Best_Score <- 0
      Gaus_Best_Model <- 0

      # Begin bin loops
      if (ZI) {  # don't do if ZI=FALSE
        for (j in tc) {   # list permutations of tree complexity allowed
          for (k in lr) {   # list permutations of learning rate allowed
            for (l in bf) {   # list permutations of bag fraction allowed

              ####4. Binomial BRT####
              assign(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l),
                     gbm.step(data = samples,
                              gbm.x = expvar,
                              gbm.y = brvcol,
                              family = fam1,
                              tree.complexity = j,
                              learning.rate = k,
                              bag.fraction = l,
                              n.trees = ntf1,
                              ...)
              )
              dev.print(file = paste0("./",names(samples[i]),"/pred_dev_bin.jpeg"), device = jpeg, width = 600)
              print(paste0("Done Bin_BRT",".tc",j,".lr",k,".bf",l))
              print(warnings())
              ####5. Select best bin model####
              if (n == 1) {
                Bin_Best_Score <- get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
                Bin_Best_Model <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)
              }  else if (get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]] > Bin_Best_Score) {
                Bin_Best_Score <- get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
                Bin_Best_Model <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)
              }

              ####6. Add bin stats to report####
              # don't do if ZI=FALSE. bin BRT stats
              if (ZI) {Report[1:6,(3 + n)] <- c(paste0("trees: ",get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$n.trees),
                                                paste0("Training Data Correlation: ",get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]),
                                                paste0("CV Mean Deviance: ",get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.mean),
                                                paste0("CV Deviance SE: ",get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.se),
                                                paste0("CV Mean Correlation: ",get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.mean),
                                                paste0("CV Correlation SE: ",get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.se))
              # bin BRT name
              colnames(Report)[3 + n] <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)
              } # close ZI if

              if (alerts) beep(2) # progress printer, right aligned
              if (gaus) {
                print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tc)*length(lr)*length(bf)) + (length(tcgaus)*length(lrgaus)*length(bfgaus)), "     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
              } else {
                print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tc)*length(lr)*length(bf)), "     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
              }
              n <- n + 1   # Add to print counter
            } # close bf
          } # close lr
        } # close tc
      } # close ZI option, making all bin BRT objects & continuing through model selection

      # Begin Gaus loops
      if (gaus) for (j in tcgaus) {   # list permutations of tree complexity allowed
        for (k in lrgaus) {   # list permutations of learning rate allowed
          for (l in bfgaus) {   # list permutations of bag fraction allowed
            ####7. Gaussian BRT####
            assign(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l),
                   gbm.step(data = grv_yes,
                            gbm.x = expvar,
                            gbm.y = grvcol,
                            family = fam2,
                            tree.complexity = j,
                            learning.rate = k,
                            bag.fraction = l,
                            n.trees = ntf2,
                            ...)
            )
            dev.print(file = paste0("./",names(samples[i]),"/pred_dev_gaus.jpeg"), device = jpeg, width = 600)
            print(paste0("Done Gaus_BRT",".tc",j,".lr",k,".bf",l))
            print(warnings())
            ####8. Select best Gaus model####
            if (m == 1)
            {Gaus_Best_Score <- get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
            Gaus_Best_Model <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)
            } else if (get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]] > Gaus_Best_Score)
            {Gaus_Best_Score <- get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
            Gaus_Best_Model <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)}

            if (alerts) beep(2) # progress printer, right aligned for visibility
            if (ZI) {print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tc)*length(lr)*length(bf)) + (length(tcgaus)*length(lrgaus)*length(bfgaus)),"     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
            } else {print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tcgaus)*length(lrgaus)*length(bfgaus)),"     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))}

            ####9. Add gaus stats to report####
            Report[1:6,(3 + n)] <- c(paste0("trees: ",get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$n.trees),
                                     paste0("Training Data Correlation: ",get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]),
                                     paste0("CV Mean Deviance: ",get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.mean),
                                     paste0("CV Deviance SE: ",get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.se),
                                     paste0("CV Mean Correlation: ",get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.mean),
                                     paste0("CV Correlation SE: ",get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.se))
            # Gaus BRT name
            colnames(Report)[3 + n] <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)

            n <- n + 1 # Add to print/loop counter for every bin or gaus BRT loop
            m <- m + 1 # Add to loop counter for Gaus best model selection

          } # close bfgaus
        } # close lrgaus
      } # close tcgaus, making all Gaus BRT objects & continuing through model selection

      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        Closed Loops         XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####10. Test simplification benefit, do so if better####
      # copy Bin/Gaus_Best_Model to Name in case not created by Simp
      Bin_Best_Name <- Bin_Best_Model
      if (gaus) Gaus_Best_Name <- Gaus_Best_Model

      # if simp TRUE & ZI=TRUE, run simplification test on best bin model
      if (simp) {
        if (ZI) {Bin_Best_Simp_Check <- gbm.simplify(get(Bin_Best_Model))
        dev.print(file = paste0("./",names(samples[i]),"/simp_drops_bin.jpeg"), device = jpeg, width = 600)
        # if best number of variables to remove isn't 0 (i.e. it's worth simplifying),
        # re-run best model (Bin_Best_Model, using gbm.call to get its values) with
        # just-calculated best number of variables to remove, removed. gbm.x asks which
        # number of drops has the minimum mean (lowest point on the line) & that calls
        # up the list of predictor variables with those removed, from $pred.list
        if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0) {
          assign("Bin_Best_Simp",
                 gbm.step(data = samples,
                          gbm.x = Bin_Best_Simp_Check$pred.list[[which.min(Bin_Best_Simp_Check$deviance.summary$mean)]],
                          gbm.y = get(Bin_Best_Model)$gbm.call$gbm.y,
                          tree.complexity = get(Bin_Best_Model)$gbm.call$tree.complexity,
                          learning.rate = get(Bin_Best_Model)$gbm.call$learning.rate,
                          family = get(Bin_Best_Model)$gbm.call$family,
                          bag.fraction = get(Bin_Best_Model)$gbm.call$bag.fraction,
                          ...))
          dev.print(file = paste0("./",names(samples[i]),"/pred_dev_bin_simp.jpeg"), device = jpeg, width = 600)}

        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Bin model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        } # close ZI

        # Same for Gaus
        if (gaus) {Gaus_Best_Simp_Check <- gbm.simplify(get(Gaus_Best_Model))
        dev.print(file = paste0("./",names(samples[i]),"/simp_drops_gaus.jpeg"), device = jpeg, width = 600)
        if (min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0)
          assign("Gaus_Best_Simp",
                 gbm.step(data = grv_yes,
                          gbm.x = Gaus_Best_Simp_Check$pred.list[[which.min(Gaus_Best_Simp_Check$deviance.summary$mean)]],
                          gbm.y = get(Gaus_Best_Model)$gbm.call$gbm.y,
                          tree.complexity = get(Gaus_Best_Model)$gbm.call$tree.complexity,
                          learning.rate = get(Gaus_Best_Model)$gbm.call$learning.rate,
                          family = get(Gaus_Best_Model)$gbm.call$family,
                          bag.fraction = get(Gaus_Best_Model)$gbm.call$bag.fraction,
                          ...))
        dev.print(file = paste0("./",names(samples[i]),"/pred_dev_gaus_simp.jpeg"), device = jpeg, width = 600)
        if (alerts) beep(2)
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Gaus model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))}

        ## Select final best models
        if (ZI) {  # don't do if ZI=FALSE. If Bin_Best has a simplified model:
          if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0)
            # & if the simplified model has better correlation than Bin_Best itself
            if (Bin_Best_Simp$self.statistics$correlation > Bin_Best_Score[1])
              # then replace Bin_Best score/model values with those from the simplified model
            {Bin_Best_Score <- Bin_Best_Simp$self.statistics$correlation
            Bin_Best_Name <- paste0(Bin_Best_Model, "_Simp")
            Bin_Best_Model <- "Bin_Best_Simp"}} # assign simp to best & close ZI

        # Same for Gaus:
        if (gaus) {if (min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0)
          if (Gaus_Best_Simp$self.statistics$correlation > Gaus_Best_Score[1])
          {Gaus_Best_Score <- Gaus_Best_Simp$self.statistics$correlation
          Gaus_Best_Name <- paste0(Gaus_Best_Model, "_Simp")
          Gaus_Best_Model <- "Gaus_Best_Simp"}}
      } # close simp optional

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Best models selected     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####11. Line plots####
      # All plots on one image for Bin & Gaus
      if (multiplot) { # don't do if multiplot=FALSE
        if (ZI) {  # don't do if ZI=FALSE
          op <- par(oma = c(5,7,1,1)) # younes
          par(mar = rep(2, 4)) # for Younes' Error in plot.new() : figure margins too large
          png(filename = paste0("./",names(samples[i]),"/Bin_Best_line.png"),
              width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
          gbm.plot(get(Bin_Best_Model),
                   n.plots = length(get(Bin_Best_Model)$contributions$var),
                   write.title = F, y.label = "Marginal Effect",
                   plot.layout = c(ceiling(sqrt(length(get(Bin_Best_Model)$contributions$var))),
                                   ifelse(sqrt(length(get(Bin_Best_Model)$contributions$var))
                                          - floor(sqrt(length(get(Bin_Best_Model)$contributions$var))) < 0.5,
                                          floor(sqrt(length(get(Bin_Best_Model)$contributions$var))),
                                          floor(sqrt(length(get(Bin_Best_Model)$contributions$var))) + 1)))
          dev.off()
          par(op)} # close ZI     # younes

        if (gaus) {png(filename = paste0("./",names(samples[i]),"/Gaus_Best_line.png"),
                       width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
          gbm.plot(get(Gaus_Best_Model),
                   n.plots = length(get(Gaus_Best_Model)$contributions$var),
                   write.title = F, y.label = "Marginal Effect",
                   plot.layout = c(ceiling(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
                                   ifelse(sqrt(length(get(Gaus_Best_Model)$contributions$var))
                                          - floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))) < 0.5,
                                          floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
                                          floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))) + 1)))
          dev.off()} #close plot device & gaus if
      } # close multiplot if

      # All plots individually, named by explanatory variable, bin & gaus
      if (ZI) {  # don't do if ZI=FALSE
        for (o in 1:length(get(Bin_Best_Model)$contributions$var)) {
          png(filename = paste0("./",names(samples[i]),"/Bin_Best_line_",as.character(get(Bin_Best_Model)$gbm.call$predictor.names[o]),".png"),
              width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = pngtype)
          par(mar = c(2.3,5,0.3,0.6), fig = c(0,1,0,1), las = 1, lwd = 8, bty = "n", mgp = c(1.25,0.5,0), xpd = NA)
          gbm.plot(get(Bin_Best_Model),
                   variable.no = o, # order of variable.no =! order of get(Bin_Best_Model)$contributions$var
                   n.plots = 1,
                   common.scale = FALSE, #added to try to get cvs values to match pngs
                   smooth = smooth,
                   rug = TRUE,
                   write.title = FALSE,
                   y.label = "",
                   x.label = NULL,
                   show.contrib = TRUE,
                   plot.layout = c(1, 1)) # ... for cex.axis, cex.lab etc
          mtext("Marginal Effect", side = 2, line = 4.05, las = 0)

          # create lines data to export to file. Need to recreate transformations from gbm.plot
          # Next 6 lines from GNG answer https://stats.stackexchange.com/a/144871/43360 which uses gbm.plot code
          if (linesfiles) {s <- match(get(Bin_Best_Model)$contributions$var[o],
                                      get(Bin_Best_Model)$gbm.call$predictor.names)

          # create dataframe
          plotgrid <- plot.gbm(get(Bin_Best_Model), s, return.grid = TRUE)

          #If factor variable
          if (is.factor(plotgrid[,1])) {
            plotgrid[,1] <- factor(plotgrid[,1], levels = levels(get(Bin_Best_Model)$gbm.call$dataframe[,get(Bin_Best_Model)$gbm.call$gbm.x[s]]))}

          # replace Y values in place with average-centred values
          plotgrid[,2] <- plotgrid[,2] - mean(plotgrid[,2])

          #Put Y values on a log scale
          plotgrid[,2] <- 1 / (1 + exp(-plotgrid[,2]))

          #Center the response to have zero mean over the data distribution
          plotgrid[,2] <- scale(plotgrid[,2], scale = FALSE)

          # write out csv
          write.csv(plotgrid, row.names = FALSE, na = "",
                    file = paste0("./", names(samples[i]), "/Bin_Best_line_",
                                  as.character(get(Bin_Best_Model)$contributions$var[o]),
                                  ".csv"))} #close linesfiles
          dev.off() }} # close ZI option

      if (gaus) {for (p in 1:length(get(Gaus_Best_Model)$contributions$var)) {
        png(filename = paste0("./",names(samples[i]),"/Gaus_Best_line_",as.character(get(Gaus_Best_Model)$gbm.call$predictor.names[p]),".png"),
            width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(2.3,5,0.3,0.6), fig = c(0,1,0,1), las = 1, lwd = 8, bty = "n", mgp = c(1.25,0.5,0), xpd = NA)
        gbm.plot(get(Gaus_Best_Model),
                 variable.no = p,
                 n.plots = 1,
                 common.scale = FALSE, #added to try to get cvs values to match pngs
                 smooth = smooth,
                 rug = TRUE,
                 write.title = FALSE,
                 y.label = "",
                 x.label = NULL,
                 show.contrib = TRUE,
                 plot.layout = c(1, 1))
        mtext("Marginal Effect", side = 2, line = 4.05, las = 0)

        if (linesfiles) {u <- match(get(Gaus_Best_Model)$contributions$var[p],
                                    get(Gaus_Best_Model)$gbm.call$predictor.names)
        plotgrid <- plot.gbm(get(Gaus_Best_Model), u, return.grid = TRUE)
        if (is.factor(plotgrid[,1])) {
          plotgrid[,1] <- factor(plotgrid[,1], levels = levels(get(Gaus_Best_Model)$gbm.call$dataframe[,get(Gaus_Best_Model)$gbm.call$gbm.x[u]]))}
        plotgrid[,2] <- plotgrid[,2] - mean(plotgrid[,2])
        plotgrid[,2] <- 1 / (1 + exp(-plotgrid[,2]))
        plotgrid[,2] <- scale(plotgrid[,2], scale = FALSE)
        write.csv(plotgrid, row.names = FALSE, na = "",
                  file = paste0("./", names(samples[i]), "/Gaus_Best_line_",
                                as.character(get(Gaus_Best_Model)$contributions$var[p]),
                                ".csv"))} #close linesfiles
        dev.off() }}

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Line plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####12. Dot plots####
      if (ZI) {  # don't do if ZI=FALSE
        png(filename = paste0("./",names(samples[i]),"/Bin_Best_dot.png"),
            width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
        gbm.plot.fits(get(Bin_Best_Model))
        dev.off()} # close ZI

      if (gaus) {png(filename = paste0("./",names(samples[i]),"/Gaus_Best_dot.png"),
                     width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
        gbm.plot.fits(get(Gaus_Best_Model))
        dev.off()}

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Dot plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####13. 3D plot TODO####
      # gbm.perspec(Bin_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)
      # gbm.perspec(Gaus_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)

      ####14. Bar plots of variable influence####
      if (ZI) {  # create tables. Don't do if ZI=FALSE
        Bin_Bars <- summary(get(Bin_Best_Model),
                            cBars = length(get(Bin_Best_Model)$var.names),
                            n.trees = get(Bin_Best_Model)$n.trees,
                            plotit = FALSE, order = TRUE, normalize = TRUE, las = 1, main = NULL)
        write.csv(Bin_Bars, file = paste0("./", names(samples[i]), "/Binary BRT Variable contributions.csv"), row.names = FALSE)} # close ZI

      if (gaus) {Gaus_Bars <- summary(get(Gaus_Best_Model),
                                      cBars = length(get(Gaus_Best_Model)$var.names),
                                      n.trees = get(Gaus_Best_Model)$n.trees,
                                      plotit = FALSE, order = TRUE, normalize = TRUE, las = 1, main = NULL)
      write.csv(Gaus_Bars, file = paste0("./", names(samples[i]), "/Gaussian BRT Variable contributions.csv"), row.names = FALSE)}

      if (alerts) beep(2)# progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Bar plot csvs created    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      if (ZI) {  # produce graphics. Don't do bin if ZI=FALSE
        pointlineseqbin <- seq(0, length(Bin_Bars[,2]) - 1, 1)
        png(filename = paste0("./",names(samples[i]),"/Bin_Bars.png"),
            width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "",
            type = pngtype)
        par(mar = c(2.5,0.3,0,0.5), fig = c(0,1,0,1), cex.lab = 0.5, mgp = c(1.5,0.5,0), cex = 1.3, lwd = 6)
        barplot(rev(Bin_Bars[,2]), cex.lab = 1.2, las = 1, # axis labs horizontal
                horiz = TRUE, # make horizontal
                xlab = "Relative Influence %", col = NA, border = NA, # no border, lwd redundant
                xlim = c(0, 2.5 + ceiling(max(Bin_Bars[,2]))),
                ylim = c(0, length(Bin_Bars[,2])), # figure height as a proportion of nBars
                beside = T) # juxtaposed not stacked
        #points(rev(Bin_Bars[,2]), pointlineseqbin, pch = 20, cex = 1.75, col = "black") #black dots at line ends
        revseq <- rev(pointlineseqbin)
        for (q in 1:length(Bin_Bars[,2])) {
          lines(c(0, Bin_Bars[q,2]), c(revseq[q], revseq[q]), col = "black", lwd = 8)}
        text(0.1, pointlineseqbin + (length(Bin_Bars[,2])/55), labels = rev(Bin_Bars[,1]), adj = 0, cex = 0.8)
        axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
        dev.off()} # close ZI

      if (gaus) {pointlineseqgaus <- seq(0, length(Gaus_Bars[,2]) - 1, 1)
      png(filename = paste0("./",names(samples[i]),"/Gaus_Bars.png"),
          width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "",
          type = pngtype)
      par(mar = c(2.5,0.3,0,0.5), fig = c(0,1,0,1), cex.lab = 0.5, mgp = c(1.5,0.5,0), cex = 1.3, lwd = 6)
      barplot(rev(Gaus_Bars[,2]), cex.lab = 1.2, las = 1, # axis labs horizontal
              horiz = TRUE, # make horizontal
              xlab = "Relative Influence %", col = NA, border = NA, # no border, lwd redundant
              xlim = c(0, 2.5 + ceiling(max(Gaus_Bars[,2]))),
              ylim = c(0, length(Gaus_Bars[,2])), # figure height as a proportion of nBars
              beside = T) # juxtaposed not stacked
      #points(rev(Gaus_Bars[,2]), pointlineseqgaus, pch = 20, cex = 1.75, col = "black")
      revseq <- rev(pointlineseqgaus)
      for (r in 1:length(Gaus_Bars[,2])) {
        lines(c(0, Gaus_Bars[r,2]), c(revseq[r], revseq[r]), col = "black", lwd = 8)}
      text(0.1, pointlineseqgaus + (length(Gaus_Bars[,2])/55), labels = rev(Gaus_Bars[,1]), adj = 0, cex = 0.8)
      axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
      dev.off()} #close PNG
      # col = rev(expvarcols[match(Bin_Bars[,1],expvarcols[,2]),1]), #in case I want to colour the bars/points later
      # elements of barplot lines+points code adapted from Jane Elith code donated to Agustín De Wysiecki & then shared with SD

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Bar plots plotted      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####15. Variable interactions####
      if (varint) {
        if (ZI) find.int_Bin <- gbm.interactions(get(Bin_Best_Model))
        if (gaus) find.int_Gaus <- gbm.interactions(get(Gaus_Best_Model))
        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Variable interactions done XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close varint if

      ####16. Save model objects####
      if (savegbm) { # Save model objects if switched on
        if (ZI) {Bin_Best_Model_Object <- get(Bin_Best_Model)
        Bin_Best_Model <<- Bin_Best_Model_Object}
        if (gaus) {Gaus_Best_Model_Object <- get(Gaus_Best_Model)
        Gaus_Best_Model <<- Gaus_Best_Model_Object
        save(Gaus_Best_Model_Object,file = paste0("./",names(samples[i]),"/Gaus_Best_Model"))}
        if (ZI) {save(Bin_Best_Model_Object,file = paste0("./",names(samples[i]),"/Bin_Best_Model"))} #only save bin if ZI=TRUE
        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Model objects saved     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      }

      ####17. Finalise & Write Report####
      if (ZI) { # only do bin bits if ZI; move 7 cols left if no gaus run
        if (gaus) {
          Report[1:5,(reportcolno - 13)] <- c(paste0("Model combo: ", Bin_Best_Name),
                                              paste0("Model CV score: ", Bin_Best_Score),
                                              paste0("Training data AUC score: ", get(Bin_Best_Model)$self.statistics$discrimination),
                                              paste0("CV AUC score: ", get(Bin_Best_Model)$cv.statistics$discrimination.mean),
                                              paste0("CV AUC se: ", get(Bin_Best_Model)$cv.statistics$discrimination.se))
        } else {
          Report[1:5,(reportcolno - 6)] <- c(paste0("Model combo: ", Bin_Best_Name),
                                             paste0("Model CV score: ", Bin_Best_Score),
                                             paste0("Training data AUC score: ", get(Bin_Best_Model)$self.statistics$discrimination),
                                             paste0("CV AUC score: ", get(Bin_Best_Model)$cv.statistics$discrimination.mean),
                                             paste0("CV AUC se: ", get(Bin_Best_Model)$cv.statistics$discrimination.se))
        } # close gaus bin report

        if (simp) { # bin & gaus simp stats (or no simp notes)
          if (gaus) {
            Report[1:dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1], (reportcolno - 12)] <- as.character(subset(Bin_Best_Simp_Check$final.drops, order > 0)$preds)
            # listing simp predictors kept: rows 1 to 'howevermany are left in simp' i.e. above 0
            # [1] is the first item of the dim list of rows & columns, i.e. rows
            Report[1:(length(Bin_Best_Simp_Check$final.drops$preds) - dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1]),(reportcolno - 11)] <-
              as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Bin_Best_Simp_Check$final.drops$preds)])
            # listing simp predictors dropped.
            if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0) {
              Report[1:6,(reportcolno - 10)] <- c(paste0("trees: ", Bin_Best_Simp$n.trees),
                                                  paste0("Training Data Correlation: ", Bin_Best_Simp$self.statistics$correlation[[1]]),
                                                  paste0("CV Mean Deviance: ", Bin_Best_Simp$cv.statistics$deviance.mean),
                                                  paste0("CV Deviance SE: ", Bin_Best_Simp$cv.statistics$deviance.se),
                                                  paste0("CV Mean Correlation: ", Bin_Best_Simp$cv.statistics$correlation.mean),
                                                  paste0("CV Correlation SE: ", Bin_Best_Simp$cv.statistics$correlation.se))
            } else {
              Report[1,(reportcolno - 10)] <- paste0("No simplification benefit")
            } # close min else
          } else { # bin predictors etc
            Report[1:dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1], (reportcolno - 5)] <- as.character(subset(Bin_Best_Simp_Check$final.drops, order > 0)$preds)
            Report[1:(length(Bin_Best_Simp_Check$final.drops$preds) - dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1]),(reportcolno - 4)] <-
              as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Bin_Best_Simp_Check$final.drops$preds)])
            if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0)
            {Report[1:6,(reportcolno - 3)] <- c(paste0("trees: ", Bin_Best_Simp$n.trees),
                                                paste0("Training Data Correlation: ", Bin_Best_Simp$self.statistics$correlation[[1]]),
                                                paste0("CV Mean Deviance: ", Bin_Best_Simp$cv.statistics$deviance.mean),
                                                paste0("CV Deviance SE: ", Bin_Best_Simp$cv.statistics$deviance.se),
                                                paste0("CV Mean Correlation: ", Bin_Best_Simp$cv.statistics$correlation.mean),
                                                paste0("CV Correlation SE: ", Bin_Best_Simp$cv.statistics$correlation.se))
            } else {Report[1,(reportcolno - 3)] <- paste0("No simplification benefit")} # close min else
          } # close bin half of bin/gaus option. Next line is 2nd half of simp option i.e. not simplified
        } else if (gaus) {
          Report[1,(reportcolno - 12):(reportcolno - 10)] <- c(paste0("simp turned off"),
                                                               paste0("simp turned off"),
                                                               paste0("simp turned off"))
        } else {# if not running gaus, report cols are changed so needs adjustment (not gaus not simp)
          Report[1,(reportcolno - 5):(reportcolno - 3)] <- c(paste0("simp turned off"),
                                                             paste0("simp turned off"),
                                                             paste0("simp turned off"))
        } # close 2nd half of simp else, i.e. nosimp bin. Closes bin & gaus simp stats (or no simp notes)

        if (gaus) {
          Report[1:(length(Bin_Bars[,1])),(reportcolno - 9)] <- as.character(Bin_Bars$var)
        } else { #bin only
          Report[1:(length(Bin_Bars[,1])),(reportcolno - 2)] <- as.character(Bin_Bars$var)
        } # close ifelse best bin variables rel inf names ordered

        if (gaus) {
          Report[1:(length(Bin_Bars[,2])),(reportcolno - 8)] <- as.character(Bin_Bars$rel.inf)
        } else {
          Report[1:(length(Bin_Bars[,2])),(reportcolno - 1)] <- as.character(Bin_Bars$rel.inf)
        } # close ifelse best bin variables rel inf scores

        if (varint) { # only do final variable interaction lines if varint=TRUE
          if (gaus) {
            Report[1:2,(reportcolno - 7)] <- c(paste0(find.int_Bin$rank.list$var1.names[1]," and ",find.int_Bin$rank.list$var2.names[1],". Size: ",find.int_Bin$rank.list$int.size[1]),
                                               paste0(find.int_Bin$rank.list$var1.names[2]," and ",find.int_Bin$rank.list$var2.names[2],". Size: ",find.int_Bin$rank.list$int.size[2]))
          } else { # close varint yes gaus yes
            Report[1:2,(reportcolno)] <- c(paste0(find.int_Bin$rank.list$var1.names[1]," and ",find.int_Bin$rank.list$var2.names[1],". Size: ",find.int_Bin$rank.list$int.size[1]),
                                           paste0(find.int_Bin$rank.list$var1.names[2]," and ",find.int_Bin$rank.list$var2.names[2],". Size: ",find.int_Bin$rank.list$int.size[2]))
          } # close varint yes gaus no
        } else { # varint no
          if (gaus) { # varint no gaus yes
            Report[1,(reportcolno - 7)] <- paste0("varint turned off")
          } else { # varint no gaus no
            Report[1,(reportcolno)] <- paste0("varint turned off")
          } # close not varint not gaus
        } # close not varint
      } # close ZI way further up start of report section

      if (gaus) {
        Report[1:2,(reportcolno - 6)] <- c(paste0("Model combo: ", Gaus_Best_Name), paste0("Model CV score: ", Gaus_Best_Score))
        if (simp) {
          Report[1:dim(subset(Gaus_Best_Simp_Check$final.drops,order > 0))[1], (reportcolno - 5)] <- as.character(subset(Gaus_Best_Simp_Check$final.drops ,order > 0)$preds)
          Report[1:(length(Gaus_Best_Simp_Check$final.drops$preds) - dim(subset(Gaus_Best_Simp_Check$final.drops, order > 0))[1]), (reportcolno - 4)] <-
            as.character(Gaus_Best_Simp_Check$final.drops$preds[((dim(subset(Gaus_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Gaus_Best_Simp_Check$final.drops$preds)])
          if (min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0) {
            Report[1:6,(reportcolno - 3)] <- c(paste0("trees: ", Gaus_Best_Simp$n.trees),
                                               paste0("Training Data Correlation: ", Gaus_Best_Simp$self.statistics$correlation[[1]]),
                                               paste0("CV Mean Deviance: ", Gaus_Best_Simp$cv.statistics$deviance.mean),
                                               paste0("CV Deviance SE: ", Gaus_Best_Simp$cv.statistics$deviance.se),
                                               paste0("CV Mean Correlation: ", Gaus_Best_Simp$cv.statistics$correlation.mean),
                                               paste0("CV Correlation SE: ", Gaus_Best_Simp$cv.statistics$correlation.se))
          } else { # close stats where simp benefit true, open note where no simp benefit
            Report[1,(reportcolno - 3)] <- paste0("No simplification benefit")
          } # close simp benefit check
        } else { # close gaus yes simp yes, do gaus yes simp no
          Report[1,(reportcolno - 5):(reportcolno - 3)] <- c(paste0("simp turned off"),
                                                             paste0("simp turned off"),
                                                             paste0("simp turned off"))
        } # close simp, still in gaus yes
        Report[1:(length(Gaus_Bars[,1])),(reportcolno - 2)] <- as.character(Gaus_Bars$var)
        Report[1:(length(Gaus_Bars[,2])),(reportcolno - 1)] <- as.character(Gaus_Bars$rel.inf)
        if (varint) { # gaus yes varint yes
          Report[1:2,(reportcolno)] <- c(paste0(find.int_Gaus$rank.list$var1.names[1]," and ",find.int_Gaus$rank.list$var2.names[1],". Size: ",find.int_Gaus$rank.list$int.size[1]),
                                         paste0(find.int_Gaus$rank.list$var1.names[2]," and ",find.int_Gaus$rank.list$var2.names[2],". Size: ",find.int_Gaus$rank.list$int.size[2]))
        } else { # gaus yes varint no
          Report[1,(reportcolno)] <- paste0("varint turned off")
        } # close varint yes no, still in gaus yes
      } # close gaus if

      write.csv(Report, row.names = FALSE, na = "", file = paste0("./", names(samples[i]), "/Report.csv"))

      #18. Machine learning evaluation metrics####
      if (MLEvaluate) { if (any(fam1 == "bernoulli", fam2 == "bernoulli")) {
        whichbin <- which(c(fam1 == "bernoulli", fam2 == "bernoulli"))
        if (whichbin == 1) getmodel <- "Bin_Best_Model" else getmodel <- "Gaus_Best_Model"
        preds <- predict.gbm(get(get(getmodel)),
                             samples,
                             n.trees = get(get(getmodel))$gbm.call$best.trees,
                             type = "response")
        #If type="response" then gbm converts back to the same scale as the outcome.
        # Currently the only effect this will have is returning probabilities for
        # bernoulli and expected counts for poisson. For the other distributions
        # "response" and "link" return the same. gbm:::predict.gbm

        # dev reported later but not used otherwise
        dev <- calc.deviance(obs = samples[, get(get(getmodel))$gbm.call$gbm.y],
                             pred = preds,
                             family = "bernoulli") # change fam if using
        # One of "binomial", "bernoulli", "poisson", "laplace", or "gaussian"
        samples <- cbind(samples, preds)
        pres <- samples[samples[, brvcol] == 1, "preds"] # check brvcol indexed properly, ditto last col is preds
        abs <- samples[samples[, brvcol] == 0, "preds"]
        e <- evaluate(p = pres,
                      a = abs)

        # Fielding, A. H. & J.F. Bell, 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24: 38-49
        # Liu, C., M. White & G. Newell, 2011. Measuring and comparing the accuracy of species distribution models with presence-absence data. Ecography 34: 232-243.
        MLEvalLength <- 31
        # Improve descriptions####
        MLEval <- data.frame(Statistic = rep(NA, MLEvalLength),
                                Description = rep(NA, MLEvalLength),
                                Value = rep(NA, MLEvalLength))
        MLEval[1,] <- c("Presence",
                           "n of presence data used",
                           round(e@np), 3)
        MLEval[2,] <- c("Absence",
                           "n of absence data used",
                        round(e@na), 3)
        MLEval[3,] <- c("AUC",
                           "Area under the receiver operator (ROC) curve",
                        round(e@auc), 3)
        if (length(e@pauc) == 0) e@pauc <- 0 # pauc may be missing, numeric(0), if so replace with 0
        MLEval[4,] <- c("pAUC",
                           "p-value for the AUC (for the Wilcoxon test W statistic)",
                        round(e@pauc), 3)
        MLEval[5,] <- c("Cor",
                           "Correlation coefficient",
                        round(e@cor[[1]]), 3)
        MLEval[6,] <- c("cor",
                           "p-value for correlation coefficient",
                        round(e@pcor), 3)
        # Steph Brodie's TSS which produces the same result as Allouche
        # -1 just makes the output range is 0:1 instead of 1:2 I think.
        # If so this means Sensitivity is e@TPR[which.max(e@TPR + e@TNR)], which doesn't include
        # (e@TPR + e@FNR) but it's a vector of 1s so is redundant. Same for Specificity
        MLEval[7,] <- c("TSS",
                           "True Skill Statistic",
                        round(max(e@TPR + e@TNR - 1)), 3)
        # sensitivity: TP/(TP+FN)
        MLEval[8,] <- c("Sensitivity",
                           "Sensitivity",
                        round(e@TPR[which.max(e@TPR + e@TNR)]), 3)
        # specificity: TN/(FP+TN)
        Specificity <- e@TNR[which.max(e@TPR + e@TNR)]
        MLEval[9,] <- c("Specificity",
                           "Specificity",
                        round(Specificity), 3)
        # Accuracy: TP+TN / TP+TN+FP+FN true false positive negative.
        # TP+TN is just TSS + 1, TP+TN+FP+FN #Sums to 2, redundant
        MLEval[10,] <- c("Accuracy",
                            "Accuracy",
                         round((e@TPR[which.max(e@TPR + e@TNR)] + e@TNR[which.max(e@TPR + e@TNR)]) / 2), 3)
        # Precision: TP/TP+FP. Ignores true negatives. “X% of the predictions are right”
        Precision <- e@TPR[which.max(e@TPR + e@TNR)] / (e@TPR[which.max(e@TPR + e@TNR)] + e@FPR[which.max(e@TPR + e@TNR)])
        MLEval[11,] <- c("Precision",
                            "X% of the predictions are right",
                         round(Precision), 3)
        # Recall: TP/TP+FN: “Y% of actually existing things are captured”.
        Recall <- e@TPR[which.max(e@TPR + e@TNR)] / (e@TPR[which.max(e@TPR + e@TNR)] + e@FNR[which.max(e@TPR + e@TNR)])
        MLEval[12,] <- c("Recall",
                            "Y% of actually existing things are captured",
                         round(Recall), 3)
        # https://www.corvil.com/kb/what-is-a-false-positive-rate
        # Allouche et al 2006:
        # overall accuracy: (TP+TN)/n
        # this seems like a weird metric since the numerator is 0:2 or 1:2 and the divisor could be tiny or huge
        MLEval[13,] <- c("OverallAccuracy",
                            "Overall Accuracy",
                         round((e@TPR[which.max(e@TPR + e@TNR)] + e@TNR[which.max(e@TPR + e@TNR)])/nrow(samples)), 3)
        # Balanced Accuracy, (Recall + Specificity) / 2
        MLEval[14,] <- c("BalancedAccuracy",
                            "Balanced Accuracy",
                         round((Recall + Specificity) / 2), 3)
        # Number of samples. Useful to include in the list
        MLEval[15,] <- c("nSamples",
                            "Number of samples",
                         round(nrow(samples)), 3)
        # Balance: precision vs recall curve. Workhorses.
        # PxR/P+R = F score (P+R = harmonic mean).
        # F1 score: P & R are equally rated. This is the most common one. F1 score importance depends on the project.
        MLEval[16,] <- c("F1score",
                            "P & R equally rated, score importance depends on project",
                         round(2 * ((Precision * Recall) / (Precision + Recall))), 3)
        # F2 score: weighted average of Precision & Recall
        MLEval[17,] <- c("F2score",
                            "weighted average of P & R",
                         round(5 * ((Precision * Recall) / (4 * Precision + Recall))), 3)
        # Threshold which produces the best combo of TPR & TNR
        # t: vector of thresholds used to compute confusion matrices
        MLEval[18,] <- c("Threshold",
                            "Threshold which produced best combo of TPR & TNR",
                         round(e@t[which.max(e@TPR + e@TNR)]), 3)
        # e@prevalence: Prevalence
        MLEval[19,] <- c("Prevalence",
                            "Prevalence",
                         round(e@prevalence[which.max(e@TPR + e@TNR)]), 3)
        # e@ODP: Overall diagnostic power
        MLEval[20,] <- c("ODP",
                            "Overall diagnostic power",
                         round(e@ODP[which.max(e@TPR + e@TNR)]), 3)
        # e@CCR: Correct classification rate
        MLEval[21,] <- c("CCR",
                            "Correct classification rate",
                         round(e@CCR[which.max(e@TPR + e@TNR)]), 3)
        # e@TPR: True positive rate
        MLEval[22,] <- c("TPR",
                            "True positive rate",
                         round(e@TPR[which.max(e@TPR + e@TNR)]), 3)
        # e@TNR: True negative rate
        MLEval[23,] <- c("TNR",
                            "True negative rate",
                         round(e@TNR[which.max(e@TPR + e@TNR)]), 3)
        # e@FPR: False positive rate
        MLEval[24,] <- c("FPR",
                            "False positive rate",
                         round(e@FPR[which.max(e@TPR + e@TNR)]), 3)
        # e@FNR: False negative rate
        MLEval[25,] <- c("FNR",
                            "False negative rate",
                         round(e@FNR[which.max(e@TPR + e@TNR)]), 3)
        # e@PPP: Positive predictive power
        MLEval[26,] <- c("PPP",
                            "Positive predictive power",
                         round(e@PPP[which.max(e@TPR + e@TNR)]), 3)
        # e@NPP: Negative predictive power
        MLEval[27,] <- c("NPP",
                            "Negative predictive power",
                         round(e@NPP[which.max(e@TPR + e@TNR)]), 3)
        # e@MCR: Misclassification rate
        MLEval[28,] <- c("MCR",
                            "Misclassification rate",
                         round(e@MCR[which.max(e@TPR + e@TNR)]), 3)
        # e@OR: Odds-ratio
        MLEval[29,] <- c("OR",
                            "Odds-ratio",
                         round(e@OR[which.max(e@TPR + e@TNR)]), 3)
        # e@kappa: Cohen's kappa
        MLEval[30,] <- c("kappa",
                            "Cohen's kappa",
                         round(e@kappa[which.max(e@TPR + e@TNR)]), 3)
        # dev from calc.deviance from dismo
        MLEval[31,] <- c("dev",
                         "deviance from 2 vecs, obs & pred vals",
                         round(dev), 3)
      } # close if any bernoulli at top of TSS section

      write.csv(MLEval, row.names = FALSE, na = "", file = paste0("./", names(samples[i]), "/MLEvalMetrics.csv"))

      evalmetrics <- c("ROC", "kappa", "prevalence", "TPR", "TNR", "FPR", "FNR", "CCR", "PPP", "NPP", "MCR", "OR")
      for (s in evalmetrics) {
        png(filename = paste0("./",names(samples[i]),"/Eval_", s, ".png"))
        plot(e, s)
        dev.off()
      }
      # can do calc.deviance for gaus also, ditto poisson
      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXX     Evaluation Metrics Processed     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close if MLEvaluate

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Report CSV written      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
    } # close loadgbm isnull

    #avoid sections 19-25 if not predicting to grids
    if (!is.null(grids)) {

      # Load model objects if loadgbm set
      if (!is.null(loadgbm)) {
        if (ZI) {  # don't do if ZI=FALSE
          load(paste0(loadgbm,"Bin_Best_Model"))
          Bin_Best_Model <- "Bin_Best_Model_Object"
        } # close ZI if
        if (gaus) {
          load(paste0(loadgbm,"Gaus_Best_Model"))
          Gaus_Best_Model <- "Gaus_Best_Model_Object"
        } # close gaus if
        dir.create(names(samples[i])) # create resvar-named directory for outputs
      } # close loadgbm optional

      ####19. Binomial predictions####
      if (ZI) {  # don't do if ZI=FALSE
        gbm.predict.grids(get(Bin_Best_Model), grids, want.grids = F, sp.name = "Bin_Preds") #with want.grids=F this is just predict.gbm
        grids$Bin_Preds <- Bin_Preds
      } # close ZI

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Binomial predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####20. Gaussian predictions####
      if (gaus) gbm.predict.grids(get(Gaus_Best_Model), grids, want.grids = F, sp.name = "Gaus_Preds")
      if (gaus) {
        if (ZI) {
          grids$Gaus_Preds <- Gaus_Preds

          if (alerts) beep(2) # progress printer, right aligned for visibility
          print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Gaussian predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

          ####21. Backtransform logged Gaus to unlogged####
          if (gaus) grids$Gaus_Preds_Unlog <- exp(Gaus_Preds + 1/2 * sd(get(Gaus_Best_Model)$residuals, na.rm = FALSE) ^ 2)

          ####22. BIN*positive abundance = final abundance####
          grids$PredAbund <- grids$Gaus_Preds_Unlog * grids$Bin_Preds
        } else { # close gaus yes zi yes run gaus yes zi no
          grids$PredAbund <- Gaus_Preds} #if ZI=TRUE, unlog gaus & multiply by bin. Else just use gaus preds.
      } else grids$PredAbund <- grids$Bin_Preds # if only doing Bin, preds are just bin preds

      predabund <- which(colnames(grids) == "PredAbund") # predicted abundance column number for writecsv

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Final abundance calculated  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####23. Final saves####
      # CSV of Predicted values at each site inc predictor variables' values.
      write.csv(grids, row.names = FALSE, file = paste0("./", names(samples[i]), "/Abundance_Preds_All.csv"))
      # CSV of Predicted values at each site without predictor variables' values.
      # coerce character gridslat/lon into numeric since predabund is given as numeric & you can't mix
      if (is.character(gridslat)) gridslat <- which(colnames(samples) == gridslat)
      if (is.character(gridslon)) gridslon <- which(colnames(samples) == gridslon)
      write.csv(grids[c(gridslat,gridslon,predabund)], row.names = FALSE, file = paste0("./", names(samples[i]), "/Abundance_Preds_only.csv"))
      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Output CSVs written     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####24. Unrepresentativeness surface builder####
      # builds doesn't plot surface. If built, plotted by map maker.
      if (RSB) {
        rsbdf_bin <- gbm.rsb(samples, grids, expvarnames, gridslat, gridslon)
        pos_samples <- subset(samples, brv > 0)
        if (gaus) {
          rsbdf_gaus <- gbm.rsb(pos_samples, grids, expvarnames, gridslat, gridslon)
          rsbdf_both <- data.frame(rsbdf_bin, "Unrep_Gaus" = rsbdf_gaus[,"Unrepresentativeness"], "Unrep_Both" = (rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"]))
          write.csv(rsbdf_both, row.names = FALSE, file = paste0("./", names(samples[i]), "/RSB.csv"))
        } else write.csv(rsbdf_bin, row.names = FALSE, file = paste0("./", names(samples[i]), "/RSB.csv")) # if not gaus
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       RSB CSV written       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close if rsbdf

      ####25. Map maker####
      if (map == TRUE) {   # generate output image & set parameters
        png(filename = paste0("./",names(samples[i]),"/PredAbundMap_",names(samples[i]),".png"),
            width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
        # run gbm.map function with generated parameters
        gbm.map(x = grids[,gridslon],
                y = grids[,gridslat],
                z = grids[,predabund],
                species = names(samples[i]),
                shape = shape, #either autogenerated or set by user so never blank
                ...)  # allows gbm.auto's optional terms to be passed to subfunctions:
        # byx, byy, mapmain, heatcol, mapback, landcol, lejback, legendloc, grdfun, zero, quantile, heatcolours, colournumber
        dev.off()

        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Reticulating splines     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour map generated     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

        if (BnW) { # if BnW=TRUE, run again in black & white for journal submission
          png(filename = paste0("./",names(samples[i]),"/PredAbundMap_BnW_",names(samples[i]),".png"),
              width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
          par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
          gbm.map(x = grids[,gridslon],
                  y = grids[,gridslat],
                  z = grids[,predabund],
                  species = names(samples[i]),
                  shape = shape, #either autogenerated or set by user so never blank
                  landcol = grey.colors(1, start = 0.8, end = 0.8), #light grey. 0=black 1=white
                  mapback = "white",
                  heatcolours = grey.colors(8, start = 1, end = 0))
          dev.off()
          if (alerts) beep(2)  # progress printer, right aligned for visibility
          print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Black & white map generated XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        } # close & save plotting device & close BnW optional

        if (RSB) { # if RSB called, plot that surface separately
          linear01seq <- seq(from = 0, to = 1, length.out = 9) #linear sequence from 0:1, 9 bins
          exp01seq <- expm1(4*linear01seq)/expm1(4) # exponentiate to change shape then scale back to 1

          png(filename = paste0("./",names(samples[i]),"/RSB_Map_Bin_",names(samples[i]),".png"),
              width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
          par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
          gbm.map(x = grids[,gridslon],
                  y = grids[,gridslat],
                  z = rsbdf_bin[,"Unrepresentativeness"],
                  mapmain = "Unrepresentativeness: ",
                  species = names(samples[i]),
                  legendtitle = "UnRep 0-1",
                  shape = shape, #either autogenerated or set by user so never blank
                  # breaks = expm1(breaks.grid(log(2000), ncol = 8, zero = TRUE))/2000) #old failing breaks
                  breaks = exp01seq)
          dev.off() #high value log breaks mean first ~5 values cluster near 0 for high
          # res there, but high values captures in the last few bins.

          if (alerts) beep(2) # progress printer, right aligned for visibility
          print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Colour RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

          if (gaus) {
            png(filename = paste0("./",names(samples[i]),"/RSB_Map_Gaus_",names(samples[i]),".png"),
                width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
            gbm.map(x = grids[,gridslon],
                    y = grids[,gridslat],
                    z = rsbdf_gaus[,"Unrepresentativeness"],
                    mapmain = "Unrepresentativeness: ",
                    species = names(samples[i]),
                    legendtitle = "UnRep 0-1",
                    shape = shape, #either autogenerated or set by user so never blank
                    breaks = exp01seq)
            dev.off()

            if (alerts) beep(2) # progress printer, right aligned for visibility
            print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

            png(filename = paste0("./",names(samples[i]),"/RSB_Map_Both_",names(samples[i]),".png"),
                width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
            gbm.map(x = grids[,gridslon],
                    y = grids[,gridslat],
                    z = rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"],
                    mapmain = "Unrepresentativeness: ",
                    species = names(samples[i]),
                    legendtitle = "UnRep 0-2",
                    shape = shape, #either autogenerated or set by user so never blank
                    breaks = exp01seq)
            dev.off()

            if (alerts) beep(2) # progress printer, right aligned for visibility
            print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))}

          if (BnW) {     # if BnW=TRUE, do again for b&w
            png(filename = paste0("./",names(samples[i]),"/RSB_Map_BnW_Bin_",names(samples[i]),".png"),
                width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
            gbm.map(x = grids[,gridslon],
                    y = grids[,gridslat],
                    z = rsbdf_bin[,"Unrepresentativeness"],
                    mapmain = "Unrepresentativeness: ",
                    mapback = "white",
                    species = names(samples[i]),
                    heatcolours = grey.colors(8, start = 1, end = 0), #default 8 greys
                    ####BUG:setting heatcolours & colournumber overrides this####
                    landcol = grey.colors(1, start = 0.8, end = 0.8), #light grey. 0=black 1=white
                    legendtitle = "UnRep 0-1",
                    shape = shape, #either autogenerated or set by user so never blank
                    breaks = exp01seq)
            dev.off()

            if (alerts) beep(2) # progress printer, right aligned for visibility
            print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     B&W RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

            if (gaus) {
              png(filename = paste0("./",names(samples[i]),"/RSB_Map_BnW_Gaus_",names(samples[i]),".png"),
                  width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
              par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
              gbm.map(x = grids[,gridslon],
                      y = grids[,gridslat],
                      z = rsbdf_gaus[,"Unrepresentativeness"],
                      mapmain = "Unrepresentativeness: ",
                      mapback = "white",
                      species = names(samples[i]),
                      heatcolours = grey.colors(8, start = 1, end = 0),
                      landcol = grey.colors(1, start = 0.8, end = 0.8),
                      legendtitle = "UnRep 0-1",
                      shape = shape, #either autogenerated or set by user so never blank
                      breaks = exp01seq)
              dev.off()

              if (alerts) beep(2) # progress printer, right aligned for visibility
              print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

              png(filename = paste0("./",names(samples[i]),"/RSB_Map_BnW_Both_",names(samples[i]),".png"),
                  width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
              par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
              gbm.map(x = grids[,gridslon],
                      y = grids[,gridslat],
                      z = rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"],
                      mapmain = "Unrepresentativeness: ",
                      mapback = "white",
                      species = names(samples[i]),
                      heatcolours = grey.colors(8, start = 1, end = 0),
                      landcol = grey.colors(1, start = 0.8, end = 0.8),
                      legendtitle = "UnRep 0-2",
                      shape = shape, #either autogenerated or set by user so never blank
                      breaks = exp01seq)
              dev.off()
              if (alerts) beep(2) # progress printer, right aligned for visibility
              print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))}
          } # close BnW RSBs
        } # close RSB mapper
      } # close Map Maker
    } #close grids option from above section 19
    print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Grids/maps/everything done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
  } # close response variable (resvar) loop
  gc() # Force R to release memory it is no longer using
  options(error = NULL) # reset error options to default
  if (alerts) beep(8)} # final user notification, then close the function
####END####
