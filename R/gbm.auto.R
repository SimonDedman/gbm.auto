#' Automated Boosted Regression Tree modelling and mapping suite
#'
#' Automates delta log normal boosted regression trees abundance prediction.
#' Loops through all permutations of parameters provided (learning
#' rate, tree complexity, bag fraction), chooses the best, then simplifies it.
#' Generates line, dot and bar plots, and outputs these and the predictions
#' and a report of all variables used, statistics for tests, variable
#' interactions, predictors used and dropped, etc. If selected, generates
#' predicted abundance maps, and Unrepresentativeness surfaces.
#' See www.GitHub.com/SimonDedman/gbm.auto for issues, feedback, and development
#' suggestions. See SimonDedman.com for links to walkthrough paper, and papers
#' and thesis published using this package.
#'
#' @param grids Explanatory data to predict to. Import with (e.g.) read.csv and
#' specify object name. Defaults to NULL (won't predict to grids).
#' @param samples Explanatory and response variables to predict from. Keep col
#' names short (~17 characters max), no odd characters, spaces, starting
#' numerals or terminal periods. Spaces may be converted to periods in directory
#' names, underscores won't. Can be a subset of a large dataset.
#' @param expvar Vector of names or column numbers of explanatory variables in
#' 'samples': c(1,3,6) or c("Temp","Sal"). No default.
#' @param resvar Name or column number(s) of response variable in samples: 12,
#' c(1,4), "Rockfish". No default. Column name is ideally species name.
#' @param randomvar Add a random variable (uniform distribution, 0-1) to the expvars, to see whether
#'  other expvars perform better or worse than random.
#' @param tc Permutations of tree complexity allowed, can be vector with
#' the largest sized number no larger than the number of explanatory variables
#' e.g. c(2,7), or a list of 2 single numbers or vectors, the first to be passed
#' to the binary BRT, the second to the Gaussian, e.g. tc = list(c(2,6), 2) or
#' list(6, c(2,6)).
#' @param lr Permutations of learning rate allowed. Can be a vector or a list of
#' 2 single numbers or vectors, the first to be passed to the binary BRT, the
#' second to the Gaussian, e.g. lr = list(c(0.01,0.02),0.0001) or
#' list(0.01,c(0.001, 0.0005)).
#' @param bf Permutations of bag fraction allowed, can be single number, vector
#' or list, per tc and lr. Defaults to 0.5.
#' @param offset Column number or quoted name in samples, containing offset values relating to the
#' samples. A numeric vector of length equal to the number of cases. Similar to weighting, see
#' https://towardsdatascience.com/offsetting-the-model-logic-to-implementation-7e333bc25798 .
#' @param n.trees From gbm.step, number of initial trees to fit. Can be
#' single or list but not vector i.e. list(fam1,fam2).
#' @param ZI Are data zero-inflated? TRUE FALSE "CHECK". Choose one. TRUE:
#' delta BRT, log-normalised Gaus, reverse log-norm and bias corrected. FALSE:
#' do Gaussian only, no log-normalisation. "CHECK": Tests data for you. Default is
#'  "CHECK". TRUE and FALSE aren't in quotes, "CHECK" is.
#' @param fam1 Probability distribution family for 1st part of delta process,
#' defaults to "bernoulli". Choose one.
#' @param fam2 Probability distribution family for 2nd part of delta process,
#' defaults to "gaussian". Choose one.
#' @param simp Try simplifying best BRTs?
#' @param gridslat Column number for latitude in 'grids'.
#' @param gridslon Column number for longitude in 'grids'.
#' @param samplesGridsAreaScaleFactor Scale up or down factor so values in the predict-to pixels of
#' 'grids' match the spatial scale sampled by rows in 'samples'. Default 1 means no change.
#' @param multiplot Create matrix plot of all line files? Default true.
#' turn off if big n of exp vars causes an error due to margin size problems.
#' @param cols Barplot colour vector. Assignment in order of explanatory
#' variables. Default 1*white: white bars black borders. '1*' repeats.
#' @param linesfiles Save individual line plots' data as csv's? Default TRUE.
#' @param smooth Apply a smoother to the line plots? Default FALSE.
#' @param savedir Save outputs to a temporary directory (default) else change to
#'  current directory e.g. "/home/me/folder". Do not use getwd() here.
#' @param savegbm Save gbm objects and make available in environment after
#' running? Open with load("Bin_Best_Model") Default TRUE.
#' @param loadgbm Relative or (very much preferably) absolute location of folder containing
#' Bin_Best_Model and Gaus_Best_Model. If set will skip BRT calculations and do
#' predicted maps and csvs. Avoids re-running BRT models again (the slow bit),
#' can run normally once with savegbm=T then multiple times with new grids &
#' loadgbm to predict to multiple grids e.g. different seasons, areas, etc.
#' Default NULL, character vector, "./" for working directory.
#' @param varint Calculate variable interactions? Default:TRUE, FALSE for error:
#' "contrasts can be applied only to factors with 2 or more levels".
#' @param map Save abundance map png files?
#' @param shape Set coast shapefile, else bounds calculated by gbm.map which
#' then calls gbm.basemap to download and auto-generate the base map. Read in
#' existing files by installing the shapefiles package then
#' DesiredMapName <- read.shapefile("ShapeFileName")
#' omitting the .shp extension.
#' @param RSB Run Unrepresentativeness surface builder? Default TRUE.
#' @param BnW Repeat maps in black and white e.g. for print journals. Default
#' TRUE.
#' @param alerts Play sounds to mark progress steps. Default TRUE but running
#' multiple small BRTs in a row (e.g. gbm.loop) can cause RStudio to crash.
#' @param pngtype Filetype for png files, alternatively try "quartz" on Mac.
#' Choose one.
#' @param gaus Do family2 (typically Gaussian) runs as well as family1
#' (typically Bin)? Default TRUE.
#' @param MLEvaluate do machine learning evaluation metrics & plots? Default
#' TRUE.
#' @param brv Dummy param for package testing for CRAN, ignore.
#' @param grv Dummy param for package testing for CRAN, ignore.
#' @param Bin_Preds Dummy param for package testing for CRAN, ignore.
#' @param Gaus_Preds Dummy param for package testing for CRAN, ignore.
#' @param ... Optional arguments for gbm.step (dismo package) arguments n.trees and
#' max.trees, both of which can be added in list(1,2) format to pass to fam1 and
#'  2; for gbm.mapsf colourscale, heatcolours, colournumber, and others.
#'
#' @return Line, dot and bar plots, a report of all variables used, statistics
#' for tests, variable interactions, predictors used and dropped, etc. If
#' selected, generates predicted abundance maps, and Unrepresentativeness surface. Biggest
#' Interactions in the report csv: see ?dismo::gbm.interactions .
#'
#' @details Errors and their origins:
#'
#' 1. install ERROR: dependencies ‘rgdal’, ‘rgeos’ are not available for package ‘gbm.auto’.
#' For Linux/*buntu systems, in terminal, type: 'sudo apt install libgeos-dev', 'sudo apt install
#' libproj-dev', 'sudo apt install libgdal-dev'.
#'
#' 2. Error in FUN(X\[\[i\]\], ...) : only defined on a data frame with all numeric
#' variables. Check your variable types are correct, e.g. numerics haven't been imported
#' as factors because there's an errant first row of text information before the
#' data. Remove NA rows from the response variable if present: convert blank
#' cells to NA on import with read.csv(x, na.strings = "") then
#' samples2 <- samples1\[-which(is.na(samples\[,resvar_column_number\])),\]
#'
#' 3. At bf=0.5, if nrows <= 42 gbm.step will crash. Use gbm.bfcheck to determine optimal viable bf
#' size
#'
#' 4. Maps/plots don't work/output. If on a Mac, try changing pngtype to "quartz".
#'
#' 5. Error in while (delta.deviance > tolerance.test AMPERSAND n.fitted <
#' max.trees): missing value where TRUE/FALSE needed. If running a zero-inflated delta model
#' (bernoulli/bin & gaussian/gaus), Data are expected to contain zeroes (lots of them in zero-
#' inflated cases), have you already filtered them out, i.e. are only testing the positive cases?
#' Or do you only have positive cases? If so only run (e.g.) Gaussian: set ZI to FALSE.
#'
#' 6. Error in round(gbm.object$cv.statistics$deviance.mean, 4) : non-numeric argument to
#' mathematical function. LR or BF probably too low in earlier BRT (normally Gaus run with highest
#' TC)
#'
#' 7. Error in if (n.trees > x$n.trees) { : argument is of length zero}. LR or BF probably too low
#' in earlier BRT (normally Gaus run with highest TC).
#'
#' 8. Error in gbm.fit(x, y, offset = offset, distribution = distribution, w = w) The dataset size
#' is too small or subsampling rate is too large: nTrain*bag.fraction <= n.minobsinnode. LR or BF
#' probably too low in earlier BRT (normally Gaus run with highest TC). It may be that you don't
#' have enough positive samples to run BRT modelling. Run gbm.bfcheck to check recommended minimum
#' BF size.
#'
#' 9. Warning message: In cor(y_i, u_i) : the standard deviation is zero. LR or BF probably too low
#' in earlier BRT (normally Gaus run with highest TC). It may be that you don't have enough positive
#'  samples to run BRT modelling. Run gbm.bfcheck to check recommended minimum BF size. Similarly:
#'  glm.fit: fitted probabilities numerically 0 or 1 occurred, and glm.fit: algorithm did not
#'  converge. Similarly: Error in if (get(paste0("Gaus_BRT", ".tc", j, ".lr", k, ".bf",
#'  l))$self.statistics$correlation\[\[1\]\]: argument is of length zero. See also: Error 15.
#'
#' 10. Anomalous values can obfuscate clarity in line plots e.g. salinity range 32:35ppm but dataset
#'  has errant 0 value: plot axis will be 0:35, and 99.99% of the data will be in the tiny bit at
#'  the right. Clean your data beforehand.
#'
#' 11. Error in plot.new() : figure margins too large: In RStudio, adjust plot frame (usually bottom
#'  right) to increase its size. Still fails? Set multiplot=FALSE.
#'
#' 12. Error in dev.print(file = paste0("./", names(samples\[i\]), "/pred_dev_bin.jpeg"): can only
#' print from a screen device. An earlier failed run (e.g. LR/BF too low) left a plotting device
#' open. Close it with: 'dev.off()'.
#'
#' 13. RStudio crashed: set alerts=F and pause cloud sync programs if outputting to a synced folder.
#'
#' 14. Error in grDevices::dev.copy(device = function (filename = "Rplot%03d.jpeg", could not open
#' file './P_PECTINATA..../pred_dev_bin.jpeg' (or similar). Your resvar column name contains an
#' illegal character e.g. /&'_. Fix with colnames(samples)[n] <- "BetterName".
#'
#' 15. Error in gbm.fit: Poisson requires the response to be a positive integer. If running Poisson
#' distributions, ensure the response variables are positive integers, but if they are, try a
#' smaller learning rate.
#'
#' 16. If lineplots of factorial variables include empty columns be sure to remove unused levels
#' with samples %<>% droplevels() before the gbm.auto run
#'
#' 17. Error in seq.default(from = min(x$var.levels\[\[i.var\[i\]\]\]), to = max(x$var.levels\[\[i.var\[i\]\]\])
#' :'from' must be a finite number. If you logged any expvars with log() and they has zeroes in them
#' , those zeroes became imaginary numbers. Use log1p() instead.
#'
#' 18. Error in loadNamespace...'dismo' 1.3-9 is being loaded, but >= 1.3.10 is required: first do
#' remotes::install_github("rspatial/dismo") then library(dismo).
#'
#' ALSO: check this section in the other functions run by gbm.auto e.g. gbm.map, gbm.basemap. Use
#' traceback() to find the source of errors.
#'
#' I strongly recommend that you download papers 1 to 5 (or just the doctoral thesis) on
#' <http://www.simondedman.com>, with emphasis on P4 (the guide) and P1 (statistical background).
#' Elith et al 2008 (<http://refhub.elsevier.com/S0304-3800(15)00207-0/sbref0085>) is also strongly
#' recommended.
#' Just because you CAN try every conceivable combination of tc, lr, bf, all, at once doesn't mean
#' you should. Try a range of lr in shrinking orders of magnitude from 0.1 to 0.000001, find the
#' best, THEN try tc c(2, n.expvars), find the best THEN bf c(0.5, 0.75, 0.9) and then in between if
#'  either outperform 0.5.
#'
#' @examples
#' \donttest{
#' # Not run. Note: grids file was heavily cropped for CRAN upload so output map
#' # predictions only cover patchy chunks of the Irish Sea, not the whole area.
#' # Full versions of these files:
#' # https://drive.google.com/file/d/1WHYpftP3roozVKwi_R_IpW7tlZIhZA7r
#' # /view?usp=sharing
#' library(gbm.auto)
#' data(grids)
#' data(samples)
#' # Set your working directory
#' gbm.auto(grids = grids, samples = samples, expvar = c(4:8, 10), resvar = 11,
#' tc = c(2,7), lr = c(0.005, 0.001), ZI = TRUE, savegbm = FALSE)}
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
#' @export
#' @import dismo
#' @importFrom beepr beep
#' @importFrom dplyr across mutate
#' @importFrom gbm plot.gbm
#' @importFrom grDevices dev.off dev.print graphics.off grey.colors jpeg png
#' @importFrom graphics axis barplot image legend lines mtext par text
#' @importFrom stats sd runif
#' @importFrom utils packageVersion read.csv write.csv
#' @importFrom stringi stri_split_fixed
#'
gbm.auto <- function(
    grids = NULL,         # explanatory data to predict to. Import with (e.g.)
    # read.csv and specify object name. Defaults to NULL (won't predict to grids)
    samples,  # explanatory and response variables to predict from.
    # Keep col names short, no odd characters, starting numerals or terminal periods
    # Spaces may be converted to periods in directory names, underscores won't.
    # Can be a subset.
    expvar,               # list of column numbers of explanatory variables in
    # 'samples', expected e.g. c(1,35,67,etc.). No default.
    resvar,               # column number(s) of response variable (e.g. CPUE) in
    # samples, e.g. 12 or c(4,5,6). No default. Column name should be species name.
    randomvar = FALSE,    # Add a random variable (uniform distribution, 0-1) to the expvars, to see
    # whether other expvars perform better or worse than random.
    tc = c(2),            # permutations of tree complexity allowed, can be a
    # vector with the largest sized number no larger than the number of
    # explanatory variables e.g. c(2,7), or a list of 2 single numbers or vectors,
    # the first to be passed to the binary BRT, the second to the Gaussian, e.g.
    # tc = list(c(2,6), 2) or list(6, c(2,6)).
    lr = c(0.01, 0.005),   # permutations of learning rate allowed. Can be a
    # vector or a list of 2 single numbers or vectors, the first to be passed to
    # the binary BRT, the second to the Gaussian, e.g.
    # lr = list(c(0.01,0.02),0.0001) or list(0.01,c(0.001, 0.0005)).
    bf = 0.5,             # permutations of bag fraction allowed, can be single
    # number, vector or list, per tc and lr.
    offset = NULL,        # column number or quoted name in samples, containing offset values
    # relating to the samples. A numeric vector of length equal to the number of cases. Similar to
    # weighting, see
    # https://towardsdatascience.com/offsetting-the-model-logic-to-implementation-7e333bc25798
    n.trees = 50,         # from gbm.step, number of initial trees to fit. Can be
    # single or list but not vector i.e. list(fam1, fam2).
    ZI = "CHECK", # are data zero-inflated? "CHECK"/FALSE/TRUE.
    # TRUE: delta BRT, log-normalised Gaus, reverse log-norm and bias corrected.
    # FALSE: do Gaussian only, no log-normalisation.
    # CHECK: Tests data for you. Default is TRUE.
    fam1 = c("bernoulli", "binomial", "poisson", "laplace", "gaussian"),
    # probability distribution family for 1st part of delta process, defaults to
    # "bernoulli",
    fam2 = c("gaussian", "bernoulli", "binomial", "poisson", "laplace"),
    # probability distribution family for 2nd part of delta process, defaults to
    # "gaussian",
    simp = TRUE,          # try simplifying best BRTs?
    gridslat = 2,         # column number for latitude in 'grids'
    gridslon = 1,         # column number for longitude in 'grids'
    samplesGridsAreaScaleFactor = 1, # Scale up or down factor so values in the predict-to pixels of
    # 'grids' match the spatial scale sampled by rows in 'samples'. Default 1 means no change.
    multiplot = TRUE,     # create matrix plot of all line files? Default true
    # turn off if large number of expvars causes an error due to margin size problems.
    cols = grey.colors(1,1,1), # bar-plot colour vector. Assignment in order of
    # explanatory variables. Default 1*white: white bars black borders. '1*' repeats
    linesfiles = TRUE,    # save individual line plots' data as CSVs?
    smooth = FALSE,       # apply a smoother to the line plots? Default FALSE
    savedir = tempdir(),  # save outputs to a temporary directory (default) else
    # change to current directory e.g. "/home/me/folder". Do not use getwd() here.
    savegbm = TRUE,       # save gbm objects and make available in environment after running? Open with load("Bin_Best_Model")
    loadgbm = NULL,       # relative or absolute location of folder containing
    # Bin_Best_Model and Gaus_Best_Model. If set will skip BRT calculations and do
    # predicted maps and CSVs. Default NULL, character vector, "./" for working directory
    varint = TRUE,        # calculate variable interactions? Default:TRUE, FALSE
    # for error "contrasts can be applied only to factors with 2 or more levels"
    map = TRUE,           # save abundance map png files?
    shape = NULL,         # set coast shapefile, else bounds calculated by gbm.map
    # which then calls gbm.basemap to download and auto-generate the base map.
    RSB = TRUE,           # run Unrepresentativeness surface builder?
    BnW = TRUE,           # repeat maps in black and white e.g. for print journals
    alerts = TRUE,        # play sounds to mark progress steps. Running many small
    # BRTs e.g. gbm.loop can cause RStudio to crash, if so set this to FALSE
    pngtype = c("cairo-png", "quartz", "Xlib"), # file-type for png files,
    # alternatively try "quartz" on Mac
    gaus = TRUE,          # do fam2 (typically Gaussian) runs as well as Bin? Default TRUE.
    MLEvaluate = TRUE,    # do machine learning evaluation metrics & plots? Default TRUE
    brv = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
    grv = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
    Bin_Preds = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
    Gaus_Preds = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
    ...)                  # Optional arguments for zero in breaks.grid in gbm.map,
# legend in legend.grid in gbm.map, mapmain in gbm.map
# (default = "Predicted CPUE (numbers per hour): ") and gbm.step (dismo package)
# arguments max.trees and others.
{
  # Generalised Boosting Model / Boosted Regression Tree process chain automater.
  # Simon Dedman, 2012-6 simondedman@gmail.com GitHub.com/SimonDedman/gbm.auto

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
  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldwd <- getwd()
  oldoptions <- options()
  on.exit(dev.off()) # close any open graphics devices to avoid issues later
  on.exit(par(oldpar))
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(options(oldoptions), add = TRUE)
  setwd(savedir)
  if (alerts) options(error = function() {
    beep(9)# give warning noise if it fails
    graphics.off()# kill all graphics devices
    setwd(oldwd) # reinstate original working directory. Probably redundant given on.exit
  } # close options subcurly
  ) # close options

  # @import utils # put at top
  # utils::globalVariables("where") # https://github.com/r-lib/tidyselect/issues/201#issuecomment-650547846
  # https://stackoverflow.com/questions/40251801/how-to-use-utilsglobalvariables
  # presence of list columns, even if not used, will break the write.table within write.csv for abundance prediction saving
  # if (any(as.data.frame(unlist(lapply(samples, class)))[,1] == "list")) {
  #   samples <- samples |> mutate(across(.cols = where(is.list), ~ sapply(.x, toString)))
  #   print("list columns converted to character columns in samples")
  # }
  # if (any(as.data.frame(unlist(lapply(grids, class)))[,1] == "list")) {
  #   grids <- grids |> mutate(across(.cols = where(is.list), ~ sapply(.x, toString)))
  #   print("list columns converted to character columns in grids")
  # }

  # ToDo: add to existing options(error) if present####
  # options(error = function() {.rs.recordTraceback(TRUE, 5, .rs.enqueueError)})
  # as.character(getOption("error"))  # (function() {.rs.recordTraceback(TRUE, 5, .rs.enqueueError)})()
  # class: call
  # "function () \n{\n    .rs.recordTraceback(TRUE, 5, .rs.enqueueError)\n}"
  # "function () \n{\n    beep(9)\n    graphics.off()\n}"
  # class(options("error")[[1]])
  # https://stackoverflow.com/questions/66003747/r-replace-optionserror-with-existing-contents-if-present-plus-additional

  fam1 <- match.arg(fam1) # populate object from function argument in proper way
  fam2 <- match.arg(fam2)
  pngtype <- match.arg(pngtype)
  if (fam1 == "binomial") fam1 <- "bernoulli" # gbm::gbm doesn't like binomial even though it's the same
  if (fam2 == "binomial") fam2 <- "bernoulli"
  if (gaus) if (fam2 == fam1) stop("attempting to run delta model with both families the same. Expects fam1==bernoulli & gaus==TRUE & fam2==somethingElse, OR fam1==anything & gaus==FALSE")

  # tibble's don't collapse into a vector, instead an X x 1 df, which breaks various functionality.
  if ("tbl" %in% class(grids)) grids <- as.data.frame(grids)
  if ("tbl" %in% class(samples)) samples <- as.data.frame(samples)

  # create basemap using gbm.basemap & these bounds, else basemap will be called for every map
  if (!is.null(grids)) if (map) { # create basemap grids not null, map requested, basemap not provided
    if (is.null(shape)) {
      if (!exists("gbm.basemap")) {stop("you need to install gbm.basemap to run this function")}
      bounds = c(range(grids[,gridslon]),range(grids[,gridslat]))
      #create standard bounds from data, and extra bounds for map aesthetic
      # xmid <- mean(bounds[1:2])
      # ymid <- mean(bounds[3:4])
      # xextramax <- ((bounds[2] - xmid) * 1.6) + xmid
      # xextramin <- xmid - ((xmid - bounds[1]) * 1.6)
      # yextramax <- ((bounds[4] - ymid) * 1.6) + ymid
      # yextramin <- ymid - ((ymid - bounds[3]) * 1.6)
      # extrabounds <- c(xextramin, xextramax, yextramin, yextramax) # identical code to what's in basemap
      shape <- gbm.basemap(bounds = bounds,
                           savedir = savedir,
                           extrabounds = TRUE)
    } # close isnull shape
  } # close isnull grids

  if (randomvar) { # add random variable if requested
    samples$randomvar <- runif(n = nrow(samples), min = 0, max = 1)  # make it then add to expvar & thus expvarnames
    if (is.numeric(expvar)) expvar <- c(expvar, which(colnames(samples) %in% "randomvar")) else expvar <- c(expvar, "randomvar")
  }
  expvarnames <- if (is.numeric(expvar)) names(samples[expvar]) else expvar # list of explanatory variable names
  if (!length(cols) == 1 & !length(cols) == length(expvarnames)) stop("length of cols is neither the same as the length of expvars (plus randomvar if selected) nor 1")
  if (length(cols) == 1) cols <- rep(cols, length(expvarnames)) # if cols is length 1, repeat it so it attaches properly next
  expvarcols <- cbind(cols[1:length(expvarnames)],expvarnames) # assign explanatory variables to colours

  if (!is.null(offset)) {
    if (is.character(offset)) offset <- which(colnames(samples) %in% offset) # if offset is the column name, change to column number
    colnames(samples)[offset] <- "offset" # then change name to "offset"
  }

  if (is.list(tc)) { # if lists entered for tc lr or bf, split them to bin and gaus
    if (length(tc) > 2) {stop("Only 2 tc list items allowed: 1 per family")}
    tcgaus <- tc[[2]]
    tc <- tc[[1]]
  } else {tcgaus <- tc} # else make the gaus object the same as the bin. close if else

  if (is.list(lr)) {
    if (length(lr) > 2) {stop("Only 2 lr list items allowed: 1 per family")}
    lrgaus <- lr[[2]]
    lr <- lr[[1]]
  } else {lrgaus <- lr} # close if else lr

  if (is.list(bf)) {
    if (length(bf) > 2) {stop("Only 2 bf list items allowed: 1 per family")}
    bfgaus <- bf[[2]]
    bf <- bf[[1]]
  } else {bfgaus <- bf} # close if else bf

  if (is.list(n.trees)) { # if list entered n.trees, split to fam1 and fam2
    if (length(n.trees) > 2) {stop("Only 2 n.trees list items allowed: 1 per family")}
    ntf1 <- n.trees[[1]]
    ntf2 <- n.trees[[2]]
  } else {
    ntf1 <- n.trees
    ntf2 <- n.trees} # else make fam1 and fam2 the same. close if else n.trees

  for (i in resvar) { # loop everything for each response variable (e.g. species)
    dir.create(names(samples[i])) # create resvar-named directory for outputs
    m = 0 # Gaus only loop counter to allow best gaus BRT choice
    n = 0 # Print counter for all loops of BRT combos & best bin BRT choice
    if (!is.null(grids)) if (!all(expvarnames %in% names(grids))) stop(print("Expvar column names in samples but missing from grids:"), print(expvarnames[which(!expvarnames %in% names(grids))]))
    if (anyNA(samples[i])) stop("Response variable range contains NA values, please filter out these rows with: mysamples <- mysamples[-which(is.na(mysamples[resvar])),]")
    if (all(samples[i] == 0)) stop("Response variable only contains zeroes")

    ####2. ZI check & log####
    # if user has asked code to check for ZI, check it & set new ZI status
    if (ZI == "CHECK") if (sum(samples[,i] == 0, na.rm = TRUE) / length(samples[,i]) >= 0.5) ZI = TRUE else ZI = FALSE
    # ensure resvar has zeroes (expects mix of successful & unsuccessful samples for bernoulli/binary runs)
    if (!ZI) if (min(samples[i]) > 0) print("No zeroes in response variable. If using a zero inflated model, Method expects unsuccessful, as well as successful, samples")

    # create binary (0/1) response variable, for bernoulli BRTs
    samples$brv <- ifelse(samples[i] > 0, 1, 0)
    brvcol <- which(colnames(samples) == "brv") # brv column number for BRT

    # create logged response variable, for Gaussian BRTs when data are zero-inflated (otherwise just use resvar directly)
    logem <- log1p(samples[,i]) # logs resvar i.e. containing zeroes
    dont  <- samples[,i]
    # log1p fam2 (gaussian response variable grv) resvar if bin only (fam1 bin, fam2 FALSE), OR if resvar is delta & ZI & NOT poisson (which can't be logged, must be positive integers)
    if (fam1 == "bernoulli" & (!gaus | (gaus & ZI & (fam2 != "poisson")))) {samples$grv <- logem} else {samples$grv <- dont}
    grvcol <- which(colnames(samples) == "grv") # grv column number for BRT

    if (ZI) {
      grv_yes <- subset(samples, grv > 0) # nonzero subset for gaussian/poisson BRTs if zero inflated
    } else {
      grv_yes <- samples # use the full dataset if not ZI
    }

    if (is.null(loadgbm)) { #if loadgbm is NULL i.e. you're running BRTs not
      # predicting from existing models. Skip to L1404

      ####3. Begin Report####
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) { # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        reportcolno = 3 + (length(tc)*length(lr)*length(bf)) + (length(tcgaus)*length(lrgaus)*length(bfgaus)) + 14
        # if only 1 permutation, = 19
      } else { # else zi
        reportcolno = 3 + (length(tcgaus)*length(lrgaus)*length(bfgaus)) + 7
        # if only 1 permutation = 11
      } # close if else ZI
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
      } else {
        # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {colnames(Report)[(reportcolno - 13):(reportcolno - 7)] <- c("Best Binary BRT",
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
                                                             "Biggest Interactions (Gaus)")} # close if else gaus
      # populate the final 14 column names
      Report[1:length(expvar),1] <- expvarnames # put expvar names in first column # names(samples[expvar])
      Report[1,2] <- names(samples[i]) # put resvar in col 2
      Report[1,3] <- ZI # ZI in col 3
      Report[2,3] <- paste0(round(sum(samples[,i] == 0, na.rm = TRUE) / length(samples[,i]), 3) * 100, "% zeroes") # add zeroes % under ZI in col3

      StatsObjectsList <- list()

      Bin_Best_Score <- 0 # create blanks for best results to use in loops
      Bin_Best_Model <- 0
      Gaus_Best_Score <- 0
      Gaus_Best_Model <- 0

      # Begin bin loops
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) { # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        for (j in tc) {   # list permutations of tree complexity allowed
          for (k in lr) {   # list permutations of learning rate allowed
            for (l in bf) {   # list permutations of bag fraction allowed
              n <- n + 1   # Add to print counter
              ####4. Binomial BRT####
              print(paste0("Running ", fam1, " BRT, tc=",j,", lr=",k,", bf=",l))
              assign(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l),
                     gbm.step.sd(data = samples,
                                 gbm.x = expvar,
                                 gbm.y = brvcol,
                                 family = fam1,
                                 tree.complexity = j,
                                 learning.rate = k,
                                 bag.fraction = l,
                                 n.trees = ntf1,
                                 {if (!is.null(offset)) offset = grv_yes$offset},
                                 ...)
              )
              if (is.null(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)))) { # test for BRT failure and skip this hyperparameter combo
                Report[1, (3 + n)] <- "Run failed, try smaller lr or step size"
                colnames(Report)[3 + n] <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)
                next
              }
              dev.print(file = paste0("./",names(samples[i]),"/pred_dev_", fam1, "_tc",j,"lr",k,"bf",l,".jpeg"), device = jpeg, width = 600)
              # dev.print(file = paste0("./",names(samples[i]),"/pred_dev_bin.jpeg"), device = jpeg, width = 600)
              print(paste0("Done Bin_BRT",".tc",j,".lr",k,".bf",l))
              print(warnings())
              # assign("last.warning", NULL, envir = baseenv()) # dumps warnings so subsequent printing doesn't reprint the existing warning
              ####5. Select best bin model####
              if (n == 1) { # if this is the first loop, best score & model name is this one by default
                Bin_Best_Score <- get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
                Bin_Best_Model <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)
                # else if this models self.statistics$correlation > the best model, make this the new best model
              }  else if (get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]] > Bin_Best_Score) {
                Bin_Best_Score <- get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
                Bin_Best_Model <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)
              } # close if else n==1

              ####6. Add bin stats to report####
              if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {Report[1:8,(3 + n)] <- c(paste0("trees: ",round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$n.trees, 3)),
                                                                                         paste0("Training Data Correlation: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]], 3)),
                                                                                         paste0("CV Mean Deviance: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.mean, 3)),
                                                                                         paste0("CV Deviance SE: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.se, 3)),
                                                                                         paste0("CV D squared: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$d.squared, 3)),
                                                                                         paste0("CV Mean Correlation: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.mean, 3)),
                                                                                         paste0("CV Correlation SE: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.se, 3)),
                                                                                         paste0("CV RMSE: ", round(get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$cv.rmse), 3))
              # bin BRT name
              colnames(Report)[3 + n] <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l)

              # Add Bin stats objects to StatsObjectsList
              StatsObjectsList[[length(StatsObjectsList) + 1]] <- get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$self.statistics # send to new position after last item
              names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l, "__self.statistics") # name it. new length now includes self.statistics
              StatsObjectsList[[length(StatsObjectsList) + 1]] <- get(paste0("Bin_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics
              names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0("Bin_BRT",".tc",j,".lr",k,".bf",l, "__cv.statistics")
              } # close ZI if

              if (alerts) beep(2) # progress printer, right aligned
              if (gaus) {
                print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tc)*length(lr)*length(bf)) + (length(tcgaus)*length(lrgaus)*length(bfgaus)), "     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
              } else { # close if else gaus
                print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tc)*length(lr)*length(bf)), "     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
              } # else if gaus
            } # close bf l
          } # close lr k
        } # close tc j
      } # close ZI option, making all bin BRT objects & continuing through model selection

      # Begin Gaus loops
      if (gaus) for (j in tcgaus) {   # list permutations of tree complexity allowed
        for (k in lrgaus) {   # list permutations of learning rate allowed
          for (l in bfgaus) {   # list permutations of bag fraction allowed
            n <- n + 1 # Add to print/loop counter for every bin or gaus BRT loop
            m <- m + 1 # Add to loop counter for Gaus best model selection
            ####7. Gaussian BRT####
            print(paste0("Running ", fam2, " BRT, tc=",j,", lr=",k,", bf=",l))
            write.csv(x = grv_yes[,grvcol], file = paste0("./",names(samples[i]),"/grv.csv"), row.names = FALSE)
            assign(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l),
                   gbm.step.sd(data = grv_yes,
                               gbm.x = expvar,
                               gbm.y = grvcol,
                               family = fam2,
                               tree.complexity = j,
                               learning.rate = k,
                               bag.fraction = l,
                               n.trees = ntf2,
                               {if (!is.null(offset)) offset = grv_yes$offset},
                               ...)
            )
            if (is.null(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)))) {
              Report[1, (3 + n)] <- "Run failed, try smaller lr or step size"
              colnames(Report)[3 + n] <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)
              next
            }
            # dev.print(file = paste0("./",names(samples[i]),"/pred_dev_gaus.jpeg"), device = jpeg, width = 600)
            dev.print(file = paste0("./",names(samples[i]),"/pred_dev_", fam2, "_tc",j,"lr",k,"bf",l,".jpeg"), device = jpeg, width = 600)
            print(paste0("Done Gaus_BRT",".tc",j,".lr",k,".bf",l))
            print(warnings())
            print("Note: previous warnings may be reprinted")
            ####8. Select best Gaus model####
            if (m == 1)
            {Gaus_Best_Score <- get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
            Gaus_Best_Model <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)
            } else if (get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]] > Gaus_Best_Score)
            {Gaus_Best_Score <- get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]]
            Gaus_Best_Model <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)} # close if else m==1

            if (alerts) beep(2) # progress printer, right aligned for visibility
            if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tc)*length(lr)*length(bf)) + (length(tcgaus)*length(lrgaus)*length(bfgaus)),"     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
            } else {print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Completed BRT ",n," of ", (length(tcgaus)*length(lrgaus)*length(bfgaus)),"     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))}

            ####9. Add gaus stats to report####
            Report[1:8,(3 + n)] <- c(paste0("trees: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$n.trees, 3)),
                                     paste0("Training Data Correlation: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics$correlation[[1]], 3)),
                                     paste0("CV Mean Deviance: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.mean, 3)),
                                     paste0("CV Deviance SE: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$deviance.se, 3)),
                                     paste0("CV D squared: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$d.squared, 3)),
                                     paste0("CV Mean Correlation: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.mean, 3)),
                                     paste0("CV Correlation SE: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$correlation.se, 3)),
                                     paste0("CV RMSE: ", round(get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics$cv.rmse, 3)))

            # Gaus BRT name
            colnames(Report)[3 + n] <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l)
            # Add Gaus stats objects to StatsObjectsList
            StatsObjectsList[[length(StatsObjectsList) + 1]] <- get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$self.statistics # send to new position after last item
            names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l, "_self.statistics") # name it. new length now includes self.statistics
            StatsObjectsList[[length(StatsObjectsList) + 1]] <- get(paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l))$cv.statistics
            names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0("Gaus_BRT",".tc",j,".lr",k,".bf",l, "_cv.statistics")
          } # close bfgaus
        } # close lrgaus
      } # close for j in tcgaus, making all Gaus BRT objects & continuing through model selection

      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        Closed Loops         XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####10. Test simplification benefit, do so if better####
      # copy Bin/Gaus_Best_Model to Name in case not created by Simp
      Bin_Best_Name <- Bin_Best_Model # both are 0 if not run
      if (gaus) Gaus_Best_Name <- Gaus_Best_Model

      # if simp TRUE & ZI=TRUE, run simplification test on best bin model
      if (simp) {
        # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {
          Bin_Best_Simp_Check <- gbm.simplify(get(Bin_Best_Model))
          dev.print(file = paste0("./",names(samples[i]),"/simp_drops_", fam1, "_tc",j,"lr",k,"bf",l,".jpeg"), device = jpeg, width = 600)
          # if best number of variables to remove isn't 0 (i.e. it's worth simplifying),
          # re-run best model (Bin_Best_Model, using gbm.call to get its values) with
          # just-calculated best number of variables to remove, removed. gbm.x asks which
          # number of drops has the minimum mean (lowest point on the line) & that calls
          # up the list of predictor variables with those removed, from $pred.list
          if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0) {
            assign("Bin_Best_Simp",
                   gbm.step.sd(data = samples,
                               # gbm.x = Bin_Best_Simp_Check$pred.list[[which.min(Bin_Best_Simp_Check$deviance.summary$mean)]],
                               gbm.x = as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Bin_Best_Simp_Check$final.drops$preds)]),
                               gbm.y = get(Bin_Best_Model)$gbm.call$gbm.y,
                               tree.complexity = get(Bin_Best_Model)$gbm.call$tree.complexity,
                               learning.rate = get(Bin_Best_Model)$gbm.call$learning.rate,
                               family = get(Bin_Best_Model)$gbm.call$family,
                               bag.fraction = get(Bin_Best_Model)$gbm.call$bag.fraction,
                               {if (!is.null(offset)) offset = grv_yes$offset},
                               ...))
            dev.print(file = paste0("./",names(samples[i]),"/pred_dev_", fam1, "_tc",j,"lr",k,"bf",l,"_simp.jpeg"), device = jpeg, width = 600)

            # Add Bin simp stats objects to StatsObjectsList
            StatsObjectsList[[length(StatsObjectsList) + 1]] <- Bin_Best_Simp$self.statistics # send to new position after last item
            names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0(Bin_Best_Model, "_Simp__self.statistics") # name it. new length now includes self.statistics
            StatsObjectsList[[length(StatsObjectsList) + 1]] <- Bin_Best_Simp$cv.statistics
            names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0(Bin_Best_Model, "_Simp__cv.statistics")

          } # close if min bin best simp

          if (alerts) beep(2) # progress printer, right aligned for visibility
          print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Bin model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        } # close if ZI

        # Same for Gaus
        if (gaus & exists("Gaus_Best_Model")) {
          Gaus_Best_Simp_Check <- gbm.simplify(get(Gaus_Best_Model))
          dev.print(file = paste0("./",names(samples[i]),"/simp_drops_", fam2, "_tc",j,"lr",k,"bf",l,".jpeg"), device = jpeg, width = 600)
          if (min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0) {
            assign("Gaus_Best_Simp",
                   gbm.step.sd(data = grv_yes,
                               # gbm.x = Gaus_Best_Simp_Check$pred.list[[which.min(Gaus_Best_Simp_Check$deviance.summary$mean)]],
                               # does the above line return ALL vars or just the simp predictors kept?####
                               # returns a subsample, but not the same as the one below. In bonnie example, top/original line is:
                               # 9 19 28 38
                               # below line is:
                               # "Latitude" "WTMP"
                               # code from report Gaus_BRT_simp predictors kept column, can potentially drop in place:
                               gbm.x = as.character(Gaus_Best_Simp_Check$final.drops$preds[((dim(subset(Gaus_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Gaus_Best_Simp_Check$final.drops$preds)]),
                               # if this works do the same for bin
                               gbm.y = get(Gaus_Best_Model)$gbm.call$gbm.y,
                               tree.complexity = get(Gaus_Best_Model)$gbm.call$tree.complexity,
                               learning.rate = get(Gaus_Best_Model)$gbm.call$learning.rate,
                               family = get(Gaus_Best_Model)$gbm.call$family,
                               bag.fraction = get(Gaus_Best_Model)$gbm.call$bag.fraction,
                               {if (!is.null(offset)) offset = grv_yes$offset},
                               ...))
            dev.print(file = paste0("./",names(samples[i]),"/pred_dev_", fam2, "_tc",j,"lr",k,"bf",l,"_simp.jpeg"), device = jpeg, width = 600)

            # Add Gaus simp stats objects to StatsObjectsList
            StatsObjectsList[[length(StatsObjectsList) + 1]] <- Gaus_Best_Simp$self.statistics # send to new position after last item
            names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0(Gaus_Best_Model, "_Simp__self.statistics") # name it. new length now includes self.statistics
            StatsObjectsList[[length(StatsObjectsList) + 1]] <- Gaus_Best_Simp$cv.statistics
            names(StatsObjectsList)[[length(StatsObjectsList)]] <- paste0(Gaus_Best_Model, "_Simp__cv.statistics")

          } # close if min gaus best simp
          if (alerts) beep(2)
          print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Simplified Gaus model    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        } # close gaus if

        ## Select final best models
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {  # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI. If Bin_Best has a simplified model:
          if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0) {
            # & if the simplified model has better correlation than Bin_Best itself
            if (Bin_Best_Simp$self.statistics$correlation > Bin_Best_Score[1]) {
              # then replace Bin_Best score/model values with those from the simplified model
              Bin_Best_Score <- Bin_Best_Simp$self.statistics$correlation
              Bin_Best_Name <- paste0(Bin_Best_Model, "_Simp")
              Bin_Best_Model <- "Bin_Best_Simp" # assign simp to best
            } # close if bin best simp
          } # close if min bin best simp
        } # close ZI

        # Same for Gaus:
        if (gaus & exists("Gaus_Best_Model")) {
          if (min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0) {
            if (Gaus_Best_Simp$self.statistics$correlation > Gaus_Best_Score[1]) {
              Gaus_Best_Score <- Gaus_Best_Simp$self.statistics$correlation
              Gaus_Best_Name <- paste0(Gaus_Best_Model, "_Simp")
              Gaus_Best_Model <- "Gaus_Best_Simp"
            } # close if gaus best
          } # close if min gaus best
        } # close if gaus
      } # close if simp

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Best models selected     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####11. Line plots####
      # All plots on one image for Bin & Gaus
      if (multiplot) { # don't do if multiplot=FALSE
        # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {  # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
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
          par(op)
        } # close if (fam1 == "bernoulli"

        if (gaus & exists("Gaus_Best_Model")) {
          png(filename = paste0("./",names(samples[i]),"/Gaus_Best_line.png"),
              width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
          gbm.plot(get(Gaus_Best_Model),
                   n.plots = length(get(Gaus_Best_Model)$contributions$var),
                   write.title = F, y.label = "Marginal Effect",
                   plot.layout = c(ceiling(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
                                   ifelse(sqrt(length(get(Gaus_Best_Model)$contributions$var))
                                          - floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))) < 0.5,
                                          floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))),
                                          floor(sqrt(length(get(Gaus_Best_Model)$contributions$var))) + 1)))
          dev.off() #close plot device
        } # close gaus if
      } # close multiplot if

      # All plots individually, named by explanatory variable, bin & gaus
      # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {
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
          # abline(h = 0, lty = 2) # https://github.com/SimonDedman/gbm.auto/issues/7 & https://github.com/rspatial/dismo/issues/41
          mtext("Marginal Effect", side = 2, line = 4.05, las = 0)
          # gbm.plot calls plot.gbm ~L47 but then centres to have 0 mean @L53
          # Asked Robert Hijmans to add a param to omit this: https://github.com/rspatial/dismo/issues/22
          dev.off()

          # create lines data to export to file. Need to recreate transformations from gbm.plot
          # Next 6 lines from GNG answer https://stats.stackexchange.com/a/144871/43360 which uses gbm.plot code

          # CHANGED
          # s <- match(get(Bin_Best_Model)$contributions$var[o], # original
          #            get(Bin_Best_Model)$gbm.call$predictor.names)

          # create dataframe
          # plotgrid <- plot.gbm(get(Bin_Best_Model), s, return.grid = TRUE)
          plotgrid <- plot.gbm(get(Bin_Best_Model), o, return.grid = TRUE) # CHANGED

          # This section centres the values around 0,
          # Inverts their position relative to 0 (top becomes bottom),
          # Exponentiates them (midrange values push towards extremes)
          # Then rescales from 0:1 to +/- values by subtracting the mean from each value
          # 2021-10-18 per https://github.com/SimonDedman/gbm.auto/issues/45
          # Why I built this code? Copied from https://stats.stackexchange.com/a/144871/43360
          # I don't think I need it nor is it helpful. Comment out.
          # # replace Y values in place with average-centred values
          # plotgrid[,2] <- plotgrid[,2] - mean(plotgrid[,2])
          # #Put Y values on a log scale
          # plotgrid[,2] <- 1 / (1 + exp(-plotgrid[,2]))
          # #Center the response to have zero mean over the data distribution
          # plotgrid[,2] <- scale(plotgrid[,2], center = TRUE, scale = FALSE)
          # 2023-02-24 this is useful for making ggplots like the plot.gbm
          # https://github.com/SimonDedman/gbm.auto/issues/81
          plotgrid$ycentred <- plotgrid$y - mean(plotgrid$y)

          #If factor variable
          if (is.factor(plotgrid[,1])) {
            plotgrid[,1] <- factor(plotgrid[,1], levels = levels(get(Bin_Best_Model)$gbm.call$dataframe[,get(Bin_Best_Model)$gbm.call$gbm.x[o]])) # CHANGED
            # needed at all? if it's a factor won't it already have levels?
          } # close if is factor

          if (linesfiles) {
            # write out csv
            write.csv(plotgrid, row.names = FALSE, na = "",
                      file = paste0("./", names(samples[i]), "/Bin_Best_line_",
                                    as.character(get(Bin_Best_Model)$gbm.call$predictor.names[o]), # CHANGED
                                    # as.character(get(Bin_Best_Model)$contributions$var[o]),
                                    ".csv"))
          } #close linesfiles

          if (is.factor(plotgrid[,1])) {
            # overwrite plot with gbm.factorplot output
            gbm.factorplot(x = plotgrid,
                           # factorplotlevels = NULL,
                           # ggplot2guideaxisangle = 0,
                           # ggplot2labsx = "",
                           # ggplot2labsy = "Marginal Effect",
                           # ggplot2axistext = 1.5,
                           # ggplot2axistitle = 2,
                           # ggplot2legendtext = 1,
                           # ggplot2legendtitle = 1.5,
                           # ggplot2legendtitlealign = 0, # otherwise effect type title centre aligned for some reason
                           # ggplot2plotbackgroundfill = "white", # white background
                           # ggplot2plotbackgroundcolour = "grey50",
                           # ggplot2striptextx = 2,
                           # ggplot2panelbordercolour = "black",
                           # ggplot2panelborderfill = NA,
                           # ggplot2panelborderlinewidth = 1,
                           # ggplot2legendspacingx = unit(0, "cm"), # compress spacing between legend items, this is min
                           # ggplot2legendbackground = ggplot2::element_blank(),
                           # ggplot2panelbackgroundfill = "white",
                           # ggplot2panelbackgroundcolour = "grey50",
                           # ggplot2panelgridcolour = "grey90",
                           # ggplot2legendkey = ggplot2::element_blank(),
                           ggsavefilename = paste0("Bin_Best_line_", as.character(get(Bin_Best_Model)$gbm.call$predictor.names[o]), "_gg.png"),
                           # ggsaveplot = last_plot(),
                           # ggsavedevice = "png",
                           ggsavepath = paste0("./", names(samples[i]), "/"),
                           # ggsavescale = 2,
                           ggsavewidth = 4*480,
                           ggsaveheight = 4*480,
                           ggsaveunits = "px",
                           # ggsavedpi = 300,
                           # ggsavelimitsize = TRUE
                           ...
            ) # close gbm.factorplot
          } # close if is factor
        } # close for o
      } # close if fam1 bernoulli / ZI option

      if (gaus & exists("Gaus_Best_Model")) {
        for (p in 1:length(get(Gaus_Best_Model)$contributions$var)) {
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
          # abline(h = 0, lty = 2) # https://github.com/SimonDedman/gbm.auto/issues/7 & https://github.com/rspatial/dismo/issues/41
          mtext("Marginal Effect", side = 2, line = 4.05, las = 0)
          dev.off()

          # u <- match(get(Gaus_Best_Model)$contributions$var[p],
          #            get(Gaus_Best_Model)$gbm.call$predictor.names)

          plotgrid <- plot.gbm(get(Gaus_Best_Model), p, return.grid = TRUE)
          # plotgrid[,2] <- plotgrid[,2] - mean(plotgrid[,2])
          # plotgrid[,2] <- 1 / (1 + exp(-plotgrid[,2]))
          # plotgrid[,2] <- scale(plotgrid[,2], scale = FALSE)
          plotgrid$ycentred <- plotgrid$y - mean(plotgrid$y)

          if (is.factor(plotgrid[,1])) {
            plotgrid[,1] <- factor(plotgrid[,1], levels = levels(get(Gaus_Best_Model)$gbm.call$dataframe[,get(Gaus_Best_Model)$gbm.call$gbm.x[p]]))
          } # close if is factor plotgrid

          if (linesfiles) {
            write.csv(plotgrid, row.names = FALSE, na = "",
                      file = paste0("./", names(samples[i]), "/Gaus_Best_line_",
                                    as.character(get(Gaus_Best_Model)$gbm.call$predictor.names[p]),
                                    ".csv"))
          } #close linesfiles

          if (is.factor(plotgrid[,1])) {
            gbm.factorplot(x = plotgrid,
                           ggsavefilename = paste0("Gaus_Best_line_", as.character(get(Gaus_Best_Model)$gbm.call$predictor.names[p]), "_gg.png"),
                           ggsavepath = paste0("./", names(samples[i]), "/"),
                           ggsavewidth = 4*480,
                           ggsaveheight = 4*480,
                           ggsaveunits = "px",
                           ...) # close gbm.factorplot
          } # close if is factor plotgrid
        } # close for p
      } # close if gaus

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Line plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####12. Dot plots####
      # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {  # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        png(filename = paste0("./",names(samples[i]),"/Bin_Best_dot.png"),
            width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
        gbm.plot.fits(get(Bin_Best_Model))
        dev.off()} # close ZI

      if (gaus & exists("Gaus_Best_Model")) {png(filename = paste0("./",names(samples[i]),"/Gaus_Best_dot.png"),
                                                 width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
        gbm.plot.fits(get(Gaus_Best_Model))
        dev.off()} # close if gaus

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Dot plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####13. 3D plot TODO####
      # gbm.perspec(Bin_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)
      # gbm.perspec(Gaus_Best,3,2, z.range=c(0,31), theta=340, phi=35,smooth="none",border="#00000025",col="#ff003310",shade = 0.95, ltheta = 80, lphi = 50)

      ####14. Bar plots of variable influence####
      # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {  # create tables
        Bin_Bars <- summary(get(Bin_Best_Model),
                            cBars = length(get(Bin_Best_Model)$var.names),
                            n.trees = get(Bin_Best_Model)$n.trees,
                            plotit = FALSE, order = TRUE, normalize = TRUE, las = 1, main = NULL)
        write.csv(Bin_Bars, file = paste0("./", names(samples[i]), "/Binary BRT Variable contributions.csv"), row.names = FALSE)} # close ZI

      if (gaus & exists("Gaus_Best_Model")) {Gaus_Bars <- summary(get(Gaus_Best_Model),
                                                                  cBars = length(get(Gaus_Best_Model)$var.names),
                                                                  n.trees = get(Gaus_Best_Model)$n.trees,
                                                                  plotit = FALSE, order = TRUE, normalize = TRUE, las = 1, main = NULL)
      write.csv(Gaus_Bars, file = paste0("./", names(samples[i]), "/Gaussian BRT Variable contributions.csv"), row.names = FALSE)} # close if gaus

      if (alerts) beep(2)# progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Bar plot csvs created    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {  # produce graphics
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

      if (gaus & exists("Gaus_Best_Model")) {
        pointlineseqgaus <- seq(0, length(Gaus_Bars[,2]) - 1, 1)
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
          lines(c(0, Gaus_Bars[r,2]), c(revseq[r], revseq[r]), col = "black", lwd = 8)
        } # close for r
        text(0.1, pointlineseqgaus + (length(Gaus_Bars[,2])/55), labels = rev(Gaus_Bars[,1]), adj = 0, cex = 0.8)
        axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
        dev.off() #close PNG
      } # close if gaus
      # col = rev(expvarcols[match(Bin_Bars[,1],expvarcols[,2]),1]), #in case I want to colour the bars/points later
      # elements of barplot lines+points code adapted from Jane Elith code donated to Agustín De Wysiecki & then shared with SD

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Bar plots plotted      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####15. Variable interactions####
      if (varint) {
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) find.int_Bin <- gbm.interactions(get(Bin_Best_Model))
        if (gaus & exists("Gaus_Best_Model")) find.int_Gaus <- gbm.interactions(get(Gaus_Best_Model))
        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Variable interactions done XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close varint if

      ####16. Save model objects####
      if (savegbm) { # Save model objects if switched on
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {
          Bin_Best_Model_Object <- get(Bin_Best_Model)
          # Bin_Best_Model <<- Bin_Best_Model_Object # this causes Bin_Best_Model to BE the model not the name of the original
        } # close if ZI
        if (gaus & exists("Gaus_Best_Model")) {
          Gaus_Best_Model_Object <- get(Gaus_Best_Model)
          # Gaus_Best_Model <<- Gaus_Best_Model_Object
          save(Gaus_Best_Model_Object, file = paste0("./",names(samples[i]),"/Gaus_Best_Model"))
        } # close if gaus
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) {
          save(Bin_Best_Model_Object, file = paste0("./",names(samples[i]),"/Bin_Best_Model")) #only save bin if ZI=TRUE
        } # close if ZI
        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Model objects saved     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close if savegbm

      ####17. Finalise & Write Report####
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) { # only do bin bits if ZI; move 7 cols left if no gaus run
        # Combine sections ZI yes, gaus ifelse was L812,873,879####
        if (gaus & exists("Gaus_Best_Model")) {
          Report[1:15,(reportcolno - 13)] <- c(paste0("Model combo: ", Bin_Best_Name),
                                               paste0("trees: ", get(Bin_Best_Model)$n.trees),
                                               paste0("Training Data Correlation: ", round(Bin_Best_Score, 3)),
                                               paste0("Training data AUC score: ", round(get(Bin_Best_Model)$self.statistics$discrimination, 3)),
                                               paste0("CV AUC score: ", round(get(Bin_Best_Model)$cv.statistics$discrimination.mean, 3)),
                                               paste0("CV AUC se: ", round(get(Bin_Best_Model)$cv.statistics$discrimination.se, 3)),
                                               paste0("Overfitting (Training data AUC - CV AUC): ", round(get(Bin_Best_Model)$self.statistics$discrimination - get(Bin_Best_Model)$cv.statistics$discrimination.mean, 3)),
                                               paste0("CV Mean Deviance: ", round(get(Bin_Best_Model)$cv.statistics$deviance.mean, 3)),
                                               paste0("CV Deviance SE: ", round(get(Bin_Best_Model)$cv.statistics$deviance.se, 3)),
                                               paste0("CV D squared: ", round(get(Bin_Best_Model)$cv.statistics$d.squared, 3)),
                                               paste0("CV Mean Correlation: ", round(get(Bin_Best_Model)$cv.statistics$correlation.mean, 3)),
                                               paste0("CV Correlation SE: ", round(get(Bin_Best_Model)$cv.statistics$correlation.se, 3)),
                                               paste0("CV RMSE: ", round(get(Bin_Best_Model)$cv.statistics$cv.rmse, 3)),
                                               paste0("Deviance% explained relative to null, training: ", round(((get(Bin_Best_Model)$self.statistics$mean.null - get(Bin_Best_Model)$self.statistics$mean.resid) / get(Bin_Best_Model)$self.statistics$mean.null)*100, 2)),
                                               paste0("Deviance% explained relative to null, CV: ", round(((get(Bin_Best_Model)$self.statistics$mean.null - get(Bin_Best_Model)$cv.statistics$deviance.mean) / get(Bin_Best_Model)$self.statistics$mean.null)*100, 2)))
        } else {
          Report[1:15,(reportcolno - 6)] <- c(paste0("Model combo: ", Bin_Best_Name),
                                              paste0("trees: ", get(Bin_Best_Model)$n.trees),
                                              paste0("Training Data Correlation: ", round(Bin_Best_Score, 3)),
                                              paste0("Training data AUC score: ", round(get(Bin_Best_Model)$self.statistics$discrimination, 3)),
                                              paste0("CV AUC score: ", round(get(Bin_Best_Model)$cv.statistics$discrimination.mean, 3)),
                                              paste0("CV AUC se: ", round(get(Bin_Best_Model)$cv.statistics$discrimination.se, 3)),
                                              paste0("Overfitting (Training data AUC - CV AUC): ", round(get(Bin_Best_Model)$self.statistics$discrimination - get(Bin_Best_Model)$cv.statistics$discrimination.mean, 3)),
                                              paste0("CV Mean Deviance: ", round(get(Bin_Best_Model)$cv.statistics$deviance.mean, 3)),
                                              paste0("CV Deviance SE: ", round(get(Bin_Best_Model)$cv.statistics$deviance.se, 3)),
                                              paste0("CV D squared: ", round(get(Bin_Best_Model)$cv.statistics$d.squared, 3)),
                                              paste0("CV Mean Correlation: ", round(get(Bin_Best_Model)$cv.statistics$correlation.mean, 3)),
                                              paste0("CV Correlation SE: ", round(get(Bin_Best_Model)$cv.statistics$correlation.se, 3)),
                                              paste0("CV RMSE: ", round(get(Bin_Best_Model)$cv.statistics$cv.rmse, 3)),
                                              paste0("Deviance% explained relative to null, training: ", round(((get(Bin_Best_Model)$self.statistics$mean.null - get(Bin_Best_Model)$self.statistics$mean.resid) / get(Bin_Best_Model)$self.statistics$mean.null)*100, 2)),
                                              paste0("Deviance% explained relative to null, CV: ", round(((get(Bin_Best_Model)$self.statistics$mean.null - get(Bin_Best_Model)$cv.statistics$deviance.mean) / get(Bin_Best_Model)$self.statistics$mean.null)*100, 2)))
        } # close if else gaus bin report

        if (simp) { # bin & gaus simp stats (or no simp notes)
          if (gaus & exists("Gaus_Best_Model")) {
            Report[1:dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1], (reportcolno - 12)] <- as.character(subset(Bin_Best_Simp_Check$final.drops, order > 0)$preds)
            # listing simp predictors kept: rows 1 to 'howevermany are left in simp' i.e. above 0
            # [1] is the first item of the dim list of rows & columns, i.e. rows
            Report[1:(length(Bin_Best_Simp_Check$final.drops$preds) - dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1]),(reportcolno - 11)] <-
              as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Bin_Best_Simp_Check$final.drops$preds)])
            # listing simp predictors dropped.
            if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0) {
              Report[1:14,(reportcolno - 10)] <- c(paste0("trees: ", Bin_Best_Simp$n.trees),
                                                   paste0("Training Data Correlation: ", round(Bin_Best_Simp$self.statistics$correlation[[1]], 3)),
                                                   paste0("Training data AUC score: ", round(Bin_Best_Simp$self.statistics$discrimination, 3)),
                                                   paste0("CV AUC score: ", round(Bin_Best_Simp$cv.statistics$discrimination.mean, 3)),
                                                   paste0("CV AUC se: ", round(Bin_Best_Simp$cv.statistics$discrimination.se, 3)),
                                                   paste0("Overfitting (Training data AUC - CV AUC): ", round(Bin_Best_Simp$self.statistics$discrimination - Bin_Best_Simp$cv.statistics$discrimination.mean, 3)),
                                                   paste0("CV Mean Deviance: ", round(Bin_Best_Simp$cv.statistics$deviance.mean, 3)),
                                                   paste0("CV Deviance SE: ", round(Bin_Best_Simp$cv.statistics$deviance.se, 3)),
                                                   paste0("CV D squared: ", round(Bin_Best_Simp$cv.statistics$d.squared, 3)),
                                                   paste0("CV Mean Correlation: ", round(Bin_Best_Simp$cv.statistics$correlation.mean, 3)),
                                                   paste0("CV Correlation SE: ", round(Bin_Best_Simp$cv.statistics$correlation.se, 3)),
                                                   paste0("CV RMSE: ", round(Bin_Best_Simp$cv.statistics$cv.rmse, 3)),
                                                   paste0("Deviance% explained relative to null, training: ", round(((Bin_Best_Simp$self.statistics$mean.null - Bin_Best_Simp$self.statistics$mean.resid) / Bin_Best_Simp$self.statistics$mean.null)*100, 2)),
                                                   paste0("Deviance% explained relative to null, CV: ", round(((Bin_Best_Simp$self.statistics$mean.null - Bin_Best_Simp$cv.statistics$deviance.mean) / Bin_Best_Simp$self.statistics$mean.null)*100, 2)))
            } else { # if min
              Report[1,(reportcolno - 10)] <- paste0("No simplification benefit")
            } # close if min else
          } else { # bin predictors etc
            Report[1:dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1], (reportcolno - 5)] <- as.character(subset(Bin_Best_Simp_Check$final.drops, order > 0)$preds)
            Report[1:(length(Bin_Best_Simp_Check$final.drops$preds) - dim(subset(Bin_Best_Simp_Check$final.drops, order > 0))[1]),(reportcolno - 4)] <-
              as.character(Bin_Best_Simp_Check$final.drops$preds[((dim(subset(Bin_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Bin_Best_Simp_Check$final.drops$preds)])
            if (min(Bin_Best_Simp_Check$deviance.summary$mean) < 0) {
              Report[1:14,(reportcolno - 3)] <- c(paste0("trees: ", Bin_Best_Simp$n.trees),
                                                  paste0("Training Data Correlation: ", round(Bin_Best_Simp$self.statistics$correlation[[1]], 3)),
                                                  paste0("Training data AUC score: ", round(Bin_Best_Simp$self.statistics$discrimination, 3)),
                                                  paste0("CV AUC score: ", round(Bin_Best_Simp$cv.statistics$discrimination.mean, 3)),
                                                  paste0("CV AUC se: ", round(Bin_Best_Simp$cv.statistics$discrimination.se, 3)),
                                                  paste0("Overfitting (Training data AUC - CV AUC): ", round(Bin_Best_Simp$self.statistics$discrimination - Bin_Best_Simp$cv.statistics$discrimination.mean, 3)),
                                                  paste0("CV Mean Deviance: ", round(Bin_Best_Simp$cv.statistics$deviance.mean, 3)),
                                                  paste0("CV Deviance SE: ", round(Bin_Best_Simp$cv.statistics$deviance.se, 3)),
                                                  paste0("CV D squared: ", round(Bin_Best_Simp$cv.statistics$d.squared, 3)),
                                                  paste0("CV Mean Correlation: ", round(Bin_Best_Simp$cv.statistics$correlation.mean, 3)),
                                                  paste0("CV Correlation SE: ", round(Bin_Best_Simp$cv.statistics$correlation.se, 3)),
                                                  paste0("CV RMSE: ", round(Bin_Best_Simp$cv.statistics$cv.rmse, 3)),
                                                  paste0("Deviance% explained relative to null, training: ", round(((Bin_Best_Simp$self.statistics$mean.null - Bin_Best_Simp$self.statistics$mean.resid) / Bin_Best_Simp$self.statistics$mean.null)*100, 2)),
                                                  paste0("Deviance% explained relative to null, CV: ", round(((Bin_Best_Simp$self.statistics$mean.null - Bin_Best_Simp$cv.statistics$deviance.mean) / Bin_Best_Simp$self.statistics$mean.null)*100, 2)))
            } else { # if min bin best simp
              Report[1,(reportcolno - 3)] <- paste0("No simplification benefit")
            } # close if min bin best simp else
          } # close bin half of bin/gaus option. Next line is 2nd half of simp option i.e. not simplified
        } else if (gaus) { # if not simp but is gaus
          Report[1,(reportcolno - 12):(reportcolno - 10)] <- c(paste0("simp turned off"),
                                                               paste0("simp turned off"),
                                                               paste0("simp turned off"))
        } else {# not gaus not simp: report cols are changed so needs adjustment
          Report[1,(reportcolno - 5):(reportcolno - 3)] <- c(paste0("simp turned off"),
                                                             paste0("simp turned off"),
                                                             paste0("simp turned off"))
        } # close 2nd half of simp else, i.e. nosimp bin. Closes bin & gaus simp stats (or no simp notes)

        if (gaus & exists("Gaus_Best_Model")) { # if zi & gaus again, second section
          Report[1:(length(Bin_Bars[,1])),(reportcolno - 9)] <- as.character(Bin_Bars$var)
        } else { # if zi & bin only
          Report[1:(length(Bin_Bars[,1])),(reportcolno - 2)] <- as.character(Bin_Bars$var)
        } # close ifelse best bin variables rel inf names ordered

        if (gaus & exists("Gaus_Best_Model")) { # if ZI & gaus again
          Report[1:(length(Bin_Bars[,2])),(reportcolno - 8)] <- as.character(round(Bin_Bars$rel.inf), 2)
        } else { # zi & not gaus
          Report[1:(length(Bin_Bars[,2])),(reportcolno - 1)] <- as.character(round(Bin_Bars$rel.inf), 2)
        } # close ifelse best bin variables rel inf scores

        if (varint) { # only do final variable interaction lines if varint=TRUE
          if (gaus & exists("Gaus_Best_Model")) { # varint gaus
            # Report[1:2,(reportcolno - 7)] <- c(paste0(find.int_Bin$rank.list$var1.names[1]," and ",find.int_Bin$rank.list$var2.names[1],". Size: ",find.int_Bin$rank.list$int.size[1]),
            #                                    paste0(find.int_Bin$rank.list$var1.names[2]," and ",find.int_Bin$rank.list$var2.names[2],". Size: ",find.int_Bin$rank.list$int.size[2]))
            # 2020-12-17 update, see below
            for (v in 1:min(length(which(find.int_Bin$rank.list$int.size > 0)),
                            nrow(Report))) {
              Report[v, (reportcolno - 7)] <- paste0(find.int_Bin$rank.list$var1.names[v]," and ",find.int_Bin$rank.list$var2.names[v],". Size: ",find.int_Bin$rank.list$int.size[v])
            }
          } else { # varint yes gaus no
            # Report[1:2,(reportcolno)] <- c(paste0(find.int_Bin$rank.list$var1.names[1]," and ",find.int_Bin$rank.list$var2.names[1],". Size: ",find.int_Bin$rank.list$int.size[1]),
            #                                paste0(find.int_Bin$rank.list$var1.names[2]," and ",find.int_Bin$rank.list$var2.names[2],". Size: ",find.int_Bin$rank.list$int.size[2]))
            # 2020-12-17 update, see below
            for (u in 1:min(length(which(find.int_Bin$rank.list$int.size > 0)),
                            nrow(Report))) {
              Report[u, (reportcolno)] <- paste0(find.int_Bin$rank.list$var1.names[u]," and ",find.int_Bin$rank.list$var2.names[u],". Size: ",find.int_Bin$rank.list$int.size[u])
            }


          } # close varint yes gaus no
        } else { # varint no
          if (gaus & exists("Gaus_Best_Model")) { # varint no gaus yes
            Report[1,(reportcolno - 7)] <- paste0("varint turned off")
          } else { # varint no gaus no
            Report[1,(reportcolno)] <- paste0("varint turned off")
          } # close not varint not gaus
        } # close not varint
      } # close ZI way further up start of report section (L810)

      if (gaus & exists("Gaus_Best_Model")) {
        Report[1:11 ,(reportcolno - 6)] <- c(paste0("Model combo: ", Gaus_Best_Name),
                                             paste0("trees: ", get(Gaus_Best_Model)$n.trees),  # new, might not work
                                             paste0("Training Data Correlation: ", round(Gaus_Best_Score, 3)),
                                             paste0("CV Mean Deviance: ", round(get(Gaus_Best_Model)$cv.statistics$deviance.mean, 3)), # new, might not work
                                             paste0("CV Deviance SE: ", round(get(Gaus_Best_Model)$cv.statistics$deviance.se, 3)), # new, might not work
                                             paste0("CV D squared: ", round(get(Gaus_Best_Model)$cv.statistics$d.squared, 3)), # new, might not work
                                             paste0("CV Mean Correlation: ", round(get(Gaus_Best_Model)$cv.statistics$correlation.mean, 3)), # new, might not work
                                             paste0("CV Correlation SE: ", round(get(Gaus_Best_Model)$cv.statistics$correlation.se, 3)), # new, might not work
                                             paste0("CV RMSE: ", round(get(Gaus_Best_Model)$cv.statistics$cv.rmse, 3)), # new, might not work
                                             paste0("Deviance% explained relative to null, training: ", round(((get(Gaus_Best_Model)$self.statistics$mean.null - get(Gaus_Best_Model)$self.statistics$mean.resid) / get(Gaus_Best_Model)$self.statistics$mean.null)*100, 2)), # new, might not work
                                             paste0("Deviance% explained relative to null, CV: ", round(((get(Gaus_Best_Model)$self.statistics$mean.null - get(Gaus_Best_Model)$cv.statistics$deviance.mean) / get(Gaus_Best_Model)$self.statistics$mean.null)*100, 2))) # new, might not work

        if (simp) {
          Report[1:dim(subset(Gaus_Best_Simp_Check$final.drops,order > 0))[1], (reportcolno - 5)] <- as.character(subset(Gaus_Best_Simp_Check$final.drops ,order > 0)$preds)
          Report[1:(length(Gaus_Best_Simp_Check$final.drops$preds) - dim(subset(Gaus_Best_Simp_Check$final.drops, order > 0))[1]), (reportcolno - 4)] <-
            as.character(Gaus_Best_Simp_Check$final.drops$preds[((dim(subset(Gaus_Best_Simp_Check$final.drops,order > 0))[1]) + 1):length(Gaus_Best_Simp_Check$final.drops$preds)])
          if (min(Gaus_Best_Simp_Check$deviance.summary$mean) < 0) {
            Report[1:10, (reportcolno - 3)] <- c(paste0("trees: ", Gaus_Best_Simp$n.trees),
                                                 paste0("Training Data Correlation: ", round(Gaus_Best_Simp$self.statistics$correlation[[1]], 3)),
                                                 paste0("CV Mean Deviance: ", round(Gaus_Best_Simp$cv.statistics$deviance.mean, 3)),
                                                 paste0("CV Deviance SE: ", round(Gaus_Best_Simp$cv.statistics$deviance.se, 3)),
                                                 paste0("CV D squared: ", round(Gaus_Best_Simp$cv.statistics$d.squared, 3)),
                                                 paste0("CV Mean Correlation: ", round(Gaus_Best_Simp$cv.statistics$correlation.mean, 3)),
                                                 paste0("CV Correlation SE: ", round(Gaus_Best_Simp$cv.statistics$correlation.se, 3)),
                                                 paste0("CV RMSE: ", round(Gaus_Best_Simp$cv.statistics$cv.rmse, 3)),
                                                 paste0("Deviance% explained relative to null, training: ", round(((Gaus_Best_Simp$self.statistics$mean.null - Gaus_Best_Simp$self.statistics$mean.resid) / Gaus_Best_Simp$self.statistics$mean.null)*100, 2)), # new, might not work
                                                 paste0("Deviance% explained relative to null, CV: ", round(((Gaus_Best_Simp$self.statistics$mean.null - Gaus_Best_Simp$cv.statistics$deviance.mean) / Gaus_Best_Simp$self.statistics$mean.null)*100, 2)))
          } else { # else if min gaus best simp, stats where simp benefit true, open note where no simp benefit
            Report[1,(reportcolno - 3)] <- paste0("No simplification benefit")
          } # close if else simp benefit check
        } else { # close gaus yes simp yes, do gaus yes simp no
          Report[1,(reportcolno - 5):(reportcolno - 3)] <- c(paste0("simp turned off"),
                                                             paste0("simp turned off"),
                                                             paste0("simp turned off"))
        } # close simp, still in gaus yes
        Report[1:(length(Gaus_Bars[,1])),(reportcolno - 2)] <- as.character(Gaus_Bars$var)
        Report[1:(length(Gaus_Bars[,2])),(reportcolno - 1)] <- as.character(round(Gaus_Bars$rel.inf), 2)
        if (varint) { # gaus yes varint yes
          # 2020.12.17 instead of top 2 interactions, do all of them that are > 0 size (unless there are more than the report nrow)
          for (t in 1:min(length(which(find.int_Gaus$rank.list$int.size > 0)),
                          nrow(Report))) {
            Report[t, (reportcolno)] <- paste0(find.int_Gaus$rank.list$var1.names[t]," and ",find.int_Gaus$rank.list$var2.names[t],". Size: ",find.int_Gaus$rank.list$int.size[t])
          }
          # Report[1:2,(reportcolno)] <- c(paste0(find.int_Gaus$rank.list$var1.names[1]," and ",find.int_Gaus$rank.list$var2.names[1],". Size: ",find.int_Gaus$rank.list$int.size[1]),
          #                                paste0(find.int_Gaus$rank.list$var1.names[2]," and ",find.int_Gaus$rank.list$var2.names[2],". Size: ",find.int_Gaus$rank.list$int.size[2]))
        } else { # gaus yes varint no
          Report[1,(reportcolno)] <- paste0("varint turned off")
        } # close varint yes no, still in gaus yes
      } # close gaus if

      Report[, "Metadata"] <- as.character(NA)
      Report$Metadata[1:4] <- c(
        paste0("gbm.auto v", packageVersion("gbm.auto")),
        paste0("dismo v", packageVersion("dismo")),
        paste0("fam1 n ", length(samples[, i])),
        ifelse(gaus, paste0("fam2 n ", length(grv_yes[, i])), as.character(NA))
      )

      # ) was a closed bracket here, I think errantly

      write.csv(Report, row.names = FALSE, na = "", file = paste0("./", names(samples[i]), "/Report.csv"))
      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Report CSV written      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      # 18. Finalise & Write Self & CV Stats csv ####
      StatsObjectsDf <- data.frame(StatsNames = names(unlist(StatsObjectsList)),
                                   Value = round(unlist(StatsObjectsList), 3),
                                   row.names = NULL)
      StatsObjectsNames <- as.data.frame(stri_split_fixed(str = StatsObjectsDf$StatsNames,
                                                          pattern = "__",
                                                          n = 2,
                                                          simplify = TRUE))
      colnames(StatsObjectsNames) <- c("Model", "Test_Statistic")
      Self_CV_Statistics <- data.frame(StatsObjectsNames, Value = StatsObjectsDf[,"Value"])
      write.csv(x = Self_CV_Statistics, row.names = FALSE, na = "",
                file = paste0("./", names(samples[i]), "/Self_CV_Statistics.csv"))
      rm(list = c("StatsObjectsList", "StatsObjectsDf", "StatsObjectsNames", "Self_CV_Statistics"))


      #19. Machine learning evaluation metrics####
      if (MLEvaluate) { # if user wants ML evaluation
        # if (any(fam1 == "bernoulli", fam2 == "bernoulli")) {
        #   whichbin <- which(c(fam1 == "bernoulli", fam2 == "bernoulli"))
        #   if (whichbin == 1) getmodel <- "Bin_Best_Model" else getmodel <- "Gaus_Best_Model"
        # } # close if any fam 1

        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI)) & exists("Bin_Best_Model")) { # only do if bin exists, previously was: exists("Bin_Best_Model")
          preds <- predict.gbm(get(Bin_Best_Model),
                               samples, # predict back to samples, not out of bag, for performance evaluation
                               n.trees = get(Bin_Best_Model)$gbm.call$best.trees,
                               type = "response")
          #If type="response" then gbm converts back to the same scale as the outcome.
          # Currently the only effect this will have is returning probabilities for
          # bernoulli and expected counts for poisson. For the other distributions
          # "response" and "link" return the same. gbm:::predict.gbm

          # calc.deviance = remaining deviance
          pct.dev.remain.samples <- calc.deviance(obs = samples[, get(Bin_Best_Model)$gbm.call$gbm.y],
                                                  pred = preds,
                                                  family = "bernoulli") # change fam if using
          pct.dev.remain.samples <- (1 - pct.dev.remain.samples) * 100 # convert to % deviance explained
          AUC.samples <- roc(obsdat = samples[, get(Bin_Best_Model)$gbm.call$gbm.y],
                             preddat = preds)

          # One of "binomial", "bernoulli", "poisson", "laplace", or "gaussian"
          samples <- cbind(samples, preds)
          pres <- samples[samples[, brvcol] == 1, "preds"] # check brvcol indexed properly, ditto last col is preds
          abs <- samples[samples[, brvcol] == 0, "preds"]
          e <- evaluate(p = pres,
                        a = abs)
          # saveRDS(e, file = paste0("./",names(samples[i]),"/e.rds")) # saveout for diagnostics, bug since fixed

          # Fielding, A. H. & J.F. Bell, 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24: 38-49
          # Liu, C., M. White & G. Newell, 2011. Measuring and comparing the accuracy of species distribution models with presence-absence data. Ecography 34: 232-243.
          MLEvalLength <- 32
          # Improve ML STATS descriptions####
          MLEval <- data.frame(Statistic = rep(NA, MLEvalLength),
                               Description = rep(NA, MLEvalLength),
                               Value = rep(NA, MLEvalLength))
          MLEval[1,] <- c("Presence",
                          "n of presence data used",
                          e@np)
          MLEval[2,] <- c("Absence",
                          "n of absence data used",
                          e@na)
          MLEval[3,] <- c("AUC",
                          "Area under the receiver operator (ROC) curve",
                          round(e@auc, 3))
          if (length(e@pauc) == 0) e@pauc <- 0 # pauc may be missing, numeric(0), if so replace with 0
          MLEval[4,] <- c("pAUC",
                          "p-value for the AUC (for the Wilcoxon test W statistic)",
                          round(e@pauc,3))
          MLEval[5,] <- c("Cor",
                          "Correlation coefficient",
                          round(e@cor[[1]],3))
          MLEval[6,] <- c("cor",
                          "p-value for correlation coefficient",
                          round(e@pcor,3))
          # Steph Brodie's TSS which produces the same result as Allouche
          # -1 just makes the output range is 0:1 instead of 1:2 I think.
          # If so this means Sensitivity is e@TPR[which.max(e@TPR + e@TNR)], which doesn't include
          # (e@TPR + e@FNR) but it's a vector of 1s so is redundant. Same for Specificity
          MLEval[7,] <- c("TSS",
                          "True Skill Statistic",
                          round(max(e@TPR + e@TNR - 1),3))
          # sensitivity: TP/(TP+FN)
          MLEval[8,] <- c("Sensitivity",
                          "Sensitivity",
                          round(e@TPR[which.max(e@TPR + e@TNR)],3))
          # specificity: TN/(FP+TN)
          Specificity <- round(e@TNR[which.max(e@TPR + e@TNR)],3)
          MLEval[9,] <- c("Specificity",
                          "Specificity",
                          Specificity)
          # Accuracy: TP+TN / TP+TN+FP+FN true false positive negative.
          # TP+TN is just TSS + 1, TP+TN+FP+FN #Sums to 2, redundant
          MLEval[10,] <- c("Accuracy",
                           "Accuracy",
                           round((e@TPR[which.max(e@TPR + e@TNR)] + e@TNR[which.max(e@TPR + e@TNR)]) / 2, 3))
          # Precision: TP/TP+FP. Ignores true negatives. “X% of the predictions are right”
          Precision <- e@TPR[which.max(e@TPR + e@TNR)] / (e@TPR[which.max(e@TPR + e@TNR)] + e@FPR[which.max(e@TPR + e@TNR)])
          MLEval[11,] <- c("Precision",
                           "X% of the predictions are right",
                           round(Precision,3))
          # Recall: TP/TP+FN: “Y% of actually existing things are captured”.
          Recall <- e@TPR[which.max(e@TPR + e@TNR)] / (e@TPR[which.max(e@TPR + e@TNR)] + e@FNR[which.max(e@TPR + e@TNR)])
          MLEval[12,] <- c("Recall",
                           "Y% of actually existing things are captured",
                           round(Recall,3))
          # https://www.corvil.com/kb/what-is-a-false-positive-rate
          # Allouche et al 2006:
          # overall accuracy: (TP+TN)/n
          # this seems like a weird metric since the numerator is 0:2 or 1:2 and the divisor could be tiny or huge
          MLEval[13,] <- c("OverallAccuracy",
                           "Overall Accuracy",
                           round((e@TPR[which.max(e@TPR + e@TNR)] + e@TNR[which.max(e@TPR + e@TNR)])/nrow(samples),3))
          # Balanced Accuracy, (Recall + Specificity) / 2
          MLEval[14,] <- c("BalancedAccuracy",
                           "Balanced Accuracy",
                           round((Recall + Specificity) / 2,3))
          # Number of samples. Useful to include in the list
          MLEval[15,] <- c("nSamples",
                           "Number of samples",
                           nrow(samples))
          # Balance: precision vs recall curve. Workhorses.
          # PxR/P+R = F score (P+R = harmonic mean).
          # F1 score: P & R are equally rated. This is the most common one. F1 score importance depends on the project.
          MLEval[16,] <- c("F1score",
                           "P & R equally rated, score importance depends on project",
                           round(2 * ((Precision * Recall) / (Precision + Recall)),3))
          # F2 score: weighted average of Precision & Recall
          MLEval[17,] <- c("F2score",
                           "weighted average of P & R",
                           round(5 * ((Precision * Recall) / (4 * Precision + Recall)),3))
          # Threshold which produces the best combo of TPR & TNR
          # t: vector of thresholds used to compute confusion matrices
          MLEval[18,] <- c("Threshold",
                           "Threshold which produced best combo of TPR & TNR",
                           round(e@t[which.max(e@TPR + e@TNR)],3))
          # e@prevalence: Prevalence
          MLEval[19,] <- c("Prevalence",
                           "Prevalence",
                           round(e@prevalence[which.max(e@TPR + e@TNR)],3))
          # e@ODP: Overall diagnostic power
          MLEval[20,] <- c("ODP",
                           "Overall diagnostic power",
                           round(e@ODP[which.max(e@TPR + e@TNR)],3))
          # e@CCR: Correct classification rate
          MLEval[21,] <- c("CCR",
                           "Correct classification rate",
                           round(e@CCR[which.max(e@TPR + e@TNR)],3))
          # e@TPR: True positive rate
          MLEval[22,] <- c("TPR",
                           "True positive rate",
                           round(e@TPR[which.max(e@TPR + e@TNR)],3))
          # e@TNR: True negative rate
          MLEval[23,] <- c("TNR",
                           "True negative rate",
                           round(e@TNR[which.max(e@TPR + e@TNR)],3))
          # e@FPR: False positive rate
          MLEval[24,] <- c("FPR",
                           "False positive rate",
                           round(e@FPR[which.max(e@TPR + e@TNR)],3))
          # e@FNR: False negative rate
          MLEval[25,] <- c("FNR",
                           "False negative rate",
                           round(e@FNR[which.max(e@TPR + e@TNR)],3))
          # e@PPP: Positive predictive power
          MLEval[26,] <- c("PPP",
                           "Positive predictive power",
                           round(e@PPP[which.max(e@TPR + e@TNR)],3))
          # e@NPP: Negative predictive power
          MLEval[27,] <- c("NPP",
                           "Negative predictive power",
                           round(e@NPP[which.max(e@TPR + e@TNR)],3))
          # e@MCR: Misclassification rate
          MLEval[28,] <- c("MCR",
                           "Misclassification rate",
                           round(e@MCR[which.max(e@TPR + e@TNR)],3))
          # e@OR: Odds-ratio
          MLEval[29,] <- c("OR",
                           "Odds-ratio",
                           round(e@OR[which.max(e@TPR + e@TNR)],3))
          # e@kappa: Cohen's kappa
          MLEval[30,] <- c("kappa",
                           "Cohen's kappa",
                           round(e@kappa[which.max(e@TPR + e@TNR)],3))
          # (1 - pct.dev.remain.samples) * 100; pct.dev.remain.samples from calc.deviance from dismo
          MLEval[31,] <- c("pct.dev.remain.samples",
                           "% deviance explained, samples data only, calc.deviance frunction in gbm.auto, Leathwick & Elith",
                           round(pct.dev.remain.samples,3))
          MLEval[32,] <- c("AUC.samples",
                           "AUC for samples data, roc function in gbm.auto, by J.Elith",
                           round(AUC.samples,3))

          # MLEval$Value <- round(MLEval$Value, digits = 5)
          write.csv(MLEval, row.names = FALSE, na = "", file = paste0("./", names(samples[i]), "/MLEvalMetricsBin.csv"))

          # 2023-03-12
          # bug where plot(e,s) wanted type=double because ModeEvaluation wasn't an exported method from dismo:
          # https://stackoverflow.com/questions/75669228/r-plot-error-in-as-doublex-cannot-coerce-type-s4-to-vector-of-type-doubl
          # Raised in bug: https://github.com/rspatial/dismo/issues/37
          # Fixed in commit: https://github.com/rspatial/dismo/commit/e4cc66f2abaec80a46d40f0afc6f7b06f16c7228
          # Included dismo github 1.3-10 into DESCRIPTION via:
          # use_dev_package(package = "dismo",type = "Imports",remote = "rspatial")
          # adds Remotes section. Manually changed dismo version.
          # Per https://r-pkgs.org/dependencies-in-practice.html#depending-on-the-development-version-of-a-package
          # After document(): Warning message: In loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]):
          #  namespace ‘dismo’ 1.3-9 is already loaded, but >= 1.3.10 is required
          # Note I need to change to the next CRAN dismo submission before I can submit to CRAN again.


          evalmetrics <- c("ROC", "kappa", "prevalence", "TPR", "TNR", "FPR", "FNR", "CCR", "PPP", "NPP", "MCR", "OR")
          for (s in evalmetrics) {
            png(filename = paste0("./",names(samples[i]),"/Bin_Eval_", s, ".png"))
            plot(e, s)
            dev.off()
          } # close for (s in evalmetrics)
        } # close if(fam1 == "bernoulli" & (!gaus | (gaus & ZI)))


        # # can do calc.deviance for gaus also, ditto poisson
        # #ToFix GAUS ML STATS####
        # #Gaus metrics won't work as is
        # #Also code this better to reduce duplication & allow for bin & gaus runs
        # if (gaus & exists("Gaus_Best_Model")) { #
        #   preds <- predict.gbm(get(Gaus_Best_Model),
        #                        samples,
        #                        n.trees = get(Gaus_Best_Model)$gbm.call$best.trees,
        #                        type = "response")
        #   #If type="response" then gbm converts back to the same scale as the outcome.
        #   # Currently the only effect this will have is returning probabilities for
        #   # bernoulli and expected counts for poisson. For the other distributions
        #   # "response" and "link" return the same. gbm:::predict.gbm
        #
        #   # dev reported later but not used otherwise
        #   dev <- calc.deviance(obs = samples[, get(Gaus_Best_Model)$gbm.call$gbm.y],
        #                        pred = preds,
        #                        family = "Gaussian") # change fam if using
        #   # One of "binomial", "bernoulli", "poisson", "laplace", or "gaussian"
        #   samples <- cbind(samples, preds)
        #   pres <- samples[samples[, brvcol] == 1, "preds"] # check brvcol indexed properly, ditto last col is preds
        #   abs <- samples[samples[, brvcol] == 0, "preds"]
        #   # GAUS ML STATS FAILS HERE ####
        #   # THERE MAY NOT BE ANY ABSENCES IN A GAUSSIAN DISTRIBUTION
        #   # means abs = numeric(0) & evaluate() fails.
        #   e <- evaluate(p = pres,
        #                 a = abs)
        #
        #   # Fielding, A. H. & J.F. Bell, 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24: 38-49
        #   # Liu, C., M. White & G. Newell, 2011. Measuring and comparing the accuracy of species distribution models with presence-absence data. Ecography 34: 232-243.
        #   MLEvalLength <- 31
        #   # Improve descriptions####
        #   MLEval <- data.frame(Statistic = rep(NA, MLEvalLength),
        #                        Description = rep(NA, MLEvalLength),
        #                        Value = rep(NA, MLEvalLength))
        #   MLEval[1,] <- c("Presence",
        #                   "n of presence data used",
        #                   e@np)
        #   MLEval[2,] <- c("Absence",
        #                   "n of absence data used",
        #                   e@na)
        #   MLEval[3,] <- c("AUC",
        #                   "Area under the receiver operator (ROC) curve",
        #                   e@auc)
        #   if (length(e@pauc) == 0) e@pauc <- 0 # pauc may be missing, numeric(0), if so replace with 0
        #   MLEval[4,] <- c("pAUC",
        #                   "p-value for the AUC (for the Wilcoxon test W statistic)",
        #                   e@pauc)
        #   MLEval[5,] <- c("Cor",
        #                   "Correlation coefficient",
        #                   e@cor[[1]])
        #   MLEval[6,] <- c("cor",
        #                   "p-value for correlation coefficient",
        #                   e@pcor)
        #   # Steph Brodie's TSS which produces the same result as Allouche
        #   # -1 just makes the output range is 0:1 instead of 1:2 I think.
        #   # If so this means Sensitivity is e@TPR[which.max(e@TPR + e@TNR)], which doesn't include
        #   # (e@TPR + e@FNR) but it's a vector of 1s so is redundant. Same for Specificity
        #   MLEval[7,] <- c("TSS",
        #                   "True Skill Statistic",
        #                   max(e@TPR + e@TNR - 1))
        #   # sensitivity: TP/(TP+FN)
        #   MLEval[8,] <- c("Sensitivity",
        #                   "Sensitivity",
        #                   e@TPR[which.max(e@TPR + e@TNR)])
        #   # specificity: TN/(FP+TN)
        #   Specificity <- e@TNR[which.max(e@TPR + e@TNR)]
        #   MLEval[9,] <- c("Specificity",
        #                   "Specificity",
        #                   Specificity)
        #   # Accuracy: TP+TN / TP+TN+FP+FN true false positive negative.
        #   # TP+TN is just TSS + 1, TP+TN+FP+FN #Sums to 2, redundant
        #   MLEval[10,] <- c("Accuracy",
        #                    "Accuracy",
        #                    (e@TPR[which.max(e@TPR + e@TNR)] + e@TNR[which.max(e@TPR + e@TNR)]) / 2)
        #   # Precision: TP/TP+FP. Ignores true negatives. “X% of the predictions are right”
        #   Precision <- e@TPR[which.max(e@TPR + e@TNR)] / (e@TPR[which.max(e@TPR + e@TNR)] + e@FPR[which.max(e@TPR + e@TNR)])
        #   MLEval[11,] <- c("Precision",
        #                    "X% of the predictions are right",
        #                    Precision)
        #   # Recall: TP/TP+FN: “Y% of actually existing things are captured”.
        #   Recall <- e@TPR[which.max(e@TPR + e@TNR)] / (e@TPR[which.max(e@TPR + e@TNR)] + e@FNR[which.max(e@TPR + e@TNR)])
        #   MLEval[12,] <- c("Recall",
        #                    "Y% of actually existing things are captured",
        #                    Recall)
        #   # https://www.corvil.com/kb/what-is-a-false-positive-rate
        #   # Allouche et al 2006:
        #   # overall accuracy: (TP+TN)/n
        #   # this seems like a weird metric since the numerator is 0:2 or 1:2 and the divisor could be tiny or huge
        #   MLEval[13,] <- c("OverallAccuracy",
        #                    "Overall Accuracy",
        #                    (e@TPR[which.max(e@TPR + e@TNR)] + e@TNR[which.max(e@TPR + e@TNR)])/nrow(samples))
        #   # Balanced Accuracy, (Recall + Specificity) / 2
        #   MLEval[14,] <- c("BalancedAccuracy",
        #                    "Balanced Accuracy",
        #                    (Recall + Specificity) / 2)
        #   # Number of samples. Useful to include in the list
        #   MLEval[15,] <- c("nSamples",
        #                    "Number of samples",
        #                    nrow(samples))
        #   # Balance: precision vs recall curve. Workhorses.
        #   # PxR/P+R = F score (P+R = harmonic mean).
        #   # F1 score: P & R are equally rated. This is the most common one. F1 score importance depends on the project.
        #   MLEval[16,] <- c("F1score",
        #                    "P & R equally rated, score importance depends on project",
        #                    2 * ((Precision * Recall) / (Precision + Recall)))
        #   # F2 score: weighted average of Precision & Recall
        #   MLEval[17,] <- c("F2score",
        #                    "weighted average of P & R",
        #                    5 * ((Precision * Recall) / (4 * Precision + Recall)))
        #   # Threshold which produces the best combo of TPR & TNR
        #   # t: vector of thresholds used to compute confusion matrices
        #   MLEval[18,] <- c("Threshold",
        #                    "Threshold which produced best combo of TPR & TNR",
        #                    e@t[which.max(e@TPR + e@TNR)])
        #   # e@prevalence: Prevalence
        #   MLEval[19,] <- c("Prevalence",
        #                    "Prevalence",
        #                    e@prevalence[which.max(e@TPR + e@TNR)])
        #   # e@ODP: Overall diagnostic power
        #   MLEval[20,] <- c("ODP",
        #                    "Overall diagnostic power",
        #                    e@ODP[which.max(e@TPR + e@TNR)])
        #   # e@CCR: Correct classification rate
        #   MLEval[21,] <- c("CCR",
        #                    "Correct classification rate",
        #                    e@CCR[which.max(e@TPR + e@TNR)])
        #   # e@TPR: True positive rate
        #   MLEval[22,] <- c("TPR",
        #                    "True positive rate",
        #                    e@TPR[which.max(e@TPR + e@TNR)])
        #   # e@TNR: True negative rate
        #   MLEval[23,] <- c("TNR",
        #                    "True negative rate",
        #                    e@TNR[which.max(e@TPR + e@TNR)])
        #   # e@FPR: False positive rate
        #   MLEval[24,] <- c("FPR",
        #                    "False positive rate",
        #                    e@FPR[which.max(e@TPR + e@TNR)])
        #   # e@FNR: False negative rate
        #   MLEval[25,] <- c("FNR",
        #                    "False negative rate",
        #                    e@FNR[which.max(e@TPR + e@TNR)])
        #   # e@PPP: Positive predictive power
        #   MLEval[26,] <- c("PPP",
        #                    "Positive predictive power",
        #                    e@PPP[which.max(e@TPR + e@TNR)])
        #   # e@NPP: Negative predictive power
        #   MLEval[27,] <- c("NPP",
        #                    "Negative predictive power",
        #                    e@NPP[which.max(e@TPR + e@TNR)])
        #   # e@MCR: Misclassification rate
        #   MLEval[28,] <- c("MCR",
        #                    "Misclassification rate",
        #                    e@MCR[which.max(e@TPR + e@TNR)])
        #   # e@OR: Odds-ratio
        #   MLEval[29,] <- c("OR",
        #                    "Odds-ratio",
        #                    e@OR[which.max(e@TPR + e@TNR)])
        #   # e@kappa: Cohen's kappa
        #   MLEval[30,] <- c("kappa",
        #                    "Cohen's kappa",
        #                    e@kappa[which.max(e@TPR + e@TNR)])
        #   # dev from calc.deviance from dismo
        #   MLEval[31,] <- c("dev",
        #                    "deviance from 2 vecs, obs & pred vals",
        #                    dev)
        #
        #   # MLEval$Value <- round(MLEval$Value, digits = 5)
        #   write.csv(MLEval, row.names = FALSE, na = "", file = paste0("./", names(samples[i]), "/MLEvalMetricsGaus.csv"))
        #
        #   evalmetrics <- c("ROC", "kappa", "prevalence", "TPR", "TNR", "FPR", "FNR", "CCR", "PPP", "NPP", "MCR", "OR")
        #   for (s in evalmetrics) {
        #     png(filename = paste0("./",names(samples[i]),"/Gaus_Eval_", s, ".png"))
        #     plot(e, s)
        #     dev.off()
        #   } # close for (s in evalmetrics)
        # } # close if (gaus)

        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Evaluation Metrics ProcessedXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close if MLEvaluate

    } # close loadgbm isnull

    #avoid sections 19-25 if not predicting to grids
    if (!is.null(grids)) {

      # Load model objects if loadgbm set
      if (!is.null(loadgbm)) {
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {  # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
          # & exists("Bin_Best_Model") # not sure why this was in the above line; if using loadgbm then Bin_)Best_model necessarily won't exist in the environment.
          if (length(list.files(path = loadgbm, pattern = "Bin_Best_Model")) != 1) stop("Bin_Best_Model not found at loadgbm path")
          load(paste0(loadgbm, "Bin_Best_Model"))
          Bin_Best_Model <- "Bin_Best_Model_Object"
        } # close ZI if
        if (gaus) {
          #  & exists("Gaus_Best_Model")
          if (length(list.files(path = loadgbm, pattern = "Gaus_Best_Model")) != 1) stop("Gaus_Best_Model not found at loadgbm path")
          print(paste0(loadgbm, "Gaus_Best_Model"))
          load(paste0(loadgbm, "Gaus_Best_Model"))
          Gaus_Best_Model <- "Gaus_Best_Model_Object"
        } # close gaus if
        dir.create(names(samples[i])) # create resvar-named directory for outputs
      } # close if isnull loadgbm

      ####20. Binomial predictions####
      if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {  # do fam1 runs if it's bin only (fam1 bin, gaus (ie fam2) false), or if it's delta & ZI
        #  & exists("Bin_Best_Model")
        grids$Bin_Preds <- predict.gbm(object = get(Bin_Best_Model),
                                       newdata = grids,
                                       n.trees = get(Bin_Best_Model)$gbm.call$best.trees,
                                       type = "response")
        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Binomial predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close if ZI

      ####21. Gaussian predictions####
      if (gaus) {
        #  & exists("Gaus_Best_Model")
        Gaus_Preds <- predict.gbm(object = get(Gaus_Best_Model),
                                  newdata = grids,
                                  n.trees = get(Gaus_Best_Model)$gbm.call$best.trees,
                                  type = "response")

        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Gaussian predictions done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {
          #  & exists("Bin_Best_Model")
          grids$Gaus_Preds <- Gaus_Preds

          ####22. Backtransform logged Gaus to unlogged####
          if (fam2 == "poisson") {
            grids$Gaus_Preds_Unlog <- Gaus_Preds
          } else {
            grids$Gaus_Preds_Unlog <- expm1(Gaus_Preds + 1/2 * sd(get(Gaus_Best_Model)$residuals, na.rm = FALSE) ^ 2)
            # exp for log, expm1 for log1p, L395
          }

          ####23. BIN*positive abundance = final abundance####
          grids$PredAbund <- grids$Gaus_Preds_Unlog * grids$Bin_Preds
        } else { # close gaus yes zi yes run gaus yes zi no
          grids$PredAbund <- Gaus_Preds #if ZI=TRUE, unlog gaus & multiply by bin. Else just use gaus preds.
        } # close ifelse zi
      } else { # if not gaus
        grids$PredAbund <- grids$Bin_Preds # if only doing Bin, preds are just bin preds
      } # close ifelse gaus

      # scale up/down predictions so values in grids pixels relate to the same area sampled in samples
      grids$PredAbund <- grids$PredAbund * samplesGridsAreaScaleFactor

      predabund <- which(colnames(grids) == "PredAbund") # predicted abundance column number for writecsv

      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Final abundance calculated  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####24. Final saves####
      # CSV of Predicted values at each site inc predictor variables' values.
      write.csv(grids, row.names = FALSE, file = paste0("./", names(samples[i]), "/Abundance_Preds_All.csv"))
      # CSV of Predicted values at each site without predictor variables' values.
      # coerce character gridslat/lon into numeric since predabund is given as numeric & you can't mix
      if (is.character(gridslat)) gridslat <- which(colnames(samples) == gridslat)
      if (is.character(gridslon)) gridslon <- which(colnames(samples) == gridslon)
      write.csv(grids[c(gridslat,gridslon,predabund)], row.names = FALSE, file = paste0("./", names(samples[i]), "/Abundance_Preds_only.csv"))
      if (alerts) beep(2) # progress printer, right aligned for visibility
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Output CSVs written     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

      ####25. Unrepresentativeness surface builder####
      # builds doesn't plot surface. If built, plotted by map maker.
      if (RSB) {
        rsbdf_bin <- gbm.rsb(samples, grids, expvarnames, gridslat, gridslon)
        pos_samples <- subset(samples, brv > 0)
        if (gaus & exists("Gaus_Best_Model")) {
          rsbdf_gaus <- gbm.rsb(pos_samples, grids, expvarnames, gridslat, gridslon)
          rsbdf_both <- data.frame(rsbdf_bin, "Unrep_Gaus" = rsbdf_gaus[,"Unrepresentativeness"], "Unrep_Both" = (rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"]))
          write.csv(rsbdf_both, row.names = FALSE, file = paste0("./", names(samples[i]), "/RSB.csv"))
        } else { # not gaus
          write.csv(rsbdf_bin, row.names = FALSE, file = paste0("./", names(samples[i]), "/RSB.csv")) # if not gaus
        } # close if else gaus
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       RSB CSV written       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      } # close if RSB

      ####26. Map maker####
      if (map) {   # generate output image & set parameters
        # png(filename = paste0("./",names(samples[i]),"/PredAbundMap_",names(samples[i]),".png"),
        #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
        # # run gbm.map function with generated parameters
        # gbm.map(x = grids[,gridslon],
        #         y = grids[,gridslat],
        #         z = grids[,predabund],
        #         species = names(samples[i]),
        #         shape = shape, #either autogenerated or set by user so never blank
        #         ...)  # allows gbm.auto's optional terms to be passed to subfunctions:
        # # byx, byy, mapmain, heatcol, mapback, landcol, lejback, legendloc, grdfun, zero, quantile, heatcolours, colournumber
        # dev.off()


        gbm.mapsf(predabund = grids[c(gridslat, gridslon, predabund)],
                  # predabundlon = 2, # Longitude column number.
                  # predabundlat = 1, # Latitude column number.
                  # predabundpreds = 3, # Predicted abundance column number.
                  # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                  # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                  # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                  # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                  # mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                  # googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion factors since
                  # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                  # maptype = "satellite",
                  # darkenproportion = 0, # amount to darken the basemap, 0-1.
                  # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                  shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                  # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                  # colourscale = "viridis", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                  # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                  # heatcolours = c("white", "yellow", "orange","red", "brown4"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                  # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                  studyspecies = names(samples[i]),
                  # plottitle = paste0("Predicted abundance of ", studyspecies),
                  # plotsubtitle = "CPUE", # data %>% distinct(ID) %>% nrow() # 13
                  # legendtitle = "CPUE",
                  # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                  # axisxlabel = "Longitude",
                  # axisylabel = "Latitude",
                  legendposition = c(0.05, 0.18),
                  # fontsize = 12,
                  # fontfamily = "Times New Roman",
                  # filesavename = paste0(lubridate::today(), "_", studyspecies, "_", legendtitle, ".png"),
                  savedir = paste0("./",names(samples[i]), "/"),
                  # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                  # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                  # receivernames = NULL, # vector of names for receivers to be plotted
                  # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                  # recpointscol = "black", # Colour of receiver centrepoint outlines.
                  # recpointsfill = "white", # Colour of receiver centrepoint fills.
                  # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                  # recpointssize = 1, # Size of receiver points.
                  # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                  # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                  # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                  # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                  # reclabcol = "black", # Receiver label text colour.
                  # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                  # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                  # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                  # reclabpad = 0, # Receiver label padding in lines.
                  # reclabrad = 0.15, # Receiver label radius in lines.
                  # reclabbord = 0 # Receiver label border in mm.
                  ...
        )

        # gbm.auto param existing shapefile format, google choice, API, stamen etc passthrough


        if (alerts) beep(2) # progress printer, right aligned for visibility
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Reticulating splines     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    Colour map generated     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

        if (BnW) { # if BnW=TRUE, run again in black & white for journal submission
          # png(filename = paste0("./",names(samples[i]),"/PredAbundMap_BnW_",names(samples[i]),".png"),
          #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
          # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
          # gbm.map(x = grids[,gridslon],
          #         y = grids[,gridslat],
          #         z = grids[,predabund],
          #         species = names(samples[i]),
          #         shape = shape, #either autogenerated or set by user so never blank
          #         landcol = grey.colors(1, start = 0.8, end = 0.8), #light grey. 0=black 1=white
          #         mapback = "white",
          #         heatcolours = grey.colors(8, start = 1, end = 0),
          #         savedir = savedir, # passes to gbm.map's ... which passes to gbm.basemap
          #         ...)  # allows gbm.auto's optional terms to be passed to subfunctions:
          # # byx, byy, mapmain, heatcol, mapback, landcol, lejback, legendloc, grdfun, zero, quantile, heatcolours, colournumber
          # dev.off()

          gbm.mapsf(predabund = grids[c(gridslat, gridslon, predabund)],
                    # predabundlon = 2, # Longitude column number.
                    # predabundlat = 1, # Latitude column number.
                    # predabundpreds = 3, # Predicted abundance column number.
                    # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                    # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                    # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                    # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                    mapsource = "stamen", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                    googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion factors since
                    # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                    maptype = "toner-lite",
                    # darkenproportion = 0, # amount to darken the basemap, 0-1.
                    # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                    shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                    # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                    colourscale = "gradient", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                    # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                    heatcolours = c("white", "grey80", "grey60","grey40", "grey20", "black"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                    # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                    studyspecies = names(samples[i]),
                    # plottitle = paste0("Predicted abundance of ", studyspecies),
                    # plotsubtitle = "CPUE", # data %>% distinct(ID) %>% nrow() # 13
                    # legendtitle = "CPUE",
                    # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                    # axisxlabel = "Longitude",
                    # axisylabel = "Latitude",
                    legendposition = c(0.05, 0.18),
                    # fontsize = 12,
                    # fontfamily = "Times New Roman",
                    filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_CPUE_BnW.png"),
                    savedir = paste0("./",names(samples[i]), "/")
                    # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                    # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                    # receivernames = NULL, # vector of names for receivers to be plotted
                    # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                    # recpointscol = "black", # Colour of receiver centrepoint outlines.
                    # recpointsfill = "white", # Colour of receiver centrepoint fills.
                    # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                    # recpointssize = 1, # Size of receiver points.
                    # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                    # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                    # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                    # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                    # reclabcol = "black", # Receiver label text colour.
                    # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                    # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                    # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                    # reclabpad = 0, # Receiver label padding in lines.
                    # reclabrad = 0.15, # Receiver label radius in lines.
                    # reclabbord = 0 # Receiver label border in mm.
          )


          if (alerts) beep(2)  # progress printer, right aligned for visibility
          print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Black & white map generated XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        } # close & save plotting device & close BnW optional

        if (RSB) { # if RSB called, plot that surface separately
          linear01seq <- seq(from = 0, to = 1, length.out = 9) #linear sequence from 0:1, 9 bins
          exp01seq <- expm1(4*linear01seq)/expm1(4) # exponentiate to change shape then scale back to 1

          if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {
            #  & exists("Bin_Best_Model")
            # png(filename = paste0("./",names(samples[i]),"/RSB_Map_Bin_",names(samples[i]),".png"),
            #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
            # gbm.map(x = grids[,gridslon],
            #         y = grids[,gridslat],
            #         z = rsbdf_bin[,"Unrepresentativeness"],
            #         mapmain = "Unrepresentativeness: ",
            #         species = names(samples[i]),
            #         legendtitle = "UnRep 0-1",
            #         shape = shape, #either autogenerated or set by user so never blank
            #         # breaks = expm1(breaks.grid(log(2000), ncol = 8, zero = TRUE))/2000) #old failing breaks
            #         breaks = exp01seq)
            # dev.off() #high value log breaks mean first ~5 values cluster near 0 for high
            # # res there, but high values captures in the last few bins.

            gbm.mapsf(predabund = data.frame(Latitude = grids[, gridslat],
                                             Longitude = grids[, gridslon],
                                             Unrepresentativeness = rsbdf_bin[,"Unrepresentativeness"]),
                      # predabundlon = 2, # Longitude column number.
                      # predabundlat = 1, # Latitude column number.
                      # predabundpreds = 3, # Predicted abundance column number.
                      # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                      # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                      # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                      # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                      # mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                      # googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion factors since
                      # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                      # maptype = "satellite",
                      # darkenproportion = 0, # amount to darken the basemap, 0-1.
                      # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                      shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                      # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                      colourscale = "gradient", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                      # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                      heatcolours = c("white", "yellow", "orange","red", "brown4"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                      colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                      studyspecies = names(samples[i]),
                      plottitle = paste0("Unrepresentativeness of samples data for ", names(samples[i])),
                      plotsubtitle = "Unrepresentativeness", # data %>% distinct(ID) %>% nrow() # 13
                      legendtitle = "UnRep 0-1",
                      # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                      # axisxlabel = "Longitude",
                      # axisylabel = "Latitude",
                      legendposition = c(0.05, 0.18),
                      # fontsize = 12,
                      # fontfamily = "Times New Roman",
                      filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_RSB_Map_Bin.png"),
                      savedir = paste0("./",names(samples[i]), "/")
                      # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                      # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                      # receivernames = NULL, # vector of names for receivers to be plotted
                      # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                      # recpointscol = "black", # Colour of receiver centrepoint outlines.
                      # recpointsfill = "white", # Colour of receiver centrepoint fills.
                      # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                      # recpointssize = 1, # Size of receiver points.
                      # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                      # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                      # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                      # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                      # reclabcol = "black", # Receiver label text colour.
                      # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                      # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                      # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                      # reclabpad = 0, # Receiver label padding in lines.
                      # reclabrad = 0.15, # Receiver label radius in lines.
                      # reclabbord = 0 # Receiver label border in mm.
            )

            if (alerts) beep(2) # progress printer, right aligned for visibility
            print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  Colour RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
          } # close if zi bin

          if (gaus) {
            #  & exists("Gaus_Best_Model")
            # png(filename = paste0("./",names(samples[i]),"/RSB_Map_Gaus_",names(samples[i]),".png"),
            #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
            # gbm.map(x = grids[,gridslon],
            #         y = grids[,gridslat],
            #         z = rsbdf_gaus[,"Unrepresentativeness"],
            #         mapmain = "Unrepresentativeness: ",
            #         species = names(samples[i]),
            #         legendtitle = "UnRep 0-1",
            #         shape = shape, #either autogenerated or set by user so never blank
            #         breaks = exp01seq)
            # dev.off()

            gbm.mapsf(predabund = data.frame(Latitude = grids[, gridslat],
                                             Longitude = grids[, gridslon],
                                             Unrepresentativeness = rsbdf_gaus[,"Unrepresentativeness"]),
                      # predabundlon = 2, # Longitude column number.
                      # predabundlat = 1, # Latitude column number.
                      # predabundpreds = 3, # Predicted abundance column number.
                      # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                      # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                      # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                      # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                      # mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                      # googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion factors since
                      # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                      # maptype = "satellite",
                      # darkenproportion = 0, # amount to darken the basemap, 0-1.
                      # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                      shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                      # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                      # colourscale = "viridis", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                      # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                      # heatcolours = c("white", "yellow", "orange","red", "brown4"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                      # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                      studyspecies = names(samples[i]),
                      plottitle = paste0("Unrepresentativeness of samples data for ", names(samples[i])),
                      plotsubtitle = "Unrepresentativeness", # data %>% distinct(ID) %>% nrow() # 13
                      legendtitle = "UnRep 0-1",
                      # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                      # axisxlabel = "Longitude",
                      # axisylabel = "Latitude",
                      legendposition = c(0.05, 0.18),
                      # fontsize = 12,
                      # fontfamily = "Times New Roman",
                      filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_RSB_Map_Gaus.png"),
                      savedir = paste0("./",names(samples[i]), "/"),
                      # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                      # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                      # receivernames = NULL, # vector of names for receivers to be plotted
                      # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                      # recpointscol = "black", # Colour of receiver centrepoint outlines.
                      # recpointsfill = "white", # Colour of receiver centrepoint fills.
                      # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                      # recpointssize = 1, # Size of receiver points.
                      # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                      # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                      # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                      # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                      # reclabcol = "black", # Receiver label text colour.
                      # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                      # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                      # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                      # reclabpad = 0, # Receiver label padding in lines.
                      # reclabrad = 0.15, # Receiver label radius in lines.
                      # reclabbord = 0 # Receiver label border in mm.
                      ...
            )

            if (alerts) beep(2) # progress printer, right aligned for visibility
            print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
          } # close gaus map

          if (ZI & gaus) {
            #  & exists("Gaus_Best_Model")
            # png(filename = paste0("./",names(samples[i]),"/RSB_Map_Both_",names(samples[i]),".png"),
            #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
            # gbm.map(x = grids[,gridslon],
            #         y = grids[,gridslat],
            #         z = rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"],
            #         mapmain = "Unrepresentativeness: ",
            #         species = names(samples[i]),
            #         legendtitle = "UnRep 0-2",
            #         shape = shape, #either autogenerated or set by user so never blank
            #         breaks = exp01seq)
            # dev.off()

            gbm.mapsf(predabund = data.frame(Latitude = grids[, gridslat],
                                             Longitude = grids[, gridslon],
                                             Unrepresentativeness = rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"]),
                      # predabundlon = 2, # Longitude column number.
                      # predabundlat = 1, # Latitude column number.
                      # predabundpreds = 3, # Predicted abundance column number.
                      # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                      # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                      # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                      # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                      # mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                      # googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion factors since
                      # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                      # maptype = "satellite",
                      # darkenproportion = 0, # amount to darken the basemap, 0-1.
                      # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                      shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                      # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                      # colourscale = "viridis", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                      # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                      # heatcolours = c("white", "yellow", "orange","red", "brown4"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                      # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                      studyspecies = names(samples[i]),
                      plottitle = paste0("Unrepresentativeness of samples data for ", names(samples[i])),
                      plotsubtitle = "Unrepresentativeness", # data %>% distinct(ID) %>% nrow() # 13
                      legendtitle = "UnRep 0-2",
                      # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                      # axisxlabel = "Longitude",
                      # axisylabel = "Latitude",
                      legendposition = c(0.05, 0.18),
                      # fontsize = 12,
                      # fontfamily = "Times New Roman",
                      filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_RSB_Map_Combo.png"),
                      savedir = paste0("./",names(samples[i]), "/"),
                      # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                      # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                      # receivernames = NULL, # vector of names for receivers to be plotted
                      # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                      # recpointscol = "black", # Colour of receiver centrepoint outlines.
                      # recpointsfill = "white", # Colour of receiver centrepoint fills.
                      # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                      # recpointssize = 1, # Size of receiver points.
                      # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                      # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                      # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                      # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                      # reclabcol = "black", # Receiver label text colour.
                      # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                      # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                      # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                      # reclabpad = 0, # Receiver label padding in lines.
                      # reclabrad = 0.15, # Receiver label radius in lines.
                      # reclabbord = 0 # Receiver label border in mm.
                      ...
            )

            if (alerts) beep(2) # progress printer, right aligned for visibility
            print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Colour RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
          } # close both map

          if (BnW) {     # if BnW=TRUE, do again for b&w
            if (fam1 == "bernoulli" & (!gaus | (gaus & ZI))) {
              #  & exists("Bin_Best_Model")
              # png(filename = paste0("./",names(samples[i]),"/RSB_Map_BnW_Bin_",names(samples[i]),".png"),
              #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
              # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
              # gbm.map(x = grids[,gridslon],
              #         y = grids[,gridslat],
              #         z = rsbdf_bin[,"Unrepresentativeness"],
              #         mapmain = "Unrepresentativeness: ",
              #         mapback = "white",
              #         species = names(samples[i]),
              #         heatcolours = grey.colors(8, start = 1, end = 0), #default 8 greys
              #         ####BUG:setting heatcolours & colournumber overrides this####
              #         landcol = grey.colors(1, start = 0.8, end = 0.8), #light grey. 0=black 1=white
              #         legendtitle = "UnRep 0-1",
              #         shape = shape, #either autogenerated or set by user so never blank
              #         breaks = exp01seq)
              # dev.off()

              gbm.mapsf(predabund = data.frame(Latitude = grids[, gridslat],
                                               Longitude = grids[, gridslon],
                                               Unrepresentativeness = rsbdf_bin[,"Unrepresentativeness"]),
                        # predabundlon = 2, # Longitude column number.
                        # predabundlat = 1, # Latitude column number.
                        # predabundpreds = 3, # Predicted abundance column number.
                        # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                        # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                        # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                        # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                        mapsource = "stamen", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                        googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion factors since
                        # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                        maptype = "toner-lite",
                        # darkenproportion = 0, # amount to darken the basemap, 0-1.
                        # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                        shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                        # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                        colourscale = "gradient", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                        # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                        heatcolours = c("white", "grey80", "grey60","grey40", "grey20", "black"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                        # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                        studyspecies = names(samples[i]),
                        plottitle = paste0("Unrepresentativeness of samples data for ", names(samples[i])),
                        plotsubtitle = "Unrepresentativeness", # data %>% distinct(ID) %>% nrow() # 13
                        legendtitle = "UnRep 0-1",
                        # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                        # axisxlabel = "Longitude",
                        # axisylabel = "Latitude",
                        legendposition = c(0.05, 0.18),
                        # fontsize = 12,
                        # fontfamily = "Times New Roman",
                        filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_RSB_Map_Bin_BnW.png"),
                        savedir = paste0("./",names(samples[i]), "/")
                        # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                        # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                        # receivernames = NULL, # vector of names for receivers to be plotted
                        # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                        # recpointscol = "black", # Colour of receiver centrepoint outlines.
                        # recpointsfill = "white", # Colour of receiver centrepoint fills.
                        # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                        # recpointssize = 1, # Size of receiver points.
                        # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                        # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                        # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                        # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                        # reclabcol = "black", # Receiver label text colour.
                        # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                        # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                        # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                        # reclabpad = 0, # Receiver label padding in lines.
                        # reclabrad = 0.15, # Receiver label radius in lines.
                        # reclabbord = 0 # Receiver label border in mm.
              )

              if (alerts) beep(2) # progress printer, right aligned for visibility
              print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     B&W RSB bin map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
            } # close bin RSB

            if (gaus) {
              #  & exists("Gaus_Best_Model")
              # png(filename = paste0("./",names(samples[i]),"/RSB_Map_BnW_Gaus_",names(samples[i]),".png"),
              #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
              # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
              # gbm.map(x = grids[,gridslon],
              #         y = grids[,gridslat],
              #         z = rsbdf_gaus[,"Unrepresentativeness"],
              #         mapmain = "Unrepresentativeness: ",
              #         mapback = "white",
              #         species = names(samples[i]),
              #         heatcolours = grey.colors(8, start = 1, end = 0),
              #         landcol = grey.colors(1, start = 0.8, end = 0.8),
              #         legendtitle = "UnRep 0-1",
              #         shape = shape, #either autogenerated or set by user so never blank
              #         breaks = exp01seq)
              # dev.off()

              gbm.mapsf(predabund = data.frame(Latitude = grids[, gridslat],
                                               Longitude = grids[, gridslon],
                                               Unrepresentativeness = rsbdf_gaus[,"Unrepresentativeness"]),
                        # predabundlon = 2, # Longitude column number.
                        # predabundlat = 1, # Latitude column number.
                        # predabundpreds = 3, # Predicted abundance column number.
                        # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                        # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                        # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                        # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                        mapsource = "stamen", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                        googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion factors since
                        # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                        maptype = "toner-lite",
                        # darkenproportion = 0, # amount to darken the basemap, 0-1.
                        # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                        shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                        # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                        colourscale = "gradient", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                        # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                        heatcolours = c("white", "grey80", "grey60","grey40", "grey20", "black"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                        # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                        studyspecies = names(samples[i]),
                        plottitle = paste0("Unrepresentativeness of samples data for ", names(samples[i])),
                        plotsubtitle = "Unrepresentativeness", # data %>% distinct(ID) %>% nrow() # 13
                        legendtitle = "UnRep 0-1",
                        # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                        # axisxlabel = "Longitude",
                        # axisylabel = "Latitude",
                        legendposition = c(0.05, 0.18),
                        # fontsize = 12,
                        # fontfamily = "Times New Roman",
                        filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_RSB_Map_Gaus_BnW.png"),
                        savedir = paste0("./",names(samples[i]), "/")
                        # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                        # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                        # receivernames = NULL, # vector of names for receivers to be plotted
                        # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                        # recpointscol = "black", # Colour of receiver centrepoint outlines.
                        # recpointsfill = "white", # Colour of receiver centrepoint fills.
                        # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                        # recpointssize = 1, # Size of receiver points.
                        # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                        # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                        # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                        # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                        # reclabcol = "black", # Receiver label text colour.
                        # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                        # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                        # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                        # reclabpad = 0, # Receiver label padding in lines.
                        # reclabrad = 0.15, # Receiver label radius in lines.
                        # reclabbord = 0 # Receiver label border in mm.
              )

              if (alerts) beep(2) # progress printer, right aligned for visibility
              print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB Gaus map done    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
            } # close gaus RSB

            if (ZI & gaus) {
              #  & exists("Gaus_Best_Model")
              # png(filename = paste0("./",names(samples[i]),"/RSB_Map_BnW_Both_",names(samples[i]),".png"),
              #     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
              # par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
              # gbm.map(x = grids[,gridslon],
              #         y = grids[,gridslat],
              #         z = rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"],
              #         mapmain = "Unrepresentativeness: ",
              #         mapback = "white",
              #         species = names(samples[i]),
              #         heatcolours = grey.colors(8, start = 1, end = 0),
              #         landcol = grey.colors(1, start = 0.8, end = 0.8),
              #         legendtitle = "UnRep 0-2",
              #         shape = shape, #either autogenerated or set by user so never blank
              #         breaks = exp01seq)
              # dev.off()

              gbm.mapsf(predabund = data.frame(Latitude = grids[, gridslat],
                                               Longitude = grids[, gridslon],
                                               Unrepresentativeness = rsbdf_bin[,"Unrepresentativeness"] + rsbdf_gaus[,"Unrepresentativeness"]),
                        # predabundlon = 2, # Longitude column number.
                        # predabundlat = 1, # Latitude column number.
                        # predabundpreds = 3, # Predicted abundance column number.
                        # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
                        # trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
                        # scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
                        # gmapsAPI = NULL, # enter your Google maps API here, quoted character string
                        mapsource = "stamen", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
                        googlemap = FALSE, # If pulling basemap from Google maps, this sets expansion factors since
                        # # Google Maps tiling zoom setup doesn't align to myLocation extents.
                        maptype = "toner-lite",
                        # darkenproportion = 0, # amount to darken the basemap, 0-1.
                        # mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
                        shape = shape, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
                        # expandfactor = 0, # extents expansion factor for basemap. default was 1.6
                        colourscale = "gradient", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
                        # colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
                        heatcolours = c("white", "grey80", "grey60","grey40", "grey20", "black"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
                        # colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
                        studyspecies = names(samples[i]),
                        plottitle = paste0("Unrepresentativeness of samples data for ", names(samples[i])),
                        plotsubtitle = "Unrepresentativeness", # data %>% distinct(ID) %>% nrow() # 13
                        legendtitle = "UnRep 0-1",
                        # plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
                        # axisxlabel = "Longitude",
                        # axisylabel = "Latitude",
                        legendposition = c(0.05, 0.18),
                        # fontsize = 12,
                        # fontfamily = "Times New Roman",
                        filesavename = paste0(lubridate::today(), "_", names(samples[i]), "_RSB_Map_Combo_BnW.png"),
                        savedir = paste0("./",names(samples[i]), "/")
                        # receiverlats = NULL, # vector of latitudes for receivers to be plotted
                        # receiverlons = NULL, # vector of longitudes for receivers to be plotted
                        # receivernames = NULL, # vector of names for receivers to be plotted
                        # receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
                        # recpointscol = "black", # Colour of receiver centrepoint outlines.
                        # recpointsfill = "white", # Colour of receiver centrepoint fills.
                        # recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
                        # recpointssize = 1, # Size of receiver points.
                        # recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
                        # recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
                        # recbuffill = "grey", # Colour of the receiver buffer circle fills.
                        # recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
                        # reclabcol = "black", # Receiver label text colour.
                        # reclabfill = NA, # Receiver label fill colour, NA for no fill.
                        # reclabnudgex = 0, # Receiver label offset nudge in X dimension.
                        # reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
                        # reclabpad = 0, # Receiver label padding in lines.
                        # reclabrad = 0.15, # Receiver label radius in lines.
                        # reclabbord = 0 # Receiver label border in mm.
              )

              if (alerts) beep(2) # progress printer, right aligned for visibility
              print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    B&W RSB combo map done   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
            } # close gaus (&combo) B&W RSB
          } # close BnW RSBs
        } # close RSB mapper
      } # close Map Maker
    } #close !isnull grids option from above section 19
    print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Grids/maps/everything done  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
  } # close for i in resvar response variable loop
  gc() # Force R to release memory it is no longer using
  options(error = NULL) # reset error options to default
  if (alerts) beep(8)  # final user notification, then close the function
}

