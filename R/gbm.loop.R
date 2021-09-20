#' Calculate Coefficient Of Variation surfaces for gbm.auto predictions
#'
#' Bagging introduces stochasticity which can result in sizeable variance in
#' output predictions by gbm.auto for small datasets. This function runs a user-
#' specified number of loops through the same gbm.auto parameter combinations
#' and calculates the Coefficient Of Variation in the predicted abundance scores
#'  for each site aka cell. This can be mapped, to spatially demonstrate the
#'  output variance range.
#'
#' @param loops The number of loops required, integer.
#' @param savedir Save outputs to a temporary directory (default) else change to
#'  current directory e.g. "/home/me/folder". Do not use getwd() here.
#' @param savecsv Save coefficients of variation in simple & extended format.
#' @param calcpreds Calculate coefficients of variation of predicted abundance?
#' @param varmap Create a map of the coefficients of variation outputs?
#' @param measure Map legend, coefficients of variation of what? Default CPUE.
#' @param cleanup Remove gbm.auto-generated directory each loop? Default FALSE.
#' @param grids See gbm.auto help for all subsequent params.
#' @param samples See gbm.auto help.
#' @param expvar See gbm.auto help.
#' @param resvar See gbm.auto help.
#' @param tc See gbm.auto help.
#' @param lr See gbm.auto help.
#' @param bf See gbm.auto help.
#' @param n.trees See gbm.auto help.
#' @param ZI See gbm.auto help. Choose one.
#' @param fam1 See gbm.auto help. Choose one.
#' @param fam2 See gbm.auto help. Choose one.
#' @param simp See gbm.auto help.
#' @param gridslat See gbm.auto help.
#' @param gridslon See gbm.auto help.
#' @param multiplot See gbm.auto help. Default False
#' @param cols See gbm.auto help.
#' @param linesfiles See gbm.auto help; TRUE or linesfiles calculations fail.
#' @param smooth See gbm.auto help.
#' @param savegbm See gbm.auto help.
#' @param loadgbm See gbm.auto help.
#' @param varint See gbm.auto help.
#' @param map See gbm.auto help.
#' @param shape See gbm.auto help.
#' @param RSB See gbm.auto help.
#' @param BnW See gbm.auto help.
#' @param alerts See gbm.auto help; default FALSE as frequent use can crash
#' RStudio.
#' @param pngtype See gbm.auto help. Choose one.
#' @param gaus See gbm.auto help.
#' @param MLEvaluate See gbm.auto help.
#' @param runautos Run gbm.autos, default TRUE, turn off to only collate
#' numbered-folder results.
#' @param Min.Inf Dummy param for package testing for CRAN, ignore.
#' @param ... Additional params for gbm.auto sub-functions including gbm.step.
#'
#' @return Returns a data frame of lat, long, 1 predicted abundance per loop,
#' and a final variance score per cell.
#'
#' @export
#' @importFrom beepr beep
#' @importFrom grDevices dev.off grey.colors png
#' @importFrom graphics axis barplot legend lines mtext par text
#' @importFrom stats var
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \donttest{
#' # Not run: downloads and saves external data.
#' library("gbm.auto")
#' data(grids) # load grids
#' data(samples) # load samples
#' gbmloopexample <- gbm.loop(loops = 2, samples = samples,
#' grids = grids, expvar = c(4:10), resvar = 11, simp = F)
#' }
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
gbm.loop <- function(loops = 10, # the number of loops required, integer
                     savedir = tempdir(), # save outputs to a temporary directory (default) else
                     # change to current directory e.g. "/home/me/folder". Do not use getwd() here.
                     savecsv = TRUE, # save the coefficients of variation in simple & extended format
                     calcpreds = TRUE, # calculate coefficients of variation of predicted abundance?
                     varmap = TRUE, # create a map of the coefficients of variation outputs?
                     measure = "CPUE", # map legend, coefficients of variation of what? Default CPUE
                     cleanup = FALSE, # remove gbm.auto-generated directory each loop?
                     grids = NULL,         # explantory data to predict to. Import with (e.g.)
                     # read.csv and specify object name.
                     samples,  # explanatory and response variables to predict from.
                     # Keep col names short, no odd characters, starting numerals or terminal periods
                     # Spaces may be converted to periods in directory names, underscores won't.
                     # Can be a subset
                     expvar,               # list of column numbers of explanatory variables in
                     # 'samples', expected e.g. c(1,35,67,etc.). No default
                     resvar,               # column number of response variable (e.g. CPUE) in
                     # samples. Expected, e.g. 12. No default. Column name should be species name
                     tc = c(2),            # permutations of tree complexity allowed, can be a
                     # vector with the largest sized number no larger than the number of
                     # explanatory variables e.g. c(2,7), or a list of 2 single numbers or vectors,
                     # the first to be passed to the binary BRT, the second to the Gaussian, e.g.
                     # tc = list(c(2,6), 2) or list(6, c(2,6))
                     lr = c(0.01),   # permutations of learning rate allowed. Can be a
                     # vector or a list of 2 single numbers or vectors, the first to be passed to
                     # the binary BRT, the second to the Gaussian, e.g.
                     # lr = list(c(0.01,0.02),0.0001) or list(0.01,c(0.001, 0.0005))
                     bf = 0.5,             # permutations of bag fraction allowed, can be single
                     # number, vector or list, per tc and lr
                     n.trees = 50,         # from gbm.step, number of initial trees to fit. Can be
                     # single or list but not vector i.e. list(fam1, fam2)
                     ZI = "CHECK", # Are data zero-inflated? "CHECK"/FALSE/TRUE.
                     # Choose one.
                     # TRUE: delta BRT, log-normalised Gaus, reverse log-norm and bias corrected.
                     # FALSE: do Gaussian only, no log-normalisation.
                     # CHECK: Tests data for you. Default is TRUE.
                     fam1 = c("bernoulli", "binomial", "poisson", "laplace", "gaussian"),
                     # probability distribution family for 1st part of delta process, defaults to
                     # "bernoulli",
                     fam2 = c("gaussian", "bernoulli", "binomial", "poisson", "laplace"),
                     # probability distribution family for 2nd part of delta process, defaults to
                     # "gaussian",
                     simp = TRUE,          # try simplfying best BRTs?
                     gridslat = 2,         # column number for latitude in 'grids'
                     gridslon = 1,         # column number for longitude in 'grids'
                     multiplot = FALSE,     # create matrix plot of all line files? Default false
                     # turn off if large number of expvars causes an error due to margin size problems.
                     cols = grey.colors(1,1,1), # barplot colour vector. Assignment in order of
                     # explanatory variables. Default 1*white: white bars black borders. '1*' repeats
                     linesfiles = TRUE,   # save individual line plots' data as csv's?
                     smooth = FALSE,       # apply a smoother to the line plots? Default FALSE
                     savegbm = FALSE,       # save gbm objects and make available in environment after running? Open with load("Bin_Best_Model")
                     loadgbm = NULL,       # relative or absolute location of folder containing
                     # Bin_Best_Model and Gaus_Best_Model. If set will skip BRT calculations and do
                     # predicted maps and CSVs. Default NULL, character vector, "./" for working directory
                     varint = FALSE,        # calculate variable interactions? Default:TRUE, FALSE
                     # for error "contrasts can be applied only to factors with 2 or more levels"
                     map = TRUE,           # save abundance map png files?
                     shape = NULL,      # set coast shapefile, else downloaded and autogenerated
                     RSB = FALSE,           # run Unrepresentativeness surface builder?
                     BnW = FALSE,           # repeat maps in black and white e.g. for print journals
                     alerts = FALSE,        # play sounds to mark progress steps
                     pngtype = c("cairo-png", "quartz", "Xlib"), # file-type for png files,
                     # alternatively try "quartz" on Mac. Choose one.
                     gaus = TRUE,          # do Gaussian runs as well as Bin? Default TRUE.
                     MLEvaluate = TRUE,    # do machine learning evaluation metrics & plots? Default TRUE
                     # brv = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
                     # grv = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
                     # Bin_Preds = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
                     # Gaus_Preds = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
                     runautos = TRUE,      # run gbm.autos, default TRUE, turn off to only collate numbered-folder results
                     Min.Inf = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
                     ...) {

  # Generalised Boosting Model / Boosted Regression Tree process chain automater
  # Simon Dedman, 2012-8 simondedman@gmail.com github.com/SimonDedman/gbm.auto

  ####TODO####
  # See how many loops until things stabilise, i.e. variance decreases, average smooths out, etc, is that even logical?
  # Yes. min max av var should be similar between e.g. 1,000,000 loops & 1,000,001, but less likely between 1 & 2.
  # But based on what though? Just do a line of x:loop# vs y: minmin/maxmax/avav/avvar?
  # when change in variance from 1:2 to 1:3 to 1:n drops below a percentage threshold?
  # Fix csvs colnames, see https://github.com/SimonDedman/gbm.auto/issues/37
  # for factorial variables, need to change from lines to bars
  # Runautos doesn't work: binbars.df & gausbars.df are created & incrementally grown within the autos loop then accessed afterwards.
  # Would need to do this separately somehow, possibly a separate loop to pull these data from a source file csv?

  # utils::globalVariables("Min.Inf") # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals

  # if (alerts) if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
  # if (alerts) require(beepr)
  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldwd <- getwd()
  oldoptions <- options()
  on.exit(par(oldpar))
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(options(oldoptions), add = TRUE)
  setwd(savedir)
  if (alerts) options(error = function() {
    beep(9)
    graphics.off()})  # give warning noise if it fails

  fam1 <- match.arg(fam1) # populate object from function argument in proper way
  fam2 <- match.arg(fam2)
  pngtype <- match.arg(pngtype)

  binbars.df <- data.frame(var = rep(NA, length(expvar)),
                           rel.inf = rep(NA, length(expvar)))
  gausbars.df <- binbars.df # blank dataframes for bin & gaus bars data
  report.df <- data.frame(BinCV = rep(NA, length(loops)),
                          AUC = rep(NA, length(loops)),
                          GausCV = rep(NA, length(loops)))
  if (calcpreds) var.df <- grids[,c(gridslon, gridslat)] # create df with just lat & longs
  if (runautos) { # run gbm.autos unless turned off
    for (i in 1:loops) { # loop through all gbm.autos
      # subsavedir <- paste0("./", i)
      subsavedir <- paste0(savedir, "/", i)
      dir.create(subsavedir) # create i'th folder
      setwd(subsavedir) # move to it
      gbm.auto(grids = grids, # run i'th gbm.auto
               samples = samples,
               expvar = expvar,
               resvar = resvar,
               tc = tc,
               lr = lr,
               bf = bf,
               n.trees = n.trees,
               ZI = ZI,
               fam1 = fam1,
               fam2 = fam2,
               simp = simp,
               gridslat = gridslat,
               gridslon = gridslon,
               multiplot = multiplot,
               cols = cols,
               linesfiles = linesfiles,
               smooth = smooth,
               savedir = subsavedir,
               savegbm = savegbm,
               loadgbm = loadgbm,
               varint = varint,
               map = map,
               shape = shape,
               RSB = RSB,
               BnW = BnW,
               alerts = alerts,
               pngtype = pngtype,
               gaus = gaus,
               MLEvaluate = MLEvaluate,
               ...) # accept other gbm.auto values than these basics

      setwd(paste0("./", colnames(samples[resvar]))) # set wd to species named subfolder

      if (file.exists("Binary BRT Variable contributions.csv")) {
        binbarstmp <- read.csv("Binary BRT Variable contributions.csv") # temp container for bin bars
        if (i == 1) {binbars.df <- binbarstmp} else {# csv file to df unless df exists
          binbars.df <- rbind(binbars.df, binbarstmp)} # if so add to bottom of existing
        bin = TRUE} else bin = FALSE

      # loop thru variables name linesfiles e.g. Bin_Best_line_[varname].csv
      # adding i'th loop's values as new column
      if (bin) for (j in colnames(samples[expvar])) {
        if (!file.exists(paste0("Bin_Best_line_", j, ".csv"))) {tmp <- data.frame(x = rep(NA,100), y = rep(NA,100))}
        #if file not created because simp, populate with 0s
        if (file.exists(paste0("Bin_Best_line_", j, ".csv"))) {tmp <- read.csv(paste0("Bin_Best_line_", j, ".csv"))}
        #else use values

        colnames(tmp)[2] <- paste0("Loop",i)
        if (i == 1) {assign(paste0("binline_", j), tmp)
        } else {
          assign(paste0("binline_", j), cbind(get(paste0("binline_", j)),
                                              tmp[,2]))
          if (is.na(get(paste0("binline_", j))[1,1])) { #if the first cell is NA (all 1st col, x, is na)
            assign(paste0("binline_", j), #rebuild same obj as df
                   data.frame(x = tmp[,1], #start with this loop's x values, hopefully not NA also
                              get(paste0("binline_", j))[,2:(i + 1)]))} #then add the remainder of the existing obj cols
        }}

      if (file.exists("Gaussian BRT Variable contributions.csv")) {
        gausbarstmp <- read.csv("Gaussian BRT Variable contributions.csv") # temp container for Gaus lines
        if (i == 1) {gausbars.df <- gausbarstmp} else {
          gausbars.df <- rbind(gausbars.df, gausbarstmp)}
        gaus = TRUE} else gaus = FALSE

      if (gaus) for (k in colnames(samples[expvar])) {
        if (!file.exists(paste0("Gaus_Best_line_", k, ".csv"))) {tmp <- data.frame(x = rep(NA,100), y = rep(NA,100))}
        #if the first loop is simplified then the first col of gausline will be NAs which should be the X for the linefiles
        #else use existing csv file, 2 columns
        if (file.exists(paste0("Gaus_Best_line_", k, ".csv"))) {tmp <- read.csv(paste0("Gaus_Best_line_", k, ".csv"))}
        colnames(tmp)[2] <- paste0("Loop",i)
        if (i == 1) {assign(paste0("gausline_", k), tmp)
        } else {
          assign(paste0("gausline_", k), cbind(get(paste0("gausline_", k)),
                                               tmp[,2]))
          if (is.na(get(paste0("gausline_", k))[1,1])) { #if the first cell is NA (all 1st col, x, is na)
            assign(paste0("gausline_", k), #rebuild same obj as df
                   data.frame(x = tmp[,1], #start with this loop's x values, hopefully not NA also
                              get(paste0("gausline_", k))[,2:(i + 1)]))} #then add the remainder of the existing obj cols
          #column cbound but not named. Can name as string "col name" = 1:10, or
          #objectname ColName = 1:10 but not formulaicly paste0("Col","Name") = 1:10
          #or anything evaluated e.g. colnames(tmp)[2] = tmp[,2]
          #colnames(paste0("gausline_", k))[i + 1] <- paste0("loop", i) #rename last column (loop# + 1)
        }}
      if (!file.exists("Abundance_Preds_only.csv")) calcpreds = FALSE
      if (calcpreds) {predtmp <- read.csv("Abundance_Preds_only.csv") # temp container for latest preds
      var.df <- cbind(var.df, predtmp[,3]) # cbind preds to existing lat/longs or other preds
      colnames(var.df)[2 + i] <- paste0("Loop_", i)} # label newly added preds column

      #Collect report CV & AUC scores
      reporttmp <- read.csv("Report.csv") # temp container for bin bars

      if ("Best.Binary.BRT" %in% colnames(reporttmp)) {
        bincvtmp <- as.character(reporttmp$Best.Binary.BRT[2])
        bincvspltmp <- strsplit(bincvtmp, "Model CV score: ")
        bincvsplnumtmp <- as.numeric(bincvspltmp[[1]][2])
        report.df[i,1] <- bincvsplnumtmp # copy BinCV score from this loop's report to allreport

        auctmp <- as.character(reporttmp$Best.Binary.BRT[3])
        aucspltmp <- strsplit(auctmp, "Training data AUC score: ")
        aucsplnumtmp <- as.numeric(aucspltmp[[1]][2])
        report.df[i,2] <- aucsplnumtmp} # copy AUC score from this loop's report to allreport

      if ("Best.Gaussian.BRT" %in% colnames(reporttmp)) {
        gauscvtmp <- as.character(reporttmp$Best.Gaussian.BRT[2])
        gauscvspltmp <- strsplit(gauscvtmp, "Model CV score: ")
        gauscvsplnumtmp <- as.numeric(gauscvspltmp[[1]][2])
        report.df[i,3] <- gauscvsplnumtmp} # copy GausCV score from this loop's report to allreport

      setwd("../../") # move back up to root folder
      if (cleanup) unlink(i, recursive = TRUE)
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Loop ",i," complete        XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
    } # close i loop & go to the next i
  } # close runautos if optional

  #runautos fix####
  #move runautos binbars.df gausbars.df code to here, in separate loops
  # binbars.df <- data.frame(var = rep(NA, length(expvar)),
  #                          rel.inf = rep(NA, length(expvar)))
  # gausbars.df <- binbars.df # blank dataframes for bin & gaus bars data
  #
  # if (file.exists("Binary BRT Variable contributions.csv")) {
  #   binbarstmp <- read.csv("Binary BRT Variable contributions.csv") # temp container for bin bars
  #   if (i == 1) {binbars.df <- binbarstmp} else {# csv file to df unless df exists
  #     binbars.df <- rbind(binbars.df, binbarstmp)} # if so add to bottom of existing
  #   bin = TRUE} else bin = FALSE
  #
  # if (file.exists("Gaussian BRT Variable contributions.csv")) {
  #   gausbarstmp <- read.csv("Gaussian BRT Variable contributions.csv") # temp container for Gaus lines
  #   if (i == 1) {gausbars.df <- gausbarstmp} else {
  #     gausbars.df <- rbind(gausbars.df, gausbarstmp)}
  #   gaus = TRUE} else gaus = FALSE

  ####loops done create dfs####
  # create bin & Gaus barplot stats data frames
  if (bin) {binbars <- data.frame(Min.Inf = with(binbars.df, tapply(rel.inf, var, min)),
                                  Av.Inf = with(binbars.df, tapply(rel.inf, var, mean)),
                                  Max.Inf = with(binbars.df, tapply(rel.inf, var, max)),
                                  Inf.variance = with(binbars.df, tapply(rel.inf, var, var)),
                                  row.names = levels.default(binbars.df$var))
  binbars <- binbars[order(-binbars[,"Av.Inf"]),]
  binbarsgood <- subset(binbars, Min.Inf > 0)} #also create barplots for nonzero vars

  if (gaus) {gausbars <- data.frame(Min.Inf = with(gausbars.df, tapply(rel.inf, var, min)),
                                    Av.Inf = with(gausbars.df, tapply(rel.inf, var, mean)),
                                    Max.Inf = with(gausbars.df, tapply(rel.inf, var, max)),
                                    Inf.variance = with(gausbars.df, tapply(rel.inf, var, var)),
                                    row.names = levels.default(gausbars.df$var))
  gausbars <- gausbars[order(-gausbars[,"Av.Inf"]),]
  gausbarsgood <- subset(gausbars, Min.Inf > 0)}

  # create linesfiles end-column stats for each variable
  if (bin) for (l in colnames(samples[expvar])) {
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "MinLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, min, na.rm = TRUE)))
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "AvLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, mean, na.rm = TRUE)))
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "MaxLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, max, na.rm = TRUE)))
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "VarLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, var, na.rm = TRUE)))}

  if (gaus) for (m in colnames(samples[expvar])) {
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "MinLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, min, na.rm = TRUE)))
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "AvLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, mean, na.rm = TRUE)))
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "MaxLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, max, na.rm = TRUE)))
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "VarLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, var, na.rm = TRUE)))}

  # apply variances to a new column at the end of var.df
  if (calcpreds) var.df[,"C of V"] <- apply(var.df[,(3:(2 + loops))], MARGIN = 1, var, na.rm = TRUE)

  # Build CV & AUC stats report by scraping individual loops' reports
  Minima <- c(min(report.df[,1], na.rm = TRUE), min(report.df[,2], na.rm = TRUE), min(report.df[,3], na.rm = TRUE))
  Averages <- c(mean(report.df[,1], na.rm = TRUE), mean(report.df[,2], na.rm = TRUE), mean(report.df[,3], na.rm = TRUE))
  Maxima <- c(max(report.df[,1], na.rm = TRUE), max(report.df[,2], na.rm = TRUE), max(report.df[,3], na.rm = TRUE))
  Variances <- c(var(report.df[,1], na.rm = TRUE), var(report.df[,2], na.rm = TRUE), var(report.df[,3], na.rm = TRUE))
  # if all loops' values are NA (e.g. BinCV & AUC when only gaus), min will be Inf & max -Inf. A bit messy but Unimportant
  report.df <- rbind(report.df, Minima, Averages, Maxima, Variances)
  rep.len <- dim(report.df)[1]
  rownames(report.df) <- c((1:(rep.len - 4)), "Minima", "Averages", "Maxima", "Variances")

  ####save csvs####
  # create resvar named subfolder & go to it
  dir.create(names(samples[resvar]))
  setwd(names(samples[resvar]))

  if (savecsv) {
    if (bin) {write.csv(binbars, file = "BinBarsLoop.csv", row.names = T)
      write.csv(binbarsgood, file = "BinBarsGoodLoop.csv", row.names = T)}
    if (gaus) {write.csv(gausbars, file = "GausBarsLoop.csv", row.names = T)
      write.csv(gausbarsgood, file = "GausBarsGoodLoop.csv", row.names = T)}

    if (bin) for (n in colnames(samples[expvar])) {
      write.csv(get(paste0("binline_", n)), file = paste0("BinLineLoop_", n, ".csv"), row.names = F)}
    if (gaus) for (o in colnames(samples[expvar])) {
      write.csv(get(paste0("gausline_", o)), file = paste0("GausLineLoop_", o, ".csv"), row.names = F)}

    if (calcpreds) {write.csv(var.df, file = "VarAll.csv", row.names = F)
      write.csv(var.df[,c(1,2,(3 + loops))], file = "VarOnly.csv", row.names = F)}
    write.csv(report.df, file = "Report.csv", row.names = TRUE)
    print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      csv files created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))}

  ####plot linesfiles####
  if (bin) for (p in colnames(samples[expvar])) {
    yrange <- c(min(get(paste0("binline_", p))[,"MinLine"]), max(get(paste0("binline_", p))[,"MaxLine"]))
    png(filename = paste0("Bin_Loop_lines_", p, ".png"),
        width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.3,5,0.3,0.4), fig = c(0,1,0,1), las = 1, lwd = 8, bty = "n", mgp = c(1.25,0.5,0), xpd = NA)
    plot(get(paste0("binline_", p))[,1],
         get(paste0("binline_", p))[,"AvLine"],
         type = "l",
         #xlab = colnames(get(paste0("binline_", p)))[1],
         xlab = paste0(p, " (", round(binbars[p, "Av.Inf"],1), "%)"),
         ylab = "",
         main = "",
         yasx = "r",
         ylim = yrange)
    mtext("Marginal Effect", side = 2, line = 4.05, las = 0)
    lines(get(paste0("binline_", p))[,1], get(paste0("binline_", p))[,"MinLine"], col = "grey66") #[,1] is 1st column, X values, always the same
    lines(get(paste0("binline_", p))[,1], get(paste0("binline_", p))[,"MaxLine"], col = "grey33")
    legend("topleft", legend = c("Max","Av.","Min"), col = c("grey33","black","grey66"),
           lty = 1, pch = "-")
    dev.off()}

  if (gaus) for (q in colnames(samples[expvar])) {
    yrange <- c(min(get(paste0("gausline_", q))[,"MinLine"]), max(get(paste0("gausline_", q))[,"MaxLine"]))
    png(filename = paste0("Gaus_Loop_lines_", q, ".png"),
        width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.3,5,0.3,0.4), fig = c(0,1,0,1), las = 1, lwd = 8, bty = "n", mgp = c(1.25,0.5,0), xpd = NA)
    plot(get(paste0("gausline_", q))[,1],
         get(paste0("gausline_", q))[,"AvLine"],
         type = "l",
         #xlab = colnames(get(paste0("gausline_", q)))[1],
         xlab = paste0(q, " (", round(gausbars[q, "Av.Inf"],1), "%)"),
         ylab = "",
         main = "",
         yasx = "r",
         ylim = yrange)
    mtext("Marginal Effect", side = 2, line = 4.05, las = 0)
    lines(get(paste0("gausline_", q))[,1], get(paste0("gausline_", q))[,"MinLine"], col = "grey66")
    lines(get(paste0("gausline_", q))[,1], get(paste0("gausline_", q))[,"MaxLine"], col = "grey33")
    legend("topleft", legend = c("Max","Av.","Min"), col = c("grey33","black","grey66"),
           lty = 1, pch = "-")
    dev.off()}

  print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Line plots created      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

  # bin & gaus barplots####
  if (bin) {
    pointlineseqbin <- seq(0, length(binbars[,2]) - 1, 1)
    revseq <- rev(pointlineseqbin)
    png(filename = "BinBarsLoop.png", width = 4*480, height = 4*480, units = "px",
        pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.5,0.3,0,0.5), fig = c(0,1,0,1), cex.lab = 0.5, mgp = c(1.5,0.5,0), cex = 1.3, lwd = 6)
    midpoints <- barplot(rev(binbars[,2]), cex.lab = 1.2, las = 1,
                         horiz = TRUE, cex.names = 0.8, xlab = "Av. Influence %",
                         col = NA, border = NA,
                         xlim = c(0,2.5 + ceiling(max(binbars[,2]))),
                         ylim = c(0, length(binbars[,2])))
    for (r in 1:length(binbars[,2])) {
      lines(c(0, binbars[r,2]), c(revseq[r], revseq[r]), col = "black", lwd = 8)} #draw lines
    text(0.1, pointlineseqbin + (length(binbars[,2])/55), labels = rev(rownames(binbars)), adj = 0, cex = 0.8)
    axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
    dev.off()

    pointlineseqbin <- seq(0, length(binbarsgood[,2]) - 1, 1)
    revseq <- rev(pointlineseqbin)
    png(filename = "BinBarsGoodLoop.png", width = 4*480, height = 4*480, units = "px",
        pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.5,0.3,0,0.5), fig = c(0,1,0,1), cex.lab = 0.5, mgp = c(1.5,0.5,0), cex = 1.3, lwd = 6)
    midpoints <- barplot(rev(binbarsgood[,2]), cex.lab = 1.2, las = 1,
                         horiz = TRUE, cex.names = 0.8, xlab = "Av. Influence %",
                         col = NA, border = NA,
                         xlim = c(0,2.5 + ceiling(max(binbarsgood[,2]))),
                         ylim = c(0, length(binbarsgood[,2])))
    for (r in 1:length(binbarsgood[,2])) {
      lines(c(0, binbarsgood[r,2]), c(revseq[r], revseq[r]), col = "black", lwd = 8)} #draw lines
    text(0.1, pointlineseqbin + (length(binbarsgood[,2])/55), labels = rev(rownames(binbarsgood)), adj = 0, cex = 0.8)
    axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
    dev.off()}

  if (gaus) {
    pointlineseqgaus <- seq(0, length(gausbars[,2]) - 1, 1)
    revseq <- rev(pointlineseqgaus)
    png(filename = "GausBarsLoop.png", width = 4*480, height = 4*480, units = "px",
        pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.5,0.3,0,0.5), fig = c(0,1,0,1), cex.lab = 0.5, mgp = c(1.5,0.5,0), cex = 1.3, lwd = 6)
    midpoints <- barplot(rev(gausbars[,2]), cex.lab = 1.2, las = 1,
                         horiz = TRUE, cex.names = 0.8, xlab = "Av. Influence %",
                         col = NA, border = NA,
                         xlim = c(0,2.5 + ceiling(max(gausbars[,2]))),
                         ylim = c(0, length(gausbars[,2])))
    for (s in 1:length(gausbars[,2])) {
      lines(c(0, gausbars[s,2]), c(revseq[s], revseq[s]), col = "black", lwd = 8)}
    text(0.1, pointlineseqgaus + (length(gausbars[,2])/55), labels = rev(rownames(gausbars)), adj = 0, cex = 0.8)
    axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
    dev.off()

    pointlineseqgaus <- seq(0, length(gausbarsgood[,2]) - 1, 1)
    revseq <- rev(pointlineseqgaus)
    png(filename = "GausBarsGoodLoop.png", width = 4*480, height = 4*480, units = "px",
        pointsize = 4*12, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.5,0.3,0,0.5), fig = c(0,1,0,1), cex.lab = 0.5, mgp = c(1.5,0.5,0), cex = 1.3, lwd = 6)
    midpoints <- barplot(rev(gausbarsgood[,2]), cex.lab = 1.2, las = 1,
                         horiz = TRUE, cex.names = 0.8, xlab = "Av. Influence %",
                         col = NA, border = NA,
                         xlim = c(0,2.5 + ceiling(max(gausbarsgood[,2]))),
                         ylim = c(0, length(gausbarsgood[,2])))
    for (s in 1:length(gausbarsgood[,2])) {
      lines(c(0, gausbarsgood[s,2]), c(revseq[s], revseq[s]), col = "black", lwd = 8)}
    text(0.1, pointlineseqgaus + (length(gausbarsgood[,2])/55), labels = rev(rownames(gausbarsgood)), adj = 0, cex = 0.8)
    axis(side = 1, lwd = 6, outer = TRUE, xpd = NA)
    dev.off()}
  print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      Bar plots plotted      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))

  ####map predabund CofVs####
  if (calcpreds) if (varmap) { # if mapping requested,
    if (is.null(shape)) { # and shape not set, check presence of basemap
      if (!exists("gbm.basemap")) {stop("you need to install gbm.basemap to run this function")}
      bounds = c(range(grids[,gridslon]),range(grids[,gridslat])) #then create bounds
      #create standard bounds from data, and extra bounds for map aesthetic
      xmid <- mean(bounds[1:2])
      ymid <- mean(bounds[3:4])
      xextramax <- ((bounds[2] - xmid) * 1.6) + xmid
      xextramin <- xmid - ((xmid - bounds[1]) * 1.6)
      yextramax <- ((bounds[4] - ymid) * 1.6) + ymid
      yextramin <- ymid - ((ymid - bounds[3]) * 1.6)
      extrabounds <- c(xextramin, xextramax, yextramin, yextramax)
      shape <- gbm.basemap(bounds = extrabounds)
    } else {shape <- shape} # if shape not null then use it.

    png(filename = "CofVMap.png", width = 4*1920, height = 4*1920, units = "px",
        pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
    gbm.map(x = var.df[,gridslon], # add Unrepresentativeness alpha surface
            y = var.df[,gridslat],
            z = var.df[,"C of V"],
            mapmain = "Coefficient of Variation: ",
            species = names(samples[resvar]),
            legendtitle = paste0("C. of V.: ", measure),
            ####change this####
            shape = shape) #either autogenerated or set by user so never blank
    #breaks = expm1(breaks.grid(log(2000), ncol = 8, zero = FALSE))/2000)
    dev.off() #high value log breaks mean first ~5 values cluster near 0 for high
    # res there, but high values captures in the last few bins.
  } # close map optional
  setwd("../") # go back up from samples-named wd to original parent
  if (alerts) beep(3)
  if (calcpreds) return(var.df) #return output
} # close function
