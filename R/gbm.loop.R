#' Calculate Coefficient Of Variation surfaces for gbm.auto predictions
#'
#' Processes a user-specified number of loops through the same gbm.auto
#' parameter combinations and calculates the Coefficient Of Variation in the
#' predicted abundance scores for each site aka cell. This can be mapped to
#' spatially demonstrate the output variance range.
#'
#' @param loops The number of loops required, integer
#' @param savecsv Save the coefficients of variation in simple & extended format
#' @param calcpreds Calculate coefficients of variation of predicted abundance?
#' @param varmap Create a map of the coefficients of variation outputs?
#' @param measure Map legend, coefficients of variation of what? Default CPUE
#' @param cleanup Remove gbm.auto-generated directory each loop? Default FALSE
#' @param grids See gbm.auto help for all subsequent params
#' @param samples See gbm.auto help
#' @param expvar See gbm.auto help
#' @param resvar See gbm.auto help
#' @param tc See gbm.auto help
#' @param lr See gbm.auto help
#' @param bf See gbm.auto help
#' @param ZI See gbm.auto help
#' @param simp See gbm.auto help
#' @param gridslat See gbm.auto help
#' @param gridslon See gbm.auto help
#' @param cols See gbm.auto help
#' @param linesfiles See gbm.auto help
#' @param savegbm See gbm.auto help
#' @param varint See gbm.auto help
#' @param map See gbm.auto help
#' @param shape See gbm.auto help
#' @param RSB See gbm.auto help
#' @param BnW See gbm.auto help
#' @param alerts See gbm.auto help
#' @param pngtype See gbm.auto help
#' @param ... Additional params for gbm.auto subfunctions inc gbm.step
#'
#' @return Returns a data frame of lat, long, 1 predicted abundance per loop,
#' and a final variance score per cell.
#'
#' @export
#' @importFrom beepr beep
#'
#' @examples
#' library("gbm.auto")
#' mygrids <- gbm.auto::grids # load grids
#' mysamples <- gbm.auto::samples # load samples
#' gbmloopexample <- gbm.loop(loops = 2, samples = mysamples,
#' grids = mygrids, expvar = c(4:10), resvar = 11, simp = F)
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
gbm.loop <- function(loops = 10, # the number of loops required, integer
                     savecsv = TRUE, # save the coefficients of variation in simple & extended format
                     calcpreds = TRUE, # calculate coefficients of variation of predicted abundance?
                     varmap = T, # create a map of the coefficients of variation outputs?
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
                     ZI = "CHECK",         # are data zero-inflated? TRUE/FALSE/"CHECK".
                     # TRUE: delta BRT, log-normalised Gaus, reverse log-norm and bias corrected.
                     # FALSE: do Gaussian only, no log-normalisation.
                     # CHECK: Tests data for you. Default is TRUE.
                     simp = TRUE,          # try simplfying best BRTs?
                     gridslat = 2,         # column number for latitude in 'grids'
                     gridslon = 1,         # column number for longitude in 'grids'
                     cols = grey.colors(1,1,1), # barplot colour vector. Assignment in order of
                     # explanatory variables. Default 1*white: white bars black borders. '1*' repeats
                     linesfiles = F,   # save individual line plots' data as csv's?
                     savegbm = F,       # save gbm objects and make available in environment after running? Open with load("Bin_Best_Model")
                     varint = F,        # calculate variable interactions? Default:TRUE, FALSE
                     # for error "contrasts can be applied only to factors with 2 or more levels"
                     map = F,           # save abundance map png files?
                     shape = NULL,      # set coast shapefile, else downloaded and autogenerated
                     RSB = F,           # run Unrepresentativeness surface builder?
                     BnW = F,           # repeat maps in black and white e.g. for print journals
                     alerts = T,        # play sounds to mark progress steps
                     pngtype = "cairo-png",# filetype for png files, alternatively try "quartz"
                     ...) {

  # Generalised Boosting Model / Boosted Regression Tree process chain automater
  # Simon Dedman, 2012-6 simondedman@gmail.com github.com/SimonDedman/gbm.auto

  if (alerts) if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
  if (alerts) require(beepr)

  binbars.df <- data.frame(var = rep(NA, length(expvar)),
                           rel.inf = rep(NA, length(expvar)))
  gausbars.df <- binbars.df # blnk dataframes for bin & gaus bars data
  if (calcpreds) var.df <- grids[,c(gridslon, gridslat)] # create df with just lat & longs

  for (i in 1:loops) { # loop through all gbm.autos
    dir.create(paste0("./", i)) # create i'th folder
    setwd(paste0("./", i)) # move to it
    gbm.auto(grids = grids, # run i'th gbm.auto
             samples = samples,
             expvar = expvar,
             resvar = resvar,
             tc = tc,
             lr = lr,
             bf = bf,
             ZI = ZI,
             simp = simp,
             gridslat = gridslat,
             gridslon = gridslon,
             cols = cols,
             linesfiles = linesfiles,
             savegbm = savegbm,
             varint = varint,
             map = map,
             shape = shape,
             RSB = RSB,
             BnW = BnW,
             alerts = alerts,
             pngtype = pngtype,
             ...) # accept other gbm.auto values than these basics
    setwd(paste0("./", colnames(samples[resvar]))) # set wd to species named subfolder

    if (file.exists("Binary BRT Variable contributions.csv")) {
    binbarstmp <- read.csv("Binary BRT Variable contributions.csv") # temp container for bin bars
    if (i == 1) {binbars.df <- binbarstmp} else {# csv file to df unless df exists
      binbars.df <- rbind(binbars.df, binbarstmp)} # if so add to bottom of existing
    bin = TRUE} else bin = FALSE

    # loop thru variables name linesfiles e.g. Bin_Best_line_[varname].csv
    # adding i'th loop's values as new column
    if (bin) for (j in colnames(samples)[expvar]) {
      assign(paste0("bintmp_", j), read.csv(paste0("Bin_Best_line_", j, ".csv"))) # temp container for bin bars
      if (i == 1) {assign(paste0("binline_", j), get(paste0("bintmp_", j)))
        #colnames(paste0("binline_", j))[1 + i] <- paste0("Loop_", i) # label newly added column
      } else {
        assign(paste0("binline_", j), cbind(get(paste0("binline_", j)), get(paste0("bintmp_", j))[,2]))
        #colnames(paste0("binline_", j))[1 + i] <- paste0("Loop_", i)
        }}

    if (file.exists("Gaussian BRT Variable contributions.csv")) {
    gausbarstmp <- read.csv("Gaussian BRT Variable contributions.csv") # temp container for Gaus bars
    if (i == 1) {gausbars.df <- gausbarstmp} else {
      gausbars.df <- rbind(gausbars.df, gausbarstmp)}
    gaus = TRUE} else gaus = FALSE

    if (gaus) for (k in colnames(samples)[expvar]) {
      # if (!file.exists(paste0("Gaus_Best_line_", k, ".csv"))) { #if variable line csv doesn't exist
      #   assign(paste0("gaustmp_", k),
      assign(paste0("gaustmp_", k), read.csv(paste0("Gaus_Best_line_", k, ".csv")))
      ##fails if variable influence is 0 and not plotted or has been removed by simp, goes to read csv but csv not present.
      if (i == 1) {assign(paste0("gausline_", k), get(paste0("gaustmp_", k)))
        #colnames(get(paste0("gausline_", k)))[1 + i] <- paste0("Loop_", i) # label newly added column
        # target of assignment expands to non-language object
        # see https://stackoverflow.com/questions/14464442/using-get-with-replacement-functions
        } else {
        assign(paste0("gausline_", k), cbind(get(paste0("gausline_", k)), get(paste0("gaustmp_", k))[,2]))
          #colnames(paste0("gausline_", k))[1 + i] <- paste0("Loop_", i)
          }}

    if (!file.exists("Abundance_Preds_only.csv")) calcpreds = FALSE
    if (calcpreds) {predtmp <- read.csv("Abundance_Preds_only.csv") # temp container for latest preds
    var.df <- cbind(var.df, predtmp[,3]) # cbind preds to existing lat/longs or other preds
    colnames(var.df)[2 + i] <- paste0("Loop_", i)} # label newly added preds column

    setwd("../../") # move back up to root folder
    if (cleanup) unlink(i, recursive = T)
  } # close i loop & go to the next i

####loops done create dfs####
  # create bin & Gaus barplot stats data frames
  if (bin) {binbars <- data.frame(Min.Inf = with(binbars.df, tapply(rel.inf, var, min)),
                         Av.Inf = with(binbars.df, tapply(rel.inf, var, mean)),
                         Max.Inf = with(binbars.df, tapply(rel.inf, var, max)),
                         Inf.variance = with(binbars.df, tapply(rel.inf, var, var)),
                         row.names = levels.default(binbars.df$var))}
  if (gaus) {gausbars <- data.frame(Min.Inf = with(gausbars.df, tapply(rel.inf, var, min)),
                         Av.Inf = with(gausbars.df, tapply(rel.inf, var, mean)),
                         Max.Inf = with(gausbars.df, tapply(rel.inf, var, max)),
                         Inf.variance = with(gausbars.df, tapply(rel.inf, var, var)),
                         row.names = levels.default(gausbars.df$var))}

  # create linesfiles end-column stats for each variable
    if (bin) for (l in colnames(samples)[expvar]) {
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "MinLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, min)))
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "AvLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, mean)))
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "MaxLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, max)))
    assign(paste0("binline_", l), cbind(get(paste0("binline_", l)),
                                        "VarLine" = apply(get(paste0("binline_", l))[, (2:(1 + loops))], MARGIN = 1, var)))}

  if (gaus) for (m in colnames(samples)[expvar]) {
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "MinLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, min)))
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "AvLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, mean)))
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "MaxLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, max)))
    assign(paste0("gausline_", m), cbind(get(paste0("gausline_", m)),
                                         "VarLine" = apply(get(paste0("gausline_", m))[, (2:(1 + loops))], MARGIN = 1, var)))}

  # apply variances to a new column at the end of var.df
  if (calcpreds) var.df[,"C of V"] <- apply(var.df[,(3:(2 + loops))], MARGIN = 1, var)

####save csvs####
  if (savecsv) {
    if (bin) write.csv(binbars, file = "BinBarsLoop.csv", row.names = T)
    if (gaus) write.csv(gausbars, file = "GausBarsLoop.csv", row.names = T)

    if (bin) for (n in colnames(samples)[expvar]) {
      write.csv(get(paste0("binline_", n)), file = paste0("BinLineLoop_", n, ".csv"), row.names = F)}
    if (gaus) for (o in colnames(samples)[expvar]) {
      write.csv(get(paste0("gausline_", o)), file = paste0("GausLineLoop_", o, ".csv"), row.names = F)}

    if (calcpreds) {write.csv(var.df, file = "VarAll.csv", row.names = F)
    write.csv(var.df[,c(1,2,(3 + loops))], file = "VarOnly.csv", row.names = F)}}

####plot linesfiles####
  if (bin) for (p in colnames(samples)[expvar]) {
    png(filename = paste0("Bin_Loop_lines_", p, ".png"),
        width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.3,5,0.3,0.4), fig = c(0,1,0,1), las = 1, lwd = 8, bty = "n", mgp = c(1.25,0.5,0), xpd = NA)
    plot(get(paste0("binline_", p))[,1],
       get(paste0("binline_", p))[,"MaxLine"],
       type = "l",
       xlab = colnames(get(paste0("binline_", p)))[1],
       ylab = "Marginal Effect",
       main = "",
       yasx = "r")
  lines(get(paste0("binline_", p))[,1], get(paste0("binline_", p))[,"MinLine"], col = "red")
  lines(get(paste0("binline_", p))[,1], get(paste0("binline_", p))[,"AvLine"], col = "blue")
  legend("topleft", legend = c("Max","Av.","Min"), col = c("black","blue","red"),
         lty = 1, pch = "-")
  dev.off()}

  if (gaus) for (q in colnames(samples)[expvar]) {
    png(filename = paste0("Gaus_Loop_lines_", q, ".png"),
        width = 4*480, height = 4*480, units = "px", pointsize = 80, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(2.3,5,0.3,0.4), fig = c(0,1,0,1), las = 1, lwd = 8, bty = "n", mgp = c(1.25,0.5,0), xpd = NA)
    plot(get(paste0("gausline_", q))[,1],
         get(paste0("gausline_", q))[,"MaxLine"],
         type = "l",
         xlab = colnames(get(paste0("gausline_", q)))[1],
         ylab = "Marginal Effect",
         main = "",
         yasx = "r")
    lines(get(paste0("gausline_", q))[,1], get(paste0("gausline_", q))[,"MinLine"], col = "red")
    lines(get(paste0("gausline_", q))[,1], get(paste0("gausline_", q))[,"AvLine"], col = "blue")
    legend("topleft", legend = c("Max","Av.","Min"), col = c("black","blue","red"),
           lty = 1, pch = "-")
    dev.off()}
  ## need to change from lines for factorial variables
  ## need to fix the y axis range lengths, between minmin & maxmax
  ## marginal effect y axis label values are raw values not the +/- from gbm.plot

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
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
    gbm.map(x = var.df[,gridslon], # add Unrepresentativeness alpha surface
            y = var.df[,gridslat],
            z = var.df[,"Variances"],
            mapmain = "Coefficient of Variation: ",
            species = names(samples[resvar]),
            legendtitle = paste0("C. of V.: ", measure),
            ####change this####
            shape = shape) #either autogenerated or set by user so never blank
    #breaks = expm1(breaks.grid(log(2000), ncol = 8, zero = FALSE))/2000)
    dev.off() #high value log breaks mean first ~5 values cluster near 0 for high
    # res there, but high values captures in the last few bins.
  } # close map optional
  if (alerts) beep(3)
  if (calcpreds) return(var.df) #return output
} # close function
