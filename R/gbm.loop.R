#' Calculate variance surfaces for gbm.auto predictions
#'
#' Processes a user-specified number of loops through the same gbm.auto
#' parameter combinations and calculates the variance in the predicted abundance
#' scores for each site aka cell. This can be mapped to spatially demonstrate
#' the output variance range.
#'
#' @param loops The number of loops required, integer
#' @param savecsv Save the variances in simple & extended format
#' @param varmap Create a map of the variance outputs
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
#' @param mapshape See gbm.auto help
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
#'
#' @examples gbmloopexample <- gbm.loop(loops = 2, samples = mysamples,
#' grids = mygrids, expvar = c(4:10), resvar = 11, simp = F, RSB = F)
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
gbm.loop <- function(loops = 10, # the number of loops required, integer
                     savecsv = T, # save the variances in simple & extended format
                     varmap = F, # create a map of the variance outputs?
                     grids,         # explantory data to predict to. Import with (e.g.)
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
                     mapshape = NULL,      # set coast shapefile, else downloaded and autogenerated
                     RSB = TRUE,           # run Unrepresentativeness surface builder?
                     BnW = F,           # repeat maps in black and white e.g. for print journals
                     alerts = T,        # play sounds to mark progress steps
                     pngtype = "cairo-png",# filetype for png files, alternatively try "quartz"
                     ...) {
  var.df <- grids[,c(gridslon, gridslat)] # create df with just lat & longs
  for (i in 1:loops) { # loop through all gbm.autos
    dir.create(paste("./", i, sep = ""))
    setwd(paste("./", i, sep = ""))
    gbm.auto(grids = grids,
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
             mapshape = mapshape,
             RSB = RSB,
             BnW = BnW,
             alerts = alerts,
             pngtype = pngtype,
             ...) # accept other gbm.auto values than these basics
    setwd(paste("./", colnames(samples[resvar]), sep = "")) # set wd to species named subfolder
    predtmp <- read.csv("Abundance_Preds_only.csv") # temp container for latest preds
    var.df <- cbind(var.df, predtmp[,3]) # cbind preds to existing lat/longs or other preds
    colnames(var.df)[2 + i] <- paste("Loop_", i, sep = "") # label newly added preds column
    setwd("../../") # move back up to root folder
  } # close i loop & go to the next one

  var.df[,"Variances"] <- apply(var.df[,(3:(2 + loops))], MARGIN = 1, var)
  # apply variances to a new column at the end of var.df

  if (savecsv) {write.csv(var.df, file = "VarAll.csv", row.names = F)
    write.csv(var.df[,c(1,2,(3 + loops))], file = "VarOnly.csv", row.names = F)}

  if (varmap) { # if mapping requested,
    if (is.null(mapshape)) { # and mapshape not set, check presence of basemap
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
    } else {shape <- mapshape} # if mapshape not null then use it.

    png(filename = "VarMap.png", width = 4*1920, height = 4*1920, units = "px",
        pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
    gbm.map(x = var.df[,gridslon], # add Unrepresentativeness alpha surface
            y = var.df[,gridslat],
            z = var.df[,"Variances"],
            mapmain = "Variance: ",
            species = names(samples[resvar]),
            legendtitle = "SCALE!!",
            ####change this####
            shape = shape) #either autogenerated or set by user so never blank
    #breaks = expm1(breaks.grid(log(2000), ncol = 8, zero = FALSE))/2000)
    dev.off() #high value log breaks mean first ~5 values cluster near 0 for high
    # res there, but high values captures in the last few bins.
  } # close map optional

  return(var.df) #return output

#### add a final alert####
} # close function
