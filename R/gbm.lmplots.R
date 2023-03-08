#' Plot linear models for all expvars against the resvar
#'
#' Loops the lmplot function, shows linear model plots for all expvars against the resvar. Good
#' practice to do this before running gbm.auto so you have a sense of the basic relationship of the
#' variables.
#'
#' @param samples Explanatory and response variables to predict from. Keep col
#' names short (~17 characters max), no odd characters, spaces, starting
#' numerals or terminal periods. Spaces may be converted to periods in directory
#' names, underscores won't. Can be a subset of a large dataset.
#' @param expvar Vector of names or column numbers of explanatory variables in
#' 'samples': c(1,3,6) or c("Temp","Sal"). No default.
#' @param resvar Name or column number(s) of response variable in samples: 12,
#' c(1,4), "Rockfish". No default. Column name is ideally species name.
#' @param expvarnames Vector of names same length as expvars, if you want nicer names.
#' @param resvarname Single character object, if you want a nicer resvar name.
#' @param savedir Save location, end with "/".
#' @param plotname Character vector of plot names else expvarnames else expvar will be used.
#' @param pngtype Filetype for png files, alternatively try "quartz" on Mac.
#' @param r2line Plot rsquared trendline, default TRUE.
#' @param pointtext Label each point? Default FALSE.
#' @param pointlabs Point labels, defaults to resvar value.
#' @param pointcol Points colour, default "black".
#' @param ... Allows controlling of text label params e.g. adj cex &.
#'
#' @return Invisibly saves png plots into savedir.
#'
#' @details Errors and their origins:
#'
#' @export
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
## Simon Dedman 2018.08.30 & 2023-03-07
## To Run:
# lmplot(x = samples[,expvar],
#        y = samples[,resvar],
#        xname = expvar,
#        yname = resvar)


gbm.lmplots <- function(samples = NULL, # dataframe
                        expvar = NULL, # character vector of colnames in samples
                        resvar = NULL, # single character colname in samples
                        expvarnames = NULL, # character vector same length as expvars
                        resvarname = NULL, # single character to overwrite resvar name if desired
                        savedir = NULL,
                        plotname = NULL,
                        r2line = TRUE, # plot rsquared trendline, default TRUE
                        pointtext = FALSE, # label each point? Default false
                        pointlabs = resvar, # point labels, defaults to resvar value
                        pointcol = "black" # points colour
){

  # removed so plot not imported which overwrites the use of dismo::plot
  #  #' @importFrom grDevices dev.off png
  #  #' @importFrom graphics legend par mtext abline text


  for (i in expvar) {
    # overwrite xname if expvarnames present
    xname <- ifelse(test = is.null(expvarnames),
                    yes = i,
                    no = expvarnames[which(expvar %in% i)])
    # overwrite plotname if plotname present
    plotname <- ifelse(test = is.null(plotname),
                       yes = xname,
                       no = plotname[which(expvar %in% i)])
    # overwrite yname if resvarname present
    yname <- ifelse(test = is.null(resvarname),
                    yes = resvar,
                    no = resvarname)
    # run plot
    gbm.auto::lmplot(
      x = samples[, i],
      y = samples[, resvar],
      xname = xname,
      yname = yname,
      savedir = savedir,
      pngtype = c("cairo-png", "quartz", "Xlib"),
      # xlab = xname, # x axis label, parsed from xname unless specified
      # ylab = yname, # y axis label, parsed from yname unless specified
      plotname = plotname, # filename for png, parsed from xname unless specified
      r2line = r2line, # plot rsquared trendline, default TRUE
      pointtext = pointtext, # label each point? Default false
      pointlabs = pointlabs, # point labels, defaults to resvar value
      pointcol = pointcol, # points colour
    )
  }
}
