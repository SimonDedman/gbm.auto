#' Plot linear model for two variables with R2 & P printed and saved
#'
#' Simple function to plot and name a linear model
#'
#' @param x Explanatory variable data.
#' @param y Response variable data.
#' @param xname Variable name for plot header.
#' @param yname Variable name for plot header.
#' @param pngtype Filetype for png files, alternatively try "quartz" on Mac.
#' @param xlab X axis label, parsed from xname unless specified.
#' @param ylab Y axis label, parsed from yname unless specified.
#' @param plotname Filename for png, parsed from xname unless specified.
#' @param r2line Plot rsquared trendline, default TRUE.
#' @param pointtext Label each point? Default FALSE.
#' @param pointlabs Point labels, defaults to resvar value.
#' @param pointcol Points colour, default "black".
#' @param savedir Save location, end with "/".
#' @param ... Allows controlling of text label params e.g. adj cex &.
#'
#' @return Invisibly saves png plot into savedir.
#'
#' @details Errors and their origins:
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend par mtext abline text
#' @export
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
## Simon Dedman 2018.08.30 & 2023-03-07
## To Run:
# lmplot(x = Data[,Expvar],
#        y = Data[,Resvar],
#        xname = Expvar,
#        yname = Resvar)

lmplot <- function(x, # explanatory variable data
                   y, # response variable data
                   xname = "X variable", # variable name for plot header
                   yname = "Y variable", # variable name for plot header
                   pngtype = c("cairo-png", "quartz", "Xlib"),
                   xlab = xname, # x axis label, parsed from xname unless specified
                   ylab = yname, # y axis label, parsed from yname unless specified
                   plotname = xname, # filename for png, parsed from xname unless specified
                   r2line = TRUE, # plot rsquared trendline, default TRUE
                   pointtext = FALSE, # label each point? Default false
                   pointlabs = x, # point labels, defaults to resvar value
                   pointcol = "black", # points colour
                   savedir = "", # save location, end with /
                   ...){ # allows controlling of text label params e.g. adj cex &

  options(scipen = 5) # prevents plot using 1e+05 godawful notation for 100000
  pngtype <- match.arg(pngtype)
  fit <- lm(y ~ x) # linear model, save as object for summary
  parmardefault <- par("mar") # save original par values to reinstate
  png(filename = paste0(savedir, plotname, ".png"),
      width = 1920,
      height = 1920,
      units = "px",
      pointsize = 48,
      bg = "white",
      res = NA,
      family = "",
      type = pngtype)
  par("mar" = c(4.6, 4.1, 1.8, 1)) # bltr
  plot(x, y, xlab = xlab, ylab = ylab, pch = 20, col = pointcol) # plot points
  if (r2line) abline(fit) # plot trendline if requested
  if (pointtext) text(x = x, y = y, labels = pointlabs, ...) # add text labels if requested, allow colouring etc
  mtext(paste0(xname,
               " (x) vs ",
               yname,
               " (y). ",
               "Rsquared: ",
               round(summary(fit)$r.squared, digits = 3),
               ", P: ",
               round(summary(fit)$coefficients[,4][[2]], digits = 3)),
        side = 3,
        line = 0.5) # add title
  dev.off() # close plotting device
  options(scipen = 0) # revert scipen option to default
  par("mar" = parmardefault) # revert par options to default: bltr: (5.1 4.1 4.1 2.1)
  return(c(r2 = round(summary(fit)$r.squared, digits = 3),
           p = round(summary(fit)$coefficients[,4][[2]], digits = 3))) # return r2 & p for outside use if desired
}
