#' Defines breakpoints for draw.grid and legend.grid; mapplots fork
#'
#' Defines breakpoints from values in grd with options to exclude outliers,
#' set number of bins, and include a dedicated zero column.
#' Forked by SD 05/01/2019 to add 'lo', else bins always begin at 0, killing
#' plotting when all data are in a tight range at high values e.g. 600:610
#'
#' @param grd An array produced by make.grid or a list produced by
#'  make.multigrid or a vector of positive values.
#' @param quantile The maximum value of the breaks will be determined by the
#' quantile given here. This can be used to deal with outlying values in grd.
#' If quantile = 1 then the maximum value of the breaks will be the same as the
#' maximum value in grd.
#' @param ncol Number of colours to be used, always one more than the number of
#' breakpoints. Defaults to 12.
#' @param zero Logical, should zero be included as a separate category? Defaults
#'  to TRUE.
#'
#' @export
#' @return A vector of breakpoints for draw.grid in mapplots
#' @importFrom graphics legend
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Hans Gerritsen
#' @examples
#' breaks.grid(100,ncol=6)
#' breaks.grid(100,ncol=5,zero=FALSE)
#'
#' # create breaks on the log scale
#' exp(breaks.grid(log(10000),ncol=4,zero=FALSE))

breaks.grid <- function(grd, quantile = 0.975, ncol = 12, zero = TRUE)
{
  if (is.list(grd) == FALSE)
    grd <- list(grd)
  qua <- max(c(unlist(lapply(grd, quantile, probs = quantile,
                             na.rm = TRUE)), 0), na.rm = T)
  if (qua > 0)
    grd <- lapply(grd, function(x) ifelse(x > qua, qua, x))
  hi <- max(unlist(lapply(grd, max, na.rm = T)))
  lo <- min(unlist(lapply(grd, min, na.rm = T))) #added by SD
  if (zero) {
    len <- ncol
    breaks <- c(0, seq(lo, hi, length = len)) #lo replaces 0
  }
  else {
    len <- ncol + 1
    breaks <- seq(lo, hi, length = len) #lo replaces 0
  }
  return(breaks)
}
