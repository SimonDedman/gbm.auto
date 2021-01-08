#' Representativeness Surface Builder
#'
#' Loops through explanatory variables comparing their histogram in 'samples' to
#' their histogram in 'grids' to see how well the explanatory variable range in
#' samples represents the range being predicted to in grids. Assigns a
#' representativeness score per variable per site in grids, and takes the
#' average score per site if there's more than 1 expvar. Saves this to a CSV;
#' it's plotted by gbm.map if called in gbm.auto. This shows you which areas
#' have the most and least representative coverage by samples, therefore where you
#' can have the most/least confidence in the predictions from gbm.predict.grids.
#' Can be called directly, and choosing a subset of expvars allows one to see
#' their individual / collective representativeness.
#'
#' @param samples Data frame with response and explanatory variables.
#' @param grids Data frame of (more/different) explanatory variables and no
#' response variable, to be predicted to by gbm.predict.grids.
#' @param expvarnames Vector of column names of explanatory variables being
#' tested. Can be length 1. Names must match in samples and grids.
#' @param gridslat Column number for latitude in 'grids'.
#' @param gridslon Column number for longitude in 'grids'.
#'
#' @return Gridded data table of representativeness values which is then mapped
#' with gbm.map and also saved as a csv
#' @export
#' @importFrom graphics hist
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples None
#'
gbm.rsb <- function(samples, grids, expvarnames, gridslat, gridslon){
# Generalised Boosting Models, Representativeness Surface Builder. Simon Dedman, 2014, simondedman@gmail.com

# Loops through explanatory variables comparing their histogram in samples to their histogram in grids to see how well the explanatory
# variable range in samples represents the range being predicted to in grids. Assigns a representativeness score per variable per site in
# grids, and takes the average score per site if there's more than 1 expvar. Saves this to a CSV; it's plotted by gbm.map if called in
# gbm.auto. This shows you which areas have the most and least representative coverage by samples, therefore where you can have the most /
# least confidence in the predictions from gbm.predict.grids. Can be called directly, and choosing a subset of expvars allows one to see
# their individual / collective representativeness.

# samples: data frame with response and explanatory variables
# grids: data frame of (more/different) explanatory variables and no response variable, to be predicted to by gbm.predict.grids
# expvarnames: vector of column names of explanatory variables being tested. Can be length 1. Names must match in samples and grids.
# gridslat: column number for latitude in 'grids'
# gridslon: column number for longitude in 'grids'

  # loop through explanatory variables
  for (q in seq(from = 1, to = length(expvarnames))) {
    # range min = lowest value per variable
    nmin <- min(grids[,expvarnames[q]],samples[,expvarnames[q]],na.rm = TRUE)
    # ditto for max
    nmax <- max(grids[,expvarnames[q]],samples[,expvarnames[q]],na.rm = TRUE)
    # bin range is the length between the two
    binrange <- nmax - nmin
    # 10 bins. Length of one bin = binrange/10
    bin <- binrange/10
    # set breaks, min to max, 10 binrange increments. 0.01 added as findInterval (later) needs x to be < nmax, and some will == nmax, causing NAs.
    binbreaks <- c(nmin, nmin + bin, nmin + (bin * 2), nmin + (bin * 3), nmin + (bin * 4), nmin + (bin * 5), nmin + (bin * 6),nmin + (bin * 7), nmin + (bin * 8), nmin + (bin*9), nmax + 0.01)
    # make object from samples histogram
    assign(paste0("hist_samples_", expvarnames[q]), hist(samples[,expvarnames[q]], breaks = binbreaks, plot = FALSE))
    # make object from grids histogram
    assign(paste0("hist_grids_", expvarnames[q]), hist(grids[,expvarnames[q]], breaks = binbreaks, plot = FALSE))
    # calculate difference between frequencies, assign to object
    assign(paste0("hist_diff_", expvarnames[q]), (get(paste0("hist_samples_", expvarnames[q]))$density*bin - get(paste0("hist_grids_",expvarnames[q]))$density * bin))
    # calculate modulus of that #could use abs() to do this
    assign(paste0("hist_diff_mod_",expvarnames[q]),sqrt(get(paste0("hist_diff_", expvarnames[q])) ^ 2))
    # create a vector for the diff lookup results: from that expvar's dataframe, get the diff value (col4) for the bin range number corresponding to the expvar value in grids
    assign(paste0(expvarnames[q],"_hist_diff"),get(paste0("hist_diff_", expvarnames[q]))[findInterval(as.numeric(unlist(grids[expvarnames[q]])), binbreaks)])
    # create a vector for the modulus lookup results: from that expvar's dataframe, get the modulus value (col5) for the bin range number corresponding to the expvar value in grids
    assign(paste0(expvarnames[q],"_hist_diff_mod"),get(paste0("hist_diff_mod_", expvarnames[q]))[findInterval(as.numeric(unlist(grids[expvarnames[q]])), binbreaks)])
    # put those 2 vectors in a dafa frame (first expvar) or add them to the existing one (latter expvars)
    ifelse(q == 1,
           rsbdf <- data.frame(get(paste0(expvarnames[q],"_hist_diff")), get(paste0(expvarnames[q], "_hist_diff_mod"))),
           rsbdf <- data.frame(rsbdf, get(paste0(expvarnames[q], "_hist_diff")), get(paste0(expvarnames[q], "_hist_diff_mod"))))
    # name those columns
    colnames(rsbdf)[(length(rsbdf) - 1):length(rsbdf)] <- c(paste0(expvarnames[q],"_hist_diff"), paste0(expvarnames[q], "_hist_diff_mod"))
  }  # close expvar loop
  # create vector of sum of mod diffs, scaled to score out of 1. Add to rsbdf. Globally assign so it's available to gbm.map as Z later. Will cause problems in loops?
  rsbdf <- data.frame("Latitude" = grids[,gridslat],
                      "Longitude" = grids[,gridslon],
                      rsbdf,
                      "Unrepresentativeness" = rowMeans(rsbdf[ls(pattern = "_hist_diff_mod")]))
  rsbdf}   # return rsbdf as the object result of this function, for use elsewhere
