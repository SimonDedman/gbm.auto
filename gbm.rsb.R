gbm.rsb <- function(samples, grids, expvarnames,gridslat,gridslon,rsbres){
# Generalised Boosting Models, Representativeness Surface Builder. Simon Dedman, 2014, simondedman@gmail.com

# Loops through explanatory variables comparing their histogram in samples to their histogram in grids to see how well the explanatory
# variable range in samples represents the range being predicted to in grids. Assigns a representativeness score per variable per site in
# grids, and takes the average score per site if there's more than 1 expvar. Saves this to a CSV; it's plotted by gbm.map if called in
# gbm.auto. This shows you which areas have the most & least representative coverage by samples, therefore where you can have the most /
# least confidence in the predictions from gbm.predict.grids. Can be called directly, and choosing a subset of expvars allows one to see
# their individual / collective representativeness.

# samples: data frame with response & explanatory variables
# grids: data frame of (more/different) explanatory variables & no response variable, to be predicted to by gbm.predict.grids
# expvarnames: vector of column names of explanatory variables being tested. Can be length 1. Names must match in samples & grids.
# gridslat: column number for latitude in 'grids'
# gridslon: column number for longitude in 'grids'
# rsbres: column number of response variable in samples [i]

  # loop through explanatory variables
  for (q in seq(from=1, to=length(expvarnames))){
    # range min = lowest value per variable
    nmin<-min(grids[,expvarnames[q]],samples[,expvarnames[q]],na.rm=TRUE)
    # ditto for max
    nmax<-max(grids[,expvarnames[q]],samples[,expvarnames[q]],na.rm=TRUE)
    # bin range is the length between the two
    binrange<-nmax-nmin
    # 10 bins. Length of one bin = binrange/10
    bin<-binrange/10
    # set breaks, min to max, 10 binrange increments. 0.01 added as findInterval (later) needs x to be < nmax, and some will == nmax, causing NAs.
    binbreaks<-c(nmin,nmin+bin,nmin+(bin*2),nmin+(bin*3),nmin+(bin*4),nmin+(bin*5),nmin+(bin*6),nmin+(bin*7),nmin+(bin*8),nmin+(bin*9),nmax+0.01)
    # make object from samples histogram
    assign(paste("hist_samples_",expvarnames[q],sep=""),hist(samples[,expvarnames[q]],breaks=binbreaks,plot=FALSE))
    # make object from grids histogram
    assign(paste("hist_grids_",expvarnames[q],sep=""),hist(grids[,expvarnames[q]],breaks=binbreaks,plot=FALSE))
    # calculate difference between frequencies, assign to object
    assign(paste("hist_diff_",expvarnames[q], sep=""), (get(paste("hist_samples_",expvarnames[q],sep=""))$density*bin - get(paste("hist_grids_",expvarnames[q],sep=""))$density*bin))
    # calculate modulus of that
    assign(paste("hist_diff_mod_",expvarnames[q], sep=""),sqrt(get(paste("hist_diff_",expvarnames[q], sep=""))^2))
    # create a vector for the diff lookup results: from that expvar's dataframe, get the diff value (col4) for the bin range number corresponding to the expvar value in grids
    assign(paste(expvarnames[q],"_hist_diff",sep=""),get(paste("hist_diff_",expvarnames[q], sep=""))[findInterval(as.numeric(unlist(grids[expvarnames[q]])),binbreaks)])
    # create a vector for the modulus lookup results: from that expvar's dataframe, get the modulus value (col5) for the bin range number corresponding to the expvar value in grids
    assign(paste(expvarnames[q],"_hist_diff_mod",sep=""),get(paste("hist_diff_mod_",expvarnames[q], sep=""))[findInterval(as.numeric(unlist(grids[expvarnames[q]])),binbreaks)])
    # put those 2 vectors in a dafa frame (first expvar) or add them to the existing one (latter expvars)
    ifelse(q==1, rsbdf<-data.frame(get(paste(expvarnames[q],"_hist_diff",sep="")),get(paste(expvarnames[q],"_hist_diff_mod",sep=""))),
                 rsbdf<-data.frame(rsbdf,get(paste(expvarnames[q],"_hist_diff",sep="")),get(paste(expvarnames[q],"_hist_diff_mod",sep=""))))
    # name those columns
    colnames(rsbdf)[(length(rsbdf)-1):length(rsbdf)] <- c(paste(expvarnames[q],"_hist_diff",sep=""),paste(expvarnames[q],"_hist_diff_mod",sep=""))
    # close expvar loop
  }
  # create vector of sum of mod diffs, scaled to score out of 1. Add to rsbdf. Globally assign so it's available to gbm.map as Z later. Will cause problems in lops?
  rsbdf<-data.frame(rsbdf,"Unrepresentativeness"=rowMeans(rsbdf[ls(pattern="_hist_diff_mod")]))
  # write out grids to a csv
  write.csv(c(grids[c(gridslat,gridslon)],rsbdf), row.names=FALSE, file=paste("./",names(samples[rsbres]),"/RSB.csv",sep=""))
  # return rsbdf as the object result of this function, for use elsewhere
  rsbdf
}