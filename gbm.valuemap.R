gbm.valuemap<-function(data,  # data.frame to load. Expects Lon, Lat & data columns: predicted abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
                       loncolno = 1, # column number in data which has longitudes
                       latcolno = 2, # column number in data which has latitudes
                       goodcols,  # which column numbers are abundances (where higher = better)?
                       badcols,  # which column numbers are 'negative' elements e.g. fishing (where higher = worse)?
                       goodname = "Species", #the response varible(s) name(s), e.g. species name(s), or collective term if agglomerating >1 response variable. Single character string, not a vector. No spaces or terminal periods.
                       badname = "Impacts",
                       bothname = "Conservation",
                       legendtitleV = "Conservation Value", # legend tile default, sends to legendtitle in gbm.map, differentially named to allow different
                       steps = 100, # number of steps for good vs bad mapping curve / Bpa
                       #scalerange,  # which column numbers to be scaled to 1 (the data columns). Colname suffix "s" added. Optional: uses raw values if blank
                       #goodweight,  # single or vector of weighting multiple(s) for goodcols array, no default
                       #badweight,  # ditto for badcols array, no default
                       #mapmain = paste(get(paste(p,"name",sep=""))," Value: ",sep=""), #default uses badcols steps value in map title
                       ...){  # optional terms for goodweight & badweight. And for gbm.map: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile species

####WRITE PROPER DESCRIPTION####
# Generalised Boosting Models, automated conservation value map generator.
# Simon Dedman, 2015, simondedman@gmail.com, https://github.com/SimonDedman/gbm.auto

# Produces bothdata, gooddata and baddata, and maps each of them.
  # what does this DO??
  # It takes lat/long/abundance & other variables
  # and performs scale and arithmetic functions on them
  # producing an output vector to be fed into gbm.map's "z"


####todo####
# multi threading, still
# produce line graph of E vs Bpa. Thicken & clean lines. Add Bpa line. Automatically get intercept?
# change my-specific terms to universal e.g. Ewgt, Bpa, E, etc.
# good/bad/both maps optional? Maybe have all default to TRUE
# not using the right data: hansF should be E!
# bothdata maps title bpa numbers going in right direction but mostly massively negative. Wait til I'm using the right datasets then look at.
  
# Check & load gbm.map
if(!exists("gbm.map")) {stop("you need to install gbm.map to run this function")}
if(!require(beepr)) {stop("you need to install the beepr package to run this function")}
require(beepr)
  
# Check goodweight & badweight are the same length as goodcols & badcols
if(exists("goodweight")) if(!length(goodweight)==length(goodcols)) stop("number of goodweights doesn't match number of goodcols")
if(exists("badweight")) if(!length(badweight)==length(badcols)) stop("number of badweights doesn't match number of badcols")

# note original number of columns for use later
datacoln <- ncol(data)

# If scalerange has been entered, scale those columns, then multiply by weighting factors and sum as objects
ifelse(exists("scalerange"),{ #test, open yes  
  # Scale values to 1 and overwrite in place, renaming
  for(i in scalerange){
    data[,i]<-data[i]/max(data[i],na.rm=TRUE)  # for each column to scale (i), bind a column to data with i/max(i)
    colnames(data)[i]<- paste(names(data)[i],"s",sep="")}

  # If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
  ifelse(exists("goodweight"),
         gooddata <- as.matrix(data[,goodcols]) %*% goodweight,
         gooddata <- rowSums(data[,goodcols,drop=FALSE])) #Drop stops length=1 cols/scales results dropping to a list & breaking rowSums
  ifelse(exists("badweight"),
         baddata <- as.matrix(data[,badcols]) %*% badweight, #is matrix
         baddata <- rowSums(data[,badcols,drop=FALSE]))  #is numeric
  
} , { # scalerange is NA: omitted by user & created as a NA object by function
  
  # If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
  ifelse(exists("goodweight"),
    gooddata <- as.matrix(data[,goodcols]) %*% goodweight,
    gooddata <- rowSums(data[,goodcols,drop=FALSE]))
  ifelse(exists("badweight"),
    baddata <- as.matrix(data[,badcols]) %*% badweight, #matrix
    baddata <- rowSums(data[,badcols,drop=FALSE])) #numeric

}) # close scale optional & scalerange ifselse


####map gooddata####
dir.create(goodname) # create directory for map plots
breaks <- seq(0, max(baddata,na.rm=TRUE), length = 12)  # set breaks for gbm.map instead of breaks.grid.
####?Use ncols somehow?####

png(filename = paste("./",goodname,"/",goodname,"_Map.png",sep=""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# run gbm.map function with generated parameters
gbm.map(x = data[,loncolno],
        y = data[,latcolno],
        z = gooddata,
        species = goodname,
        legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
        ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
dev.off()

####map baddata####
png(filename = paste("./",goodname,"/",badname,"_Map.png",sep=""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# run gbm.map function with generated parameters
gbm.map(x = data[,loncolno],
        y = data[,latcolno],
        z = baddata,
        mapmain = "Fishing Effort: ",
        species = badname,
        legendtitle="Effort Level", # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
        ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
dev.off()

####loop & map good-bad####
# generate n steps between the highest effort value, and 0.
Esteps <- rev(seq(from=0, to=max(baddata,na.rm=TRUE), length.out=steps+1))
# what would the effort be? from no effort to full effort.

for(q in Esteps){
  # add a new column with baddata-step (minimums to 0), named based on step as a %
  data[, paste(Mod((q/max(Esteps))-1)*100,"%E",sep="")] <- pmax(0, baddata - q)

  png(filename = paste("./",goodname,"/",bothname,"_Map",round(Mod((q/max(Esteps))-1)*100,digits=2),".png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
  # run gbm.map function with generated parameters
  gbm.map(x = data[,loncolno],
          y = data[,latcolno],
          z = gooddata - data[,ncol(data)], #probably fine if they're both on the same scale
          mapmain = paste(round(Mod((q/max(Esteps))-1)*100,digits=2),"%E, Total Bpa: ", round(sum(gooddata)-sum(data[,ncol(data)],na.rm=TRUE),digits=2), sep=""), # Print bothdata sum on map
          species = "",
          legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
          ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
  dev.off()
 } #close q loop

  ####make & save good vs bad line plot####
png(filename = paste("./",goodname,"/",goodname," vs ",badname," plot.png",sep=""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE, lwd=6)

plot(x=round(Mod((Esteps/max(Esteps))-1)*100,digits=2),
     y=sum(gooddata)-colSums(data[,(datacoln+1):ncol(data)],na.rm=TRUE),
     type="o",
     main= "Fishing E weighting Vs Abundance",
     xlab= "Fishing E weighting",
     ylab= "Ray Abundance")

dev.off()
  ####write csvs####
  write.csv(data,row.names=FALSE, file = paste("./",goodname,"/","Scaled Data.csv",sep=""))
  
  # make & write map report: E weighting as row names, totals as row values.
  bpa_map_report<-data.frame(Totals=colSums(data[,(datacoln+1):ncol(data)],na.rm=TRUE))
  write.csv(bpa_map_report, file=paste("./",goodname,"/","Bpa Map Report.csv",sep=""),row.names=TRUE)
  beep(8) # notify the user with a noise, since this process can take a long time.
} #close function