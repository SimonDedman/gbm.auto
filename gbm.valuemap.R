# Load gbm.map
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")
# load data
matf<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F/Abundance_Preds_plusF.csv", header = TRUE)
juve<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/Abundance_Preds_plusF.csv", header = TRUE)

#automate this as a function notes#
# x/y/z data source as object to set at start e.g. matf

#automate this as a function notes#


# set wd
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Analysis")
# load data
data<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F/Abundance_Preds_plusF.csv", header = TRUE)

#Data input format: LATITUDE/LONGITUDE/{columns to be considered e.g. F, PredAbunds}
#Need a system for including all predabunds and F.
#Start with adding all predabunds
data$PredAbundAll <- data$PredAbund_C + data$PredAbund_T + data$PredAbund_B + data$PredAbund_S # add column in grids for predicted abundance

####23. Map maker####
if (!require(mapplots)) {stop("you need to install the mapplots package to run this function")}
require(mapplots)
  # generate output image & set parameters
  png(filename = "Conservation Value Mature Females.png",
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
  # run gbm.map function with generated parameters
  gbm.map(x = data[,"LONGITUDE"],
          y = data[,"LATITUDE"],
          z = data[,"PredAbundAll"],
          mapmain = "Conservation Value: ",
          species = "Mature Females",
          shape = coast,
          landcol = "darkgreen",
          legendloc = "bottomright",
          legendtitle = "Conservation Value",
          grids=data,
          gridslon=data[,"LONGITUDE"],
          gridslat=data[,"LATITUDE"],
          predabund=predabund)  # hopefully parses grids dataset to gbm.map to use
  dev.off()

mean(data[,"PredAbund_C"])
mean(data[,"PredAbund_T"])
mean(data[,"PredAbund_B"])
mean(data[,"PredAbund_S"])

# Scale Predicted abundances to 1 and add as columns to data
data$PredAbund_Cs <- data$PredAbund_C / max(data$PredAbund_C,na.rm=TRUE)
data$PredAbund_Ts <- data$PredAbund_T / max(data$PredAbund_T,na.rm=TRUE)
data$PredAbund_Bs <- data$PredAbund_B / max(data$PredAbund_B,na.rm=TRUE)
data$PredAbund_Ss <- data$PredAbund_S / max(data$PredAbund_S,na.rm=TRUE)
# then add them together as column
data$PredAbundAlls <- data$PredAbund_Cs + data$PredAbund_Ts + data$PredAbund_Bs + data$PredAbund_Ss
# then map that

png(filename = "Conservation Value Mature Females.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA
# run gbm.map function with generated parameters
gbm.map(x = data[,"LONGITUDE"],
        y = data[,"LATITUDE"],
        z = data[,"PredAbundAlls"],
        mapmain = "Conservation Value: ",
        species = "Mature Females",
        shape = coast,
        landcol = "darkgreen",
        legendloc = "bottomright",
        legendtitle = "Conservation Value",
        grids=data,
        gridslon=data[,"LONGITUDE"],
        gridslat=data[,"LATITUDE"],
        predabund=predabund)  # hopefully parses grids dataset to gbm.map to use
dev.off()

# manually run core of gbm.map in case
byx=NULL
byy=NULL
  # work out cell size for uniform square gridded data: Create blank vector for grid length calcs
  bydist<-rep(NA,length(data[,"LONGITUDE"]))
  # and attach it to data
  data<-cbind(data,bydist)
  # fill it: if [next longitude minus current longitude] equals [current longitude minus previous longitude], that's a uniform cell.
  # data rounded to prevent tiny fluctuations counting as differences. Need to set that tolerance.
  # Could do 10% of average distance between uniform points, but you don't know that til the end!
  data[2:(length(data[,"LONGITUDE"])-1),"bydist"] <-
    ifelse(round(data[2:(length(data[,"LONGITUDE"])-1),"LONGITUDE"]-data[1:(length(data[,"LONGITUDE"])-2),"LONGITUDE"],digits=5)
           ==
             round(data[3:length(data[,"LONGITUDE"]),"LONGITUDE"]-data[2:(length(data[,"LONGITUDE"])-1),"LONGITUDE"],digits=5),
           round(data[2:(length(data[,"LONGITUDE"])-1),"LONGITUDE"]-data[1:(length(data[,"LONGITUDE"])-2),"LONGITUDE"],digits=5),NA)
  # Take an average of those distances, they should all be identical anyway. Apply it to byx & byy.
  byx<-mean(data$bydist,na.rm=TRUE)
  byy<-byx

png(filename = "Conservation Value Mature Females Scaled.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA

grd <- make.grid(data[,"LONGITUDE"], data[,"LATITUDE"], data[,"PredAbundAlls"], byx, byy, xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]),fun=mean)
breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Conservation Value: ","Mature Females",sep=""))
draw.grid(grd,breaks) # plot grd data w/ breaks for colour breakpoints
draw.shape(coast, col="darkgreen") # add coastline
legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value")
dev.off()

# Scale Fishery catch to 1 and add as columns to data
data$Fs <- data$HansF_LPUE / max(data$HansF_LPUE,na.rm=TRUE)
# Subtract it from Predictions & add as new column
data$PredsFs <- data$PredAbundAlls - data$Fs


png(filename = "Cons.Val. Mat.F Scaled w F.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA

grd <- make.grid(data[,"LONGITUDE"], data[,"LATITUDE"], data[,"PredsFs"], byx, byy, xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]),fun=mean)
breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Conservation Value: ","Mature Females",sep=""))
draw.grid(grd,breaks) # plot grd data w/ breaks for colour breakpoints
draw.shape(coast, col="darkgreen") # add coastline
legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value")
dev.off()

# Makes NO difference!
summary(data$HansF_LPUE)
# ah. Crazy distribution:
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#    0.00     0.00     0.00    14.07     0.00 16980.00        1
# insane peak with low average means everything scales to 0 except a few samples

# plot fishing CPUE to test

png(filename = "Fishery LPUE.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)  #mgp:c:2,0.5,0, xpd=NA

grd <- make.grid(data[,"LONGITUDE"], data[,"LATITUDE"], data[,"HansF_LPUE"], byx, byy, xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]),fun=mean)
breaks <- breaks.grid(grd,zero=TRUE,quantile=1,ncol=6) #ncol to try to make 6 colours (inc 0)
basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Conservation Value: ","Mature Females",sep=""))
draw.grid(grd,breaks)
draw.shape(coast, col="darkgreen")
legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Fishery LPUE")
dev.off()