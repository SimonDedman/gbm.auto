###automate this as a function notes#
# in maps, possible to have an option to set 0 as alpha?
###automate this as a function notes#


# Load gbm.map
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")

####MatF####
# set wd
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F/ConsValMaps")
# load data
data<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Mature Females plus Hans' F/MatF_Abundance_Preds_plusE.csv", header = TRUE)
# note original number of columns for use later
datacoln <- ncol(data)

####Scale columns to 1####
# Set column numbers to scale to 1
scalerange <- c(3:7)
# which column numbers are abundances (where higher = better)?
goodcols <- c(4:7)
# which column numbers are 'negative' elements e.g. fishing (where higher = worse)?
badcols <- c(3)
# weighting multiple for goodcols array, no default
goodweight <- c(1.5,1,3,1) #cuckoo 1.5, blonde 3
# weighting multiple for goodcols array, no default
#badweight <- rep(1,length(badcols))
badweight <- 4  #fishing 4 times as important, for example.
  
# Scale values to 1 and add as columns to data
for(i in scalerange){
  data<-cbind(data,data[i]/max(data[i],na.rm=TRUE))
  colnames(data)[ncol(data)]<- paste(names(data)[i],"s",sep="")
  }

# If weighting factors given, multiply scaled values & overwrite
if(exists("goodweight")) data[,match(goodcols,scalerange)+datacoln] <- goodweight * data[,match(goodcols,scalerange)+datacoln]
if(exists("badweight")) data[,match(badcols,scalerange)+datacoln] <- badweight * data[,match(badcols,scalerange)+datacoln]

# # column numbers for scaled variables in data:
# data[,match(goodcols,scalerange)+datacoln]
# data[,match(badcols,scalerange)+datacoln]

# Sum the good & bad data as objects. Create object for good-bad
ifelse(length(goodcols)>1,
       gooddata<-rowSums(data[,match(goodcols,scalerange)+datacoln]),
       gooddata<-data[,match(goodcols,scalerange)+datacoln])
ifelse(length(badcols)>1,
       baddata<-rowSums(data[,match(badcols,scalerange)+datacoln]),
       baddata<-data[,match(badcols,scalerange)+datacoln])
bothdata<-gooddata-baddata
# summary(bothdata)


#Map maker
 if (!require(mapplots)) {stop("you need to install the mapplots package to run this function")}
 require(mapplots)
 library(mapplots)
 data(coast)

# manually run core of gbm.map
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

####map1####
png(filename = "Conservation Value All Mature Females Scaled.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)

grd <- make.grid(data[,"LONGITUDE"],
                 data[,"LATITUDE"],
                 gooddata,
                 byx,
                 byy,
                 xlim=range(data[,"LONGITUDE"]),
                 ylim=range(data[,"LATITUDE"]),
                 fun=mean)
breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Conservation Value All Mature Females Scaled",sep=""))
draw.grid(grd,breaks) # plot grd data w/ breaks for colour breakpoints
draw.shape(coast, col="darkgreen") # add coastline
legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value")
dev.off()

#####map2####
png(filename = "Fishing Effort Scaled.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)

grd <- make.grid(data[,"LONGITUDE"],
                 data[,"LATITUDE"],
                 baddata,
                 byx,
                 byy,
                 xlim=range(data[,"LONGITUDE"]),
                 ylim=range(data[,"LATITUDE"]),
                 fun=mean)
breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Fishing Effort",sep=""))
draw.grid(grd,breaks) # plot grd data w/ breaks for colour breakpoints
draw.shape(coast, col="darkgreen") # add coastline
legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value")
dev.off()

#####map3####
png(filename = "All Mature Females Minus Fishing Effort Scaled.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)

grd <- make.grid(data[,"LONGITUDE"],
                 data[,"LATITUDE"],
                 bothdata,
                 byx,
                 byy,
                 xlim=range(data[,"LONGITUDE"]),
                 ylim=range(data[,"LATITUDE"]),
                 fun=mean)
breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("All Mature Females Minus Fishing Effort Scaled.png",sep=""))
draw.grid(grd,breaks) # plot grd data w/ breaks for colour breakpoints
draw.shape(coast, col="darkgreen") # add coastline
legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value")
dev.off()

####write csv####
write.csv(data,row.names=FALSE, file = "Scaled Data.csv")