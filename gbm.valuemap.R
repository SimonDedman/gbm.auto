# what does this DO??
# It takes lat/long/abundance & other variables
# and performs scale and arithmetic functions on them
# producing an output vector to be fed into gbm.map's "z"

###automate this as a function notes#
# option to set alpha yes/no, zero yes/no, data, scalerange, goodcols, badcols, goodweight, badweight, coast shape, quantile,
###automate this as a function notes#

####running this (in gbm.auto to begin with)####
source("C:/Users/Simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")  #SD specific, remove @ end
####MatF####
# set wd
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/ConsValMaps")
# load data. Expecting Longitude, Latitude, data columns: predicted abundances, fishing effort etc.. E.g.: Abundance_Preds_All.csv from gbm.auto
data<-read.csv("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Both_Abundance_Preds_plusE.csv", header = TRUE)

gbm.valuemap(data,
             loncolno = 1,
             latcolno = 2,
             scalerange = c(3:11),
             goodcols = c(4:11),
             badcols = c(3),
             goodweight = c(1.5,1,3,1), #cuckoo 1.5, blonde
             badweight = 4,  #fishing 4 times as important, for example.
             species = names(samples[i]),
             ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile


####manual run & dave loop####


#Esteps <- seq(from=no effort, to=full effort, by=100)  # what would the effort be?
# trying to calculate z, which is what? bothdata. which is what?

# for(q in Esteps){
#   png(filename = paste("./",species,"/",bothname,"_Map",".png",sep=""),
#       width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
#   par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
#   # run gbm.map function with generated parameters
#   gbm.map(x = data[,loncolno],
#           y = data[,latcolno],
#           z = bothdata,
#           mapmain = paste(bothname," Value: ",sep=""),
#           species = bothname,
#           legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
#           ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
#   dev.off()}

setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Juveniles/Individual Predators/Thornback")
fullpreddata <- read.csv(file="Abundance_Preds_All.csv")
preddata <- fullpreddata[,c(1,2,9,47)]
badcols <- 3 #set badcol for later but also for esteps

# generate 100 steps between the highest effort value, and 0.
Esteps <- rev(seq(from=0, to=max(preddata[,badcols],na.rm=TRUE), length.out=10))  # what would the effort be? from no effort to full effort. Change to length.out=100

for(q in Esteps){
  preddata[,5] <- preddata[,badcols] - Esteps[q]  # due badcols to final column & reduce by increment q
      gbm.valuemap(data=preddata,
               loncolno = 1,
               latcolno = 2,
               scalerange = c(3,5),  # need to scale for the moment since I don't have bpa and therefore the values are all over the place
               goodcols = c(5), # only one: resvar
               badcols = c(3), # only one: E
               #goodweight = 1, #not weighting resvar in proper run?
               #badweight = q,  #.
               species = names(samples[i]), #or set manually
               loop = "both", # only map bothdata,
               mapmain = paste("FishingE weighting: ",badweight,", Total Bpa: ",bpasum=sum(bothdata),sep=""), # Print bothdata sum on map
               ...)
}




####function####
gbm.valuemap<-function(data,  # data.frame to load. Expects Lon, Lat & data columns: predicted abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
                       loncolno = 1, # column number in data which has longitudes
                       latcolno = 2, # column number in data which has latitudes
                       scalerange,  # which column numbers to be scaled to 1 (the data columns). Colname suffix "s" added. Optional: uses raw values if blank
                       goodcols,  # which column numbers are abundances (where higher = better)?
                       badcols,  # which column numbers are 'negative' elements e.g. fishing (where higher = worse)?
                       goodweight,  # single or vector of weighting multiple(s) for goodcols array, no default
                       badweight,  # ditto for badcols array, no default
                       species, # the response varible(s) name(s), e.g. species name(s), or collective term if agglomerating >1 response variable. Single character string, not a vector. No spaces or terminal periods.
                       goodname = "Species",
                       badname = "Impacts",
                       bothname = "Conservation",
                       legendtitleV="Conservation Value", # legend tile default, sends to legendtitle in gbm.map, differentially named to allow different
                       loop = c("good","bad","both"), # vector of good/bad/both to loop through for maps for loop+"data" & loop+"name"
                       ...)  # optional terms for gbm.map: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile

####WRITE PROPER DESCRIPTION####
#(this is gbm.map's). add note in paper that "species" in gbm.map defaults to "Response Variable" so only gbm.map necessaries are x,y,x (currently assuming BrIreland; ideally will be fixed)

# Produces bothdata, gooddata and baddata, and maps each of them.

# Generalised Boosting Models, automated map generator. Simon Dedman, 2014, simondedman@gmail.com, https://github.com/SimonDedman/gbm.auto
  
# Generates maps from the outputs of gbm.step then gbm.predict.grids, handled automatically within gbm.auto but can be run alone, and
# generates representativeness surfaces from the output of gbm.rsb (suggest: z = rsbdf[,"Unrepresentativeness"],
# mapmain = "Unrepresentativeness: ",legendtitle = "UnRep 0-1"). Suggested code for outputting to e.g. png:
# png(...); par(...); gbm.map(...); dev.off()

####todo####
# option to run only one?
# better to have it run only one anyway?
# balance between "precise & modular control for R people" and "simple and intuitive for non R people"
# change it to make a vector of goodname,badname,bothname, and loop through it making maps rather than writing 3 maps

# when called in gbm.auto, will legendtitle (as set here, conservation value) be overruled by something else? gbm.map's default is "CPUE" but this is ABOVE gbm.map, CALLING gbm.map
# nonetheless, if the user sets legendtitle in gbm.auto, it'll propogate to gbm.map>legend.grid AND gbm.valuemap>gbm.map>legend.grid.
# call it something other than legendtitle in gbm.valuemap e.g. legendtitleV, so the gbm.map call below is gbm.map(..., legendtitle=legendtitleV)
# Do I need to do this for other things?

# Esteps: work out what this should be.
# Add beepr at the end if the looping takes a long time.
# organise this better: get the function working w/ w/o scales & weights, test changing names etc, then
# delete all unused crap so it's a tight function
# then make the dave loop another function?


# Check & load gbm.map
if(!exists("gbm.map")) {stop("you need to install gbm.map to run this function")}

# note original number of columns for use later
datacoln <- ncol(data)

#Scale columns to 1#
# Set which column numbers to be scaled to 1
#scalerange <- c(3:11)

# which column numbers are abundances (where higher = better)?
#goodcols <- c(4:11)

# which column numbers are 'negative' elements e.g. fishing (where higher = worse)?
#badcols <- c(3)

# weighting multiple for goodcols array, no default
  #goodweight <- c(1.5,1,3,1) #cuckoo 1.5, blonde 3

# weighting multiple for badcols array, no default
  #badweight <- 4  #fishing 4 times as important, for example.
  
# If scalerange has been entered, scale those columns, then multiply by weighting factors and sum as objects
if(!is.na(scalerange)){   #Error: object 'scalerange' not found

# Scale values to 1 and add as columns to data
for(i in scalerange){
  data<-cbind(data,data[i]/max(data[i],na.rm=TRUE))  # for each column to scale (i), bind a column to data with i/max(i)
  colnames(data)[ncol(data)]<- paste(names(data)[i],"s",sep="") # then add column name "(i)s"
  }

# If weighting factors given, multiply scaled values & overwrite
if(exists("goodweight")) data[,match(goodcols,scalerange)+datacoln] <- goodweight * data[,match(goodcols,scalerange)+datacoln]
if(exists("badweight")) data[,match(badcols,scalerange)+datacoln] <- badweight * data[,match(badcols,scalerange)+datacoln]

# # column numbers for scaled variables in data:
# data[,match(goodcols,scalerange)+datacoln]
# data[,match(badcols,scalerange)+datacoln]

# Sum the good & bad data as objects. Create object for good minus bad
ifelse(length(goodcols)>1,
       gooddata<-rowSums(data[,match(goodcols,scalerange)+datacoln]),
       gooddata<-data[,match(goodcols,scalerange)+datacoln])
ifelse(length(badcols)>1,
       baddata<-rowSums(data[,match(badcols,scalerange)+datacoln]),
       baddata<-data[,match(badcols,scalerange)+datacoln])
} else { # scale optional: if scale omitted & therefore is creared as a NA object once function is run
gooddata <- sumprod(goodcols,goodweight)  #need the right term to multiply a vector of >=1 columns by the same no of numerics, to return a vector of results
baddata <- sumprod(badcols,badweight)  #ditto
  } # close scale optional

bothdata<-gooddata-baddata  #Error: object 'gooddata' not found

#####Map maker####
#  if (!require(mapplots)) {stop("you need to install the mapplots package to run this function")}
#  require(mapplots)
#  library(mapplots)
#  data(coast)

####manually run core of gbm.map####
 # incorporate into gbm.auto, have it auto-produce standard map (weight=1)? 
 # Gbm.auto = 1 per species so standard map already exists; produce bad map (fishing) and both map (CPUE-fishing)

# byx=NULL
# byy=NULL
#   # work out cell size for uniform square gridded data: Create blank vector for grid length calcs
#   bydist<-rep(NA,length(data[,"LONGITUDE"]))
#   # and attach it to data
#   data<-cbind(data,bydist)
#   # fill it: if [next longitude minus current longitude] equals [current longitude minus previous longitude], that's a uniform cell.
#   # data rounded to prevent tiny fluctuations counting as differences. Need to set that tolerance.
#   # Could do 10% of average distance between uniform points, but you don't know that til the end!
#   data[2:(length(data[,"LONGITUDE"])-1),"bydist"] <-
#     ifelse(round(data[2:(length(data[,"LONGITUDE"])-1),"LONGITUDE"]-data[1:(length(data[,"LONGITUDE"])-2),"LONGITUDE"],digits=5)
#            ==
#              round(data[3:length(data[,"LONGITUDE"]),"LONGITUDE"]-data[2:(length(data[,"LONGITUDE"])-1),"LONGITUDE"],digits=5),
#            round(data[2:(length(data[,"LONGITUDE"])-1),"LONGITUDE"]-data[1:(length(data[,"LONGITUDE"])-2),"LONGITUDE"],digits=5),NA)
#   # Take an average of those distances, they should all be identical anyway. Apply it to byx & byy.
#   byx<-mean(data$bydist,na.rm=TRUE)
#   byy<-byx

####map1####
# png(filename = "Conservation Value All Mature Females Scaled.png",
#     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)

# grd <- make.grid(data[,"LONGITUDE"],
#                  data[,"LATITUDE"],
#                  gooddata,
#                  byx,
#                  byy,
#                  xlim=range(data[,"LONGITUDE"]),
#                  ylim=range(data[,"LATITUDE"]),
#                  fun=mean)
# breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
# basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Conservation Value All Mature Females Scaled",sep=""))
  
# if(breaks[2]==0) { # if zero=true in breaks
#   draw.grid(grd,breaks,col=c("#00000000",colorRampPalette(c("lightyellow", "yellow", "orange", "red", "brown4"))(length(breaks)-2))) #add alpha as 1st colour in palette
# } else {
#   draw.grid(grd,breaks) # else don't
# }
# 
# draw.shape(coast, col="darkgreen") # add coastline
# 
# if(breaks[2]==0) { # if zero=true in breaks
#   legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value", col=c("#00000000",colorRampPalette(c("lightyellow", "yellow", "orange", "red", "brown4"))(length(breaks)-2))) #add alpha as 1st colour in palette
# } else {
#   legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value") # else don't
# }
#dev.off()

####newmap1####
# png(filename = paste("./",species,"/",goodname,"_Map",".png",sep=""),
#     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# # run gbm.map function with generated parameters
# gbm.map(x = data[,loncolno],
#         y = data[,latcolno],
#         z = gooddata,
#         mapmain = paste(goodname," Value: ",sep=""), # manually set this for valuemap 1 2 & 3
#         species = goodname,
#         legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
#         ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
# dev.off()

####map2####
# png(filename = "Fishing Effort Scaled.png",
#     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# 
# grd <- make.grid(data[,"LONGITUDE"],
#                  data[,"LATITUDE"],
#                  baddata,
#                  byx,
#                  byy,
#                  xlim=range(data[,"LONGITUDE"]),
#                  ylim=range(data[,"LATITUDE"]),
#                  fun=mean)
# breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
# basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("Fishing Effort",sep=""))
# 
# if(breaks[2]==0) { # if zero=true in breaks
#   draw.grid(grd,breaks,col=c("#00000000",colorRampPalette(c("lightyellow", "yellow", "orange", "red", "brown4"))(length(breaks)-2))) #add alpha as 1st colour in palette
# } else {
#   draw.grid(grd,breaks) # else don't
# }
# 
# draw.shape(coast, col="darkgreen") # add coastline
# 
# if(breaks[2]==0) { # if zero=true in breaks
#   legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value", col=c("#00000000",colorRampPalette(c("lightyellow", "yellow", "orange", "red", "brown4"))(length(breaks)-2))) #add alpha as 1st colour in palette
# } else {
#   legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value") # else don't
# }
# 
# dev.off()


####newmap2####
# png(filename = paste("./",species,"/",badname,"_Map",".png",sep=""),
#     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# # run gbm.map function with generated parameters
# gbm.map(x = data[,loncolno],
#         y = data[,latcolno],
#         z = baddata,
#         mapmain = paste(badname," Value: ",sep=""),
#         species = badname,
#         legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
#         ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
# dev.off()



####map3####
# png(filename = "All Mature Females Minus Fishing Effort Scaled.png",
#     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# 
# grd <- make.grid(data[,"LONGITUDE"],
#                  data[,"LATITUDE"],
#                  bothdata,
#                  byx,
#                  byy,
#                  xlim=range(data[,"LONGITUDE"]),
#                  ylim=range(data[,"LATITUDE"]),
#                  fun=mean)
# breaks <- breaks.grid(grd,zero=TRUE,quantile=1) # define breakpoints from grd, allow 0 category, max=max Z from grd
# basemap(xlim=range(data[,"LONGITUDE"]), ylim=range(data[,"LATITUDE"]), main=paste("All Mature Females Minus Fishing Effort Scaled.png",sep=""))
# 
# if(breaks[2]==0) { # if zero=true in breaks
#   draw.grid(grd,breaks,col=c("#00000000",colorRampPalette(c("lightyellow", "yellow", "orange", "red", "brown4"))(length(breaks)-2))) #add alpha as 1st colour in palette
# } else {
#   draw.grid(grd,breaks) # else don't
# }
# 
# draw.shape(coast, col="darkgreen") # add coastline
# 
# if(breaks[2]==0) { # if zero=true in breaks
#   legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value", col=c("#00000000",colorRampPalette(c("lightyellow", "yellow", "orange", "red", "brown4"))(length(breaks)-2))) #add alpha as 1st colour in palette
# } else {
#   legend.grid("bottomright", breaks=breaks, type=2, inset=0, bg="white", title="Conservation Value") # else don't
# }
# 
# dev.off()

####newmap3####
# png(filename = paste("./",species,"/",bothname,"_Map",".png",sep=""),
#     width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# # run gbm.map function with generated parameters
# gbm.map(x = data[,loncolno],
#         y = data[,latcolno],
#         z = bothdata,
#         mapmain = paste(bothname," Value: ",sep=""),
#         species = bothname,
#         legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
#         ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
# dev.off()


####data loop####
for(p in loop){   # Error: object 'loop' not found
  # for bothdata runs only, append sum of bothdata (i.e. resulting CPUE/Bpa) to an object
  if(p=="both"){
  if(!exists("bpa_map_report")) {bpa_map_report <- data.frame(badweight=badweight,bpasum=sum(bothdata),stringsAsFactors=FALSE)
  } else {
    bpa_map_report<-rbind(bpa_map_report,c(badweight,sum(bothdata)))
  }}
  
png(filename = paste("./",species,"/",get(paste(p,"name",sep="")),"_Map",".png",sep=""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# run gbm.map function with generated parameters
gbm.map(x = data[,loncolno],
        y = data[,latcolno],
        z = get(paste(p,"data",sep="")),
        mapmain = paste(get(paste(p,"name",sep=""))," Value: ",sep=""),
        species = get(paste(p,"name",sep="")),
        legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
        ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
dev.off()}

####write csvs####
write.csv(data,row.names=FALSE, file = "Scaled Data.csv")   #Error in as.data.frame.default(x[[i]], optional = TRUE):cannot coerce class ""function"" to a data.frame
write.csv(bpa_map_report, file="Bpa Map Report.csv",row.names=FALSE)  # Error in is.data.frame(x) : object 'bpa_map_report' not found