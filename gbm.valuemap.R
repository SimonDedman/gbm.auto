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
                       plotthis = c("good","bad","both","line","fish","fishermen"), # what to plot, defaults to everything, can delete any, to delete all set to NULL
                       savethis = c("data","report"), #which csvs to export, defaults to everything, can delete any, to delete all set to NULL
                       limitline = NULL, # minimum threshold for gooddata, to be added as flat line to plot
                       scalerange = NULL,  # which column numbers to be scaled to 1 (the data columns). Colname suffix "s" added. Optional: uses raw values if blank
                       goodweight = NULL,  # single or vector of weighting multiple(s) for goodcols array, no default
                       badweight = NULL,  # ditto for badcols array, no default
                       m = 1, # multiplication factor for Bpa units, default 1. 1000 to convert tonnes to kilos, 0.001 kilos to tonnes
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
# Automatically get intercept
# change my-specific terms to universal e.g. Ewgt, Bpa, E, etc.

# CPUEmap(gooddata) * ?         = Bmap (biodata) #? is the CPUE to Biomass conversion
 # CPUEmap(gooddata) * Emap(baddata)  = CatchMaps(bothdata, kinda. catchdata) #do this in (e.g.) 100 percentage steps of E
 # In which case why wouldn't we use Fishery F which is already a CatchMap but direct from source?
 # assumedly we have the CPUE and the E for the fishery? I already do right?
 # But the CPUE would have been derived by Hans from C & E, so we already have C! 
 # Just use C directly
# for 100 steps of E, CPUE*E[1:100] = C[1:100]
# biodata - C[1:100] = bioleftdata_C[1:100]
 
# So I still need that CPUE to B conversion from Sam.
# Am I expecting THAT, or just a final Bpa tonnage?
# BOTH. Chased 8/9/15
  
  # simple case: take highest ray abundance & close it. then next, then next. Done that, bothdatadf
  # or close areas with highest E: minEmap is wrong. Need to take highest E, close that and sum RAY CPUE at each point
  # keep eliminating E from max to 0. Find intercept where E=Bpa
  # wanna develop an elimination procedure which combines both the ray abundance and E.
  # eg. scale both Ray CPUE & E to 1. FOr each cell, add cpue & e. Make new map. High cpue & high E = 2.
  # already done, subtract.
  # looking to close area which would then encompass Bpa value
  # this is currently an effort management regime: it lowers the effort surface but doesn't actually close areas
  # it just says how much effort should be for Bpa
  
  # 0-1 rays and 1-0 effort (1= no effort), add them to get max of 2. "2" cells will have 0 E & max CPUE.
# close those. Close your way down this mountain adding B as you go til Bpa. That's your closed area.
   ##done but doesn't go on the same map, asked Hans
  
  # remaining thing: do you give them equal weight? Only weight if you scale. Only scale if you don't use raw values
  # weighting equally minimises displacement while maximising results for fish.
  # scale to 1 equal weighting also models the situation 'as is'
  # whatever the weighting is, you'll still keep eliminating rectangles until you reach Bpa.
  # what will weighting change. Dave thinks weighting E more will close low E zones preferentially.
  # As you weiht E down you close more good E and cause more displacement
  # line plot curve thingy for this as well: print the loop stages:
  # Bpa sum @ stage (already does) & percent open areas
  #
  # Closure map for most abundant species this. Blonde exclusion map/mask from above method.
  # Start with this as zero effort input integrated into the actual effort map i.e. delete this out of the effor tlayer
  # then re-run the closed effort loop for species 2 e.g. cuckoo.
  # at each loop run it should recalculate B for EVERY SPECIES.
  
  # Then expore what happens with any 1 species as you change effort.#
  
  # THEN: use this as input into P5 map 
  # Since loop1 blonde B will be Bpa
  # loop 2 cuckoo B will be Bpa BUT blonde B will be > than Bpa
  # by loop4 blonde B will be >>> than Bpa.
  # So when you put this as an input into P5, the closed area will be pre-populated
  # & the bars will already be filled
  
# Check & load gbm.map
if(!exists("gbm.map")) {stop("you need to install gbm.map to run this function")}
if(!require(beepr)) {stop("you need to install the beepr package to run this function")}
require(beepr)
  
# Check goodweight & badweight are the same length as goodcols & badcols
if(!is.null(goodweight)) if(!length(goodweight)==length(goodcols)) stop("number of goodweights doesn't match number of goodcols")
if(!is.null(badweight)) if(!length(badweight)==length(badcols)) stop("number of badweights doesn't match number of badcols")

# note original number of columns for use later
datacoln <- ncol(data)

# If scalerange has been entered, scale those columns, then multiply by weighting factors and sum as objects
ifelse(!is.null(scalerange),{ #test, open yes  
  # Scale values to 1 and overwrite in place, renaming
  for(i in scalerange){
    data[,i]<-data[i]/max(data[i],na.rm=TRUE)  # for each column to scale (i), bind a column to data with i/max(i)
    colnames(data)[i]<- paste(names(data)[i],"s",sep="")}

  # If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
  ifelse(!is.null(goodweight),
         gooddata <- as.matrix(data[,goodcols]) %*% goodweight,
         gooddata <- rowSums(data[,goodcols,drop=FALSE])) #Drop stops length=1 cols/scales results dropping to a list & breaking rowSums
  ifelse(!is.null(badweight),
         baddata <- as.matrix(data[,badcols]) %*% badweight, #is matrix
         baddata <- rowSums(data[,badcols,drop=FALSE]))  #is numeric
  
} , { # scalerange is NA: omitted by user & created as a NA object by function
  
  # If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
  ifelse(!is.null(goodweight),
    gooddata <- as.matrix(data[,goodcols]) %*% goodweight,
    gooddata <- rowSums(data[,goodcols,drop=FALSE]))
  ifelse(!is.null(badweight),
    baddata <- as.matrix(data[,badcols]) %*% badweight, #matrix
    baddata <- rowSums(data[,badcols,drop=FALSE])) #numeric

}) # close scale optional & scalerange ifselse

####change gooddata to biodata!####
####map gooddata####
dir.create(goodname) # create directory for map plots
breaks <- seq(0, max(baddata,na.rm=TRUE), length = 12)  # set breaks for gbm.map instead of breaks.grid.
####?Use ncols somehow?####
if("good" %in% plotthis){
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
dev.off()}

####map baddata####
if("bad" %in% plotthis){
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
dev.off()}

# ####loop & map good+bad####
# # generate n steps between the highest effort value, and 0.
 Esteps <- rev(seq(from=0, to=max(baddata,na.rm=TRUE), length.out=steps+1))
# # what would the effort be? from no effort to full effort.
# 
# for(q in Esteps){
#   # add a new column with baddata-step (minimums to 0), named based on step as a %
#   # this is n steps of e.g. effort, from 0 to 100%
#   data[, paste(round(Mod((q/max(Esteps))-1)*100,digits=1),"%E",sep="")] <- pmax(0, baddata - q)
# 
#   if("both" %in% plotthis){
#   png(filename = paste("./",goodname,"/",bothname,"_Map",round(Mod((q/max(Esteps))-1)*100,digits=2),".png",sep=""),
#       width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
#   par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
#   # run gbm.map function with generated parameters
#   gbm.map(x = data[,loncolno],
#           y = data[,latcolno],
#           #z = (gooddata - data[,ncol(data)])*m, #probably fine if they're both on the same scale
#           #z = (gooddata - baddata)*m, #no need to use ncol: that's the final column, which is baddata-q when q=0
#           z = (gooddata + baddata)*m, #updated post dave reid train chat
#           mapmain = paste(round(Mod((q/max(Esteps))-1)*100,digits=2),"%E, Total Bpa: ", round(sum(gooddata)-sum(data[,ncol(data)],na.rm=TRUE),digits=2),"t", sep=""), # Print bothdata sum on map
#           species = "",
#           legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
#           ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
#   dev.off()
#  }} #close q loop & plot option

# invert baddata
badmax<-max(baddata)
ifelse((baddata-badmax<=0),{baddata= -(baddata-badmax)},{baddata= (baddata-badmax)})
# create bothdata: good + (inverted) bad
bothdata <- gooddata + baddata


if("both" %in% plotthis){
  png(filename = paste("./",goodname,"/",goodname," plus ",badname,"_Map.png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
  # run gbm.map function with generated parameters
  gbm.map(x = data[,loncolno],
          y = data[,latcolno],
          z = bothdata*m,
          mapmain = "Fishing Effort: ",
          species = "",
          legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
          ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
  dev.off()}


# ####make & save good vs bad line plot####
# if("line" %in% plotthis){
# png(filename = paste("./",goodname,"/",goodname," vs ",badname," plot.png",sep=""),
#     width = 1920, height = 1920, units = "px", pointsize = 48, bg = "white", res = NA, family = "", type = "cairo-png")
# par(mar=c(3.2,4.8,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE, lwd=6)
# 
# plot(x=round(Mod((Esteps/max(Esteps))-1)*100,digits=2),
#      y=(sum(gooddata)-colSums(data[,(datacoln+1):ncol(data)],na.rm=TRUE))*m,
#      type="o",
#      main= "Fishing E weighting Vs Abundance",
#      xlab= "Fishing E weighting, %",
#      ylab= "")
# if(!is.null(limitline)) abline(h=limitline)
# axis(2)
# mtext("Ray Biomass, CPUE, tonnes",side=2, line=3.7, las=0)
# 
# dev.off()} #close line optional
####error####
#Error in `[.data.frame`(data, , (datacoln + 1):ncol(data)) :  undefined columns selected 

####map min Bpa & min E####
if("fish" %in% plotthis){  # min area for Bpa, starting with best fish cells
 
    png(filename = paste("./",goodname,"/","fish conservation combo map.png",sep=""),   # when it gets above bpa plot it
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
    gbm.map(x = data[,loncolno],
            y = data[,latcolno],
            z = bothdata*m,
            mapmain = "Bpa threshold starting with best CPUE",
            species = "",
            legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
            ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
    # dev.off()
    
        #minBmap <- data.frame(data[,loncolno], data[,latcolno], pmax(0,gooddata - baddata)*m)
        bothdatadf <- data.frame(data[,loncolno], data[,latcolno], data[,goodcols], bothdata*m) #goodcols is only single, this is a hack. *m? check
        # ####hack:goodcols. fix####
        bothdatadf <- bothdatadf[ order(-bothdatadf[,4]),] #sort largest to smallest bothdata value, i.e. highest cpue lowest e combo
        n <- 1 # set current row at start to 1
        for(i in bothdatadf[,3]){ #run down gooddata from best BOTHDATA to worst
          if(sum(bothdatadf[1:n,3]) < limitline) {n = n+1} else { #if sum of rows to this point < Bpa add another row and sum again
            bothdatadf <- bothdatadf[1:n,] # subset df for only the pertinent rows once reached Bpa

#             png(filename = paste("./",goodname,"/","closed areas.png",sep=""),   # when it gets above bpa plot it
#                 width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
#             par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
            gbm.map(x = c(bothdatadf[,1],min(data[,loncolno]),max(data[,loncolno])), #map those sites, i.e. closed area. Overlap map onto same plot. Works?
                    y = c(bothdatadf[,2],min(data[,latcolno]),max(data[,latcolno])), #use c() to add min & max from main dataset to expand range
                    z = c(bothdatadf[,4],0,0), #black closed areas overlaid. Two blanks added for the corner positions
                mapmain = "",
                species = "",
                heatcol=rep("black",12),
                legendtitle="")
                dev.off()
       break}}
    }

if("fishermen" %in% plotthis){  # min area for Bpa, starting with best fisherman effort cells
  minEmap <- data.frame(data[,loncolno], data[,latcolno], pmax(0,gooddata - data[,ncol(data)])*m)
  minEmap <- minEmap[ order(minEmap[,3]),] # smallest to largest B
  n <- 1 # reset current row at start to 1
 for(i in minEmap[,3]){
  if(sum(minEmap[1:n,3]) < limitline) {n = n+1} else { #add another row and sum again
    minEmap <- minEmap[1:n,]# subset df for only the pertinent rows
    png(filename = paste("./",goodname,"/","minEmap.png",sep=""),   # when it gets above bpa plot it
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
    par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
    gbm.map(x = minEmap[,1],
            y = minEmap[,2],
            z = minEmap[,3],
            mapmain = "Bpa threshold starting with best E",
            species = "",
            legendtitle=legendtitleV, # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
            ...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
    dev.off()
    break}}}

  ####write csvs####
if("data" %in% savethis) write.csv(data,row.names=FALSE, file = paste("./",goodname,"/","Scaled Data.csv",sep=""))
  # make & write map report: E weighting as row names, totals as row values.
if("report" %in% savethis) {
  bpa_map_report<-data.frame(Totals=(sum(gooddata)-colSums(data[,(datacoln+1):ncol(data)],na.rm=TRUE))*m)
  write.csv(bpa_map_report, file=paste("./",goodname,"/","Bpa Map Report.csv",sep=""),row.names=TRUE)}
  
  beep(8) # notify the user with a noise, since this process can take a long time.
} #close function