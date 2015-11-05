gbm.valuemap<-function(dbase,  # data.frame to load. Expects Lon, Lat & data columns: predicted abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
                       loncolno = 1, # column number in x which has longitudes
                       latcolno = 2, # column number in x which has latitudes
                       goodcols,  # which column numbers are abundances (where higher = better)? List them in order of highest conservation importance first e.g. c(3,1,2,4)
                       badcols,  # which column numbers are 'negative' elements e.g. fishing (where higher = worse)?
                       #legendtitleV = "Conservation Value", # legend tile default, sends to legendtitle in gbm.map, differentially named to allow different
                       plotthis = c("good","bad","both","close"), # what to plot, defaults to everything, can delete any, to delete all set to NULL
                       savethis = c("data","close"), #which csvs to export, defaults to everything, can delete any, to delete all set to NULL
                       HRMSY = NULL, # maximum percentage of each goodcols stock which can be removed each year, as decimal e.g. 0.15 = 15%. Single number or vector. Same order as goodcols
                       goodweight = NULL,  # single or vector of weighting multiple(s) for goodcols array, no default
                       badweight = NULL,  # ditto for badcols array, no default
                       m = 1, # multiplication factor for Bpa units, default 1. 1000 to convert tonnes to kilos, 0.001 kilos to tonnes. Assumedly the same for all goodcols.
                       #mapmain = paste(get(paste(p,"name",sep=""))," Value: ",sep=""), #default uses badcols steps value in map title
                       ...){  # optional terms for goodweight & badweight. And for gbm.map: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile species
# Check & load gbm.map
if(!exists("gbm.map")) {stop("you need to install gbm.map to run this function")}
if(!require(beepr)) {stop("you need to install the beepr package to run this function")}
require(beepr)
require(mapplots)
library(mapplots)
data(coast,package="mapplots")  # get Britain & Ireland coast data. I'm looking to make this global but am having a problem w/ the maps packge

goodname = colnames(dbase)[goodcols] #the response varible(s) name(s), e.g. species name(s), or collective term if agglomerating >1 response variable. Single character string, not a vector. No spaces or terminal periods.
badname = colnames(dbase)[badcols] #ditto for badcols. Both of these moved out of function parameters list to foce user to specify in colnames

# Check goodweight & badweight are the same length as goodcols & badcols
if(!is.null(goodweight)) if(!length(goodweight)==length(goodcols)) stop("number of goodweights doesn't match number of goodcols")
if(!is.null(badweight)) if(!length(badweight)==length(badcols)) stop("number of badweights doesn't match number of badcols")

# note original number of columns for use later
dbasecoln <- ncol(dbase)

# Scale values to 1 and add as columns at end of 'data' then name them
  for(i in 1:length(c(goodcols,badcols))){
    dbase[,dbasecoln+i]<-dbase[,(c(goodcols,badcols))[i]]/max(dbase[,(c(goodcols,badcols))[i]],na.rm=TRUE)  # for each column to scale (i), bind a column to data with i/max(i)
    colnames(dbase)[dbasecoln+i]<- paste(names(dbase)[(c(goodcols,badcols))[i]],"_s",sep="")}

# If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
ifelse(!is.null(goodweight),
         gooddata <- as.matrix(dbase[,dbasecoln+seq(1,length(goodcols))]) %*% diag(goodweight), #diag converts vector to matrix
         gooddata <- dbase[,dbasecoln+seq(1,length(goodcols))])
ifelse(!is.null(badweight), #we're currently assuming there's only going to be one baddata column. But this still works, leave it.
         baddata <- as.matrix(dbase[,dbasecoln+length(goodcols)+seq(1,length(badcols))]) %*% diag(badweight), #is matrix
         baddata <- dbase[,dbasecoln+length(goodcols)+seq(1,length(badcols))])  #is numeric

####map baddata####
if("bad" %in% plotthis){
  png(filename = paste("./",badname,"_Map.png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
  gbm.map(x = dbase[,loncolno],  # run gbm.map function with generated parameters
          y = dbase[,latcolno],
          z = baddata,
          mapmain = "Fishing Effort",
          species = "",
          legendtitle="Effort Level") # passes the lgenedtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
          #...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
  dev.off()}

# invert baddata
badmax<-max(baddata)
ifelse((baddata-badmax<=0),{baddata= -(baddata-badmax)},{baddata= (baddata-badmax)})
bothdata <- gooddata + baddata # create bothdata: good + (inverted) bad. ncols=gooddata+baddata
dbase <- cbind(dbase,bothdata) # add bothdata to data
for(i in 1:length(goodcols)){ #name bothdata cols: append w (Weighted) c (Combined)
  colnames(dbase)[ncol(dbase)-length(goodcols)+i]<- paste(names(dbase)[goodcols[i]],"_swc",sep="")}

####create DF for HRMSY masking maps in loop####
bothdatarange <- (ncol(dbase)-length(goodcols)+1):ncol(dbase)
dbase[,bothdatarange] <- bothdata*m

####Start Species Loop 1####
for(j in 1:length(goodcols)){  #loop through gooddata columns

####map gooddata####
if("good" %in% plotthis){
png(filename = paste("./",goodname[j],"_Map.png",sep=""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
# run gbm.map function with generated parameters
gbm.map(x = dbase[,loncolno],
        y = dbase[,latcolno],
        z = dbase[,goodcols[j]],
        species = goodname[j],
        legendtitle="CPUE") # passes the legendtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
        #...)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
dev.off()}

# map bothdata i.e. predabund scaled to 1 + effort scaled to 1 inverted.
if("both" %in% plotthis){
  png(filename = paste("./",goodname[j]," plus ",badname,"_Map.png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)
  # run gbm.map function with generated parameters
  gbm.map(x = dbase[,loncolno],
          y = dbase[,latcolno],
          z = dbase[,bothdatarange[j]], #bothdata[,j]*m,
          mapmain = "Predicted Abundance + Fishing Effort: ",
          species = goodname[j],
          #heatcolours = c("white", "yellow", "orange","red", "brown4"),
          heatcolours = c("red", "lightyellow","green"),
          legendtitle="0 - 2", # passes the legendtitleV set by user in gbm.valuemap call to legendtitle in gbm.map
          byxout = TRUE)  # optional terms: byx byy mapmain heatcol shape mapback landcol legendtitle lejback legendloc grdfun zero quantile
  dev.off()
  byx<-byxport
  byy<-byx}
  } # close species loop 1

maploopnames <- c("Combo","Biomass","Effort")
maploopcodes <- c("dbase[order(-dbase[,bothdatarange[1]-1+j]),]","dbase[order(-dbase[,2+j]),]","dbase[order(dbase[,badcols]),]")
for(o in 1:3){ # start o loop through combination, goodcols & badcols
j <- 1 # set / reset J so it restarts the loops at 1 rather than max
####Start Species Loop 2####
for(j in 1:length(goodcols)){  # j loop through gooddata columns
  n <- 1 # set current row at start to 1, for HRMSY subloop: preserves previous value for next species.
  ####HRMSY limit map####
if("close" %in% plotthis){  # min area for MRMSY, starting with best fish cells
  # Sort by bothdata descending then map that then overlay 15% biomass by highest bothdata
  CPUEMSY <- (sum(dbase[,2+j]) * HRMSY[j])
  dbase <- eval(parse(text = maploopcodes[o])) #sort largest to smallest bothdata value for that species, i.e. highest cpue lowest e combo.
   for(k in 1:nrow(dbase)){ # run the loop for every row, if required
      print(paste(n,", ",round((sum(dbase[1:n,2+j])/CPUEMSY)*100,2),"%",sep=""))
          if(sum(dbase[1:n,2+j]) < CPUEMSY) {n = n+1} else { #if sum of rows to this point < CPUEMSY add another row and sum again, ELSE:
            assign(paste("sort",maploopnames[o],j,sep=""),rep(0,nrow(dbase))) # create a vector of zeroes called sort[j number]
            dbase <- cbind(dbase,get(paste("sort",maploopnames[o],j,sep=""))) # bind it to x
            colnames(dbase)[ncol(dbase)] <- paste("sort",maploopnames[o],j,sep="") # reinstate its name (lost because bound with get())
            dbase[1:n-1,ncol(dbase)] <- rep(1,n-1) # populate first n-1 rows with 1 [k:n]
            badcut <- sum(dbase[1:n,badcols]) # sum badcols displaced
            badall <- sum(dbase[,badcols]) # total badcols
            badpct <- round((badcut/badall)*100,1) # badcols percent

            # INDIVIDUAL MAPS: Map bothdata. Then Overlay map onto bothdata map: ncol=1 zeroes=TRUE. Heatcol="black". Hopefully this makes zeroes invisible and closed area black.
            png(filename = paste("./ClosedValueMap",maploopnames[o],"_",goodname[j],".png",sep=""), # map species j's bothdata with black closed areas overlaid
                width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
            par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)

            # run gbm.map function's internal code. Set parameters
            x = dbase[,loncolno] #vector of longitudes, from make.grid in mapplots; x
            y = dbase[,latcolno] #vector of latitudes, from make.grid in mapplots; grids[,gridslat]
            z = dbase[,2*(length(goodcols)+length(badcols))+2+j] #vector of abundances generated by gbm.predict.grids, from make.grid in mapplots; grids[,predabund]
            heatcolours = c("red", "lightyellow","green")
            #heatcolours = c("white", "yellow", "orange","red", "brown4") #abundance colour scale, defaults to the heatcol from legend.grid & draw.grid in mapplots.
            colournumber = 8   #number of colours to spread heatcol over, default:8
            shape = coast   #basemap shape to draw, from draw.shape in mapplots. Defaults: 'coast': UK & Ire
            landcol = "grey80" #colour for 'null' area of map, if appropriate, from draw.shape in mapplots. Was "darkgreen" changed to light grey
            mapback = "lightblue" #basemap background
            legendloc = "bottomright" #location on map of legend box, from legend.grid in mapplots
            legendtitle = "0 - 2" #the metric of abundance, e.g. CPUE for fisheries, from legend.grid in mapplots
            lejback = "white"  #backgroud colour of legend, from legend.grid in mapplots
            zero = TRUE # allow 0 category in breaks.grid & thus legend? Default TRUE
            quantile = 1 # set max breakpoint; lower this to cutoff outliers
            #...) # breaks. vector of breakpoints for colour scales; default blank, generated automatically

              if(!exists("byx")) {    # if users hasn't entered byx or byy values, generate them from the data
                bydist<-rep(NA,length(x))   # work out cell size for uniform square gridded data: Create blank vector for grid length calcs
                cells<-data.frame(LONGITUDE=x,bydist=bydist,stringsAsFactors=FALSE)   # and attach it to grids
                cells[2:(length(x)-1),"bydist"] <- ifelse(round(cells[2:(length(x)-1),1]-cells[1:(length(x)-2),1],digits=5) ==
                  round(cells[3:length(x),1]-cells[2:(length(x)-1),1],digits=5),round(cells[2:(length(x)-1),1]-cells[1:(length(x)-2),1],digits=5),NA)
                byx<-mean(cells$bydist,na.rm=TRUE)  # Take an average of those distances, they should all be identical anyway. Apply it to byx & byy.
                byy<-byx}

              # Plot first map: same as bothdata
              grd <- make.grid(x, y, z, byx, byy, xlim=range(x), ylim=range(y),fun=mean) #create gridded data. fun defaults to sum which is bad
              heatcol = colorRampPalette(heatcolours)(colournumber) #create heatcol from component parts
              breaks <- breaks.grid(grd,zero=zero,quantile=quantile,ncol=length(heatcol))  #if breaks specified, do nothing (it'll be used later). Else generate it.
              if(zero){heatcol=c("#00000000",colorRampPalette(heatcol)(length(heatcol)-1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
              basemap(xlim=range(x), ylim=range(y), main=paste("Predicted Abundance + Fishing Effort: ",goodname[j],sep=""), bg=mapback, xlab = "Longitude (°W)", ylab = "Latitude (°N)")
              #remove xlab & ylab above for general code
              draw.grid(grd,breaks,col=heatcol) # plot grd data w/ breaks for colour breakpoints

              # Plot second map: closed area overlay only
              x = dbase[,loncolno] #set parameters for closed area map
              y = dbase[,latcolno]
              z = dbase[,ncol(dbase)]
              heatcolours2 = c("black","black")
              colournumber2 = 2
              grd2 <- make.grid(x, y, z, byx, byy, xlim=range(x), ylim=range(y),fun=mean) #create gridded data. fun defaults to sum which is bad
              heatcol2 = colorRampPalette(heatcolours2)(colournumber2) #create heatcol from component parts
              breaks2 <- breaks.grid(grd,zero=zero,quantile=quantile,ncol=length(heatcol2))  #if breaks specified, do nothing (it'll be used later). Else generate it.
              if(zero){heatcol2=c("#00000000",colorRampPalette(heatcol2)(length(heatcol2)-1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
              #basemap(xlim=range(x), ylim=range(y), main=paste("Predicted Abundance + Fishing Effort: ",goodname[j],sep=""), bg=mapback, xlab = "Longitude (°W)", ylab = "Latitude (°N)")
              #remove xlab & ylab above for general code
              draw.grid(grd2,breaks2,col=heatcol2) # plot grd data w/ breaks for colour breakpoints

              draw.shape(shape=shape, col=landcol) # add coastline
              legend.grid(legendloc, breaks=breaks, type=2, inset=0, bg=lejback, title=paste(badpct,"% E closed",sep=""), col=heatcol) #breaks=breaks/1000 was causing odd legend? From make.grid help, Hans using to convert kg to t?
              
            dev.off()

         break} # end of ELSE section
          } # end of FOR loop k (data rows)
      } # end of "close" optional calculation & mapping section
   } # end of 2nd FOR loop j (species)

if("close" %in% plotthis){ #re-open a "close" optional section now j loop finished, to process dbase
# sort largest to smallest species 1 value, then 2, 3 4 etc   [k:n]
  SortCol0 <- ncol(dbase)-length(goodcols)
  dbase <- dbase[order(-dbase[,SortCol0+1],-dbase[,SortCol0+2],-dbase[,SortCol0+3],-dbase[,SortCol0+4]),]
####how to do this algorithmically so I don't have to know how many species there are?####

# then add a column which is max(k:n). This will be 1s and 0s and will be the full extent (try append) [p]
assign(paste("AllClosed_",maploopnames[o],sep=""),pmax(dbase[,SortCol0+1],dbase[,SortCol0+2],dbase[,SortCol0+3],dbase[,SortCol0+4]))
dbase <- cbind(dbase,get(paste("AllClosed_",maploopnames[o],sep="")))
colnames(dbase)[ncol(dbase)] <- paste("AllClosed_",maploopnames[o],sep="") # reinstate its name (lost because bound with get())

# then add a column which is sum(k:n). This will be 0:length(goodcols) and will be the full extent. Both are similar [q]
assign(paste("SumClosed_",maploopnames[o],sep=""),rowSums(dbase[,(SortCol0+1):(SortCol0+4)]))
dbase <- cbind(dbase,get(paste("SumClosed_",maploopnames[o],sep="")))
colnames(dbase)[ncol(dbase)] <- paste("SumClosed_",maploopnames[o],sep="") # reinstate its name (lost because bound with get())

# then create a new df with a column of zeroes called Zeroes, length = nrow(x)
#Zeroes <- rep(0,nrow(dbase))
assign(paste("Zeroes_",maploopnames[o],sep=""),rep(0,nrow(dbase)))
MPAgrow <- as.data.frame(get(paste("Zeroes_",maploopnames[o],sep="")))
colnames(MPAgrow)[ncol(MPAgrow)] <- paste("Zeroes_",maploopnames[o],j,sep="") # reinstate its name (lost because bound with get())
assign(paste("Closure1_",maploopnames[o],sep=""),pmax(dbase[,SortCol0+1],MPAgrow[,1]))
MPAgrow <- cbind(MPAgrow,get(paste("Closure1_",maploopnames[o],sep="")))
colnames(MPAgrow)[ncol(MPAgrow)] <- paste("Closure1_",maploopnames[o],j,sep="") # reinstate its name (lost because bound with get())
assign(paste("Closure2_",maploopnames[o],sep=""),pmax(dbase[,SortCol0+2],MPAgrow[,2]))
MPAgrow <- cbind(MPAgrow,get(paste("Closure2_",maploopnames[o],sep="")))
colnames(MPAgrow)[ncol(MPAgrow)] <- paste("Closure2_",maploopnames[o],j,sep="") # reinstate its name (lost because bound with get())
assign(paste("Closure3_",maploopnames[o],sep=""),pmax(dbase[,SortCol0+3],MPAgrow[,3]))
MPAgrow <- cbind(MPAgrow,get(paste("Closure3_",maploopnames[o],sep="")))
colnames(MPAgrow)[ncol(MPAgrow)] <- paste("Closure3_",maploopnames[o],j,sep="") # reinstate its name (lost because bound with get())
assign(paste("Closure4_",maploopnames[o],sep=""),pmax(dbase[,SortCol0+4],MPAgrow[,4]))
MPAgrow <- cbind(MPAgrow,get(paste("Closure4_",maploopnames[o],sep="")))
colnames(MPAgrow)[ncol(MPAgrow)] <- paste("Closure4_",maploopnames[o],j,sep="") # reinstate its name (lost because bound with get())
# bind & save as global object (not global as now at top level, change to global if we have to move this inside the loop)
dbase <- cbind(dbase,MPAgrow)

# loop through final 4 cols of x i.e. [t:w] & map them: cumulative CPUEMSY limit maps
for (l in 1:4){

  # Generate badcol reduction percentage
  badcoldata <- dbase[,badcols] #vector of badcol data
  closecol <- dbase[,ncol(dbase)-4+l] #closures column
  badcut2 <- sum(badcoldata[closecol==1]) #sum of badcol data in closed cells
  badpct2 <- round((badcut2/badall)*100,1) # badcols percent
  
  png(filename = paste("./CumulativeClosedArea",maploopnames[o],"Map_",goodname[l],".png",sep=""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar=c(3.2,3,1.3,0), las=1, mgp=c(2.1,0.5,0),xpd=FALSE)

  gbm.map(x = dbase[,loncolno],
          y = dbase[,latcolno],
          z = dbase[,ncol(dbase)-4+l],
          mapmain = "Cumulative Closed Area: ",
          species = paste(maploopnames[o]," ",goodname[l],sep=""),
          heatcolours = c("black","black"),
          colournumber = 2,
          shape = coast,
          landcol = "grey80",
          mapback = "white",
          legendloc = "bottomright",
          legendtitle = paste(badpct2,"% E closed",sep=""),
          lejback = "white",
          zero = TRUE,
          quantile = 1,
          byx = byx,
          byy = byy)
    dev.off()
  } # end of l loop (column numbers of growing MPA in x)
 } # end of second "close" optional section
} # end o loop through combination, goodcols & badcols

####save csvs####
if("data" %in% savethis) write.csv(dbase,row.names=FALSE, file = paste("./ProcessedData.csv",sep=""))
  beep(8) # notify the user with a noise, since this process can take a long time.
} #close function