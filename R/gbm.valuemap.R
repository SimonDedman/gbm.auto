#' Decision Support Tool that generates (Marine) Protected Area options using species predicted abundance maps
#'
#' Scales response variable data, maps a user-defined explanatory variable to be
#' avoided, e.g. fishing effort, combines them into a map showing areas to
#' preferentially close. Bpa – the precautionary biomass required to protect the
#' spawning stock – is used to calculate MPA size. MPA is then grown to add
#' subsequent species starting from the most conservationally at-risk species,
#' resulting in one MPA map per species, and a multicolour MPA map of all.
#' All maps list the percentage of the avoid-variable’s total that is overlapped
#' by the MPA in the map legend.
#'
#' @param dbase data.frame to load. Expects Lon, Lat & data columns: predicted abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
#' @param loncolno Column number in dbase which has longitudes
#' @param latcolno Column number in dbase which has latitudes
#' @param goodcols Which column numbers are abundances (where higher = better)? List them in order of highest conservation importance first e.g. c(3,1,2,4)
#' @param badcols Which col no.s are 'negative' e.g. fishing (where higher = worse)?
#' @param conservecol Conservation column, from gbm.conserve
#' @param plotthis To plot? delete any,or all w/ NULL
#' @param maploops Sort loops to run
#' @param savethis Export all data as csv?
#' @param HRMSY Maximum % of each goodcols stock which can be removed yearly, as decimal (0.15 = 15%). Must protect remainder: 1-HRMSY. Single number or vector. Same order as goodcols
#' @param goodweight Single/vector weighting multiple(s) for goodcols array
#' @param badweight Ditto for badcols array
#' @param m Multiplication factor for Bpa units. 1000 to convert tonnes to kilos, 0.001 kilos to tonnes. Assumedly the same for all goodcols.
#' @param alerts Play sounds to mark progress steps
#' @param BnW Also produce greyscale images for print publications
#' @param mapshape Set coastline shapefile, else uses British Isles. Generate your own with gbm.basemap
#' @param pngtype Filetype for png files, alternatively try "quartz"
#' @param ... Optional terms for gbm.map
#'
#' @return Species abundance, abundance vs 'avoid variable', and MPA maps per
#' species and sort type, in b&w if set. CSVs of all maps if set.
#'
#' @export
#'
#' @examples None
gbm.valuemap <- function(
  dbase,  # data.frame to load. Expects Lon, Lat & data columns: predicted
# abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
  loncolno = 1, # column number in dbase which has longitudes
  latcolno = 2, # column number in dbase which has latitudes
  goodcols,  # which column numbers are abundances (where higher = better)? List
# them in order of highest conservation importance first e.g. c(3,1,2,4)
  badcols,  # which col no.s are 'negative' e.g. fishing (where higher = worse)?
  conservecol = NULL, # conservation column, from gbm.conserve
  plotthis = c("good","bad","both","close"), #to plot? delete any,or all w/ NULL
  maploops = c("Combo","Biomass","Effort","Conservation"), # sort loops to run
  savethis = TRUE, # export all data as csv?
  HRMSY = NULL, # maximum % of each goodcols stock which can be removed yearly,
# as decimal (0.15 = 15%). Must protect remainder: 1-HRMSY.
# Single number or vector. Same order as goodcols
  goodweight = NULL,  # single/vector weighting multiple(s) for goodcols array
  badweight = NULL,  # ditto for badcols array
  m = 1, # multiplication factor for Bpa units. 1000 to convert tonnes to kilos,
# 0.001 kilos to tonnes. Assumedly the same for all goodcols.
  alerts = TRUE,  # play sounds to mark progress steps
  BnW = TRUE,  # also produce greyscale images for print publications
  mapshape = NULL, #  set coastline shapefile, else
# uses British Isles. Generate your own with gbm.basemap
  pngtype = "cairo-png", # filetype for png files, alternatively try "quartz"
  ...) {  # optional terms for gbm.map

# Check & load gbm.map
if (!exists("gbm.map")) {stop("you need to install gbm.map to run this function")}
  #require(gbm.map) #for mapping; can't use require on non-CRAN?
if (alerts) if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
  if (alerts) require(beepr) #for progress noises
if (alerts) options(error = function() {beep(9)})  # warn for fails

if (is.null(mapshape)) {
  if (!exists("gbm.basemap")) {stop("you need to install gbm.basemap to run this function")}
  mapshape <- gbm.basemap(bounds = c(min(dbase[,loncolno]),
                                     max(dbase[,loncolno]),
                                     min(dbase[,latcolno]),
                                     max(dbase[,latcolno])))
  #data(coast,package = "mapplots")
  #mapshape <- coast
  }

goodname = colnames(dbase)[goodcols] #the response varible(s) name(s), e.g. species name(s), or collective term if agglomerating >1 response variable. Single character string, not a vector. No spaces or terminal periods.
badname = colnames(dbase)[badcols] #ditto for badcols. Both of these moved out of function parameters list to foce user to specify in colnames

# Check goodweight & badweight are the same length as goodcols & badcols
if (!is.null(goodweight)) if (!length(goodweight) == length(goodcols)) stop("number of goodweights doesn't match number of goodcols")
if (!is.null(badweight)) if (!length(badweight) == length(badcols)) stop("number of badweights doesn't match number of badcols")

dbasecoln <- ncol(dbase)  # note original number of columns for use later

# Scale values to 1 and add as columns at end of 'data' then name them, in goodcols order NOT original order.
  for (i in 1:length(c(goodcols,badcols))) {
    dbase[,dbasecoln + i] <- dbase[,(c(goodcols,badcols))[i]]/max(dbase[,(c(goodcols,badcols))[i]],na.rm = TRUE)  # for each column to scale (i), bind a column to data with i/max(i)
    colnames(dbase)[dbasecoln + i] <- paste(names(dbase)[(c(goodcols,badcols))[i]], "_s", sep = "")} #name them

# If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
ifelse(!is.null(goodweight),
         gooddata <- as.matrix(dbase[,dbasecoln + seq(1,length(goodcols))]) %*% diag(goodweight,ncol = length(goodcols)), #diag converts vector to matrix
         gooddata <- dbase[,dbasecoln + seq(1,length(goodcols))])
ifelse(!is.null(badweight), # assuming there's only going to be one baddata column.
         baddata <- as.matrix(dbase[,dbasecoln + length(goodcols) + seq(1,length(badcols))]) %*% diag(badweight, ncol = length(badcols)), #is matrix
         baddata <- dbase[,dbasecoln + length(goodcols) + seq(1, length(badcols))])  #is numeric

####map baddata####
if ("bad" %in% plotthis) {
  print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       Mapping Baddata     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
  png(filename = paste("./",badname,"_Map.png",sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
  gbm.map(x = dbase[,loncolno],  # run gbm.map function with generated parameters
          y = dbase[,latcolno], z = baddata, mapmain = "Fishing Effort",
          species = "", legendtitle = "Effort Level")
  dev.off()

  if (BnW) {
    png(filename = paste("./",badname,"_Map_BnW.png",sep = ""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
    gbm.map(x = dbase[,loncolno],
            y = dbase[,latcolno], z = baddata,
            mapmain = "Fishing Effort", species = "",
            legendtitle = "Effort Level",
            landcol = grey.colors(1, start = 0.8, end = 0.8), #light grey. 0=black 1=white
            mapback = "white", heatcolours = grey.colors(8, start = 1, end = 0))
    dev.off()
    if (alerts) beep(2) # alert user of success
  }} # close badplot & BnW

badmax <- max(baddata) # invert baddata
ifelse((baddata - badmax <= 0), {baddata = -(baddata - badmax)}, {baddata = (baddata - badmax)})
bothdata <- gooddata + baddata # create bothdata: good+(inverted) bad. ncols=gooddata+baddata
dbase <- cbind(dbase,bothdata) # add bothdata to data
for (i in 1:length(goodcols)) { #name bothdata cols: append w (Weighted) c (Combined)
  colnames(dbase)[ncol(dbase) - length(goodcols) + i] <- paste(names(dbase)[goodcols[i]],"_swc",sep = "")}

####create DF for HRMSY masking maps in loop####
bothdatarange <- (ncol(dbase) - length(goodcols) + 1):ncol(dbase)
dbase[,bothdatarange] <- bothdata * m

####Start Species Loop 1####
for (j in 1:length(goodcols)) {  #loop through gooddata columns

####map gooddata####
if ("good" %in% plotthis) {
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Mapping Gooddata ",j,"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
png(filename = paste("./", goodname[j], "_Map.png", sep = ""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
gbm.map(x = dbase[,loncolno], y = dbase[,latcolno], z = dbase[,goodcols[j]],
        species = goodname[j])
dev.off()

if (BnW) {
  png(filename = paste("./",goodname[j],"_Map_BnW.png",sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
  gbm.map(x = dbase[,loncolno], y = dbase[,latcolno], z = dbase[,goodcols[j]],
          species = goodname[j],
          landcol = grey.colors(1, start = 0.8, end = 0.8),
          mapback = "white",
          heatcolours = grey.colors(8, start = 1, end = 0))
  dev.off()
  if (alerts) beep(2) # alert user of success
}} # close goodplot & BnW

####map bothdata####
  #i.e. predabund scaled to 1 + effort scaled to 1 inverted.
if ("both" %in% plotthis) {
 print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Mapping Bothdata ",j,"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
  png(filename = paste("./",goodname[j]," plus ", badname, "_Map.png", sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
  gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
          z = dbase[,bothdatarange[j]],
          mapmain = "Predicted Abundance + Fishing Effort: ",
          species = goodname[j],
          heatcolours = c("red", "lightyellow","blue"),
          legendtitle = "0 - 2", byxout = TRUE)
  dev.off()
  byx <- byxport
  byy <- byx

  if (BnW) {
    png(filename = paste("./",goodname[j]," plus ",badname,"_Map_BnW.png",sep = ""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
    gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
            z = dbase[,bothdatarange[j]],
            mapmain = "Predicted Abundance + Fishing Effort: ",
            species = goodname[j], legendtitle = "0 - 2",
            byxout = TRUE,
            landcol = grey.colors(1, start = 0.8, end = 0.8),
            mapback = "white",
            heatcolours = grey.colors(3, start = 1, end = 0))
    dev.off()} # close BnW
  if (alerts) beep(2) } # alert user, close bothplot
  } # close species loop 1

if ("close" %in% plotthis) {
# set parameters for the sorting loops
maploopnames <- c("Combo","Biomass","Effort","Conservation")
# 1: max to min bothdata value for that species: highest cpue lowest e combo
# 2: max to min gooddata value for that species, i.e. highest biomass.
# 3: min to max baddata value (universal) THEN max to min gooddata value for that species, i.e. lowest effort THEN highest biomass
# 4: max to min conserve value (universal), i.e. max conservation value (all species)
maploopcodes <- c("dbase[order(-dbase[,bothdatarange[1]-1+j]),]","dbase[order(-dbase[,goodcols[j]]),]","dbase[order(dbase[,badcols],-dbase[,goodcols[j]]),]","dbase[order(-dbase[,conservecol]),]")
loopindex <- which(maploopnames == maploops) # which loops to run? set by user
maploopnames <- maploopnames[loopindex] # update names for only those loops
maploopcodes <- maploopcodes[loopindex] # update codes for only those loops

for (o in 1:length(maploopnames)) { # start o loop through maploops
j <- 1 # set / reset J so it restarts the loops at 1 rather than max

####Start Species Loop 2####
for (j in 1:length(goodcols)) {  # j loop through gooddata (species) columns
  n <- nrow(dbase) # set current row at start @ final row for HRMSY subloop: resets value for next species.
  ####HRMSY limit map####
  # min area for HRMSY, starting with best fish cells
  # Sort by bothdata descending then map that then overlay 15% biomass by highest bothdata
  CPUEMSY <- (sum(dbase[,goodcols[j]]) * HRMSY[j]) # HRMSY% * total biomass for J'th species = biomass to protect i.e. 'sum-to threshold'
  # CPUEMSY won't exist if you don't run
  dbase <- eval(parse(text = maploopcodes[o])) #sort according to O'th maploopcode
   for (k in 1:nrow(dbase)) { # run the loop for every row, if required
      print(paste("Run ",((o - 1) * length(maploopcodes)) + j," of ", length(goodcols) * length(maploopcodes),"; ",n,", ",round((sum(dbase[n:nrow(dbase),goodcols[j]])/CPUEMSY) * 100, 3),"%",sep = ""))  # progress printer
          #if (sum(dbase[n:nrow(dbase),goodcols[j]]) < CPUEMSY) {n = n - 1} else {#if sum of rows from end to this point < CPUEMSY aka HRMSY% move 1 row up & resum, ELSE:
            n <- which.min((cumsum(dbase[nrow(dbase):1,goodcols[j]]) - CPUEMSY) ^ 2)
            # coilins <- which.min((cumsum(dbase[nrow(dbase):1,"CPUE"]) - CPUEMSY) ^ 2)
            print(n)
            # cumsum from end (n) upwards towards 1 of CPUE until the row X
            # where end:X = CPUEMSY i.e. 'open' cells; close rest i.e. n-X
            assign(paste("sort",maploopnames[o],"_",names(dbase)[goodcols[j]],sep = ""),rep(0,nrow(dbase))) # create a vector of zeroes called sort[j name]
            dbase <- cbind(dbase,get(paste("sort",maploopnames[o],"_",names(dbase)[goodcols[j]],sep = ""))) # bind it to dbase
            colnames(dbase)[ncol(dbase)] <- paste("sort",maploopnames[o],"_",names(dbase)[goodcols[j]],sep = "") # reinstate its name (lost because bound with get())
            dbase[1:n + 1,ncol(dbase)] <- rep(1,n) # populate first n+1 rows with 1 [k:n]
            #dbase[1:n, ncol(dbase)] <- rep(1,n) # populate first n+1 rows with 1 [k:n]
            badcut <- sum(dbase[1:n + 1,badcols]) # sum badcols values (i.e. total badcol value in closed area)
            #badcut <- sum(dbase[1:n, badcols]) # sum badcols values (i.e. total badcol value in closed area)
            badall <- sum(dbase[,badcols])    # total badcols values
            badpct <- round((badcut/badall)*100,1) # percent of badcols values in closed area
            # this counts UP from the WORST row (last) until reaching HRMSY% (e.g. 8%) then closes the inverse
            # This is massively quicker than counting DOWN from the best row to 1-HRMSY, e.g. counting through 92% of the data

            # INDIVIDUAL MAPS: Map bothdata. Then Overlay map onto bothdata map: ncol=1 zeroes=TRUE. Heatcol="black".
            png(filename = paste("./ClosedValueMap_",maploopnames[o],"_",goodname[j],".png",sep = ""), # map species j's bothdata with black closed areas overlaid
                width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
            par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)

            # run gbm.map function's internal code. Set parameters
            x = dbase[,loncolno] #vector of longitudes, from make.grid in mapplots
            y = dbase[,latcolno] #vector of latitudes, from make.grid in mapplots; grids[,gridslat]
            z = dbase[,bothdatarange[j]] # scaled & weighted (bothdata) *m
            heatcolours = c("red", "lightyellow","blue")
            colournumber = 8   #number of colours to spread heatcol over, default:8
            shape = mapshape   #basemap shape to draw, from draw.shape in mapplots
            landcol = "grey80" #colour for 'null' (land) area of map, from draw.shape in mapplots
            mapback = "lightblue" # basemap background (sea) colour
            legendloc = "bottomright" #location on map of legend box, from legend.grid in mapplots
            legendtitle = "0 - 2" #the metric of abundance, e.g. CPUE for fisheries, from legend.grid in mapplots
            lejback = "white"  #backgroud colour of legend, from legend.grid in mapplots
            zero = TRUE # allow 0 category in breaks.grid & thus legend?
            quantile = 1 # set max breakpoint; lower this to cutoff outliers

print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Overlay Map ",((o - 1)*length(maploopcodes)) + j," of ",length(goodcols)*length(maploopcodes),"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
              if (!exists("byx")) {    # if users hasn't entered byx or byy values, generate them from the data
                bydist <- rep(NA,length(x))   # work out cell size for uniform square gridded data: Create blank vector for grid length calcs
                cells <- data.frame(LONGITUDE = x, bydist = bydist, stringsAsFactors = FALSE)   # and attach it to grids
                cells[2:(length(x) - 1),"bydist"] <- ifelse(round(cells[2:(length(x) - 1),1] - cells[1:(length(x) - 2),1],digits = 5) ==
                  round(cells[3:length(x),1] - cells[2:(length(x) - 1),1], digits = 5),round(cells[2:(length(x) - 1),1] - cells[1:(length(x) - 2),1], digits = 5), NA)
                byx <- mean(cells$bydist, na.rm = TRUE)  # Take an average of those distances, they should all be identical anyway. Apply it to byx & byy.
                byy <- byx}

              # Plot first map: same as bothdata
              grd <- make.grid(x, y, z, byx, byy, xlim = range(x), ylim = range(y),fun = mean) #create gridded data. fun defaults to sum which is bad
              heatcol = colorRampPalette(heatcolours)(colournumber) #create heatcol from component parts
              breaks <- breaks.grid(grd, zero = zero, quantile = quantile, ncol = length(heatcol))  #if breaks specified, do nothing (it'll be used later). Else generate it.
              if (zero) {heatcol = c("#00000000", colorRampPalette(heatcol)(length(heatcol) - 1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
              basemap(xlim = range(x), ylim = range(y), main = paste(maploopnames[o], "-Sorted Closed Area: ", goodname[j], sep = ""), bg = mapback, xlab = "Longitude (°W)", ylab = "Latitude (°N)")
####remove xlab & ylab above for general code####
              draw.grid(grd, breaks, col = heatcol) # plot grd data w/ breaks for colour breakpoints

              # Plot second map: closed area overlay only
              x = dbase[,loncolno] #set parameters for closed area map
              y = dbase[,latcolno]
              z = dbase[,ncol(dbase)] #last column, cbound @ L187, sort column, zeroes populated with 1s.
              heatcolours2 = c("black","black")
              colournumber2 = 2
              grd2 <- make.grid(x, y, z, byx, byy, xlim = range(x), ylim = range(y), fun = mean) #create gridded data. fun defaults to sum which is bad
              heatcol2 = colorRampPalette(heatcolours2)(colournumber2) #create heatcol from component parts
              breaks2 <- breaks.grid(grd, zero = zero, quantile = quantile, ncol = length(heatcol2))  #if breaks specified, do nothing (it'll be used later). Else generate it.
              if (zero) {heatcol2 = c("#00000000", colorRampPalette(heatcol2)(length(heatcol2) - 1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
              draw.grid(grd2, breaks2, col = heatcol2) # plot grd data w/ breaks for colour breakpoints
              draw.shape(shape = shape, col = landcol) # add coastline
              legend.grid(legendloc, breaks = breaks, type = 2, inset = 0, bg = lejback, title = paste(badpct, "% E closed", sep = ""), col = heatcol)
            dev.off()

            # create csvs for closed area comparisons
            ca_x <- dbase[,loncolno]
            ca_y <- dbase[,latcolno]
            ca_z <- dbase[,ncol(dbase)]
            ca_df <- data.frame(lon = ca_x, lat = ca_y, closed = ca_z)
            write.csv(ca_df, row.names = FALSE, file = paste("./ClosedArea_", maploopnames[o], "_", goodname[j], ".csv", sep = ""))

            if (BnW) {
              png(filename = paste("./ClosedValueMap_",maploopnames[o],"_",goodname[j],"_BnW.png",sep = ""),
                  width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
              par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)

              # run gbm.map function's internal code. Set parameters
              x = dbase[,loncolno]
              y = dbase[,latcolno]
              z = dbase[,bothdatarange[j]] # scaled & weighted (bothdata) *m
              heatcolours = grey.colors(3, start = 1, end = 0)
              colournumber = 8
              shape = mapshape
              landcol = grey.colors(1, start = 0.8, end = 0.8)
              mapback = "white"
              legendloc = "bottomright"
              legendtitle = "0 - 2"
              lejback = "white"
              zero = TRUE
              quantile = 1

              if (!exists("byx")) {    # if users hasn't entered byx or byy values, generate them from the data
                bydist <- rep(NA, length(x))   # work out cell size for uniform square gridded data: Create blank vector for grid length calcs
                cells <- data.frame(LONGITUDE = x, bydist = bydist, stringsAsFactors = FALSE)   # and attach it to grids
                cells[2:(length(x) - 1), "bydist"] <- ifelse(round(cells[2:(length(x) - 1), 1] - cells[1:(length(x) - 2), 1], digits = 5) ==
                                                            round(cells[3:length(x),1] - cells[2:(length(x) - 1), 1], digits = 5), round(cells[2:(length(x) - 1), 1] - cells[1:(length(x) - 2), 1], digits = 5), NA)
                byx <- mean(cells$bydist, na.rm = TRUE)  # Take an average of those distances, they should all be identical anyway. Apply it to byx & byy.
                byy <- byx}

              # Plot first map: same as bothdata
              grd <- make.grid(x, y, z, byx, byy, xlim = range(x), ylim = range(y), fun = mean) #create gridded data. fun defaults to sum which is bad
              heatcol = colorRampPalette(heatcolours)(colournumber) #create heatcol from component parts
              breaks <- breaks.grid(grd, zero = zero, quantile = quantile, ncol = length(heatcol))  #if breaks specified, do nothing (it'll be used later). Else generate it.
              if (zero) {heatcol = c("#00000000", colorRampPalette(heatcol)(length(heatcol) - 1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
              basemap(xlim = range(x), ylim = range(y), main = paste(maploopnames[o], "-Sorted Closed Area: ", goodname[j], sep = ""), bg = mapback, xlab = "Longitude (°W)", ylab = "Latitude (°N)")
####remove xlab & ylab above for general code####
              draw.grid(grd, breaks, col = heatcol)

              # Plot second map: closed area overlay only
              x = dbase[,loncolno]
              y = dbase[,latcolno]
              z = dbase[,ncol(dbase)]
              heatcolours2 = c("black","black")
              colournumber2 = 2
              grd2 <- make.grid(x, y, z, byx, byy, xlim = range(x), ylim = range(y), fun = mean)
              heatcol2 = colorRampPalette(heatcolours2)(colournumber2)
              breaks2 <- breaks.grid(grd, zero = zero, quantile = quantile, ncol = length(heatcol2))
              if (zero) {heatcol2 = c("#00000000",colorRampPalette(heatcol2)(length(heatcol2) - 1))}
              draw.grid(grd2, breaks2, col = heatcol2)
              draw.shape(shape = shape, col = landcol)
              legend.grid(legendloc, breaks = breaks, type = 2, inset = 0, bg = lejback, title = paste(badpct, "% E closed", sep = ""), col = heatcol)
              dev.off()
            } # close BnW

         #break} # end of ELSE section
          } # end of FOR loop k (data rows)
  if (alerts) beep(2)} # alert user & end of 2nd FOR loop j (species)

####Build closed areas####
  SortCol0 <- ncol(dbase) - length(goodcols) # last col before goodcols
  # sort largest to smallest species 1 value, then 2, 3 4
  # loop through the columns in reverse, sort them 4, 3, 2, 1
for (q in length(goodcols):1) {dbase <- dbase[order(-dbase[,SortCol0 + q]),]}

# then add a column which is max(k:n). This will be 1s and 0s and will be the full extent (try append)
  assign(paste("AllClosed_", maploopnames[o], sep = ""),
         do.call(pmax, dbase[,(SortCol0 + 1):(SortCol0 + length(goodcols))]))
         # do.call allows df format for pmax
 dbase <- cbind(dbase, get(paste("AllClosed_", maploopnames[o],sep = "")))
 # reinstate its name (lost because bound with get())
 colnames(dbase)[ncol(dbase)] <- paste("AllClosed_", maploopnames[o],sep = "")

# then add a column which is sum(k:n). This will be 0:length(goodcols) and will be the full extent. Both are similar
assign(paste("SumClosed_", maploopnames[o], sep = ""), rowSums(dbase[,(SortCol0 + 1):(SortCol0 + length(goodcols))]))
dbase <- cbind(dbase,get(paste("SumClosed_", maploopnames[o], sep = "")))
colnames(dbase)[ncol(dbase)] <- paste("SumClosed_", maploopnames[o], sep = "")

# then create a new df with a column of zeroes, length = nrow(dbase)
assign(paste("Zeroes_", maploopnames[o], sep = ""), rep(0, nrow(dbase)))
MPAgrow <- as.data.frame(get(paste("Zeroes_", maploopnames[o], sep = "")))

# then create a vec of the max of zeroes & species 1 value
for (r in 1:length(goodcols)) { # loop through 1:length(goodcols)
  assign(paste("Closure", r, "_", maploopnames[o], sep = ""), pmax(dbase[,SortCol0 + r], MPAgrow[,r]))
  # bind that to the zeroes & name it
  MPAgrow <- cbind(MPAgrow, get(paste("Closure", r, "_", maploopnames[o], sep = "")))
  colnames(MPAgrow)[ncol(MPAgrow)] <- paste("Closure", r, "_", maploopnames[o], sep = "")}

dbase <- cbind(dbase, MPAgrow[,2:ncol(MPAgrow)]) # bind to dbase, removing zeroes

# loop through final length(goodcols) cols of dbase & map them: cumulative CPUEMSY limit maps
for (l in 1:length(goodcols)) {
  # Generate badcol reduction percentage
  badcoldata <- dbase[,badcols] #vector of badcol data
  closecol <- dbase[, ncol(dbase) - length(goodcols) + l] #closures column
  badcut2 <- sum(badcoldata[closecol == 1]) #sum of badcol data in closed cells
  badpct2 <- round((badcut2/badall) * 100, 1) # badcols percent

####Cumulative closed area maps####
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Cumulative Closed Area Map ",((o - 1) * length(maploopcodes)) + l," of ",length(goodcols)*length(maploopcodes)," XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
  png(filename = paste("./CumulativeClosedArea",maploopnames[o],"Map_",goodname[l],".png",sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)

  gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
          z = dbase[,ncol(dbase) - length(goodcols) + l],
          mapmain = "Cumulative Closed Area: ",
          species = paste(maploopnames[o]," ",goodname[l],sep = ""),
          heatcolours = c("black","black"), colournumber = 2,
          shape = mapshape, mapback = "white",
          legendtitle = paste(badpct2,"% E closed",sep = ""),
          byx = byx, byy = byy)
    dev.off()

    if (alerts) beep(2)} # alert user & end l loop (col no.s MPAgrow species)

#SpeciesGrow: like MPAgrow but records the species number only when it grows the closed area (instead of max to 1 or sum)
SpeciesGrow <- rep(0, nrow(dbase))
MPAgrow2 <- cbind(SpeciesGrow,dbase[,(SortCol0 + 1):(SortCol0 + length(goodcols))])
for (p in 1:length(goodcols)) {
  MPAgrow2[,p + 1] <- ifelse(MPAgrow2[,p + 1] == 1, rep(p, length(MPAgrow2[,p + 1])), MPAgrow2[,p + 1]) # replace 1s in sortcols with species numbers
  MPAgrow2[,1] <- ifelse(MPAgrow2[,1] == 0, MPAgrow2[,p + 1], MPAgrow2[,1])} # for areas in col1 not already closed, that are closed in this species' col, put the species no.

dbase <- cbind(dbase, MPAgrow2[,1]) # add completed SpeciesGrow column (MPAgrow2 col1 to dbase. Is single column.
colnames(dbase)[ncol(dbase)] <- paste("SpeciesGrow_", maploopnames[o], sep = "") # Name the column

####Per species closed area maps####
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Per Species Closed Area Map ",o," of ",length(goodcols)," XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
png(filename = paste("./PerSpeciesClosedArea",maploopnames[o],"Map.png",sep = ""),
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)

gbm.map(x = dbase[,loncolno], y = dbase[,latcolno], z = dbase[,ncol(dbase)],
        mapmain = "Per Species Closed Area: ", species = maploopnames[o],
        heatcolours = c("black",rainbow(length(goodcols) - 1)), #black then rainbow colours. Total n is length(goodcols): blank, black, then goodcols-1 other colours
        colournumber = length(goodcols) + 1, # Colours are set as the breaks, but painted from the midpoints
        shape = mapshape, mapback = "white",
        legendtitle = paste(badpct2,"% E closed",sep = ""),
        byx = byx, byy = byy, breaks = c(0,0:length(goodcols)))
dev.off()

if (BnW) {
  png(filename = paste("./PerSpeciesClosedArea_BnW", maploopnames[o], "Map.png", sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)

  gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
          z = dbase[,ncol(dbase)], mapmain = "Per Species Closed Area: ",
          species = maploopnames[o],
          heatcolours = grey.colors(8, start = 0.7, end = 0),
          colournumber = length(goodcols) + 1,
          shape = mapshape, mapback = "white",
          legendtitle = paste(badpct2, "% E closed", sep = ""),
          byx = byx, byy = byy, breaks = c(0,0:length(goodcols)))
    dev.off()} # close BnW
if (alerts) beep(2)} # alert user & end of "close" optional section
} # end o loop through combination, goodcols & badcols

####Save csvs####
print(paste("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         Saving CSV        XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",sep = ""))
if (savethis) write.csv(dbase,row.names = FALSE, file = paste("./ProcessedData.csv", sep = ""))
beep(8)} # notify user & close function
