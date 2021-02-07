#' Decision Support Tool that generates (Marine) Protected Area options using
#' species predicted abundance maps
#'
#' Scales response variable data, maps a user-defined explanatory variable to be
#' avoided, e.g. fishing effort, combines them into a map showing areas to
#' preferentially close. Bpa, the precautionary biomass required to protect the
#' spawning stock, is used to calculate MPA size. MPA is then grown to add
#' subsequent species starting from the most conservationally at-risk species,
#' resulting in one MPA map per species, and a multicolour MPA map of all.
#' All maps list the percentage of the avoid-variables total that is overlapped
#' by the MPA in the map legend.
#'
#' @param dbase Data.frame to load. Expects Lon, Lat & data columns: predicted
#' abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto.
#' @param loncolno Column number in dbase which has longitudes.
#' @param latcolno Column number in dbase which has latitudes.
#' @param goodcols Which column numbers are abundances (where higher = better)?
#' List them in order of highest conservation importance first e.g. c(3,1,2,4).
#' @param badcols Which col no.s are 'negative' e.g. fishing (where higher =
#' worse)?
#' @param conservecol Conservation column, from gbm.cons.
#' @param plotthis To plot? delete any,or all w/ NULL.
#' @param maploops Sort loops to run.
#' @param savethis Export all data as csv?
#' @param HRMSY Maximum percent of each goodcols stock which can be removed
#' yearly, as decimal (0.15 = 15 pct). Must protect remainder: 1-HRMSY. Single
#' number or vector. If vector, same order as goodcols. Required.
#' @param goodweight Single/vector weighting multiple(s) for goodcols array.
#' @param badweight Ditto for badcols array.
#' @param m Multiplication factor for Bpa units. 1000 to convert tonnes to
#' kilos, 0.001 kilos to tonnes. Assumedly the same for all goodcols.
#' @param alerts Play sounds to mark progress steps.
#' @param BnW Also produce greyscale images for print publications.
#' @param shape Set coastline shapefile, else uses British Isles. Generate your
#' own with gbm.basemap.
#' @param pngtype Filetype for png files, alternatively try "quartz".
#' @param byxport Dummy param for package testing for CRAN, ignore.
#' @param ... Optional terms for gbm.map.
#'
#' @return Species abundance, abundance vs avoid variable, and MPA maps per
#' species and sort type, in b&w if set. CSVs of all maps if set.
#'
#' @export
#' @import mapplots
#' @importFrom beepr beep
#' @importFrom grDevices colorRampPalette dev.off grey.colors png rainbow
#' @importFrom graphics image legend par
#' @importFrom utils write.csv
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
gbm.valuemap <- function(
  dbase,  # data.frame to load. Expects Lon, Lat & data columns: predicted
  # abundances, fishing effort etc. E.g.: Abundance_Preds_All.csv from gbm.auto
  loncolno = 1, # column number in database which has longitudes
  latcolno = 2, # column number in database which has latitudes
  goodcols,  # which column numbers are abundances (where higher = better)? List
  # them in order of highest conservation importance first e.g. c(3,1,2,4)
  badcols,  # which col no.s are 'negative' e.g. fishing (where higher = worse)?
  conservecol = NULL, # conservation column, from gbm.cons
  plotthis = c("good","bad","both","close"), #to plot? delete any,or all w/ NULL
  maploops = c("Combo","Biomass","Effort","Conservation"), # sort loops to run
  savethis = TRUE, # export all data as csv?
  HRMSY = 0.15, # maximum % of each goodcols stock which can be removed yearly,
  # as decimal (0.15 = 15%). Must protect remainder: 1-HRMSY. Single number or
  # vector. If vector, same order as goodcols. Required.
  goodweight = NULL,  # single/vector weighting multiple(s) for goodcols array
  badweight = NULL,  # ditto for badcols array
  m = 1, # multiplication factor for Bpa units. 1000 to convert tonnes to kilos,
  # 0.001 kilos to tonnes. Assumedly the same for all goodcols.
  alerts = TRUE,  # play sounds to mark progress steps
  BnW = TRUE,  # also produce greyscale images for print publications
  shape = NULL, #  set coastline shapefile, else
  # uses British Isles. Generate your own with gbm.basemap
  pngtype = "cairo-png", # filetype for png files, alternatively try "quartz"
  byxport = NULL, # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
  ...) {  # optional terms for gbm.map

  # utils::globalVariables("byxport") # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals

  # Check & load gbm.map
  # if (!exists("gbm.map")) {stop("you need to install gbm.map to run this function")}
  #require(gbm.map) #for mapping; can't use require on non-CRAN?
  # if (alerts) if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
  # if (alerts) require(beepr) #for progress noises
  oldpar <- par(no.readonly = TRUE) # defensive block, thanks to Gregor Sayer
  oldoptions <- options()
  on.exit(par(oldpar))
  on.exit(options(oldoptions), add = TRUE)
  if (alerts) options(error = function() {
    beep(9)
    graphics.off()})  # give warning noise if it fails

  if ("Conservation" %in% maploops & is.null(conservecol)) stop("conservecol must be specified if Conservation presesent in maploops")

  if (is.null(shape)) {
    # if (!exists("gbm.basemap")) {stop("you need to install gbm.basemap to run this function")}
    bounds = c(range(dbase[,loncolno]),range(dbase[,latcolno]))
    #create standard bounds from data, and extra bounds for map aesthetic
    shape <- gbm.basemap(bounds = bounds, extrabounds = TRUE)
  } # close isnull shape

  goodname = colnames(dbase)[goodcols] #the response variable(s) name(s), e.g. species name(s), or collective term if agglomerating >1 response variable. Single character string, not a vector. No spaces or terminal periods.
  badname = colnames(dbase)[badcols] #ditto for badcols. Both of these moved out of function parameters list to foce user to specify in colnames

  # Check goodweight & badweight are the same length as goodcols & badcols
  if (!is.null(goodweight)) if (!length(goodweight) == length(goodcols)) stop("number of goodweights doesn't match number of goodcols")
  if (!is.null(badweight)) if (!length(badweight) == length(badcols)) stop("number of badweights doesn't match number of badcols")

  dbasecoln <- ncol(dbase)  # note original number of columns for use later

  # Scale values to 1 and add as columns at end of 'data' then name them, in goodcols order NOT original order.
  for (i in 1:length(c(goodcols,badcols))) {
    dbase[,dbasecoln + i] <- dbase[,(c(goodcols,badcols))[i]]/max(dbase[,(c(goodcols,badcols))[i]],na.rm = TRUE)  # for each column to scale (i), bind a column to data with i/max(i)
    colnames(dbase)[dbasecoln + i] <- paste0(names(dbase)[(c(goodcols,badcols))[i]], "_s")} #name them

  # If weighting factors given, multiply then sum scaled values & create objects, else sum & create objects
  if (!is.null(goodweight)) gooddata <- as.matrix(dbase[,dbasecoln + seq(1, length(goodcols))]) %*% diag(goodweight, ncol = length(goodcols)) #diag converts vector to matrix
  if (is.null(goodweight)) gooddata = dbase[,dbasecoln + seq(1,length(goodcols))]
  # assuming there's only going to be one baddata column.
  if (!is.null(badweight)) baddata = as.matrix(dbase[,dbasecoln + length(goodcols) + seq(1,length(badcols))]) %*% diag(badweight, ncol = length(badcols)) #is matrix
  if (is.null(badweight)) baddata = dbase[,dbasecoln + length(goodcols) + seq(1, length(badcols))]  #is numeric

  ####map baddata####
  if ("bad" %in% plotthis) {
    print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       Mapping Baddata     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
    png(filename = paste0("./",badname,"_Map.png"),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
    gbm.map(x = dbase[,loncolno],  # run gbm.map function with generated parameters
            y = dbase[,latcolno], z = baddata, mapmain = "Fishing Effort",
            species = "", legendtitle = "Effort Level", shape = shape)
    dev.off()

    if (BnW) {
      png(filename = paste0("./",badname,"_Map_BnW.png"),
          width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
      par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
      gbm.map(x = dbase[,loncolno],
              y = dbase[,latcolno], z = baddata,
              mapmain = "Fishing Effort", species = "",
              legendtitle = "Effort Level",
              shape = shape,
              landcol = grey.colors(1, start = 0.8, end = 0.8), #light grey. 0=black 1=white
              mapback = "white", heatcolours = grey.colors(8, start = 1, end = 0))
      dev.off()
      if (alerts) beep(2) # alert user of success
    } # close BnW
  } # close badplot

  badmax <- max(baddata) # invert baddata
  ifelse((baddata - badmax <= 0), {baddata = -(baddata - badmax)}, {baddata = (baddata - badmax)})
  bothdata <- gooddata + baddata # create bothdata: good+(inverted) bad. ncols=gooddata+baddata
  dbase <- cbind(dbase,bothdata) # add bothdata to data
  for (i in 1:length(goodcols)) { #name bothdata cols: append w (Weighted) c (Combined)
    colnames(dbase)[ncol(dbase) - length(goodcols) + i] <- paste0(names(dbase)[goodcols[i]],"_swc")}

  ####create DF for HRMSY masking maps in loop####
  bothdatarange <- (ncol(dbase) - length(goodcols) + 1):ncol(dbase)
  dbase[,bothdatarange] <- bothdata * m

  ####Start Species Loop 1####
  for (j in 1:length(goodcols)) {  #loop through gooddata columns

    ####map gooddata####
    if ("good" %in% plotthis) {
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Mapping Gooddata ",j,"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      png(filename = paste0("./", goodname[j], "_Map.png"),
          width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
      par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
      gbm.map(x = dbase[,loncolno], y = dbase[,latcolno], z = dbase[,goodcols[j]],
              species = goodname[j], shape = shape)
      dev.off()

      if (BnW) {
        png(filename = paste0("./",goodname[j],"_Map_BnW.png"),
            width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
        gbm.map(x = dbase[,loncolno], y = dbase[,latcolno], z = dbase[,goodcols[j]],
                species = goodname[j],
                landcol = grey.colors(1, start = 0.8, end = 0.8),
                mapback = "white",
                heatcolours = grey.colors(8, start = 1, end = 0),
                shape = shape)
        dev.off()
        if (alerts) beep(2) # alert user of success
      } # close bnw
    } # close goodplot

    ####map bothdata####
    #i.e. predabund scaled to 1 + effort scaled to 1 inverted.
    if ("both" %in% plotthis) {
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Mapping Bothdata ",j,"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      png(filename = paste0("./",goodname[j]," plus ", badname, "_Map.png"),
          width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
      par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
      gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
              z = dbase[,bothdatarange[j]],
              mapmain = "Predicted Abundance + Fishing Effort: ",
              species = goodname[j],
              heatcolours = c("red", "lightyellow","blue"),
              legendtitle = "0 - 2", byxout = TRUE,
              shape = shape)
      dev.off()
      # byx <- byxport # byxport exported <<- from gbm.map
      # byy <- byx
      # if byxport is null, these create null byx byy causing issues later. Code later is copied from gbm.map to do this anyway

      if (BnW) {
        png(filename = paste0("./",goodname[j]," plus ",badname,"_Map_BnW.png"),
            width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
        gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
                z = dbase[,bothdatarange[j]],
                mapmain = "Predicted Abundance + Fishing Effort: ",
                species = goodname[j], legendtitle = "0 - 2",
                byxout = TRUE,
                landcol = grey.colors(1, start = 0.8, end = 0.8),
                mapback = "white",
                heatcolours = grey.colors(3, start = 1, end = 0),
                shape = shape)
        dev.off()
      } # close BnW
      if (alerts) beep(2) # alert user
    } # close bothplot
  } # close species loop 1, for j

  if ("close" %in% plotthis) {
    # set parameters for the sorting loops
    maploopnames <- c("Combo","Biomass","Effort","Conservation")
    # 1: max to min bothdata value for that species: highest cpue lowest e combo
    # 2: max to min gooddata value for that species, i.e. highest biomass.
    # 3: min to max baddata value (universal) THEN max to min gooddata value for that species, i.e. lowest effort THEN highest biomass
    # 4: max to min conserve value (universal), i.e. max conservation value (all species)
    maploopcodes <- c("dbase[order(-dbase[,bothdatarange[1]-1+j]),]",
                      "dbase[order(-dbase[,goodcols[j]]),]",
                      "dbase[order(dbase[,badcols],-dbase[,goodcols[j]]),]",
                      "dbase[order(-dbase[,conservecol]),]")
    loopindex <- which(maploopnames %in% maploops) # which loops to run? set by user
    maploopnames <- maploopnames[loopindex] # update names for only those loops
    maploopcodes <- maploopcodes[loopindex] # update codes for only those loops
    counterA <- 1 # counter for overlay map plots
    counterB <- 1 # counter for cumulative closed area map plots
    for (o in 1:length(maploopnames)) { # start o loop through maploops
      j <- 1 # set / reset J so it restarts the loops at 1 rather than max(length(goodcols))

      ####Start Species Loop 2####
      for (j in 1:length(goodcols)) {  # j loop through gooddata (species) columns
        n <- nrow(dbase) # set current row at start @ final row for HRMSY subloop: resets value for next species.
        ####HRMSY limit map####
        # min area for HRMSY, starting with best fish cells
        # Sort by bothdata descending then map that then overlay 15% biomass by highest bothdata
        CPUEMSY <- (sum(dbase[,goodcols[j]]) * HRMSY[j]) # HRMSY% * total biomass for J'th species = biomass to protect i.e. 'sum-to threshold'
        dbase <- eval(parse(text = maploopcodes[o])) #sort according to O'th maploopcode
        n <- which.min((cumsum(dbase[nrow(dbase):1,goodcols[j]]) - CPUEMSY) ^ 2)
        # cumsum from end (n) upwards towards 1 of CPUE until the row X
        # where end:X = CPUEMSY i.e. 'open' cells; close rest i.e. n-X
        assign(paste0("sort",maploopnames[o],"_",names(dbase)[goodcols[j]]),rep(0,nrow(dbase))) # create a vector of zeroes called sort[j name]
        dbase <- cbind(dbase,get(paste0("sort",maploopnames[o],"_",names(dbase)[goodcols[j]]))) # bind it to dbase
        colnames(dbase)[ncol(dbase)] <- paste0("sort",maploopnames[o],"_",names(dbase)[goodcols[j]]) # reinstate its name (lost because bound with get())
        dbase[1:n + 1,ncol(dbase)] <- rep(1,n) # populate first n+1 rows with 1 [k:n]
        badcut <- sum(dbase[1:n + 1,badcols]) # sum badcols values (i.e. total badcol value in closed area)
        badall <- sum(dbase[,badcols])    # total badcols values
        badpct <- round((badcut/badall)*100,1) # percent of badcols values in closed area
        # this counts UP from the WORST row (last) until reaching HRMSY% (e.g. 8%) then closes the inverse
        # This is massively quicker than counting DOWN from the best row to 1-HRMSY, e.g. counting through 92% of the data
        # INDIVIDUAL MAPS: Map bothdata. Then Overlay map onto bothdata map: ncol=1 zeroes=TRUE. Heatcol="black".
        png(filename = paste0("./ClosedValueMap_",maploopnames[o],"_",goodname[j],".png"), # map species j's bothdata with black closed areas overlaid
            width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)

        # run gbm.map function's internal code. Set parameters
        x = dbase[,loncolno] #vector of longitudes, from make.grid in mapplots
        y = dbase[,latcolno] #vector of latitudes, from make.grid in mapplots; grids[,gridslat]
        z = dbase[,bothdatarange[j]] # scaled & weighted (bothdata) *m
        heatcolours = c("red", "lightyellow","blue")
        colournumber = 8   #number of colours to spread heatcol over, default:8
        shape = shape   #basemap shape to draw, from draw.shape in mapplots
        landcol = "grey80" #colour for 'null' (land) area of map, from draw.shape in mapplots
        mapback = "lightblue" # basemap background (sea) colour
        legendloc = "bottomright" #location on map of legend box, from legend.grid in mapplots
        legendtitle = "0 - 2" #the metric of abundance, e.g. CPUE for fisheries, from legend.grid in mapplots
        lejback = "white"  #backgroud colour of legend, from legend.grid in mapplots
        zero = TRUE # allow 0 category in breaks.grid & thus legend?
        quantile = 1 # set max breakpoint; lower this to cutoff outliers

        # print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Overlay Map ",((o - 1)*length(maploopcodes)) + j," of ",length(goodcols)*length(maploopcodes),"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     Overlay Map ", counterA," of ",length(goodcols)*length(maploopcodes),"    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        counterA <- counterA + 1 # increment counter
        if (!exists("byx")) {    # if users hasn't entered byx or byy values, generate them from the data
          bydist <- rep(NA,length(x))   # work out cell size for uniform square gridded data: Create blank vector for grid length calcs
          cells <- data.frame(LONGITUDE = x, bydist = bydist, stringsAsFactors = FALSE)   # and attach it to grids
          cells[2:(length(x) - 1),"bydist"] <- ifelse(round(cells[2:(length(x) - 1),1] - cells[1:(length(x) - 2),1],digits = 5) ==
                                                        round(cells[3:length(x),1] - cells[2:(length(x) - 1),1], digits = 5),round(cells[2:(length(x) - 1),1] - cells[1:(length(x) - 2),1], digits = 5), NA)
          byx <- mean(cells$bydist, na.rm = TRUE)  # Take an average of those distances, they should all be identical anyway. Apply it to byx & byy.
          byy <- byx
        } # close if (!exists("byx"))

        # Plot first map: same as bothdata
        grd <- make.grid(x, y, z, byx, byy, xlim = range(x), ylim = range(y),fun = mean) #create gridded data. fun defaults to sum which is bad
        heatcol = colorRampPalette(heatcolours)(colournumber) #create heatcol from component parts
        breaks <- breaks.grid(grd, zero = zero, quantile = quantile, ncol = length(heatcol))  #if breaks specified, do nothing (it'll be used later). Else generate it.
        if (zero) {heatcol = c("#00000000", colorRampPalette(heatcol)(length(heatcol) - 1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
        basemap(xlim = range(x), ylim = range(y), main = paste0(maploopnames[o], "-Sorted Closed Area: ", goodname[j]), bg = mapback, xlab = "Longitude", ylab = "Latitude")
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
        legend.grid(legendloc, breaks = breaks, type = 2, inset = 0, bg = lejback, title = paste0(badpct, "% E closed"), col = heatcol)
        dev.off()

        # create csvs for closed area comparisons
        ca_x <- dbase[,loncolno]
        ca_y <- dbase[,latcolno]
        ca_z <- dbase[,ncol(dbase)]
        ca_df <- data.frame(lon = ca_x, lat = ca_y, closed = ca_z)
        write.csv(ca_df, row.names = FALSE, file = paste0("./ClosedArea_", maploopnames[o], "_", goodname[j], ".csv"))

        if (BnW) {
          png(filename = paste0("./ClosedValueMap_",maploopnames[o],"_",goodname[j],"_BnW.png"),
              width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
          par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
          x = dbase[,loncolno]
          y = dbase[,latcolno]
          z = dbase[,bothdatarange[j]] # scaled & weighted (bothdata) *m
          heatcolours = grey.colors(3, start = 1, end = 0)
          colournumber = 8
          shape = shape
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
            byy <- byx
          } # close if (!exists("byx"))

          # Plot first map: same as bothdata
          grd <- make.grid(x, y, z, byx, byy, xlim = range(x), ylim = range(y), fun = mean) #create gridded data. fun defaults to sum which is bad
          heatcol = colorRampPalette(heatcolours)(colournumber) #create heatcol from component parts
          breaks <- breaks.grid(grd, zero = zero, quantile = quantile, ncol = length(heatcol))  #if breaks specified, do nothing (it'll be used later). Else generate it.
          if (zero) {heatcol = c("#00000000", colorRampPalette(heatcol)(length(heatcol) - 1))} #if zero=TRUE add alpha as 1st colour (1st 2 breakpoints)
          basemap(xlim = range(x), ylim = range(y), main = paste0(maploopnames[o], "-Sorted Closed Area: ", goodname[j]), bg = mapback, xlab = "Longitude", ylab = "Latitude")
          draw.grid(grd, breaks, col = heatcol)
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
          legend.grid(legendloc, breaks = breaks, type = 2, inset = 0, bg = lejback, title = paste0(badpct, "% E closed"), col = heatcol)
          dev.off()
        } # close BnW
        if (alerts) beep(2) # alert user
      } # end of 2nd FOR loop j (species)

      ####Build closed areas####
      SortCol0 <- ncol(dbase) - length(goodcols) # last col before goodcols
      # sort largest to smallest species 1 value, then 2, 3 4
      # loop through the columns in reverse, sort them 4, 3, 2, 1
      for (q in length(goodcols):1) {dbase <- dbase[order(-dbase[,SortCol0 + q]),]}

      # then add a column which is max(k:n). This will be 1s and 0s and will be the full extent (try append)
      if (length(goodcols) == 1) { # if only 1 goodcol then do.call fails since it wants a list
        assign(paste0("AllClosed_", maploopnames[o]), dbase[,(SortCol0 + q)])
      } else {
        assign(paste0("AllClosed_", maploopnames[o]),
               do.call(pmax, dbase[,(SortCol0 + 1):(SortCol0 + length(goodcols))]))
      }
      # do.call allows df format for pmax
      dbase <- cbind(dbase, get(paste0("AllClosed_", maploopnames[o])))
      # reinstate its name (lost because bound with get())
      colnames(dbase)[ncol(dbase)] <- paste0("AllClosed_", maploopnames[o])

      # then add a column which is sum(k:n). This will be 0:length(goodcols) and will be the full extent. Both are similar
      if (length(goodcols) == 1) { # if only 1 goodcol then rowSums fails since it wants a list
        assign(paste0("SumClosed_", maploopnames[o]), dbase[,(SortCol0 + q)])
      } else {
        assign(paste0("SumClosed_", maploopnames[o]), rowSums(dbase[,(SortCol0 + 1):(SortCol0 + length(goodcols))]))
      }
      dbase <- cbind(dbase,get(paste0("SumClosed_", maploopnames[o])))
      colnames(dbase)[ncol(dbase)] <- paste0("SumClosed_", maploopnames[o])

      # then create a new df with a column of zeroes, length = nrow(dbase)
      assign(paste0("Zeroes_", maploopnames[o]), rep(0, nrow(dbase)))
      MPAgrow <- as.data.frame(get(paste0("Zeroes_", maploopnames[o])))

      # then create a vec of the max of zeroes & species 1 value
      for (r in 1:length(goodcols)) { # loop through 1:length(goodcols)
        assign(paste0("Closure", r, "_", maploopnames[o]), pmax(dbase[,SortCol0 + r], MPAgrow[,r]))
        # bind that to the zeroes & name it
        MPAgrow <- cbind(MPAgrow, get(paste0("Closure", r, "_", maploopnames[o])))
        colnames(MPAgrow)[ncol(MPAgrow)] <- paste0("Closure", r, "_", maploopnames[o])
      } # lose for (r in 1:length(goodcols))

      dbase <- cbind(dbase, MPAgrow[,2:ncol(MPAgrow)]) # bind to dbase, removing zeroes
      for (r in 1:length(goodcols)) {colnames(dbase)[ncol(dbase) - length(goodcols) + r] <- paste0("Closure_", goodname[r], "_", maploopnames[o])} # Name the column(s)

      # loop through final length(goodcols) cols of dbase & map them: cumulative CPUEMSY limit maps
      for (l in 1:length(goodcols)) {
        # Generate badcol reduction percentage
        badcoldata <- dbase[,badcols] #vector of badcol data
        closecol <- dbase[, ncol(dbase) - length(goodcols) + l] #closures column
        badcut2 <- sum(badcoldata[closecol == 1]) #sum of badcol data in closed cells
        badpct2 <- round((badcut2/badall) * 100, 1) # badcols percent

        ####Cumulative closed area maps####
        # print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Cumulative Closed Area Map ",((o - 1) * length(maploopcodes)) + l," of ",length(goodcols)*length(maploopcodes)," XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Cumulative Closed Area Map ", counterB," of ",length(goodcols)*length(maploopcodes)," XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
        counterB <- counterB + 1 # increment counter

        png(filename = paste0("./CumulativeClosedArea",maploopnames[o],"Map_",goodname[l],".png"),
            width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)

        gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
                z = dbase[,ncol(dbase) - length(goodcols) + l],
                mapmain = "Cumulative Closed Area: ",
                species = paste0(maploopnames[o]," ",goodname[l]),
                heatcolours = c("black","black"), colournumber = 2,
                shape = shape, mapback = "white",
                legendtitle = paste0(badpct2,"% E closed"),
                byx = byx, byy = byy)
        dev.off()

        if (alerts) beep(2) # alert user
      } # end l loop (col no.s MPAgrow species)

      #SpeciesGrow: like MPAgrow but records the species number only when it grows the closed area (instead of max to 1 or sum)
      SpeciesGrow <- rep(0, nrow(dbase))
      MPAgrow2 <- cbind(SpeciesGrow,dbase[,(SortCol0 + 1):(SortCol0 + length(goodcols))])
      for (p in 1:length(goodcols)) {
        MPAgrow2[,p + 1] <- ifelse(MPAgrow2[,p + 1] == 1, rep(p, length(MPAgrow2[,p + 1])), MPAgrow2[,p + 1]) # replace 1s in sortcols with species numbers
        MPAgrow2[,1] <- ifelse(MPAgrow2[,1] == 0, MPAgrow2[,p + 1], MPAgrow2[,1]) # for areas in col1 not already closed, that are closed in this species' col, put the species no.
      } # close for (p in 1:length(goodcols))

      dbase <- cbind(dbase, MPAgrow2[,1]) # add completed SpeciesGrow column (MPAgrow2 col1 to dbase. Is single column.
      colnames(dbase)[ncol(dbase)] <- paste0("SpeciesGrow_", maploopnames[o]) # Name the column

      ####Per species closed area maps####
      print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Per Species Closed Area Map ",o," of ",length(maploopnames)," XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
      png(filename = paste0("./PerSpeciesClosedArea",maploopnames[o],"Map.png"),
          width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
      par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)

      gbm.map(x = dbase[,loncolno], y = dbase[,latcolno], z = dbase[,ncol(dbase)],
              mapmain = "Per Species Closed Area: ", species = maploopnames[o],
              heatcolours = c("black",rainbow(length(goodcols) - 1)), #black then rainbow colours. Total n is length(goodcols): blank, black, then goodcols-1 other colours
              colournumber = length(goodcols) + 1, # Colours are set as the breaks, but painted from the midpoints
              shape = shape, mapback = "white",
              legendtitle = paste0(badpct2,"% E closed"),
              byx = byx, byy = byy, breaks = c(0,0:length(goodcols)))
      dev.off()

      if (BnW) {
        png(filename = paste0("./PerSpeciesClosedArea_BnW", maploopnames[o], "Map.png"),
            width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48, bg = "white", res = NA, family = "", type = pngtype)
        par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)

        gbm.map(x = dbase[,loncolno], y = dbase[,latcolno],
                z = dbase[,ncol(dbase)], mapmain = "Per Species Closed Area: ",
                species = maploopnames[o],
                heatcolours = grey.colors(8, start = 0.7, end = 0),
                colournumber = length(goodcols) + 1,
                shape = shape, mapback = "white",
                legendtitle = paste0(badpct2, "% E closed"),
                byx = byx, byy = byy, breaks = c(0,0:length(goodcols)))
        dev.off()
      } # close BnW
      if (alerts) beep(2) # alert user
    } # close "for (o in 1:length(maploopnames)"
  } # close "if ("close" %in% plotthis)" end o loop through combination, goodcols & badcols

  ####Save csvs####
  print(paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         Saving CSV        XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
  if (savethis) write.csv(dbase,row.names = FALSE, file = paste0("./ProcessedData.csv"))
  beep(8) # notify user & close function
}
