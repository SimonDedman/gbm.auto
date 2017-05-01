#' Conservation Area Mapping
#'
#' Runs gbm.auto for multiple subsets of the same overall dataset and
#' scales the combined results, leading to maps which highlight areas of high
#' conservation importance for multiple species in the same study area e.g.
#' using juvenile and adult female subsets to locate candidate nursery grounds
#' and spawning areas respectively.
#'
#' @param mygrids Gridded lat+long+data object to predict to
#' @param subsets Subset name(s): character; single or vector, corresponding to
#' matching-named dataset objects e.g. read in by read.csv()
#' @param alerts Play sounds to mark progress steps
#' @param map Produce maps
#' @param BnW Also produce B&W maps?
#' @param resvars Vector of resvars cols from dataset objects for gbm.autos, length(subsets)*species, no default
#' @param gbmautos Do gbm.auto runs for species?
#' @param expvars List object of expvar vectors for gbm.autos, length = no. of subsets * no. of species. No default
#' @param tcs Gbm.auto parameters, autocalculated below if not provided by user
#' @param lrs Gbm.auto parameter, uses defaults if not provided by user
#' @param bfs Gbm.auto parameter, uses defaults if not provided by user
#' @param ZIs Gbm.auto parameter, autocalculated below if not provided by user
#' @param colss Gbm.auto parameter, uses defaults if not provided by user
#' @param linesfiless Gbm.auto parameter, uses defaults if not provided by user
#' @param savegbms Gbm.auto parameter, uses defaults if not provided by user
#' @param varints Gbm.auto parameter, uses defaults if not provided by user
#' @param maps Gbm.auto parameter, uses defaults if not provided by user
#' @param RSBs Gbm.auto parameter, uses defaults if not provided by user
#' @param BnWs Gbm.auto parameter, uses defaults if not provided by user
#' @param zeroes For breaks.grid, include zero-only category in colour
#' breakpoints and subsequent legend. Defaults to TRUE
#' @param shape Coastline file for gbm.map
#' @param pngtype Filetype for png files, alternatively try "quartz"
#' @param gridslat per Gbm.auto defaults to 2
#' @param gridslon per Gbm.auto defaults to 1
#'
#' @return  Maps via gbm.map & saved data as csv file
#' @export
#' @importFrom beepr beep
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples
#' gbm.cons(grids = mygrids, subsets = c("Juveniles","Adult_Females"),
#'          resvars = c(44:47,11:14),
#'          expvars = list(c(4:11,15,17,21,25,29,37),
#'                         c(4:11,15,18,22,26,30,38),
#'                         c(4:11,15,19,23,27,31),
#'                         c(4:11,15,20,24,28,32,39),
#'                         4:10, 4:10, 4:10, 4:10),
#'          tcs = list(c(2,14), c(2,14), 13, c(2,14), c(2,6), c(2,6), 6, c(2,6)),
#'          lrs = list(c(0.01,0.005), c(0.01,0.005), 0.005, c(0.01,0.005),
#'                0.005, 0.005, 0.001, 0.005),
#'          ZIs = rep(TRUE, 8),
#'          savegbms = rep(FALSE, 8),
#'          varints = rep(FALSE, 8),
#'          RSBs = rep(FALSE, 8),
#'          BnWs = rep(FALSE, 8),
#'          zeroes = rep(FALSE,8))
#'
gbm.cons <- function(mygrids,       # gridded lat+long+data object to predict to
                     subsets,       # Subset name(s): character; single or vector
                     # corresponding to matching-named dataset objects e.g. read in by read.csv()
                     alerts = TRUE, # play sounds to mark progress steps
                     map = TRUE,    # produce maps
                     BnW = TRUE,    # also produce B&W maps
                     resvars,  # vector of resvars cols from dataset objects for
                     # gbm.autos, length(subsets)*species, no default
                     gbmautos = TRUE, # do gbm.auto runs for species?
                     expvars,  # list object of expvar vectors for gbm.autos,
                     # length = no. of subsets * no. of species. No default
                     tcs = NULL, # autocalculated below if not provided by user
                     lrs = rep(list(c(0.01,0.005)),length(resvars)),
                     bfs = rep(0.5, length(resvars)),
                     ZIs = rep("CHECK", length(resvars)),
                     colss = rep(list(grey.colors(1,1,1)),length(resvars)),
                     linesfiless = rep(FALSE, length(resvars)),
                     savegbms = rep(TRUE, length(resvars)),
                     varints = rep(TRUE, length(resvars)),
                     maps = rep(TRUE, length(resvars)),
                     RSBs = rep(TRUE, length(resvars)),
                     BnWs = rep(TRUE, length(resvars)),
                     zeroes = rep(TRUE, length(resvars)),
                     shape = NULL, # coastline file for gbm.map
                     pngtype = "cairo-png", # filetype for png files, alternatively try "quartz"
                     gridslat = 2, #per Gbm.auto defaults to 2
                     gridslon = 1) #per Gbm.auto defaults to 1
{
  ####todo: make running gbm.auto optional####
  # if they've already been run.
  # Have to have subset folders & species folder names correct.
  # test this. Changes default requirement of grids. And samples? And loads of stuff.

####Load functions & data####
if (map) if (!exists("gbm.map")) {stop("you need to install the gbm.map function to run this function")}
if (is.null(shape)) {if (!exists("gbm.basemap")) {stop("you need to install gbm.basemap to run this function")}}
if (alerts) if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
  if (alerts) require(beepr)
if (alerts) options(error = function() {beep(9)})  # give warning noise if it fails
if (gbmautos) {if (is.null(tcs)) {tcs = list() #make blank then loop populate w/ 2 & expvar length
for (g in 1:length(resvars)) {tcs[[g]] <- c(2,length(expvars[[g]]))}}}

# load saved models if re-running aspects from a previous run
# load("Bin_Best_Model")
# load("Gaus_Best_Model")

# create basemap if not provided
  if (map) {
    if (is.null(shape)) {
      bounds = c(range(grids[,gridslon]),range(grids[,gridslat]))
      #create standard bounds from data, and extra bounds for map aesthetic
      xmid <- mean(bounds[1:2])
      ymid <- mean(bounds[3:4])
      xextramax <- ((bounds[2] - xmid) * 1.6) + xmid
      xextramin <- xmid - ((xmid - bounds[1]) * 1.6)
      yextramax <- ((bounds[4] - ymid) * 1.6) + ymid
      yextramin <- ymid - ((ymid - bounds[3]) * 1.6)
      extrabounds <- c(xextramin, xextramax, yextramin, yextramax)
      shape <- gbm.basemap(bounds = extrabounds)
    } else {shape <- shape}} # use user-defined otherwise; close map section

# create a list of response variables for name ranges
GS <- length(resvars)/length(subsets) # calculate group size, e.g. 8/2
resvarrange = list() # create blank list to populate with subset values
for (h in 1:length(subsets)) {  #e.g. 1:2
  resvarrange[[h]] <- (1 + (GS * (h - 1))):(GS * h)} # e.g. 1:4, 5:8

####gbm.auto loops subsets & species####
# Loop through subsets
for (i in 1:length(subsets)) {  #currently 2
  if (gbmautos) {dir.create(paste0("./", subsets[i]))} # Create WD for subset[i] name
  setwd(paste0("./", subsets[i]))  # go there

  for (j in 1:GS) {  #loop through all species in group e.g. 4

    # below: ((i-1)*GS)+j is the subset-loop-disregarding number
    # e.g. CTBS,CTBS = 1:8. Allows (e.g.) length=8 lists to be entered by user in
    # function terms & used by gbm.auto calls across subset loops

    mysamples <- get(subsets[i]) # for gbm.auto (default) & later
    # gbm.auto pulls relevant group-ignoring variable from user entries or defaults
    if (gbmautos) {gbm.auto(grids = mygrids,
                            samples = mysamples,
                            expvar = expvars[[((i - 1) * GS) + j]],
                            resvar = resvars[[((i - 1) * GS) + j]],
                            tc = tcs[[((i - 1) * GS) + j]],
                            lr = lrs[[((i - 1) * GS) + j]],
                            bf = bfs[[((i - 1) * GS) + j]],
                            ZI = ZIs[[((i - 1) * GS) + j]],
                            gridslat = gridslat,
                            gridslon = gridslon,
                            cols = colss[[((i - 1) * GS) + j]],
                            linesfiles = linesfiless[[((i - 1) * GS) + j]],
                            savegbm = savegbms[[((i - 1) * GS) + j]],
                            varint = varints[[((i - 1) * GS) + j]],
                            map = maps[[((i - 1) * GS) + j]],
                            RSB = RSBs[[((i - 1) * GS) + j]],
                            BnW = BnWs[[((i - 1) * GS) + j]],
                            zero = zeroes[[((i - 1) * GS) + j]],
                            shape = shape)
      if (alerts) beep(2)} # ping for each completion

    # create object for resulting abundance preds csv, e.g. Juveniles_Cuckoo
    assign(paste0(subsets[i], "_", names(mysamples)[(resvars)[resvarrange[[i]]]][j]),
           read.csv(paste0("./", names(mysamples)[(resvars)[resvarrange[[i]]]][j], "/Abundance_Preds_only.csv"), header = TRUE))
    # names(mysamples)[(resvars)[resvarrange[[i]]]][j] is the (species) name for the
    # column no. in samples, for the j'th response variable in this subsets' group
    print(paste0("XXXXXXXXXXXXXXXXXXXXXX           Species ", j, " of ", GS, ", Subset ", i, " of ", length(subsets), "           XXXXXXXXXXXXXXXXXXXXXXXXXX"))
  } # reloop/end j loop of species
  setwd("../") # go back up to /Maps root folder for correct placement @ restart
} # reloop/end i loop of subsets

####Conservation maps####
dir.create("ConservationMaps") # create conservation maps directory
# Loop through subsets then add them together
for (k in names(mysamples)[(resvars)[resvarrange[[length(subsets)]]]]) {
  # name list from last mysamples & last subset resvarnames list, e.g. CTBS
  # make grids-length blank objects e.g. All_Cuckoo, Scaled_Cuckoo & allscaled
  assign(paste0("All_", k),rep(0,dim(mygrids)[1]))
  assign(paste0("Scaled_", k),rep(0,dim(mygrids)[1]))
  allscaled <- rep(0,dim(mygrids)[1])

  # loop subsets
  for (i in 1:length(subsets)) {
    # replace All_Cuckoo (starts blank) w/ All_Cuckoo + e.g. Juveniles_Cuckoo
    assign(paste0("All_", k),get(paste0("All_", k)) + get(paste0(subsets[i], "_", k))[,3])
    # scale subsets' values to 1 for species k & add to blanks
    assign(paste0("Scaled_", k),
           get(paste0("Scaled_", k))
           + (get(paste0(subsets[i], "_", k))[,3]
              / max(get(paste0(subsets[i], "_", k))[,3], na.rm = TRUE)))
    ###could have option here to normalise another way than scaling to 1

    xtmp <- get(paste0(subsets[i], "_", k))[,2] # LONG for later
    ytmp <- get(paste0(subsets[i], "_", k))[,1] # LAT for later
  } # end/reloop i & add next subset of same species e.g. Adult Females_Cuckoo

  print(paste0("XXXXXXXXXXXXXXXXXXXXXX          Both subsets scaled for ", k, "          XXXXXXXXXXXXXXXXXXXXXXXXXX"))
  # create simple temp object name for e.g. All_Cuckoo
  ztmp <- get(paste0("All_", k)) # already includes the [,3]
  dir.create(paste0("./ConservationMaps/", k))
  setwd(paste0("./ConservationMaps/", k))

  ####Unscaled conservation maps####
  # map all subsets' abundance for species k
  if (map) {
    png(filename = paste0("./Conservation_Map_", k, ".png"),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
        bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
    gbm.map(x = xtmp,
            y = ytmp,
            z = ztmp,
            mapmain = "Predicted CPUE (numbers per hour): ",
            species = k,
            zero = FALSE,
            legendtitle = "CPUE",
            shape = shape) # set coast shapefile, else downloaded and autogenerated
    dev.off()

    if (BnW) {  # again in B&W
      png(filename = paste0("./Conservation_Map_BnW_", k, ".png"),
          width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
          bg = "white", res = NA, family = "", type = pngtype)
      par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
      gbm.map(x = xtmp,
              y = ytmp,
              z = ztmp,
              mapmain = "Predicted CPUE (numbers per hour): ",
              species = k,
              zero = FALSE,
              colournumber = 5,
              heatcolours = grey.colors(5, start = 1, end = 0),
              mapback = "white",
              legendtitle = "CPUE",
              shape = shape)
      dev.off()
    } # close BnW
    if (alerts) beep(2)
    print(paste0("XXXXXXXXXXXXXXXXXXXXXX   Unscaled ", k, " conservation map(s) generated  XXXXXXXXXXXXXXXXXXXXXXXXXXX"))
  } # close mapping IF

  ####Scaled-to-1 conservation maps####
  if (map) {
    png(filename = paste0("./Scale1-1_Conservation_Map_",k,".png",sep = ""),
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
        bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
    gbm.map(x = xtmp,
            y = ytmp,
            z = (get(paste0("Scaled_", k))) * (100/length(subsets)),
            mapmain = "Predicted CPUE (numbers per hour): ",
            species = k,
            zero = FALSE,
            breaks = c(0,20,40,60,80,100),
            colournumber = 5,
            legendtitle = "CPUE (Scaled %)",
            shape = shape)
    dev.off()

    if (BnW) {   # again in B&W
      png(filename = paste0("./Scale1-1_Conservation_Map_BnW_",k,".png"),
          width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
          bg = "white", res = NA, family = "", type = pngtype)
      par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
      gbm.map(x = xtmp,
              y = ytmp,
              z = (get(paste0("Scaled_", k))) * (100/length(subsets)),
              mapmain = "Predicted CPUE (numbers per hour): ",
              species = k,
              zero = FALSE,
              breaks = c(0,20,40,60,80,100),
              colournumber = 5,
              heatcolours = grey.colors(5, start = 1, end = 0),
              mapback = "white",
              legendtitle = "CPUE (Scaled %)",
              shape = shape)
      dev.off()
    } # close BnW
    if (alerts) beep(2) # ping on completion
    print(paste0("XXXXXXXXXXXXXXXXXXXXXX    Scaled ", k, " conservation map(s) generated  XXXXXXXXXXXXXXXXXXXXXXXXXXX"))
  } # close mapping IF

  setwd("../") # go back up to ConservationMaps
  setwd("../")
} # go back up to Maps & end/reloop k for next species

####Add scaled outputs, all species####
dir.create(paste0("./ConservationMaps/Combo/"))
setwd(paste0("./ConservationMaps/Combo/"))

# loop & add each species' combined scaled values
for (l in names(mysamples)[(resvars)[resvarrange[[length(subsets)]]]]) {
  allscaled <- allscaled + get(paste0("Scaled_", l))}

#save as csv
allscaleddf <- data.frame(LATITUDE = ytmp, LONGITUDE = xtmp, allscaled)
write.csv(allscaleddf, row.names = FALSE, file = paste0("./AllScaledData.csv"))

if (map) {
  png(filename = "Scaled_Conservation_Map.png",
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
      bg = "white", res = NA, family = "", type = pngtype)
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = xtmp,
          y = ytmp,
          z = allscaled * (100/length(resvars)), # multiplier raises to 100
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = "All Species",
          zero = FALSE,
          legendtitle = "CPUE (Scaled %)",
          shape = shape)
  dev.off()

  if (BnW) { # again in B&W
    png(filename = "Scaled_Conservation_Map_BnW.png",
        width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
        bg = "white", res = NA, family = "", type = pngtype)
    par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
    gbm.map(x = xtmp,
            y = ytmp,
            z = allscaled * (100/length(resvars)),
            mapmain = "Predicted CPUE (numbers per hour): ",
            species = "All Species",
            zero = FALSE,
            heatcolours = grey.colors(5, start = 1, end = 0),
            mapback = "white",
            legendtitle = "CPUE (Scaled %)",
            shape = shape)
    dev.off()
  } # close BnW
  if (alerts) beep(2) # ping on completion
  print(paste0("XXXXXXXXXXXXXXXXXXXXXX     All-scaled conservation map(s) generated      XXXXXXXXXXXXXXXXXXXXXXXXXX"))
} # close mapping IF

setwd("../../") # return to original directory
print(paste0("XXXXXXXXXXXXXXXXXXXXXX                Everything complete                XXXXXXXXXXXXXXXXXXXXXXXXXX"))
if (alerts) beep(8) # complete sound & close function
####END####
}
