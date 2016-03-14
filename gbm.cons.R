gbm.cons <- function(grids,         # csv file (inc relative location) of gridded lat+long+data to predict to
                     subsets,       # Subset name(s): character; single or vector
                     conssamples,   # single or vector of samples data csv files
                     # (inc relative location) corresponding to subsets
                     alerts = TRUE, # play sounds to mark progress steps
                     map = TRUE,    # produce maps
                     BnW = TRUE,    # also produce B&W maps
                     expvars,  # list object of expvar vectors for gbm.autos,
                     # length = no. of subsets * no. of species. No default
                     resvars,  # vector of resvars cols from conssamples for gbm.autos, length(subsets)*species, no default
                     tcs = NULL, # autocalculated below if not provided by user
                     lrs = rep(list(c(0.01,0.005)),length(resvars)),
                     bfs = rep(0.5, length(resvars)),
                     ZIs = rep("CHECK", length(resvars)),
                     gridslats = rep(2, length(resvars)),
                     gridslons = rep(1, length(resvars)),
                     colss = rep(list(grey.colors(1,1,1)),length(resvars)),
                     linesfiless = rep(FALSE, length(resvars)),
                     savegbms = rep(TRUE, length(resvars)),
                     varints = rep(TRUE, length(resvars)),
                     maps = rep(TRUE, length(resvars)),
                     RSBs = rep(TRUE, length(resvars)),
                     BnWs = rep(TRUE, length(resvars)),
                     zeroes = rep(TRUE, length(resvars))){

  ####todo: make running gbm.auto optional####
  # if they've already been run.
  # Have to have subset folders & species folder names correct.
  # test this. Changes default requirement of grids. And samples? And loads of stuff.

  ## gbm.todo.P2
  # working script for paper 2; 1/6/2015

  # Generalised Boosting Model / Boosted Regression Tree process chain automater.
  # Simon Dedman, 2014, simondedman@gmail.com, https://github.com/SimonDedman/gbm.auto

  # Function to automate the many steps required to use boosted regression trees to predict abundances in a delta process,
  # i.e. binary (0/1) proportion prediction coupled with presence-only abundance prediction to give total prediction.
  # Loops through all permutations of parameters provided (learning rate, tree complexity, bag fraction), chooses the best,
  # then tries to simplify that. Generates line, dot & bar plots, and outputs these and the predictions and a report of all
  # variables used, statistics for tests, variable interactions, predictors used & dropped, etc.. If selected, generates
  # predicted abundance maps, and Unrepresentativeness surfaces.
  #
  # Underlying functions are from packages gbm and dismo, functions from Elith et al. 2008 (bundled as gbm.utils.R), mapplots,
  # and my own functions gbm.map, gbm.rsb, gbm.valuemap


####Load functions & data####
# load gmb.utils, containing Elith's packages not in dismo, & SD's gbm functions
#### require check these, have user load em####
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.auto.R") # loads beepr (but needed earlier), gbm.map gbm.rsb, gbm.utils
  # maybe just try gbm.auto check and that will then launch the other package checks when it gets there? except beepr.
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.utils.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.rsb.R")
source("/home/simon/Dropbox/Galway/Analysis/R/gbm.auto/gbm.map.R")

if (map) if (!exists("gbm.map")) {stop("you need to install the gbm.map function to run this function")}
if (!require(beepr)) {stop("you need to install the beepr package to run this function")}
library(beepr)
if (alerts) options(error = function() {beep(9)})  # give warning noise if it fails

if (is.null(tcs)) {tcs = list() #make blank then loop populate w/ 2 & expvar length
for (g in 1:length(resvars)) {tcs[[g]] <- c(2,length(expvars[[g]]))}}

# load saved models if re-running aspects from a previous run
# load("Bin_Best_Model")
# load("Gaus_Best_Model")

mygrids <<- read.csv(grids, header = TRUE)  # load grids, <<- bad but reqd
dir.create("ConservationMaps") # create conservation maps directory

# create a list of response variables for name ranges
GS <- length(resvars)/length(subsets) # calculate group size, e.g. 8/2
resvarrange = list() # create blank list
for (h in 1:length(subsets)) {  #e.g. 1:2
  resvarrange[[h]] <- (1 + (GS * (h - 1))):(GS * h)} # e.g. 1:4, 5:8

####gbm.auto loops subsets & species####
# Loop through subsets
for (i in 1:length(subsets)) {  #currently 2

# name samples by subset name & load
assign(subsets[i], read.csv(conssamples[i],header = TRUE, row.names = NULL))

dir.create(paste("./", subsets[i], sep = "")) # Create WD for subset[i] name
setwd(paste("./", subsets[i], sep = ""))  # go there

for (j in 1:GS) {  #loop through all species in group e.g. 4

# below: ((i-1)*GS)+j is the subset-loop-disregarding number
# e.g. CTBS,CTBS = 1:8. Allows (e.g.) length=8 lists to be entered by user in
# function terms & used by gbm.auto calls across subset loops

mysamples <<- get(subsets[i]) # for gbm.auto (default) & later, <<- bad but reqd
# gbm.auto pulls relevant group-ignoring variable from user entries or defaults
gbm.auto(grids = mygrids,
         expvar = expvars[[((i - 1) * GS) + j]],
         resvar = resvars[[((i - 1) * GS) + j]],
         tc = tcs[[((i - 1) * GS) + j]],
         lr = lrs[[((i - 1) * GS) + j]],
         bf = bfs[[((i - 1) * GS) + j]],
         ZI = ZIs[[((i - 1) * GS) + j]],
         gridslat = gridslats[[((i - 1) * GS) + j]],
         gridslon = gridslons[[((i - 1) * GS) + j]],
         cols = colss[[((i - 1) * GS) + j]],
         linesfiles = linesfiless[[((i - 1) * GS) + j]],
         savegbm = savegbms[[((i - 1) * GS) + j]],
         varint = varints[[((i - 1) * GS) + j]],
         map = maps[[((i - 1) * GS) + j]],
         RSB = RSBs[[((i - 1) * GS) + j]],
         BnW = BnWs[[((i - 1) * GS) + j]],
         zero = zeroes[[((i - 1) * GS) + j]])
if (alerts) beep(2) # ping for each completion

# create object for resulting abundance preds csv, e.g. Juveniles_Cuckoo
assign(paste(subsets[i], "_", names(mysamples)[(resvars)[resvarrange[[i]]]][j], sep = ""),
       read.csv(paste("./", names(mysamples)[(resvars)[resvarrange[[i]]]][j], "/Abundance_Preds_only.csv", sep = ""), header = TRUE))
# names(mysamples)[(resvars)[resvarrange[[i]]]][j] is the (species) name for the
# column no. in samples, for the j'th response variable in this subsets' group

} # reloop/end j loop of species
print(paste("XXXXXXXXXXXXXXXXXXXXXX           Species ", j, " of ", GS, ", Subset ", i, " of ", length(subsets), "           XXXXXXXXXXXXXXXXXXXXXXXXXX", sep = ""))
setwd("../") # go back up to /Maps root folder for correct placement @ restart
} # reloop/end i loop of subsets

####Conservation maps####
# Loop through subsets then add them together

for (k in names(mysamples)[(resvars)[resvarrange[[length(subsets)]]]]) {
  # name list from last mysamples & last subset resvarnames list, e.g. CTBS
  # make grids-length blank objects e.g. All_Cuckoo, Scaled_Cuckoo & allscaled
  assign(paste("All_", k, sep = ""),rep(0,dim(mygrids)[1]))
  assign(paste("Scaled_", k, sep = ""),rep(0,dim(mygrids)[1]))
  allscaled <- rep(0,dim(mygrids)[1])

  # loop subsets
  for (i in 1:length(subsets)) {
    # replace All_Cuckoo (starts blank) w/ All_Cuckoo + e.g. Juveniles_Cuckoo
    assign(paste("All_", k, sep = ""),get(paste("All_", k, sep = "")) + get(paste(subsets[i], "_", k, sep = ""))[,3])

    # scale subsets' values to 1 for species k & add to blanks
    assign(paste("Scaled_", k, sep = ""),
           get(paste("Scaled_", k, sep = ""))
           +
           (get(paste(subsets[i], "_", k, sep = ""))[,3]
           /
           max(get(paste(subsets[i], "_", k, sep = ""))[,3], na.rm = TRUE)))

    xtmp <- get(paste(subsets[i], "_", k, sep = ""))[,2] # LONG for later
    ytmp <- get(paste(subsets[i], "_", k, sep = ""))[,1] # LAT for later
    } # end/reloop i & add next subset of same species e.g. Adult Females_Cuckoo
print(paste("XXXXXXXXXXXXXXXXXXXXXX          Both subsets scaled for ", k, "          XXXXXXXXXXXXXXXXXXXXXXXXXX", sep = ""))

  # create simple temp object name for e.g. All_Cuckoo
  ztmp <- get(paste("All_", k, sep = "")) # already includes the [,3]

  dir.create(paste("./ConservationMaps/", k, sep = ""))
  setwd(paste("./ConservationMaps/", k, sep = ""))

####Unscaled conservation maps####
  # map all subsets' abundance for species k
  if (map) {
  png(filename = paste("./Conservation_Map_", k, ".png", sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
      bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0),xpd = FALSE)
  gbm.map(x = xtmp,
          y = ytmp,
          z = ztmp,
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = k,
          zero = FALSE,
          legendtitle = "CPUE")
  dev.off()

  if (BnW) {  # again in B&W
  png(filename = paste("./Conservation_Map_BnW", k, ".png", sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
      bg = "white", res = NA, family = "", type = "cairo-png")
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
          legendtitle = "CPUE")
  dev.off()}} # close BnW & mapping IFs
  if (alerts) beep(2)
print(paste("XXXXXXXXXXXXXXXXXXXXXX      Unscaled conservation map(s) generated       XXXXXXXXXXXXXXXXXXXXXXXXXX", sep = ""))

####Scaled-to-1 conservation maps####
  if (map) {
  png(filename = paste("./Scale1-1_Conservation_Map_",k,".png",sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
      bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = xtmp,
          y = ytmp,
          z = (get(paste("Scaled_", k, sep = ""))) * (100/length(subsets)),
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = k,
          zero = FALSE,
          breaks = c(0,20,40,60,80,100),
          colournumber = 5,
          legendtitle = "CPUE (Scaled %)")
  dev.off()

  if (BnW) {   # again in B&W
  png(filename = paste("./Scale1-1_Conservation_Map_BnW_",k,".png", sep = ""),
      width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
      bg = "white", res = NA, family = "", type = "cairo-png")
  par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
  gbm.map(x = xtmp,
          y = ytmp,
          z = (get(paste("Scaled_", k, sep = ""))) * (100/length(subsets)),
          mapmain = "Predicted CPUE (numbers per hour): ",
          species = k,
          zero = FALSE,
          breaks = c(0,20,40,60,80,100),
          colournumber = 5,
          heatcolours = grey.colors(5, start = 1, end = 0),
          mapback = "white",
          legendtitle = "CPUE (Scaled %)")
  dev.off()}} # close BnW & mapping IF
  setwd("../") # go back up to ConservationMaps
  setwd("../")} # go back up to Maps & end/reloop k for next species
if (alerts) beep(2) # ping on completion
print(paste("XXXXXXXXXXXXXXXXXXXXXX       Scaled conservation map(s) generated        XXXXXXXXXXXXXXXXXXXXXXXXXX", sep = ""))

####Add scaled outputs, all species####
dir.create(paste("./ConservationMaps/Combo/", sep = ""))
setwd(paste("./ConservationMaps/Combo/", sep = ""))

# loop & add each species' combined scaled values
for (l in names(mysamples)[(resvars)[resvarrange[[length(subsets)]]]]) {
  allscaled <- allscaled + get(paste("Scaled_", l, sep = ""))}

#save as csv
allscaleddf <- data.frame(LATITUDE = ytmp, LONGITUDE = xtmp, allscaled)
write.csv(allscaleddf, row.names = FALSE, file = paste("./AllScaledData.csv", sep = ""))

if (map) {
png(filename = "Scaled_Conservation_Map.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
    bg = "white", res = NA, family = "", type = "cairo-png")
par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
gbm.map(x = xtmp,
        y = ytmp,
        z = allscaled * (100/length(resvars)), # multiplier raises to 100
        mapmain = "Predicted CPUE (numbers per hour): ",
        species = "All Species",
        zero = FALSE,
        legendtitle = "CPUE (Scaled %)")
dev.off()

if (BnW) { # again in B&W
png(filename = "Scaled_Conservation_Map_BnW.png",
    width = 4*1920, height = 4*1920, units = "px", pointsize = 4*48,
    bg = "white", res = NA, family = "", type = "cairo-png")
par(mar = c(3.2,3,1.3,0), las = 1, mgp = c(2.1,0.5,0), xpd = FALSE)
gbm.map(x = xtmp,
        y = ytmp,
        z = allscaled * (100/length(resvars)),
        mapmain = "Predicted CPUE (numbers per hour): ",
        species = "All Species",
        zero = FALSE,
        heatcolours = grey.colors(5, start = 1, end = 0),
        mapback = "white",
        legendtitle = "CPUE (Scaled %)")
dev.off()}} # close BnW & mapping IF
print(paste("XXXXXXXXXXXXXXXXXXXXXX     All-scaled conservation map(s) generated      XXXXXXXXXXXXXXXXXXXXXXXXXX", sep = ""))
print(paste("XXXXXXXXXXXXXXXXXXXXXX                Everything complete                XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX", sep = ""))
if (alerts) beep(8)} # complete sound & close function
####END####