#' Creates Basemaps for Gbm.auto mapping from your data range
#'
#' Downloads unzips crops & saves NOAAs global coastline shapefiles to user-set
#' box. Use for 'shape' in gbm.map. If downloading in RStudio uncheck
#' "Use secure download method for HTTP" in Tools > Global Options > Packages.
#' Simon Dedman, 2015/6 simondedman@gmail.com github.com/SimonDedman/gbm.auto
#'
#' @param bounds Region to crop to: c(xmin,xmax,ymin,ymax)
#' @param grids if bounds unspecified, name your grids database here
#' @param gridslat if bounds unspecified, specify which column in grids is latitude
#' @param gridslon if bounds unspecified, specify which column in grids is longitude
#' @param getzip Download & unpack GSHHS data to WD? "TRUE" else
#' absolute/relative reference to GSHHS_shp folder, inc that folder
#' @param zipvers GSHHS version, in case it updates. Please email developer (SD)
#'  if this is incorrect
#' @param savename Shapefile savename, no extension, default is "Crop_Map"
#' @param res Resolution, 1:5 (low:high) OR c,l,i,h,f (coarse, low,
#' intermediate, high, full) or "CALC" to calculate based on bounds
#' @param extrabounds Grow bounds 16pct each direction to expand rectangular
#' datasets basemaps over the entire square area created by basemap in mapplots

#'
#' @return basemap coastline file for gbm.map in gbm.auto. "cropshp"
#' SpatialPolygonsDataFramein in local environment & user-named files in
#' "CroppedMap" folder. Load later with maptools function:
#' MyMap <- readShapePoly("./CroppedMap/Crop_Map")
#'
#' @export
#' @import rgeos
#' @importFrom rgdal readOGR
#' @importFrom maptools writeSpatialShape
#' @importFrom raster crop
#' @importFrom shapefiles read.shapefile
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples
#' mybounds <- c(range(samples[,3]),range(samples[,2]))
#' gbm.basemap(bounds = mybounds, getzip = "./GSHHS_shp/", savename =
#' "My_Crop_Map", res = "f")
#' In this example GSHHS folder already downloaded to the working directory
#' hence I pointed getzip at that rather than having it download the zip again.
#'
#' @details errors and their origins:
#'
#' 1. Error in setwd(getzip) : cannot change working directory
#' If youve specified the location of the local GSHHS_shp folder, ensure youre
#' in the correct directory relative to it. This error means it looked for the
#' folder and couldnt find it.
#'
#' 2. Error in writeSpatialShape(cropshp, savename) x is aNULLobject, not a
#' compatible Spatial*DataFrame.
#' Ensure that your lats and longs are the the right way around
#'
gbm.basemap <- function(bounds = NULL, # region to crop to: c(xmin,xmax,ymin,ymax)
                        grids = NULL, # if bounds unspecified, name your grids
# database here
                        gridslat = NULL, # if bounds unspecified, specify which
# column in grids is latitude
                        gridslon = NULL, # if bounds unspecified, specify which
# column in grids is longitude
                        getzip = TRUE, # download & unpack GSHHS data to WD?
# "TRUE" else absolute/relative reference to GSHHS_shp folder, inc that folder
                        zipvers = "2.3.7", # GSHHS version, in case it updates
# Please email developer if this is incorrect
                        savename = "Crop_Map", #shapefile savename, no extension
                        res = "CALC", # resolution, 1:5 (low:high) OR c,l,i,h,f
# (coarse, low, intermediate, high, full) or "CALC" to calculate based on bounds
                        extrabounds = FALSE) { # grow bounds 16pct each direction
# to expand rectangular datasets basemaps over the entire square area created by
# basemap in mapplots

print(paste("if rgdal install fails in linux try: sudo apt-get install libgdal-dev && sudo apt-get install libproj-dev"))
if (!require(rgdal)) install.packages("rgdal")
  require(rgdal) # for readOGR
if (!require(rgeos)) install.packages("rgeos")
  require(rgeos) # subfunctions for rgdal & others
if (!require(raster)) install.packages("raster")
  require(raster) # for crop
if (!require(maptools)) install.packages("maptools")
  require(maptools) # for WriteSpatialShape
if (!require(shapefiles)) install.packages("shapefiles")
  require(shapefiles) # for read.shapefile
  ###improve these: check if installed, install if not else library####
startdir <- getwd() # record original directory

# if bounds is entered it's user below, else check grids & gridslat & gridslon
if(is.null(bounds)) {
  #check none of grids & gridslat & gridslon is null, if any are print message
  if(is.null(grids)) stop("if bounds is NULL grids needs to be specified")
  if(is.null(gridslat)) stop("if bounds is NULL gridslat needs to be specified")
  if(is.null(gridslon)) stop("if bounds is NULL gridslon needs to be specified")
  #check they're all the correct format
  if(!is.data.frame(grids)) stop("grids needs to be a data frame")
  if(!is.numeric(gridslat)) stop("gridslat needs to be a number")
  if(!is.numeric(gridslon)) stop("gridslon needs to be a number")
  # construct bounds from gridslat & gridslon ranges in grids
  bounds <- c(range(grids[,gridslon]), range(grids[,gridslat]))}

if (res == 1) res <- "c" # If res provided as number convert to letter
if (res == 2) res <- "l"
if (res == 3) res <- "i"
if (res == 4) res <- "h"
if (res == 5) res <- "f"

if (res == "CALC") { # Calculate res based on size of bounds
  scope <- max((bounds[2] - bounds[1]), (bounds[4] - bounds[3])) # distance of largest dimension, x or y
  if (scope >= 160) res <- "c" # bigger diff = larger map = lower res
  if (160 > scope & scope >= 70) res <- "l"
  if (70 > scope & scope >= 29) res <- "i"
  if (29 > scope & scope >= 9) res <- "h"
  if (9 > scope) res <- "f"}

ifelse(getzip == TRUE, { # download & unzip GSHGG if getzip = TRUE
  download.file(paste0("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-", zipvers, ".zip"), "GSHHG.zip")
  unzip("GSHHG.zip")
  setwd("./GSHHS_shp")}
  , setwd(getzip)) # else just setwd to there

setwd(paste("./", res, sep = "")) #setwd to res subfolder

if (extrabounds) { # grow bounds extents if requested
  bounds = c(range(grids[,gridslon]),range(grids[,gridslat]))
  xmid <- mean(bounds[1:2])
  ymid <- mean(bounds[3:4])
  xextramax <- ((bounds[2] - xmid) * 1.6) + xmid
  xextramin <- xmid - ((xmid - bounds[1]) * 1.6)
  yextramax <- ((bounds[4] - ymid) * 1.6) + ymid
  yextramin <- ymid - ((ymid - bounds[3]) * 1.6)
  bounds <- c(xextramin, xextramax, yextramin, yextramax)
}

# read in worldmap
world <- readOGR(dsn = paste0("GSHHS_", res, "_L1.shp"), layer = paste0("GSHHS_", res, "_L1"))
cropshp <- crop(world, bounds) # crop to extents
setwd(startdir)
dir.create("CroppedMap") # create conservation maps directory
setwd("CroppedMap")
writeSpatialShape(cropshp, savename)
print(paste("World map cropped and saved successfully"))
cropshp <- read.shapefile(savename) #reads back into env in correct format
setwd("../")
return(cropshp)}
