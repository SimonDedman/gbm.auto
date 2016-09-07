gbm.basemap <- function(bounds, # region to crop to: c(xmin,xmax,ymin,ymax)
                        getzip = TRUE, # download & unpack GSHHS data to WD?
  # "TRUE" else absolute/relative reference to GSHHS_shp folder, inc that folder
                        zipvers = "2.3.5-1", # GSHHS version, in case it updates
  # Please email developer if this is incorrect
                        savename = "Crop_Map", #shapefile savename, no extension
                        res = "CALC") { # resolution, 1:5 (low:high) OR c,l,i,h,f
# (coarse, low, intermediate, high, full) or "CALC" to calculate based on bounds

# Downloads unzips crops & saves NOAAs global coastline shapefiles to user-set
# box. Use for 'shape' in gbm.map. If downloading in RStudio uncheck
# "Use secure download method for HTTP" in Tools > Global Options > # Packages
# Simon Dedman, 2015/6 simondedman@gmail.com github.com/SimonDedman/gbm.auto

print(paste("if rgdal install fails in linux try: sudo apt-get install libgdal-dev && sudo apt-get install libproj-dev"))
if (!require(rgdal)) install.packages("rgdal")
library(rgdal) # for readOGR
if (!require(rgeos)) install.packages("rgeos")
library(rgeos) # subfunctions for rgdal & others
if (!require(raster)) install.packages("raster")
library(raster) # for crop
if (!require(maptools)) install.packages("maptools")
library(maptools) # for WriteSpatialShape
  ###improve these: check if installed, install if not else library####
startdir <- getwd() # record original directory

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
  download.file(paste("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-", zipvers, ".zip", sep = ""), "GSHHG.zip")
  unzip("GSHHG.zip")
  setwd("./GSHHS_shp")}
  , setwd(getzip)) # else just setwd to there

setwd(paste("./", res, sep = "")) #setwd to res subfolder

# read in worldmap
world <- readOGR(dsn = paste("GSHHS_", res, "_L1.shp", sep = ""), layer = paste("GSHHS_", res, "_L1", sep = ""))
cropshp <<- crop(world, bounds) # crop to extents
setwd(startdir)
dir.create("CroppedMap") # create conservation maps directory
setwd("CroppedMap")
writeSpatialShape(cropshp, savename)
print(paste("World map cropped and saved successfully"))
setwd("../")
return(cropshp)}