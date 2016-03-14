gbm.basemap <- function(getzip = TRUE, # download & unpack GSHHS data to WD?
                        bounds, # region to crop to: c(xmin,xmax,ymin,ymax)
                        res = 4, # resolution, high to low, 1:5 or c,l,i,h,f for
                        #coarse, low, intermediate, high, full. default 4
                        zipvers = "2.3.4", # GSHHS version, in case it updates
                        savename = "Crop_Map"){ #shapefile save name, no extension

# Downloads unzips crops & saves NOAAs global coastline shapefiles to user-set
# box. Use for 'shape' in gbm.map. If already downloaded & unzipped, setwd to
# "GSHHS_shp" folder. If downloading in RStudio uncheck "Use secure download
# method for HTTP" in Tools > Global Options > # Packages
# Simon Dedman, 2015 & 16. simondedman@gmail.com simondedman.com

print(paste("if rgdal install fails in linux try: sudo apt-get install libgdal-dev && sudo apt-get install libproj-dev"))
if (!require(rgdal)) install.packages("rgdal")
library(rgdal) # for readOGR
if (!require(raster)) install.packages("raster")
library(raster) # for crop
if (!require(maptools)) install.packages("maptools")
library(maptools) # for crop
  ###improve these: check if installed, install if not else library####

if (getzip) {  # download & unzip GSHGG if asked
  download.file(paste("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-", zipvers, ".zip", sep = ""), "GSHHG.zip")
  unzip("GSHHG.zip")
  setwd("./GSHHS_shp")}

if (res == 1) res <- "c" # change res number to folder letter if number used
if (res == 2) res <- "l"
if (res == 3) res <- "i"
if (res == 4) res <- "h"
if (res == 5) res <- "f"
setwd(paste("./", res, sep = "")) #setwd to that folder

# read in worldmap
world <- readOGR(dsn = "./", layer = paste("GSHHS_", res, "_L1", sep = ""))
cropshp <<- crop(world, bounds) # crop to extents
dir.create("CroppedMap") # create conservation maps directory
setwd("./CroppedMap")
writeSpatialShape(cropshp, savename)
print(paste("World map cropped and saved successfully"))}