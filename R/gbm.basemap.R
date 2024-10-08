#' Creates Basemaps for Gbm.auto mapping from your data range
#'
#' Downloads unzips crops & saves NOAAs global coastline shapefiles to user-set
#' box. Use for 'shape' in gbm.map. If downloading in RStudio uncheck
#' "Use secure download method for HTTP" in Tools > Global Options > Packages.
#' Simon Dedman, 2015/6 simondedman@gmail.com GitHub.com/SimonDedman/gbm.auto
#'
#' @param bounds Region to crop to: c(xmin,xmax,ymin,ymax).
#' @param grids If bounds unspecified, name your grids database here.
#' @param gridslat If bounds unspecified, specify which column in grids is
#' latitude.
#' @param gridslon If bounds unspecified, specify which column in grids is
#' longitude.
#' @param getzip Download & unpack GSHHS data to WD? "TRUE" else
#' absolute/relative reference to GSHHS_shp folder, including that folder.
#' @param zipvers GSHHS version, in case it updates. Please email developer (SD)
#'  if this is incorrect.
#' @param savedir Save outputs to a temporary directory (default) else change to
#'  current directory e.g. "/home/me/folder". Do not use getwd() here.
#' @param savename Shapefile save-name, no shp extension, default is "Crop_Map"
#' @param res Resolution, 1:5 (low:high) OR c,l,i,h,f (coarse, low,
#' intermediate, high, full) or "CALC" to calculate based on bounds. Choose one.
#' @param extrabounds Grow bounds 16pct each direction to expand rectangular
#' datasets basemaps over the entire square area created by basemap in mapplots.
# #' @param returnsf Return object as simple features object? Default FALSE,
# #'  returns as list format for draw.shape in mapplots, used in gbm.map.
#'
#' @return basemap coastline file for gbm.map in gbm.auto. "cropshp"
#' SpatialPolygonsDataFrame in in local environment & user-named files in
#' "CroppedMap" folder. Load later with maptools function:
#' MyMap <- sf::st_read(dsn = "./CroppedMap/Crop_Map.shp", layer = "Crop_Map, quiet = TRUE)
#'
#' @export
#' @importFrom graphics lines par
#' @importFrom utils download.file unzip
#' @importFrom sf st_crop st_read st_write sf_use_s2
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples
#' \donttest{
#' # Not run: downloads and saves external data.
#' data(samples)
#' mybounds <- c(range(samples[,3]),range(samples[,2]))
#' gbm.basemap(bounds = mybounds, getzip = "./GSHHS_shp/",
#' savename = "My_Crop_Map", res = "f")
#' # In this example GSHHS folder already downloaded to the working directory
#' # hence I pointed getzip at that rather than having it download the zip again
#' }
#'
#' @details errors and their origins:
#' 1. Error in setwd(getzip) : cannot change working directory
#' If you've specified the location of the local GSHHS_shp folder, ensure you're
#' in the correct directory relative to it. This error means it looked for the
#' folder and couldn't find it.
#'
#' 2. subscript out of bounds: can't crop world map to your bounds.
#' Check lat/lon are the right way around: check gridslat and gridslon point to the correct columns
#' for lat and lon in grids, and those columns named (something like) lat and lon, ARE ACTUALLY the
#' latitudes and longitudes, and not the wrong way around.
#'
#' 3. If your download is timing out use options(timeout = 240).
#'
#' 4. Error in if (scope >= 160) res <- "c" : missing value where TRUE/FALSE needed. Check gridslat
#'  and gridslon are indexing the correct columns in grids.
#'
gbm.basemap <- function(
    bounds = NULL, # region to crop to: c(xmin,xmax,ymin,ymax)
    grids = NULL, # if bounds unspecified, name your grids database here
    gridslat = NULL, # if bounds unspecified, specify which column in grids is latitude
    gridslon = NULL, # if bounds unspecified, specify which column in grids is longitude
    getzip = TRUE, # download & unpack GSHHS data to WD? "TRUE" else absolute/relative reference to GSHHS_shp folder, including that folder
    zipvers = "2.3.7", # GSHHS version, in case it updates. Please email developer if this is incorrect
    savedir = tempdir(), # save outputs to a temporary directory (default) else
    # change to current directory e.g. "/home/me/folder". Do not use getwd() here.
    savename = "Crop_Map", #shapefile save-name without the .shp
    res = "CALC", # Resolution, 1:5 (low:high) OR c,l,i,h,f (coarse, low, intermediate, high, full) or "CALC" to calculate based on bounds. Choose one.
    extrabounds = FALSE # grow bounds 16pct each direction to expand rectangular datasets basemaps over the entire square area created by basemap in mapplots
    # returnsf = FALSE # obviated by gbm.mapsf
) { # Return object as simple features object? Default FALSE, returns as list format for draw.shape in mapplots, used in gbm.map

  oldwd <- getwd() # record original directory
  on.exit(setwd(oldwd), add = TRUE) # defensive block, thanks to Gregor Sayer
  sf::sf_use_s2(FALSE) # 2021 addition of s2 code to sf often causes: Error in s2_geography_from_wkb(x, oriented = oriented, check = check):
  # Evaluation error: Found 1 feature with invalid spherical geometry. Loop 0 is not valid: Edge n has duplicate vertex with edge n2.
  # if bounds is entered it's user below, else check grids & gridslat & gridslon
  if (is.null(bounds)) {
    #check none of grids & gridslat & gridslon is null, if any are print message
    if (is.null(grids)) stop("if bounds is NULL grids needs to be specified")
    if (is.null(gridslat)) stop("if bounds is NULL gridslat needs to be specified")
    if (is.null(gridslon)) stop("if bounds is NULL gridslon needs to be specified")
    #check they're all the correct format
    if (!is.data.frame(grids)) stop("grids needs to be a data frame")
    if (!is.numeric(gridslat)) stop("gridslat needs to be a number")
    if (!is.numeric(gridslon)) stop("gridslon needs to be a number")
    # construct bounds from gridslat & gridslon ranges in grids
    bounds <- c(range(grids[,gridslon]), range(grids[,gridslat])) #still required later despite sf/st update
    xmin = min(grids[,gridslon]) #for sf/st upgrade
    xmax = max(grids[,gridslon])
    ymin = min(grids[,gridslat])
    ymax = max(grids[,gridslat])
  } else { # if bounds not null
    xmin = bounds[1] #for sf/st upgrade
    xmax = bounds[2]
    ymin = bounds[3]
    ymax = bounds[4]
  }

  if (extrabounds) { # grow bounds extents if requested
    xmid <- mean(bounds[1:2])
    ymid <- mean(bounds[3:4])
    xmax <- ((bounds[2] - xmid) * 1.6) + xmid #updated for sf/st
    xmin <- xmid - ((xmid - bounds[1]) * 1.6)
    ymax <- ((bounds[4] - ymid) * 1.6) + ymid
    ymin <- ymid - ((ymid - bounds[3]) * 1.6)
  }

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

  # If savedir has a terminal slash, remove it, it's added later
  if (substr(x = savedir, start = nchar(savedir), stop = nchar(savedir)) == "/") savedir = substr(x = savedir, start = 1, stop = nchar(savedir) - 1)

  ifelse(getzip == TRUE, {# download & unzip GSHGG if getzip = TRUE
    download.file(url = paste0("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-", zipvers, ".zip"),
                  destfile = paste0(savedir, "/GSHHG.zip")) # need to ensure 1 / between savedir & GSHHG.zip. Done above
    setwd(savedir)
    unzip("GSHHG.zip")
    setwd("GSHHS_shp")}
    , setwd(file.path(getzip, "GSHHS_shp")) # else just setwd to there, including GSHHS_shp
  )

  setwd(res) #setwd to res subfolder
  world <- sf::st_read(dsn = paste0("GSHHS_", res, "_L1.shp"), layer = paste0("GSHHS_", res, "_L1"), quiet = TRUE) # read in worldmap
  cropshp <- sf::st_crop(world, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) # crop to extents

  setwd(savedir) # setwd to savedir else saves CroppedMap folder in res folder
  dir.create("CroppedMap") # create conservation maps directory
  setwd("CroppedMap")
  sf::st_write(cropshp, dsn = paste0(savename, ".shp"), append = FALSE) # append FALSE overwrites existing files
  # cropshp <- st_read(dsn = paste0(savename, ".shp"), layer = savename, quiet = TRUE) # read in worldmap
  print(paste("World map cropped and saved successfully"))
  return(cropshp)}
