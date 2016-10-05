####TODO####
# Have users enter their lat & lon columns instead of bounds = c(xmin,xmax,ymin,ymax)?
# More in-keeping with standard practice for this package

# Can I unpack only the relevant folder? Or at least only GSHGG?
setwd("/home/simon/Desktop/")
ifelse(getzip == TRUE, { # download & unzip GSHGG if getzip = TRUE
  download.file(paste("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-", zipvers, ".zip", sep = ""), "GSHHG.zip")
  unzip("GSHHG.zip")
  unzip("GSHHG.zip", files = "GSHHS_shp") #fails
  unzip("GSHHG.zip", files = "LICENSE.TXT") #works
  setwd("./GSHHS_shp")}
  , setwd(getzip)) # else just setwd to there
####end todo####
