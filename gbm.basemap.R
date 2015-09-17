gbm.basemap <- function(fileloc,
                        xlim=range(grids[,1]), #
                        ylim=range(grids[,2]),
                        savename){
install.packages("shapefiles") #ask before doing this/ do if shapefiles not installed
library(shapefiles)
#world <- read.shapefile("C:/Users/Simon/Desktop/gshhg-shp-2.3.4/GSHHS_shp/h/GSHHS_h_L1")
world <- read.shapefile(fileloc)
xlim=xlim # do I need this if specified at top?
ylim=ylim # do I need this if specified at top?
basemap(xlim, ylim)
draw.shape(world, col="darkgreen") #should be good so far. But then need to cut the shape & save somehow?
# mybasemap <- cut.shape(world,xlim=xlim,ylim=ylim) #if this works, do I need basemap & draw.shape? Not seen...
write.shapefile(mybasemap,savename,T)
}

gbm.basemap(fileloc="C:/Users/Simon/Desktop/gshhg-shp-2.3.4/GSHHS_shp/h/GSHHS_h_L1",savename="OIreland")

library(maptools)
#data(world) # ensure world is an object
library(rgeos)
# create  a polygon that defines the boundary
bnds <- cbind(x=c(min(grids[,1]), min(grids[,1]), max(grids[,1]), max(grids[,1]), min(grids[,1])), y=c(min(grids[,2]), max(grids[,2]), max(grids[,2]), min(grids[,2]), min(grids[,2])))
# convert to a spatial polygons object with the same CRS
SP <- SpatialPolygons(list(Polygons(list(Polygon(bnds)), "1")), proj4string=CRS(proj4string(world)))
# find the intersection with the original SPDF
gI <- gIntersects(world, SP, byid=TRUE)
# create the new spatial polygons object.
out <- vector(mode="list", length=length(which(gI)))
ii <- 1
for (i in seq(along=gI)) if (gI[i]) {
  out[[ii]] <- gIntersection(world[i,], SP)
  row.names(out[[ii]]) <- row.names(world)[i]; ii <- ii+1
}
# use rbind.SpatialPolygons method to combine into a new object.
out1 <- do.call("rbind", out)
# Plot it. Probably don't need to do: just save as shapefile.
plot(out1, col = "khaki", bg = "azure2")