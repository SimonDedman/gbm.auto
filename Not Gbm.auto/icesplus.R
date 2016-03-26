## Script to convert ICES Rectangle codes to lat-longs
# set working directory
getwd()
setwd("C:/Users/Simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Scallop & Whelk MMO Sarah")

# Load function
icesplus<-function(file){
  require(mapplots)
  # lat & lon columns created
  file[,2:3]<-ices.rect(file[,1])
  # export the file, creating WKT polygon column for rectangle coordinates based on lat/lon midpoints
  write.csv(x=data.frame(file,WKT=paste("POLYGON ((",file[,2]-0.5," ",file[,3]+0.25,", ",file[,2]+0.5," ",file[,3]+0.25,", ",file[,2]+0.5," ",file[,3]-0.25,", ",file[,2]-0.5," ",file[,3]-0.25,", ",file[,2]-0.5," ",file[,3]+0.25,"))",sep="")),file=paste("ICES Rect to LatLonWKT.csv",sep=""),row.names=FALSE)
}

# Read the files in. Expected: csv files, 1 column of ICES rect codes.
# filename set so it can be run automatically
file<-read.csv(file="Whelk_Rects.csv")
file<-read.csv(file="Scallop_Rects.csv")
# Run function. File is auto-output - change name if using again
icesplus(file)

# remove everything
rm(list = ls())