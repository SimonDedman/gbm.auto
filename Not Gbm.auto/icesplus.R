## Script to convert ICES Rectangle codes to lat-longs
icesplus <- function(x, # data object 
                     datacol = 1, # Set for ICES rect codes
                     filename = "ICES Rect to LatLonWKT.csv"){ 
  require(mapplots)
  ncol <- length(x) # get number of columns
  col1 <- ncol + 1 # name & add 2 columns for lat/longs
  col2 <- ncol + 2
  x[,col1:col2] <- ices.rect(x[,datacol]) # lat & lon columns created
  # export the x, creating WKT polygon column for rectangle coordinates based on lat/lon midpoints
  write.csv(x = data.frame(x,
                           WKT = paste("POLYGON ((",
                                 x[,col1] - 0.5, " ", x[,col2] + 0.25, ", ",
                                 x[,col1] + 0.5, " ", x[,col2] + 0.25, ", ",
                                 x[,col1] + 0.5, " ", x[,col2] - 0.25, ", ",
                                 x[,col1] - 0.5, " ", x[,col2] - 0.25, ", ",
                                 x[,col1] - 0.5, " ", x[,col2] + 0.25, "))",
                                 sep = "")),
                           file = filename,
                           row.names = FALSE)
  return(x)
}