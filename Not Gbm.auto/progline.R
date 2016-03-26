# Progress bar for presentations
# Simon Dedman, simondedman@gmail.com, 13/9/2015

progline<-function(x,filetype="png",bg="transparent"){
  # x: number of slides
  # filetype: either "png", "bmp", "jpeg" or "tiff". Default: "png"
  # bg: defaults to transparent for png else white. Can set e.g. "blue"
  # function outputs a pic per slide, for x slides
  
  # option to have user-inputtable black dot, e.g. an image?
  # Make the default width = powerpoint slide width default. widescreen: 33.867 cm. Get all option sin design > slide size
  
  if(filetype=="png") cairotype <- "cairo-png" else cairotype <- "cairo"
  if(filetype!="png") bg="white"
  
for(i in 1:x){
  get(filetype)(filename = paste("Progress_",i,".png",sep=""),
      width = 90*x, height = 70, units = "px", pointsize = 48, bg = bg, res = NA, family = "", type = cairotype)
  
  par(mar=c(0,0,0,0), xpd=FALSE,cex=3)
  
  plot(1:x,rep(1,x),type="l",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
  
  points(1:x,rep(1,x),type="p",xlab="",ylab="",xaxt="n",yaxt="n",bty="n") #,lty="94", ljoin=0,lwd=2) # plots white dots & line

  points(i,1,pch=16, type="p",xlab="",ylab="",xaxt="n",yaxt="n",bty="n") #adds a black dot

dev.off()}}