# This code is licensed under the MIT License

# This is a code to replicate Figure 1 of the paper.

# Input: Rainfall rate data (downloaded from https://trmm.gsfc.nasa.gov/data_dir/data.html) and 
# a map of Eastern part of Indian subcontinent

# Outout: satellite pass of Figure 1 of the paper 
# execution time: few seconds (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)

# set directory 
 # setwd("/home/")

# Input: rainfall rate for the month June of year 2012
input_file="data/satellite_passes_2012/6.csv"

# read input data 
xy=read.csv(file=input_file)   
# read the map
map <- read.table(file = "data/map.txt", header = TRUE)
map[map$long==-999,1] = NA
map[map$lat==-999,2] = NA

# Time of satellite pass  
xy11=xy$second+60*(xy$minute+60*(xy$hour+24*xy$day))
len_1=length(xy11)
xy12=c(xy11[-1]-xy11[-len_1],0)
xy=cbind(xy,xy12)
xy$lon = as.numeric(as.character(xy$lon))
xy$lat = as.numeric(as.character(xy$lat))
rf_rate_1=list()

# Detection of dictinct pass of satellite  
k11=c(0,which(xy[,9]>100))
month_data=list()
month_data1=list()

# plot configuration
m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights= c(2,2,0.5))
par(mar = c(0.1, 0.7, 0.8, 0.1), oma = c(4,4,0,0), xpd = NA)           

# select one pass, for example 59th pass of the above month
# top left panel of Figure 1
m=59
k1=k11[m]+1
k2=k11[m+1]
xy1=xy[(k1):(k2),1:9]

rf_non=which(xy1[,3]>0)   # Which points have positive rainfall status  
if (length(rf_non>0)) {
  dat=xy1[,1:2]
  ch <- chull(dat)
  coords <- dat[c(ch, ch[1]), ]  # closed polygon
  
  plot(map$long, map$lat, type="l", xlim=c(81,93), ylim=c(20,30),xaxs="i", yaxs="i", xlab="", ylab="latitude (degrees)",xaxt='n')
  text(82,28.5,"Nepal")
  text(91,27.5,"Bhutan")
  text(91.5,29,"China")
  text(91,24.5,"Bangladesh")
  text(83,22,"India")
  lines(c(81,93,93,81,81),c(20,20,30,30,20))
  lines(c(84,90,90,84,84),c(21,21,30,30,21),lty=2,col=2)
  
  points(dat[,1],dat[,2],col="white") 
  lines(coords, col="black")
  points(xy1[rf_non,1],xy1[rf_non,2],cex=0.3,pch=0.01,col="blue")    # Plot of positive rainfal point in green colour
  legend_order <- matrix(1:2,ncol=1,byrow = TRUE)
  legend(84.5,30,c("Time: 20:24:33", "Date: 3-6-2012")[legend_order],pch=c(0),pt.cex=c(0),cex=1,
         col=c(1),bty = "n")
}

# top right panel of Figure 1 (second consecutive pass of the satellite)
m=m+2
k1=k11[m]+1
k2=k11[m+1]
xy1=xy[(k1):(k2),1:9]
rf_non=which(xy1[,3]>0)   
if (length(rf_non>0)) {
  dat=xy1[,1:2]
  ch <- chull(dat)
  coords <- dat[c(ch, ch[1]), ]  
  
  # Get the plot information 
  plot(map$long, map$lat, type="l", xaxs="i", yaxs="i",xlim=c(81,93), ylim=c(20,30), xlab="", ylab="",xaxt='n',yaxt='n')
  text(82,28.5,"Nepal")
  text(91,27.5,"Bhutan")
  text(91.5,29,"China")
  text(91,24.5,"Bangladesh")
  text(83,22,"India")
  lines(c(81,93,93,81,81),c(20,20,30,30,20))
  lines(c(84,90,90,84,84),c(21,21,30,30,21),lty=2,col=2)
  
  points(dat[,1],dat[,2],col="white") 
  lines(coords, col="black")
  points(xy1[rf_non,1],xy1[rf_non,2],cex=0.3,pch=0.01,col="blue")    # Plot of positive rainfal point in green colour
  legend_order <- matrix(1:2,ncol=1,byrow = TRUE)
  legend(84.5,30,c("Time: 02:56:53","Date: 4-6-2012")[legend_order],pch=c(0),pt.cex=c(0),cex=1,
         col=c(1),bty = "n")
}

# bottom left panel (third consecutive pass)
m=m+1
k1=k11[m]+1
k2=k11[m+1]
xy1=xy[(k1):(k2),1:9]
rf_non=which(xy1[,3]>0)   
if (length(rf_non>0)) {
  dat=xy1[,1:2]
  ch <- chull(dat)
  coords <- dat[c(ch, ch[1]), ]  
  plot(map$long, map$lat, type="l", xaxs="i", yaxs="i",xlim=c(81,93), ylim=c(20,30), xlab="longitude (degrees)", ylab="latitude(degrees)")
  text(82,28.5,"Nepal")
  text(91,27.5,"Bhutan")
  text(91.5,29,"China")
  text(90,22.5,"Bangladesh")
  text(83,22,"India")
  lines(c(81,93,93,81,81),c(20,20,30,30,20))
  lines(c(84,90,90,84,84),c(21,21,30,30,21),lty=2,col=2)
  
  points(dat[,1],dat[,2],col="white") 
  lines(coords, col="black")
  points(xy1[rf_non,1],xy1[rf_non,2],cex=0.3,pch=0.01,col="blue")    # Plot of positive rainfal point in green colour
  legend_order <- matrix(1:2,ncol=1,byrow = TRUE)
  legend(84.5,30,c("Time: 19:28:41", "Date: 4-6-2012")[legend_order],pch=c(0),pt.cex=c(0),
         cex=1,col=c(1),bty = "n")
}

# bottom right panel (fourth consecutive pass)
m=m+1
k1=k11[m]+1
k2=k11[m+1]
xy1=xy[(k1):(k2),1:9]
rf_non=which(xy1[,3]>0)   
if (length(rf_non>0)) {
  dat=xy1[,1:2]
  ch <- chull(dat)
  coords <- dat[c(ch, ch[1]), ]  
  plot(map$long, map$lat, type="l", xaxs="i", yaxs="i",xlim=c(81,93), ylim=c(20,30), xlab="longitude (degrees)", ylab="", yaxt='n')
  text(82,28.5,"Nepal")
  text(91,27.5,"Bhutan")
  text(91.5,29,"China")
  text(91,24.8,"Bangladesh")
  text(83,22,"India")
  lines(c(81,93,93,81,81),c(20,20,30,30,20))
  lines(c(84,90,90,84,84),c(21,21,30,30,21),lty=2,col=2)
  
  points(dat[,1],dat[,2],col="white") 
  lines(coords, col="black")
  points(xy1[rf_non,1],xy1[rf_non,2],cex=0.3,pch=0.01,col="blue")    # Plot of positive rainfal point in green colour
  legend_order <- matrix(1:2,ncol=1,byrow = TRUE)
  legend(84.5,30,c("Time: 02:01:04", "Date: 5-6-2012")[legend_order],pch=c(0),pt.cex=c(0),cex=1,col=1,bty = "n")
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 1,legend =c("Rainfall","Area covered in one pass of satellite","Target area"), 
       text.width=c(0.2,0.3,.5),pch=c(15,22,0),horiz = TRUE,cex=1,col=c(4,1,2),pt.cex=1,bty = "n")

