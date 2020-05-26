# This code is licensed under the MIT License

# Input: Boundary contours of few precipitation areas during a monsoon month
# Output: Figure 2 of the paper

# load libraries
library(sp)
library(geosphere)
library(plyr)

# set directory
# setwd("/home")

# few examples of rainfall region during June 2012. 
input_file1="data/rf_data/2012/M6_cluster_points"
input_file2="data/rf_data/2012/M6_cluster_lengths"
input_file3="data/rf_data/2012/M6_cen_spell"

# plot configuration 
mat<- matrix(c(1,2,3),nrow = 1,ncol = 3,byrow = TRUE)
layout(mat = mat,heights=c(9,2))
par(mar = c(2,2,2,2), oma = c(4,3,0,0), xpd = NA)            

# function to compute starhull 
func1=function(xx,yy,zz){
  p1=as.data.frame(xx)  # cluster of points
  p2=as.data.frame(yy)  # length of the cluster
  p11=as.data.frame(zz) # rainfall status of the points 
  p3=p1[-which(p1[,1]==0),]            
  p5=p2[-which(p2[,1]==0),]
  p111=p11[-which(p11[,1]==0),]
  p4=cumsum(p5)                      
  N=length(p4)
  p4=c(0,p4)
  contour_1=list()
  contour_prj=list()
  area1=array(0,dim=N)
  
  q=list()
  ll=NULL
  for (i in 1:N)
  {
    obj_1=NULL
    obj_1=p3[(p4[i]+1):p4[i+1],]      # Shape objects_1
    if (nrow(obj_1)>6){    
      m00<-data.frame(x=obj_1[,1], y=obj_1[,2])
      area1[i]=areaPolygon(m00)/(10^6)   # Compute the area of a polygon incoordinates (longitude/latitude)
      m1=rbind(m00,m00[1,])
      m11=c(min(m1[,1]),min(m1[,2]))
      m12=c(max(m1[,1]),min(m1[,2]))
      m21=c(min(m1[,1]),max(m1[,2]))
      m22=c(max(m1[,1]),max(m1[,2]))
      bounding_box=rbind(m11,m12,m22,m21,m11)  # Vertices of the bounding box
      cm1=(min(m1[,1])+max(m1[,1]))/2    # Long-Coordinate of the center of the bounding box
      cm2=(min(m1[,2])+max(m1[,2]))/2    # Lat-Coordinate of the center of the bounding box
      box_center=c(cm1,cm2)              # Center of the bounding box
      R=6371
      x=NULL
      y=NULL
      x=R*((cos(m1[,2]*2*pi/360)+cos(cm2*2*pi/360))/2)*(m1[,1]-cm1)*(2*pi/360)
      y=R*(m1[,2]-cm2)*(2*pi/360)
      contour_prj[[i]]=cbind(x,y)
    }
  }
  n1=intersect(intersect(which(area1>=200),which(area1<=13500)),which(p111==2))
  m=length(n1)
  
# Convex hull construction 
  s1=1000  # Number of equispaced points from the contour
  radius_dist=array(0,dim=c(s1,m))
  theta1=seq(0,2*pi,length.out=(s1+1))[-(s1+1)]
  f1=function(x,y){ifelse(x < 0, atan(y / x) + pi, ifelse(y < 0 , atan(y / x) + 2*pi, atan(y / x)))}
  dd=NULL
  for(j in 1: 3){  
    X=contour_prj[[n1[j]]] # Main data point
    hpts <- chull(X)   # Convex hull construction
    hpts <- c(hpts, hpts[1])
    xy1=X[hpts, ]       # Convex hull points
    l <- SpatialLines(list(Lines(Line(rbind(xy1, xy1[1, ])), ID=1))) # Now create a SpatialLines object from the vertices.
    pts <- spsample(l, 500, type="regular")
    pts_1=coordinates(pts) # Equispaced point from the convex hull
    mm11=colMeans(pts_1) # Center of the convex hull (same as above)
    # Below is the arrangements of the data points according to the centered  convex hull
    XY=cbind(X[,1]-mm11[1],X[,2]-mm11[2])   # This is the polygon whose convex hull is centered
    plot(XY,type="l",ann=FALSE,xlim=range(XY[,1]),ylim=range(XY[,2]),col="blue") # Correct plot with shift(for centering)
    points(0,0,col=4, pch=1,cex=.5)
    
# Star-hull construction 
    theta22=NULL
    r_theta22=NULL
    for (i1 in 1:(nrow(XY)-1)){
      xy1=XY[i1,]          # one vertex of a line segment of a polygon
      xy2=XY[i1+1,]        # other vertex
      t1=f1(xy1[1],xy1[2])         # theta of polar representation of the first vertex
      r1=sqrt(xy1[1]^2+xy1[2]^2)   # length
      t2=f1(xy2[1],xy2[2])
      r2=sqrt(xy2[1]^2+xy2[2]^2)
      theta11=NULL
      r_theta=NULL
      if  (abs(t2-t1)<pi) {theta11=theta1[intersect(which(theta1>=min(t1,t2)),which(theta1<=max(t1,t2)))]} else
      {theta11=theta1[c(which(theta1>=max(t1,t2)),which(theta1<=min(t1,t2)))]}
      r_theta=(r1*r2*sin(t1-t2))/(r1*sin(t1-theta11)+r2*sin(theta11-t2))
      theta22=c(theta22,theta11)
      r_theta22=c(r_theta22,r_theta)
      }
# To detect whethere the centoid is inside or outside the contour
    l1=array(0,dim=s1)
    for (i2 in 1:s1){
      l1[i2]=length(r_theta22[which(theta22==theta1[i2])])
    }
    if (min(l1)>=1){
      for (i2 in 1:s1){
        r_list=NULL
        r_list=r_theta22[which(theta22==theta1[i2])]
        radius_dist[i2,j]=max(r_list)    # Star-hull contour
      }
    }else{  # When centroid outside the contour
      for (i2 in 1:s1){
        r_list=NULL
        r_list=r_theta22[which(theta22==theta1[i2])]
        if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
        radius_dist[i2,j]=max(r_list)    # Star-hull contour
      }
    }
    theta=seq(0,2*pi,length=s1+1)[-(s1+1)]
    x11=radius_dist[,j]*cos(theta)
    y11=radius_dist[,j]*sin(theta)
    lines(x11,y11, col="red",lty=2)
  }
  return(radius_dist)
}

# read input
points1=read.table(input_file1)
lengths1=read.table(input_file2)
spell1=read.table(input_file3)

# generate results
func1(points1,lengths1,spell1)

title(xlab = "West-East (in km)",
      ylab = "South-North (in km)",
      outer = TRUE, line = 1,cex.lab=1)

