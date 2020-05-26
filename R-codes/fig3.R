# This code is licensed under the MIT License

# Input: Boundary contours of complete regions under rainfall over the years 2002-2012
# Output: Figure 3 of the paper

# Execution time: 5 minutes 37 seconds (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)


# load libraries
library(sp)
library(geosphere)
library(plyr)

# set directory
# setwd("/home")

# function to compute starhull approximation (relative error)
func2 <- function(xx,yy,zz){
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
   if (nrow(obj_1)>=4){     
    m00<-data.frame(x=obj_1[,1], y=obj_1[,2])
    area1[i]=areaPolygon(m00)/(10^6)   # The area of a polygon incoordinates (longitude/latitude)
    m1=rbind(m00,m00[1,])
    m11=c(min(m1[,1]),min(m1[,2]))
    m12=c(max(m1[,1]),min(m1[,2]))
    m21=c(min(m1[,1]),max(m1[,2]))
    m22=c(max(m1[,1]),max(m1[,2]))
    bounding_box=rbind(m11,m12,m22,m21,m11)  # Vertices of the bounding box
    cm1=(min(m1[,1])+max(m1[,1]))/2    # Long-Coordinate of the center of the bounding box
    cm2=(min(m1[,2])+max(m1[,2]))/2    # Lat-Coordinate of the center of the bounding box
    box_center=c(cm1,cm2)              # Center of the bounding box
    R=6371                             # Radius of Earth in kilometer
    x=NULL
    y=NULL
    x=R*((cos(m1[,2]*2*pi/360)+cos(cm2*2*pi/360))/2)*(m1[,1]-cm1)*(2*pi/360)
    y=R*(m1[,2]-cm2)*(2*pi/360)
    contour_prj[[i]]=cbind(x,y)
  }
  }
  n1=intersect(intersect(which(area1>200),which(area1<=13500)),which(p111==2))
  m=length(n1)

# Convex hull construction
  s1=1000                           # Number of equispaced points from the contour
  radius_dist=array(0,dim=c(s1,m))
  radius_dist0=array(0,dim=c(s1,m))
  ratio_d=array(0,dim=c(s1,m))
  c1=array(0,dim=m)
  theta1=seq(0,2*pi,length.out=(s1+1))[-(s1+1)]
  f1 <- function(x,y){ifelse(x < 0, atan(y / x) + pi, ifelse(y < 0 , atan(y / x) + 2*pi, atan(y / x)))}
  for(j in 1: m){
    X=contour_prj[[n1[j]]]         # Main data point
    hpts <- chull(X)               # Convex hull construction
    hpts <- c(hpts, hpts[1])
    xy1=X[hpts, ]                   # Convex hull points
    l <- SpatialLines(list(Lines(Line(rbind(xy1, xy1[1, ])), ID=1))) # Now create a SpatialLines object from the vertices.
    pts <- spsample(l, 500, type="regular")
    pts_1=coordinates(pts)           # Equispaced point from the convex hull
    mm11=colMeans(pts_1)             # Center of the convex hull
    XY=cbind(X[,1]-mm11[1],X[,2]-mm11[2])   # This is the polygon whose convex hull is centered

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
        if (length(r_list)==1) {radius_dist0[i2,j]=r_list}else{
          r_list=sort(r_list)
          d11=array(0,dim=length(r_list))
          d11[1]=r_list[1]
          for (k0 in 2:length(r_list)) {d11[k0]=r_list[k0]-r_list[k0-1]}
          radius_dist0[i2,j]=sum(d11[which(1:length(r_list)%%2!=0)])
        }
      }
    }else{  # When centroid outside the contour
      for (i2 in 1:s1){
        r_list=NULL
        r_list=r_theta22[which(theta22==theta1[i2])]
          if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
           radius_dist[i2,j]=max(r_list)    # Star-hull contour
        if (length(r_list)==1) {radius_dist0[i2,j]=r_list}else{
          r_list=sort(r_list)
          d11=array(0,dim=length(r_list))
          for (k0 in 2:length(r_list)) {d11[k0]=r_list[k0]-r_list[k0-1]}
          radius_dist0[i2,j]=sum(d11[which(1:length(r_list)%%2==0)])
          }
      }
    }
    ratio_d[,j]=radius_dist0[,j]/radius_dist[,j]
  }
  return(ratio_d)
}

# Computation of starhull approximation error 
# for monsoon months of the years 2002-2012 using the above function

h2=h3=h4=h5=h6=h7=h8=h9=h10=h11=h12=NULL
# Year 2002
h_6=func2(read.table("data/rf_data/2002/M67_cluster_points"),
               read.table("data/rf_data/2002/M67_cluster_lengths"),
               read.table("data/rf_data/2002/M67_cen_spell"))
h_7=func2(read.table("data/rf_data/2002/M89_cluster_points"),
               read.table("data/rf_data/2002/M89_cluster_lengths"),
               read.table("data/rf_data/2002/M89_cen_spell"))
h2=cbind(h_6,h_7)
# Year 2003
h_6=func2(read.table("data/rf_data/2003/M67_cluster_points"),
               read.table("data/rf_data/2003/M67_cluster_lengths"),
               read.table("data/rf_data/2003/M67_cen_spell"))
h_7=func2(read.table("data/rf_data/2003/M89_cluster_points"),
               read.table("data/rf_data/2003/M89_cluster_lengths"),
               read.table("data/rf_data/2003/M89_cen_spell"))
h3=cbind(h_6,h_7)
# Year 2004
h_6=func2(read.table("data/rf_data/2004/M6_cluster_points"),
               read.table("data/rf_data/2004/M6_cluster_lengths"),
               read.table("data/rf_data/2004/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2004/M7_cluster_points"),
               read.table("data/rf_data/2004/M7_cluster_lengths"),
               read.table("data/rf_data/2004/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2004/M8_cluster_points"),
               read.table("data/rf_data/2004/M8_cluster_lengths"),
               read.table("data/rf_data/2004/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2004/M9_cluster_points"),
               read.table("data/rf_data/2004/M9_cluster_lengths"),
               read.table("data/rf_data/2004/M9_cen_spell"))
h4=cbind(h_6,h_7,h_8,h_9)
# Year 2005
h_6=func2(read.table("data/rf_data/2005/M6_cluster_points"),
               read.table("data/rf_data/2005/M6_cluster_lengths"),
               read.table("data/rf_data/2005/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2005/M7_cluster_points"),
               read.table("data/rf_data/2005/M7_cluster_lengths"),
               read.table("data/rf_data/2005/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2005/M8_cluster_points"),
               read.table("data/rf_data/2005/M8_cluster_lengths"),
               read.table("data/rf_data/2005/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2005/M9_cluster_points"),
               read.table("data/rf_data/2005/M9_cluster_lengths"),
               read.table("data/rf_data/2005/M9_cen_spell"))
h5=cbind(h_6,h_7,h_8,h_9)
# Year 2006
h_6=func2(read.table("data/rf_data/2006/M6_cluster_points"),
               read.table("data/rf_data/2006/M6_cluster_lengths"),
               read.table("data/rf_data/2006/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2006/M7_cluster_points"),
               read.table("data/rf_data/2006/M7_cluster_lengths"),
               read.table("data/rf_data/2006/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2006/M8_cluster_points"),
               read.table("data/rf_data/2006/M8_cluster_lengths"),
               read.table("data/rf_data/2006/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2006/M9_cluster_points"),
               read.table("data/rf_data/2006/M9_cluster_lengths"),
               read.table("data/rf_data/2006/M9_cen_spell"))
h6=cbind(h_6,h_7,h_8,h_9)
# Year 2007
h_6=func2(read.table("data/rf_data/2007/M6_cluster_points"),
               read.table("data/rf_data/2007/M6_cluster_lengths"),
               read.table("data/rf_data/2007/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2007/M7_cluster_points"),
               read.table("data/rf_data/2007/M7_cluster_lengths"),
               read.table("data/rf_data/2007/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2007/M8_cluster_points"),
               read.table("data/rf_data/2007/M8_cluster_lengths"),
               read.table("data/rf_data/2007/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2007/M9_cluster_points"),
               read.table("data/rf_data/2007/M9_cluster_lengths"),
               read.table("data/rf_data/2007/M9_cen_spell"))
h7=cbind(h_6,h_7,h_8,h_9)
# Year 2008
h_6=func2(read.table("data/rf_data/2008/M6_cluster_points"),
               read.table("data/rf_data/2008/M6_cluster_lengths"),
               read.table("data/rf_data/2008/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2008/M7_cluster_points"),
               read.table("data/rf_data/2008/M7_cluster_lengths"),
               read.table("data/rf_data/2008/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2008/M8_cluster_points"),
               read.table("data/rf_data/2008/M8_cluster_lengths"),
               read.table("data/rf_data/2008/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2008/M9_cluster_points"),
               read.table("data/rf_data/2008/M9_cluster_lengths"),
               read.table("data/rf_data/2008/M9_cen_spell"))
h8=cbind(h_6,h_7,h_8,h_9)
# Year 2009
h_6=func2(read.table("data/rf_data/2009/M6_cluster_points"),
               read.table("data/rf_data/2009/M6_cluster_lengths"),
               read.table("data/rf_data/2009/M6_cen_spell"))

h_7=func2(read.table("data/rf_data/2009/M7_cluster_points"),
               read.table("data/rf_data/2009/M7_cluster_lengths"),
               read.table("data/rf_data/2009/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2009/M8_cluster_points"),
               read.table("data/rf_data/2009/M8_cluster_lengths"),
               read.table("data/rf_data/2009/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2009/M9_cluster_points"),
               read.table("data/rf_data/2009/M9_cluster_lengths"),
               read.table("data/rf_data/2009/M9_cen_spell"))
h9=cbind(h_6,h_7,h_8,h_9)
# Year 2010
h_6=func2(read.table("data/rf_data/2010/M6_cluster_points"),
               read.table("data/rf_data/2010/M6_cluster_lengths"),
               read.table("data/rf_data/2010/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2010/M7_cluster_points"),
               read.table("data/rf_data/2010/M7_cluster_lengths"),
               read.table("data/rf_data/2010/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2010/M8_cluster_points"),
               read.table("data/rf_data/2010/M8_cluster_lengths"),
               read.table("data/rf_data/2010/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2010/M9_cluster_points"),
               read.table("data/rf_data/2010/M9_cluster_lengths"),
               read.table("data/rf_data/2010/M9_cen_spell"))
h10=cbind(h_6,h_7,h_8,h_9)
# Year 2011
h_6=func2(read.table("data/rf_data/2011/M6_cluster_points"),
               read.table("data/rf_data/2011/M6_cluster_lengths"),
               read.table("data/rf_data/2011/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2011/M7_cluster_points"),
               read.table("data/rf_data/2011/M7_cluster_lengths"),
               read.table("data/rf_data/2011/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2011/M8_cluster_points"),
               read.table("data/rf_data/2011/M8_cluster_lengths"),
               read.table("data/rf_data/2011/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2011/M9_cluster_points"),
               read.table("data/rf_data/2011/M9_cluster_lengths"),
               read.table("data/rf_data/2011/M9_cen_spell"))
h11=cbind(h_6,h_7,h_8,h_9)
# Year 2012
h_6=func2(read.table("data/rf_data/2012/M6_cluster_points"),
               read.table("data/rf_data/2012/M6_cluster_lengths"),
               read.table("data/rf_data/2012/M6_cen_spell"))
h_7=func2(read.table("data/rf_data/2012/M7_cluster_points"),
               read.table("data/rf_data/2012/M7_cluster_lengths"),
               read.table("data/rf_data/2012/M7_cen_spell"))
h_8=func2(read.table("data/rf_data/2012/M8_cluster_points"),
               read.table("data/rf_data/2012/M8_cluster_lengths"),
               read.table("data/rf_data/2012/M8_cen_spell"))
h_9=func2(read.table("data/rf_data/2012/M9_cluster_points"),
               read.table("data/rf_data/2012/M9_cluster_lengths"),
               read.table("data/rf_data/2012/M9_cen_spell"))
h12=cbind(h_6,h_7,h_8,h_9)

h0=cbind(h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12) # combining errors from different years
h=replace(h0, is.na(h0), 1) # Repalce the NaN values: 0/0

h1=1-h   # defingin error

l2=c(0.85,0.9,0.95)  # Few percentiles of the error distribution
m=length(l2)
s11=nrow(h0) # Discretization length 

h10=array(0,dim=c(s11,m))
for (i in 1:s11){
  h10[i,]=quantile(h1[i,], probs=l2)*100
}
par(mfrow=c(1,1))
theta1=seq(0,2*pi,length.out = s11)
plot(theta1,h10[,1],col=1,lty=4,lwd=.1,type="l",xlim=range(theta1),
     xaxt="n",ylim=range(-3.5,h10),cex.lab=0.8,
     ylab="Percentage of approximation error",
     xlab="theta",cex.axis=.8)
axis(side=1, at=seq(0,2*pi,length.out=5),cex.axis=.9,
     labels=c(expression(paste(0,~(East))),expression(paste(pi/2, ~(North))),
              expression(paste(pi, ~(West))),expression(paste(3*pi/2, ~(South))),
              expression(paste(2*pi, ~(East)))))

lines(theta1,h10[,2],col=2,lty=3,lwd=1.5)
lines(theta1,h10[,3],col=3,lty=2,lwd=1.5)

legend(0.6,-0,c("85 Percentile","90 Percentile","95 Percentile"),col=c(1,2,3),cex=0.7,
       lty=c(4,3,2), lwd=c(1.5,1.5,1.5),box.lwd = 0,box.col = "white",bg = "white",horiz = TRUE)
lines(theta1,h10[,1],lty=4,lwd=1.5,type="l",cex.axis=.8)


