# This code is licensed under the MIT License

# Input: Boundary contours of all precipitation areas with censoring status during a monsoon month
# Output: Star-hull of complete contours of regions under rainfall for a monsoon month
# Execution time: approx 6 minutes (i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)


# load libraries
library(sp)
library(geosphere)
library(plyr)

# function computes star-hull of a given contour
func_star <- function(xx,yy,zz){
 p1=as.data.frame(xx)
 p2=as.data.frame(yy)
 p11=as.data.frame(zz)
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
   obj_1=p3[(p4[i]+1):p4[i+1],]      
   if (nrow(obj_1)>6){    
   m00<-data.frame(x=obj_1[,1], y=obj_1[,2])
   area1[i]=areaPolygon(m00)/(10^6)          # Compute the area of a polygon incoordinates
   m1=rbind(m00,m00[1,])
   m11=c(min(m1[,1]),min(m1[,2]))
   m12=c(max(m1[,1]),min(m1[,2]))
   m21=c(min(m1[,1]),max(m1[,2]))
   m22=c(max(m1[,1]),max(m1[,2]))
   bounding_box=rbind(m11,m12,m22,m21,m11)    # Vertices of the bounding box
   cm1=(min(m1[,1])+max(m1[,1]))/2            # Long-Coordinate of the center of the bounding box
   cm2=(min(m1[,2])+max(m1[,2]))/2            # Lat-Coordinate of the center of the bounding box
  box_center=c(cm1,cm2)                       # Center of the bounding box
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
  s1=1000                            # Number of equispaced points from the contour
  radius_dist=array(0,dim=c(s1,m))
  theta1=seq(0,2*pi,length.out=(s1+1))[-(s1+1)]
  f1=function(x,y){ifelse(x < 0, atan(y / x) + pi, ifelse(y < 0 , atan(y / x) + 2*pi, atan(y / x)))}
  dd=NULL
  for(j in 1: m){
    X=contour_prj[[n1[j]]]                 # Main data point
    hpts <- chull(X)                       # Convex hull construction
    hpts <- c(hpts, hpts[1])
    xy1=X[hpts, ]                         # Convex hull points
    l <- SpatialLines(list(Lines(Line(rbind(xy1, xy1[1, ])), ID=1))) # Now create a SpatialLines object from the vertices.
    pts <- spsample(l, 500, type="regular")
    pts_1=coordinates(pts)                    # Equispaced point from the convex hull
    mm11=colMeans(pts_1)                      # Center of the convex hull (same as above)
    # Below is the arrangements of the data points according to the centered  convex hull
    XY=cbind(X[,1]-mm11[1],X[,2]-mm11[2])   # This is the polygon whose convex hull is centered
    # plot(XY,type="l",xlab="West-East",ylab="South-North",col="red") # Correct plot with shift(for centering)
    
# Star-hull construction 
    theta22=NULL
    r_theta22=NULL
    for (i1 in 1:(nrow(XY)-1)){
      xy1=XY[i1,]                     # one vertex of a line segment of a polygon
      xy2=XY[i1+1,]                   # other vertex
      
      t1=f1(xy1[1],xy1[2])           # theta of polar representation of the first vertex
      r1=sqrt(xy1[1]^2+xy1[2]^2)   
      t2=f1(xy2[1],xy2[2])
      r2=sqrt(xy2[1]^2+xy2[2]^2)
      theta11=NULL
      r_theta=NULL
      
      if  (abs(t2-t1)<pi) {theta11=theta1[intersect(which(theta1>=min(t1,t2)),which(theta1<=max(t1,t2)))]} else
      {theta11=theta1[c(which(theta1>=max(t1,t2)),which(theta1<=min(t1,t2)))]}
      
      r_theta=(r1*cos(t1)*(r1*sin(t1)-r2*sin(t2))-r1*sin(t1)*(r1*cos(t1)-r2*cos(t2)))/
        (cos(theta11)*(r1*sin(t1)-r2*sin(t2))-sin(theta11)*(r1*cos(t1)-r2*cos(t2)))
      
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
    }

  return(radius_dist)
}


h2=h3=h4=h5=h6=h7=h8=h9=h10=h11=h12=NULL

# Each block corresponds computation of starhull for a year
# year 2002
h_6=func_star(read.table("data/rf_data/2002/M67_cluster_points"),
               read.table("data/rf_data/2002/M67_cluster_lengths"),
               read.table("data/rf_data/2002/M67_cen_spell"))
h_7=func_star(read.table("data/rf_data/2002/M89_cluster_points"),
               read.table("data/rf_data/2002/M89_cluster_lengths"),
               read.table("data/rf_data/2002/M89_cen_spell"))
h2=cbind(h_6,h_7)

# year 2003
h_6=func_star(read.table("data/rf_data/2003/M67_cluster_points"),
               read.table("data/rf_data/2003/M67_cluster_lengths"),
               read.table("data/rf_data/2003/M67_cen_spell"))
h_7=func_star(read.table("data/rf_data/2003/M89_cluster_points"),
               read.table("data/rf_data/2003/M89_cluster_lengths"),
               read.table("data/rf_data/2003/M89_cen_spell"))
h3=cbind(h_6,h_7)

# year 2004
h_6=func_star(read.table("data/rf_data/2004/M6_cluster_points"),
               read.table("data/rf_data/2004/M6_cluster_lengths"),
               read.table("data/rf_data/2004/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2004/M7_cluster_points"),
               read.table("data/rf_data/2004/M7_cluster_lengths"),
               read.table("data/rf_data/2004/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2004/M8_cluster_points"),
               read.table("data/rf_data/2004/M8_cluster_lengths"),
               read.table("data/rf_data/2004/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2004/M9_cluster_points"),
               read.table("data/rf_data/2004/M9_cluster_lengths"),
               read.table("data/rf_data/2004/M9_cen_spell"))
h4=cbind(h_6,h_7,h_8,h_9)

# year 2005
h_6=func_star(read.table("data/rf_data/2005/M6_cluster_points"),
               read.table("data/rf_data/2005/M6_cluster_lengths"),
               read.table("data/rf_data/2005/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2005/M7_cluster_points"),
               read.table("data/rf_data/2005/M7_cluster_lengths"),
               read.table("data/rf_data/2005/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2005/M8_cluster_points"),
               read.table("data/rf_data/2005/M8_cluster_lengths"),
               read.table("data/rf_data/2005/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2005/M9_cluster_points"),
               read.table("data/rf_data/2005/M9_cluster_lengths"),
               read.table("data/rf_data/2005/M9_cen_spell"))
h5=cbind(h_6,h_7,h_8,h_9)

# year 2006
h_6=func_star(read.table("data/rf_data/2006/M6_cluster_points"),
               read.table("data/rf_data/2006/M6_cluster_lengths"),
               read.table("data/rf_data/2006/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2006/M7_cluster_points"),
               read.table("data/rf_data/2006/M7_cluster_lengths"),
               read.table("data/rf_data/2006/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2006/M8_cluster_points"),
               read.table("data/rf_data/2006/M8_cluster_lengths"),
               read.table("data/rf_data/2006/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2006/M9_cluster_points"),
               read.table("data/rf_data/2006/M9_cluster_lengths"),
               read.table("data/rf_data/2006/M9_cen_spell"))
h6=cbind(h_6,h_7,h_8,h_9)

# year 2007
h_6=func_star(read.table("data/rf_data/2007/M6_cluster_points"),
               read.table("data/rf_data/2007/M6_cluster_lengths"),
               read.table("data/rf_data/2007/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2007/M7_cluster_points"),
               read.table("data/rf_data/2007/M7_cluster_lengths"),
               read.table("data/rf_data/2007/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2007/M8_cluster_points"),
               read.table("data/rf_data/2007/M8_cluster_lengths"),
               read.table("data/rf_data/2007/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2007/M9_cluster_points"),
               read.table("data/rf_data/2007/M9_cluster_lengths"),
               read.table("data/rf_data/2007/M9_cen_spell"))
h7=cbind(h_6,h_7,h_8,h_9)

# year 2008
h_6=func_star(read.table("data/rf_data/2008/M6_cluster_points"),
               read.table("data/rf_data/2008/M6_cluster_lengths"),
               read.table("data/rf_data/2008/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2008/M7_cluster_points"),
               read.table("data/rf_data/2008/M7_cluster_lengths"),
               read.table("data/rf_data/2008/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2008/M8_cluster_points"),
               read.table("data/rf_data/2008/M8_cluster_lengths"),
               read.table("data/rf_data/2008/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2008/M9_cluster_points"),
               read.table("data/rf_data/2008/M9_cluster_lengths"),
               read.table("data/rf_data/2008/M9_cen_spell"))
h8=cbind(h_6,h_7,h_8,h_9)

# year 2009
h_6=func_star(read.table("data/rf_data/2009/M6_cluster_points"),
               read.table("data/rf_data/2009/M6_cluster_lengths"),
               read.table("data/rf_data/2009/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2009/M7_cluster_points"),
               read.table("data/rf_data/2009/M7_cluster_lengths"),
               read.table("data/rf_data/2009/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2009/M8_cluster_points"),
               read.table("data/rf_data/2009/M8_cluster_lengths"),
               read.table("data/rf_data/2009/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2009/M9_cluster_points"),
               read.table("data/rf_data/2009/M9_cluster_lengths"),
               read.table("data/rf_data/2009/M9_cen_spell"))
h9=cbind(h_6,h_7,h_8,h_9)

# year 2010
h_6=func_star(read.table("data/rf_data/2010/M6_cluster_points"),
               read.table("data/rf_data/2010/M6_cluster_lengths"),
               read.table("data/rf_data/2010/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2010/M7_cluster_points"),
               read.table("data/rf_data/2010/M7_cluster_lengths"),
               read.table("data/rf_data/2010/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2010/M8_cluster_points"),
               read.table("data/rf_data/2010/M8_cluster_lengths"),
               read.table("data/rf_data/2010/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2010/M9_cluster_points"),
               read.table("data/rf_data/2010/M9_cluster_lengths"),
               read.table("data/rf_data/2010/M9_cen_spell"))
h10=cbind(h_6,h_7,h_8,h_9)

# year 2011
h_6=func_star(read.table("data/rf_data/2011/M6_cluster_points"),
               read.table("data/rf_data/2011/M6_cluster_lengths"),
               read.table("data/rf_data/2011/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2011/M7_cluster_points"),
               read.table("data/rf_data/2011/M7_cluster_lengths"),
               read.table("data/rf_data/2011/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2011/M8_cluster_points"),
               read.table("data/rf_data/2011/M8_cluster_lengths"),
               read.table("data/rf_data/2011/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2011/M9_cluster_points"),
               read.table("data/rf_data/2011/M9_cluster_lengths"),
               read.table("data/rf_data/2011/M9_cen_spell"))
h11=cbind(h_6,h_7,h_8,h_9)

# year 2012
h_6=func_star(read.table("data/rf_data/2012/M6_cluster_points"),
               read.table("data/rf_data/2012/M6_cluster_lengths"),
               read.table("data/rf_data/2012/M6_cen_spell"))
h_7=func_star(read.table("data/rf_data/2012/M7_cluster_points"),
               read.table("data/rf_data/2012/M7_cluster_lengths"),
               read.table("data/rf_data/2012/M7_cen_spell"))
h_8=func_star(read.table("data/rf_data/2012/M8_cluster_points"),
               read.table("data/rf_data/2012/M8_cluster_lengths"),
               read.table("data/rf_data/2012/M8_cen_spell"))
h_9=func_star(read.table("data/rf_data/2012/M9_cluster_points"),
               read.table("data/rf_data/2012/M9_cluster_lengths"),
               read.table("data/rf_data/2012/M9_cen_spell"))
h12=cbind(h_6,h_7,h_8,h_9)


############# Discretiziation of size 1000
# saveRDS(h2,file="data/starhull_complete_contours/2002.Rda") # Including censored contour
# saveRDS(h3,file="data/starhull_complete_contours/2003.Rda")
# saveRDS(h4,file="data/starhull_complete_contours/2004.Rda")
# saveRDS(h5,file="data/starhull_complete_contours/2005.Rda")
# saveRDS(h6,file="data/starhull_complete_contours/2006.Rda")
# saveRDS(h7,file="data/starhull_complete_contours/2007.Rda")
# saveRDS(h8,file="data/starhull_complete_contours/2008.Rda")
# saveRDS(h9,file="data/starhull_complete_contours/2009.Rda")
# saveRDS(h10,file="data/starhull_complete_contours/2010.Rda")
# saveRDS(h11,file="data/starhull_complete_contours/2011.Rda")
# saveRDS(h12,file="data/starhull_complete_contours/2012.Rda")
# 
# h0=cbind(h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12)



