# This code is licensed under the MIT License

# Input: Boundary contours of complete regions under rainfall over the years 2002-2012
# Output: Starhull error reported in Section 2.2 of the paper
# Execution time: approx 28 minutes (i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)

# load libraries
library(sp)
library(geosphere)
library(plyr)

print(date())
# Functions needed to measure error in starhull representation
sgnf=function(x,y){
  m1=x%*%t(rep(1,length(y)))
  m2=rep(1,length(x))%*%t(y)
  return(sign(m1-m2))
}

lray=function(c,cs){
  if (length(c)==1){len=0}else{
    len=sum((c[-1]-c[-length(c)])*(cs+1)/2)}
  return(len)
}

aray=function(c,cs){
  if(length(c)==1) {area=0}else{area=sum((c[-1]^2-c[-length(c)]^2)*(cs+1)/2)}
  return(0.5*area)
}

symmdist=function(x,y){
  xs=sign(x)
  ys=sign(y)
  x=c(0,abs(x))
  y=c(0,abs(y))
  lx=length(unique(x))
  ly=length(unique(y))
  if ((lx==1)|| (ly==1)) {h=(1+2)}
  if (lx>1 && ly>1){
    sd=sort(unique(c(x,y)))
    sd=sd[-1]
    jx=ceiling(rowSums((sgnf(sd,x)+1)/2)+1/2)-1
    jy=ceiling(rowSums((sgnf(sd,y)+1)/2)+1/2)-1
    xss=xs[jx]
    yss=ys[jy]
    xss[is.na(xss)]= -1
    yss[is.na(yss)]= -1
    sds=-xss[jx]*yss[jy]
    sd=c(0,sd)
  }
  if (lx==1 && ly>1){
    sd=y
    sds=ys  
  }
  if (lx>1 && ly==1){
    sd=x
    sds=xs
  }
  if (lx==1 && ly==1){
    sd=0
    sds=NULL
  }
  
  sds=c(sds,0)  # to make it same elength
  return(cbind(sd,sds))
}


function_2=function(xx,yy,zz){
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
      R=6371
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
  s1=s=1000                            # Number of equispaced points from the contour
  radius_dist0=list()
  radius_dist=array(0,dim=c(s1,m))
  star_error=actual_area=array(0,dim=m)
  ratio_d=array(0,dim=m)
  c1=array(0,dim=m)
  theta1=theta=seq(0,2*pi,length.out=(s1+1))[-(s1+1)]
  f1=function(x,y){ifelse(x < 0, atan(y / x) + pi, ifelse(y < 0 , atan(y / x) + 2*pi, atan(y / x)))}
  
  for(j in 1: m){
    X=contour_prj[[n1[j]]]               # Main data point
    hpts <- chull(X)                    # Convex hull construction
    hpts <- c(hpts, hpts[1])
    xy1=X[hpts, ]                       # Convex hull points
    l <- SpatialLines(list(Lines(Line(rbind(xy1, xy1[1, ])), ID=1))) # Now create a SpatialLines object from the vertices.
    pts <- spsample(l, 500, type="regular")
    pts_1=coordinates(pts)              # Equispaced point from the convex hull
    mm11=colMeans(pts_1)           # Center of the convex hull (same as above)
    # Below is the arrangements of the data points according to the centered  convex hull
    XY=cbind(X[,1]-mm11[1],X[,2]-mm11[2])   # This is the polygon whose convex hull is centered
    # plot(XY,type="l",xlab="West-East",ylab="South-North",col="red") # Correct plot with shift(for centering)
    
  # Star-hull construction 
    theta22=NULL
    r_theta22=NULL
    theta220=NULL
    r_theta220=NULL
    for (i1 in 1:(nrow(XY)-1)){
      xy1=XY[i1,]                                 # one vertex of a line segment of a polygon
      xy2=XY[i1+1,]                              # other vertex
      
      t1=f1(xy1[1],xy1[2])                     # theta of polar representation of the first vertex
      r1=sqrt(xy1[1]^2+xy1[2]^2)              # length
      
      t2=f1(xy2[1],xy2[2])
      r2=sqrt(xy2[1]^2+xy2[2]^2)
      
      theta11=NULL                               # theta variable for actual contour
      r_theta=NULL                                   # corresponding radius
      
      if  (abs(t2-t1)<pi) {theta11=theta1[intersect(which(theta1>=min(t1,t2)),which(theta1<=max(t1,t2)))]} else
      {theta11=theta1[c(which(theta1>=max(t1,t2)),which(theta1<=min(t1,t2)))]}
      
      r_theta=(r1*r2*sin(t1-t2))/(r1*sin(t1-theta11)+r2*sin(theta11-t2))
      theta22=c(theta22,theta11)
      r_theta22=c(r_theta22,r_theta)
      }
  # To detect whether the centoid is inside or outside the contour
    l1=array(0,dim=s)  # Number of cut points for each directions
    for (i2 in 1:s){
      l1[i2]=length(r_theta22[which(theta22==theta1[i2])])
    }
    if (min(l1)>=1){ #  Atleast one intersect in all directions
      for (i2 in 1:s){
        r_list=NULL
        r_list=sort(r_theta22[which(theta22==theta1[i2])])
        radius_dist[i2,j]=max(r_list)    # Star-hull contour
        sn=array(1,dim=length(r_list))
        sn[which(1:length(r_list)%%2==0)]=-1
        radius_dist0[[i2]]=r_list*sn
      }
      if (min(l1)>=2){ # a typical situation as in example 564 of out of sample
        for (i2 in 1:s){
          r_list=NULL
          r_list=sort(r_theta22[which(theta22==theta1[i2])])
          radius_dist[i2,j]=max(r_list)    # Star-hull contour
          sn=array(1,dim=length(r_list))
          sn[which(1:length(r_list)%%2!=0)]=-1
          radius_dist0[[i2]]=r_list*sn
        }
      }
    }else{  # When centroid outside the contour and few directions has no intersection
      for (i2 in 1:s){
        r_list=NULL
        r_list=sort(r_theta22[which(theta22==theta1[i2])])
        if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
        radius_dist[i2,j]=max(r_list)    # Star-hull contour
        
        if (length(r_list)==1) {radius_dist0[[i2]]=r_list}else{
          if (length(r_list)%%2==0){
            sn=array(1,dim=length(r_list))
            sn[which(1:length(r_list)%%2!=0)]=-1
            radius_dist0[[i2]]=r_list*sn
          }else{
            sn=array(1,dim=length(r_list))
            sn[which(1:length(r_list)%%2==0)]=-1
            radius_dist0[[i2]]=r_list*sn
          }
        }
      }
    }
    
    # Total actual area 
    a0=astar=array(0,dim=length(theta))
    for (i in 1:length(theta)){
    c=c(0,unlist(radius_dist0[[i]]))
    cs=sign(unlist(radius_dist0[[i]]))
      a0[i]=aray(c,cs)
      dstar=symmdist(unlist(radius_dist0[[i]]),radius_dist[i,j]) # With parametric method
      astar[i]=aray(dstar[,1],dstar[-nrow(dstar),2])
      }
    actual_area[j]=(2*pi/s)*sum(a0)
    star_error[j]=(2*pi/s)*sum(astar)
     # print(c(j,star_error[j]/actual_area[j],date()))
  ratio_d[j]=star_error[j]/actual_area[j]
  }
return(ratio_d)
}  
    
  
h2=h3=h4=h5=h6=h7=h8=h9=h10=h11=h12=NULL

# Each of the block corresponds a year
h_6=function_2(read.table("data/rf_data/2002/M67_cluster_points"),
               read.table("data/rf_data/2002/M67_cluster_lengths"),
               read.table("data/rf_data/2002/M67_cen_spell"))
h_7=function_2(read.table("data/rf_data/2002/M89_cluster_points"),
               read.table("data/rf_data/2002/M89_cluster_lengths"),
               read.table("data/rf_data/2002/M89_cen_spell"))
h2=c(h_6,h_7)

h_6=function_2(read.table("data/rf_data/2003/M67_cluster_points"),
               read.table("data/rf_data/2003/M67_cluster_lengths"),
               read.table("data/rf_data/2003/M67_cen_spell"))
h_7=function_2(read.table("data/rf_data/2003/M89_cluster_points"),
               read.table("data/rf_data/2003/M89_cluster_lengths"),
               read.table("data/rf_data/2003/M89_cen_spell"))
h3=c(h_6,h_7)


h_6=function_2(read.table("data/rf_data/2004/M6_cluster_points"),
               read.table("data/rf_data/2004/M6_cluster_lengths"),
               read.table("data/rf_data/2004/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2004/M7_cluster_points"),
               read.table("data/rf_data/2004/M7_cluster_lengths"),
               read.table("data/rf_data/2004/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2004/M8_cluster_points"),
               read.table("data/rf_data/2004/M8_cluster_lengths"),
               read.table("data/rf_data/2004/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2004/M9_cluster_points"),
               read.table("data/rf_data/2004/M9_cluster_lengths"),
               read.table("data/rf_data/2004/M9_cen_spell"))
h4=c(h_6,h_7,h_8,h_9)


h_6=function_2(read.table("data/rf_data/2005/M6_cluster_points"),
               read.table("data/rf_data/2005/M6_cluster_lengths"),
               read.table("data/rf_data/2005/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2005/M7_cluster_points"),
               read.table("data/rf_data/2005/M7_cluster_lengths"),
               read.table("data/rf_data/2005/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2005/M8_cluster_points"),
               read.table("data/rf_data/2005/M8_cluster_lengths"),
               read.table("data/rf_data/2005/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2005/M9_cluster_points"),
               read.table("data/rf_data/2005/M9_cluster_lengths"),
               read.table("data/rf_data/2005/M9_cen_spell"))
h5=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2006/M6_cluster_points"),
               read.table("data/rf_data/2006/M6_cluster_lengths"),
               read.table("data/rf_data/2006/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2006/M7_cluster_points"),
               read.table("data/rf_data/2006/M7_cluster_lengths"),
               read.table("data/rf_data/2006/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2006/M8_cluster_points"),
               read.table("data/rf_data/2006/M8_cluster_lengths"),
               read.table("data/rf_data/2006/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2006/M9_cluster_points"),
               read.table("data/rf_data/2006/M9_cluster_lengths"),
               read.table("data/rf_data/2006/M9_cen_spell"))
h6=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2007/M6_cluster_points"),
               read.table("data/rf_data/2007/M6_cluster_lengths"),
               read.table("data/rf_data/2007/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2007/M7_cluster_points"),
               read.table("data/rf_data/2007/M7_cluster_lengths"),
               read.table("data/rf_data/2007/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2007/M8_cluster_points"),
               read.table("data/rf_data/2007/M8_cluster_lengths"),
               read.table("data/rf_data/2007/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2007/M9_cluster_points"),
               read.table("data/rf_data/2007/M9_cluster_lengths"),
               read.table("data/rf_data/2007/M9_cen_spell"))
h7=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2008/M6_cluster_points"),
               read.table("data/rf_data/2008/M6_cluster_lengths"),
               read.table("data/rf_data/2008/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2008/M7_cluster_points"),
               read.table("data/rf_data/2008/M7_cluster_lengths"),
               read.table("data/rf_data/2008/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2008/M8_cluster_points"),
               read.table("data/rf_data/2008/M8_cluster_lengths"),
               read.table("data/rf_data/2008/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2008/M9_cluster_points"),
               read.table("data/rf_data/2008/M9_cluster_lengths"),
               read.table("data/rf_data/2008/M9_cen_spell"))
h8=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2009/M6_cluster_points"),
               read.table("data/rf_data/2009/M6_cluster_lengths"),
               read.table("data/rf_data/2009/M6_cen_spell"))

h_7=function_2(read.table("data/rf_data/2009/M7_cluster_points"),
               read.table("data/rf_data/2009/M7_cluster_lengths"),
               read.table("data/rf_data/2009/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2009/M8_cluster_points"),
               read.table("data/rf_data/2009/M8_cluster_lengths"),
               read.table("data/rf_data/2009/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2009/M9_cluster_points"),
               read.table("data/rf_data/2009/M9_cluster_lengths"),
               read.table("data/rf_data/2009/M9_cen_spell"))
h9=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2010/M6_cluster_points"),
               read.table("data/rf_data/2010/M6_cluster_lengths"),
               read.table("data/rf_data/2010/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2010/M7_cluster_points"),
               read.table("data/rf_data/2010/M7_cluster_lengths"),
               read.table("data/rf_data/2010/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2010/M8_cluster_points"),
               read.table("data/rf_data/2010/M8_cluster_lengths"),
               read.table("data/rf_data/2010/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2010/M9_cluster_points"),
               read.table("data/rf_data/2010/M9_cluster_lengths"),
               read.table("data/rf_data/2010/M9_cen_spell"))
h10=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2011/M6_cluster_points"),
               read.table("data/rf_data/2011/M6_cluster_lengths"),
               read.table("data/rf_data/2011/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2011/M7_cluster_points"),
               read.table("data/rf_data/2011/M7_cluster_lengths"),
               read.table("data/rf_data/2011/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2011/M8_cluster_points"),
               read.table("data/rf_data/2011/M8_cluster_lengths"),
               read.table("data/rf_data/2011/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2011/M9_cluster_points"),
               read.table("data/rf_data/2011/M9_cluster_lengths"),
               read.table("data/rf_data/2011/M9_cen_spell"))
h11=c(h_6,h_7,h_8,h_9)

h_6=function_2(read.table("data/rf_data/2012/M6_cluster_points"),
               read.table("data/rf_data/2012/M6_cluster_lengths"),
               read.table("data/rf_data/2012/M6_cen_spell"))
h_7=function_2(read.table("data/rf_data/2012/M7_cluster_points"),
               read.table("data/rf_data/2012/M7_cluster_lengths"),
               read.table("data/rf_data/2012/M7_cen_spell"))
h_8=function_2(read.table("data/rf_data/2012/M8_cluster_points"),
               read.table("data/rf_data/2012/M8_cluster_lengths"),
               read.table("data/rf_data/2012/M8_cen_spell"))
h_9=function_2(read.table("data/rf_data/2012/M9_cluster_points"),
               read.table("data/rf_data/2012/M9_cluster_lengths"),
               read.table("data/rf_data/2012/M9_cen_spell"))
h12=c(h_6,h_7,h_8,h_9)

h0=c(h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12)
print(date())