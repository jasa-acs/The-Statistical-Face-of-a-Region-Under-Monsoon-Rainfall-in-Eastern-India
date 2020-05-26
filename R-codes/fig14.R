# This code is licensed under the MIT License

# Input: Actual and Starhull contours over the years 2002-2012
# Output: Figure 14 of the paper
# Execution time: approx 70 hours (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)

# load libraries
library(fdasrvf)
library(sp)
library(geosphere)
library(plyr)
library(abind)
library(survival)

# set directory
# setwd("/home")

# function to compute grided contour
function_2 <- function(xx,yy,zz){
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
    if (nrow(obj_1)>=4){     
      m00<-data.frame(x=obj_1[,1], y=obj_1[,2])
      area1[i]=areaPolygon(m00)/(10^6)            # The area of a polygon incoordinates 
      m1=rbind(m00,m00[1,])
      m11=c(min(m1[,1]),min(m1[,2]))
      m12=c(max(m1[,1]),min(m1[,2]))
      m21=c(min(m1[,1]),max(m1[,2]))
      m22=c(max(m1[,1]),max(m1[,2]))
      bounding_box=rbind(m11,m12,m22,m21,m11)    # Vertices of the bounding box
      cm1=(min(m1[,1])+max(m1[,1]))/2            # Long-Coordinate of the center of the bounding box
      cm2=(min(m1[,2])+max(m1[,2]))/2            # Lat-Coordinate of the center of the bounding box
      box_center=c(cm1,cm2)                      # Center of the bounding box
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
  XY=list()
  shape_data1=array(0,dim=c(2,1001,m))
  for(j in 1: m){
    X=contour_prj[[n1[j]]]                      # Main data point
    hpts <- chull(X)                            # Convex hull construction
    hpts <- c(hpts, hpts[1])
    xy1=X[hpts, ]                               # Convex hull points
    l <- SpatialLines(list(Lines(Line(rbind(xy1, xy1[1, ])), ID=1))) # Now create a SpatialLines object from the vertices.
    pts <- spsample(l, 1000, type="regular")
    pts_1=coordinates(pts)                       # Equispaced point from the convex hull
    mm11=colMeans(pts_1)                         # Center of the convex hull (same as above)
    xy11=cbind(X[,1]-mm11[1],X[,2]-mm11[2])      # This is the polygon whose convex hull is centered
    l1 <- SpatialLines(list(Lines(Line(rbind(xy11, xy11[1, ])), ID=1))) # Now create a SpatialLines object from the vertices.
    pts11 <- spsample(l1, 1000, type="regular")
    pts=data.frame(pts11)
    pts<-rbind(pts,pts[1,])
    shape_data1[,,j]=t(coordinates(pts))  
  }
  shape_data1
}

# Computation of the actual rainfall regions using the above function 
# Each of the following block corresponds to a year

# year 2002
h_6=function_2(read.table("data/rf_data/2002/M67_cluster_points"),
               read.table("data/rf_data/2002/M67_cluster_lengths"),
               read.table("data/rf_data/2002/M67_cen_spell"))
h_7=function_2(read.table("data/rf_data/2002/M89_cluster_points"),
               read.table("data/rf_data/2002/M89_cluster_lengths"),
               read.table("data/rf_data/2002/M89_cen_spell"))
h2=abind(h_6,h_7,along=3)

# year 2003
h_6=function_2(read.table("data/rf_data/2003/M67_cluster_points"),
               read.table("data/rf_data/2003/M67_cluster_lengths"),
               read.table("data/rf_data/2003/M67_cen_spell"))
h_7=function_2(read.table("data/rf_data/2003/M89_cluster_points"),
               read.table("data/rf_data/2003/M89_cluster_lengths"),
               read.table("data/rf_data/2003/M89_cen_spell"))
h33=abind(h_6,h_7,along=3)

# year 2004
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
h44=abind(h_6,h_7,h_8,h_9,along=3)

# year 2005
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
h5=abind(h_6,h_7,h_8,h_9,along=3)

# year 2006
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
h6=abind(h_6,h_7,h_8,h_9,along=3)

# year 2007
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
h7=abind(h_6,h_7,h_8,h_9,along=3)

# year 2008
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
h8=abind(h_6,h_7,h_8,h_9,along=3)

# year 2009
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
h9=abind(h_6,h_7,h_8,h_9,along=3)

# year 2010
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
h10=abind(h_6,h_7,h_8,h_9,along=3)

# year 2011
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
h11=abind(h_6,h_7,h_8,h_9,along=3)

# year 2012
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
h12=abind(h_6,h_7,h_8,h_9,along=3)

shape_data=abind(h2,h33,h44,h5,h6,h7,h8,h9,h10,h11,h12,along=3)


# Star-hull of all complete bounded (>200 and <13500) data

files3<- list.files("data/starhull_complete_contours", pattern = "\\.Rda$",full.names = TRUE) 
data1=NULL
n1=NULL
for (j in 1:11){
  dt3=NULL
  dt3=readRDS(file=files3[j])           # all contours of a year 
  n1[j]=ncol(dt3)
  data1=cbind(data1,dt3)               # all complete larger contour function
}
data1=data1   
data_total=data1  
n11=cumsum(n1)
s=1000
N11=ncol(data1)                        # Number of complete contours
r_u=(rank(data1,ties.method = c("average")))/(s*N11) # ECDF
r_u=replace(r_u,which(r_u==1),0.99999) 
data2=matrix(qnorm(r_u),nrow=s)         # Normal (random variable) transformation
data_star=data2

dt3=readRDS(file="data/area_status_combined/area.Rda")
st= readRDS(file="data/area_status_combined/status.Rda")

dt4=intersect(which(dt3>=200),which(dt3<=13500))  # All contours in the bound
dt1=which(st==1)                                  # censored regions
dt2=which(st==2)                                  # complete regions
dt11=intersect(dt1,dt4)                           # range (> 200 sq km && <13500) censored spell
dt22=intersect(dt2,dt4)[-80]                    
dt22=intersect(dt2,dt4)                           # area (> 200 sq km && <13500) complete spell (except one or two)
N1=length(dt22)                                   # Number of complete observations
N2=length(dt11)                                   # Number of censored observations
time1=NULL
censor1=NULL
dt44=c(dt11,dt22)
time1=dt3[dt44]                                  # Region (Area) which lies in the given threshold
censor1=st[dt44]
censor1[censor1 < 1.5]=0                         # Censored: 0
censor1[censor1 > 1.5]=1                         # Complete: 1
time_ind=order(time1)
area_order=time1[time_ind]
censor_order=censor1[time_ind]
area_status=cbind(area_order,censor_order)
c_max=max(which(area_status[,2]==1))
area_status1=area_status[1:c_max,]
time=NULL
censor=NULL
time=rank(area_status1[,1])
censor=area_status1[,2]
mini.surv <- survfit(Surv(time, censor)~ 1)   # Kaplan-Meier computation
T1=mini.surv$time                             # Increasing time period 
F1=mini.surv$surv                            # Estimate of survial function
F2=c(1,F1[-length(F1)])-F1
F3=cbind(T1,F2)                               # (index,prob. mass)
F4=c(F3[which(F3[,2]>0),2])                   # Weights
order_com=order(dt3[dt22])
x0=order_com                                   # Order of the complete contours 
N1=length(x0)

# Out of sample 
l1=length(x0)
l2=length(F4)
l3=8921                                       # eighty percent of the complete observations
set.seed(10)
h1=sample(c(1:length(x0)), l3, replace = FALSE, prob = NULL) # eighty percent of random sample
h10=setdiff(c(1:length(x0)),h1)                              # the remaining out-sample
order1=x0[h1]                                                # indices for in-sample
wts0=F4[h1]                                                  # corresponding weights of in-sample
wts=wts0/sum(wts0)                                          # weights normalizations
mean_x=c(data_star[,order1]%*%wts)                          # Weighted mean
cov1=(tcrossprod((data_star[,order1]-mean_x)*(matrix(sqrt(wts), nrow = s, ncol = l3, byrow = TRUE))))  
h1=eigen(cov1)$values 
hv=eigen(cov1)$vectors
MM=12   
pc_1=((h1[1:MM])/sum(h1))*100;pc_1 # Percentage of the variation explained by the PCS
sum(pc_1)



# SRFV Computation
# Prerequisite function (obtained from FDASRVF package)
pvecnorm <-function(v,p=2){
  sum(abs(v)^p)^(1/p)
}

#
innerprod_q2 <- function(q1, q2){
  T1 = ncol(q1)
  val = sum(sum(q1*q2))/T1
  return(val)
}

#
curve_to_q1 <- function(beta){
  n = nrow(beta)
  T1 = ncol(beta)
  v = apply(beta,1,gradient, 1.0/(T1-1))
  v = t(v)
  q = matrix(0,n,T1)
  for (i in 1:T1){
    L = sqrt(pvecnorm(v[,i],2))
    if (L>0.0001){
      q[,i] = v[,i]/L
    } else {
      q[,i] = v[,i]*0.0001
    }
  }
  # q = q/sqrt(innerprod_q2(q, q)) # Scaling is being omitted
  return(q)
}

#
ndims <- function(x){
  return(length(dim(x)))
}

#
repmat1 <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  if (is.null(mx)){
    mx = 1
    nx = length(X)
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }else {
    nx = dim(X)[2]
    mat = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  return(mat)
}

#
trapz <- function(x,y,dims=1){
  if ((dims-1)>0){
    perm = c(dims:max(ndims(y),dims), 1:(dims-1))
  } else {
    perm = c(dims:max(ndims(y),dims))
  }
  if (ndims(y) == 0){
    m = 1
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
  }
  if (m==1){
    M = length(y)
    out = sum(diff(x)*(y[-M]+y[-1])/2)
  } else {
    slice1 = y[as.vector(outer(1:(m-1), dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+')) ] 
    dim(slice1) = c(m-1, length(slice1)/(m-1)) 
    slice2 = y[as.vector(outer(2:m, dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+'))] 
    dim(slice2) = c(m-1, length(slice2)/(m-1)) 
    out = t(diff(x)) %*% (slice1+slice2)/2.
    siz = dim(y)
    siz[1] = 1
    out = array(out, siz)
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    out = aperm(out, perm2)
    ind = which(dim(out) != 1)
    out = array(out, dim(out)[ind])
  }
  return(out)
}

#
calculatecentroid1 <- function(beta){
  n = nrow(beta)
  T1 = ncol(beta)
  betadot = apply(beta,1,gradient,1.0/(T1-1))
  betadot = t(betadot)
  normbetadot = apply(betadot,2,pvecnorm,2)
  integrand = matrix(0, n, T1)
  for (i in 1:T1){
    integrand[,i] = beta[,i] * normbetadot[i]
  }
  scale = trapz(seq(0,1,length.out=T1), normbetadot)
  centroid = trapz(seq(0,1,length.out=T1), integrand, 2)/scale
  return(centroid)
}
reparam_curve1 <- function(beta1,beta2,lambda=0,method="DP",w=0.01,rotated=F,
                           isclosed=T, mode="C"){
  n1 = nrow(beta2)
  M = ncol(beta2)
  timet = seq(0,1,length.out=M)
  skipm = 4
  auto = 2
  tau = 0
  if (method=="DPo"){
    # Optimize over SO(n) x Gamma
    q1 = curve_to_q1(beta1)
    # Optimize over SO(n)
    if (rotated){
      out = find_rotation_seed_coord(beta1, beta2, mode)
      beta2 = out$beta2
      R = out$O_hat
      tau = out$tau
    } else{
      R = diag(n1)
      tau = 0
    }
    q2 = curve_to_q1(beta2)
    # Optimize over Gamma
    q1i = q1
    dim(q1i) = c(M*n1)
    q2i = q2
    dim(q2i) = c(M*n1)
    G = rep(0,M)
    T1 = rep(0,M)
    size = 0
    ret = .Call('DPQ2', PACKAGE = 'fdasrvf', q1i, timet, q2i, timet, n1, M, M, timet, timet, M, M, G, T1, size, lambda);
    G = ret$G[1:ret$size]
    Tf = ret$T[1:ret$size]
    gam0 = approx(Tf,G,xout=timet)$y
  } else if (method=="DP") {
    # Optimize over SO(n) x Gamma
    q1 = curve_to_q1(beta1)
    # Optimize over SO(n)
    if (rotated){
      out = find_rotation_seed_coord(beta1, beta2);
      beta2 = out$beta2
      R = out$O_hat
      tau = out$tau
    } else{
      R = diag(n1)
      tau = 0
    }
    q2 = curve_to_q1(beta2)
    # Optimize over Gamma
    q1i = q1
    dim(q1i) = c(M*n1)
    q2i = q2
    dim(q2i) = c(M*n1)
    gam0 = .Call('DPQ', PACKAGE = 'fdasrvf', q1i, q2i, n1, M, lambda, 0, rep(0,M))
  } else if (method=="DP2") {
    c1 = t(beta1)
    dim(c1) = c(M*n1)
    c2 = t(beta2)
    dim(c2) = c(M*n1)
    opt = rep(0,M+n1*n1+1)
    swap = FALSE
    fopts = rep(0,5)
    comtime = rep(0,5)
    out = .Call('opt_reparam', PACKAGE = 'fdasrvf', c1,c2,M,n1,0.0,TRUE,
                rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)
    tmp = length(out$opt)
    gam0 = out$opt[1:(tmp-5)]
    R = matrix(out$opt[(tmp-4):(tmp-1)],nrow=2)
    if (out$swap){
      gam0 = invertGamma(gam0)
      R = t(R)
    }
  } else {
    c1 = t(beta1)
    dim(c1) = c(M*n1)
    c2 = t(beta2)
    dim(c2) = c(M*n1)
    opt = rep(0,M+n1*n1+1)
    swap = FALSE
    fopts = rep(0,5)
    comtime = rep(0,5)
    out = .Call('opt_reparam', PACKAGE = 'fdasrvf', c1,c2,M,n1,w,FALSE,
                rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)
    if (out$fopts[1] == 1000){
      out = .Call('opt_reparam', PACKAGE = 'fdasrvf', c1,c2,M,n1,0.0,TRUE,
                  rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)
    }
    tmp = length(out$opt)
    gam0 = out$opt[1:(tmp-5)]
    R = matrix(out$opt[(tmp-4):(tmp-1)],nrow=2)
    if (out$swap){
      gam0 = invertGamma(gam0);
      R = t(R)
    }
  }
  gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale
  return(list(gam=gam,R=R,tau=tau))
}

# Inversion of gamma 
invertGamma1 <- function(gam){
  N = length(gam)
  x = (0:(N-1))/(N-1)
  gamI = approx(gam,x,xout=x)$y
  gamI[N] = 1
  gamI = gamI/gamI[N]
  return(gamI)
}

### Group action by gamma ###
group_action_by_gamma_coord1 <- function(f, gamma){
  n = nrow(f)
  T1 = ncol(f)
  fn = matrix(0, n, T1)
  timet = seq(0, 1, length.out = T1)
  for (j in 1:n){
    fn[j,] = spline(timet, f[j,], xout=gamma)$y
  }
  return(fn)
}

#
find_basis_normal1<- function(q){
  n = nrow(q)
  T1 = ncol(q)
  f1 = matrix(0,n,T1)
  f2 = matrix(0,n,T1)
  for (i in 1:T1){
    f1[,i] = q[1,i]*q[,i]/pvecnorm(q[,i])+c(pvecnorm(q[,i]),0)
    f2[,i] = q[2,i]*q[,i]/pvecnorm(q[,i])+c(0,pvecnorm(q[,i]))
  }
  h3 = f1
  h4 = f2
  integrandb3 = rep(0,T1)
  integrandb4 = rep(0,T1)
  for (i in 1:T1){
    integrandb3[i] = t(q[,i])%*%h3[,i]
    integrandb4[i] = t(q[,i])%*%h4[,i]
  }
  b3 = h3 - q*trapz(seq(0,1,length.out=T1),integrandb3)
  b4 = h4 - q*trapz(seq(0,1,length.out=T1),integrandb4)
  
  basis = list(b3,b4)
  return(basis)
}

#
project_curve1 <- function(q){
  T1 = ncol(q)
  n = nrow(q)
  if(n==2){
    dt = 0.35
  }
  if(n==3){
    dt = 0.2
  }
  epsilon =- 1e-6
  e = diag(1,n)
  iter = 1
  res = rep(1,n)
  J = matrix(0,n,n)
  s = seq(0,1,length.out=T1)
  qnorm = rep(0,T1)
  G = rep(0,n)
  C = rep(0,301)
  qnew = q
  
  while (pvecnorm(res) > epsilon ){
    if (iter > 300 ){ 
      break
    }
# Compute Jacobian
    for (i in 1:n){
      for (j in 1:n){
        J[i,j]  = 3 * trapz(s, qnew[i,]*qnew[j,])
      }
    }
    J = J + diag(1,n)
    for (i in 1:T1){
      qnorm[i] = pvecnorm(qnew[,i])
    }
# Compute the residue
    for (i in 1:n){
      G[i] = trapz(s,qnew[i,]*qnorm)
    }
    res = -1*G
    
    if (pvecnorm(res)<epsilon){
      break
    }
    x = solve(J,res)
    C[iter] = pvecnorm(res)
    delG = find_basis_normal1(qnew)
    tmp = 0
    for (i in 1:n){
      tmp = tmp + x[i]*delG[[i]]*dt
    }
    qnew = qnew + tmp
    iter = iter + 1
  }
  qnew
  return(qnew)
}

#
inverse_exp_coord1 <- function(beta1, beta2, mode="C", rotated=F){
  T1 = ncol(beta1)
  centroid1 = calculatecentroid1(beta1)
  dim(centroid1) = c(length(centroid1),1)
  beta1 = beta1 - repmat1(centroid1, 1, T1)
  centroid2 = calculatecentroid1(beta2)
  dim(centroid2) = c(length(centroid2),1)
  beta2 = beta2 - repmat1(centroid2, 1, T1)
  q1 = curve_to_q1(beta1)
  if (mode=="C"){
    isclosed = TRUE
  }
  #Iteratively optimize over SO(n) x Gamma using old DP
  out = reparam_curve1(beta1, beta2, rotated=rotated, isclosed=isclosed, mode=mode)
  if (mode=="C")
    gamI = invertGamma1(out$gam)
  beta2 = group_action_by_gamma_coord1(beta2, gamI)
  q2n = curve_to_q1(beta2)
  if (mode=="C"){   
    q2n = project_curve1(q2n)
  }
  dist=sqrt(sum((q1-q2n)^2)*(1/T1))  # Can be improved 
  v= q1-q2n
  return(list(v=v,dist=dist))
}

#
gram_schmidt1 <- function(basis){
  b1 = basis[[1]]
  b2 = basis[[2]]
  basis1 = b1 / sqrt(innerprod_q2(b1,b1))
  b2 = b2 - innerprod_q2(basis1,b2)*basis1
  basis2 = b2 / sqrt(innerprod_q2(b2,b2))
  basis_o = list(basis1, basis2)
  return(basis_o)
}

#
project_tangent1 <- function(w, q, basis){
  w = w - innerprod_q2(w,q)*q
  bo = gram_schmidt1(basis)
  wproj = w - innerprod_q2(w, bo[[1]])*bo[[1]] - innerprod_q2(w,bo[[2]])*bo[[2]]
  return(wproj)
}

#
karcher_calc1 <- function(beta, q, betamean, mu, rotated=F, mode="C"){
  if (mode=="C"){
    basis = find_basis_normal1(mu)
  }
  # Compute shooting vector form mu to q_i
  out = inverse_exp_coord1(betamean, beta, mode, rotated)
  # Project to tangent space of manifold to obtain v_i
  if (mode=="O"){
    v = out$v
  } else {
    v = project_tangent1(out$v, q, basis)
  }
  return(list(v=v,d=out$dist))
}
# Karcher mean
curve_karcher_mean1 <- function(beta, mode="C", rotated=F, maxit=maxit0){
  tmp = dim(beta)
  n = tmp[1]
  T1 = tmp[2]
  N = tmp[3]
  q = array(0, c(n,T1,N))
  for (ii in 1:N){
    centroid1 = calculatecentroid1(beta[,,ii])
    dim(centroid1) = c(length(centroid1),1)
    beta[,,ii] = beta[,,ii] - repmat1(centroid1, 1, T1)
    q[,,ii] = curve_to_q1(beta[,,ii])
  }
  ##Initialize mu as one of the shapes
  mnq = rowMeans(q[1,,])
  dqq = sqrt(colSums((q[1,,] - matrix(mnq,ncol=N,nrow=T1))^2))
  min_ind = which.min(dqq)
  mu = q[,,min_ind]
  betamean = beta[,,min_ind]
  # compute average direction of tangent vectors v_i
  delta = 0.3
  tolv = 0.00001
  told = 5*1e-3
  itr = 1
  sumd = rep(0,maxit+1)
  v = array(0,c(n,T1,N))
  normvbar = rep(0,maxit+1)
  while (itr<maxit){
    cat(sprintf("Iteration: %d\n",itr))
    sumv = matrix(0,2,T1)
    sumd[itr] = 0
    for (i in 1:N){
      v[,,i] = q[,,i]-mu
    }
    sumv = rowSums(v,dims=2)
    # compute average direction of tangent vectors v_i
    vbar = sumv/N
    normvbar[itr] = sqrt(innerprod_q2(vbar,vbar))
    normv = normvbar[itr]
    if (normv>tolv){
      mu = mu + delta*vbar # Changed for L2
      # For closed curve 
      if (mode=="C"){    
        mu=mu
        mu = project_curve1(mu)
      }
      x = q_to_curve(mu)
      a = -1 * calculatecentroid1(x)
      dim(a) = c(length(a),1)
      betamean = x + repmat1(a,1,T1)
    } else {
      break
    }
    itr = itr + 1
  }
  return(list(mu=mu,betamean=betamean,v=v,q=q))
}

## Align curves to Karcher mean  
curve_srvf_align1 <- function(beta, mode="C", rotated=F, maxit=maxit0){
  if (mode=="C"){isclosed = TRUE}
  tmp = dim(beta)
  n = tmp[1]
  T1 = tmp[2]
  N = tmp[3]
  out = curve_karcher_mean1(beta,mode, rotated, maxit)
  mu = out$mu
  betamean = out$betamean
  v = out$v
  q = out$q
  qn = array(0, c(n,T1,N))
  betan = array(0, c(n,T1,N))
  centroid2 = calculatecentroid1(betamean)
  dim(centroid2) = c(length(centroid2),1)
  betamean = betamean - repmat1(centroid2, 1, T1)
  q_mu = curve_to_q1(betamean)
  # align to mean
  for (ii in 1:N){
    print(c(2,ii,date()))
    beta1 = beta[,,ii]
    centroid1 = calculatecentroid1(beta1)
    dim(centroid1) = c(length(centroid1),1)
    beta1 = beta1 - repmat1(centroid1,1,T1)
    # Iteratively optimize over SO(n) x Gamma
    out = reparam_curve1(betamean, beta1, isclosed = isclosed, rotated = F)
    gamI = invertGamma1(out$gam)
    if (mode=="C")
      # Apply optimal re-parameterization to the second curve
      beta1 = group_action_by_gamma_coord1(beta1, gamI)
    # Optimize over SO(n)
    if (rotated){
      out = find_rotation_seed_coord(betamean, beta1)
      qn[,,ii] = curve_to_q(out$beta2new)
      betan[,,ii] = out$beta2new
    } else {
      qn[,,ii] = curve_to_q1(beta1)
      betan[,,ii] = beta1
    }
  }
  return(list(betan=betan, qn=qn, betamean=betamean, q_mu=q_mu))
}

##  Karcher covariance for closed curve 
curve_karcher_cov1 <- function(beta_in,beta_out,order1,wts,mode="C",rotated=F,maxit=maxit0){
  tmp = dim(beta_in)
  n = tmp[1]
  T1 = tmp[2]
  N = tmp[3]
  out=curve_srvf_align1(beta_in, mode="C", rotated=F, maxit=maxit0)
  # Compute Karcher covariance of uniformly sampled mean
  betamean = resamplecurve(out$betamean, T1)
  mu = out$q_mu
  if (mode=="C"){
    mu = project_curve1(mu)
    basis = find_basis_normal1(mu)
  }
  q=out$qn
  v = array(0, c(n,T1,N))
  for (i in 1:N){
    # Project to the tangent space of manifold to obtain v_i
    v[,,i] = q[,,i]-mu
  }
  K = matrix(0, 2*T1, 2*T1)
  for (i in 1:N){
    w = v[,,i]
    w = c(w[1,], w[2,])
    K = K + (w %*% t(w))*wts[i]
  }
  ######### Alignment of the out-sample data #######
  tmp = dim(beta_out)
  n = tmp[1]
  T1 = tmp[2]
  N = tmp[3]
  qn_out = array(0, c(n,T1,N))
  #align to previpus mean: beatmean
  for (ii in 1:N){
    print(c(3,ii,date()))
    beta1 = beta_out[,,ii]
    centroid1 = calculatecentroid1(beta1)
    dim(centroid1) = c(length(centroid1),1)
    beta1 = beta1 - repmat1(centroid1,1,T1)
    #Iteratively optimize over SO(n) x Gamma
    out1 = reparam_curve1(betamean, beta1, isclosed = isclosed, rotated = F)
    gamI = invertGamma1(out1$gam)
    if (mode=="C")
      # Apply optimal re-parameterization to the second curve
      beta1 = group_action_by_gamma_coord1(beta1, gamI)
    # Optimize over SO(n)
    if (rotated){
      out2 = find_rotation_seed_coord(betamean, beta1)
      qn_out[,,ii] = curve_to_q(out2$beta2new)} else {
        qn_out[,,ii] = curve_to_q1(beta1)}
  }
  return(list(cov=K,mu=mu, qn=out$qn,qn_out=qn_out))
}
shape_data00=shape_data
shape_data_used=shape_data

maxit0=8                      
l1=length(x0)
l2=length(F4)
l3=8921                          # eighty percent of the complete observations
set.seed(10)
h1=sample(c(1:length(x0)), l3, replace = FALSE, prob = NULL) # eighty percent of random sample
h00=setdiff(c(1:length(x0)),h1)                        # the remaining out observations
order1=x0[h1]                                         # indices for in-observations
wts0=F4[h1]                                           # corresponding weights of in-observations
wts=wts0/sum(wts0)                                   # weights normalizations
shape_data_in=shape_data_used[,,h1]
shape_data_out=shape_data_used[,,h00]
K0=curve_karcher_cov1(shape_data_in,shape_data_out,order1,wts,mode="C",maxit=maxit0)

# SRVF FPCA 
K=K0$cov
h11=eigen(K)$values  # Eigen values
hv1=eigen(K)$vectors  # Eigen vectors
MM1=290 # Taking out first 290 PCs
print(sum(h11[1:MM1])/sum(h11))

## Star-hull of all complete bounded (>200 and <13500) data

k1=2002
m=dim(shape_data_out)[3]
data200=data_star[,h00]
data300=data1[,h00]

# FDA
theta=seq(0,2*pi,length=s+1)[-(s+1)]
x1=seq(1,1,length=1000)
x2=cos(theta)
x3=sin(theta)
x4=cos(2*theta)
x5=sin(2*theta)
x6=cos(3*theta)
x7=sin(3*theta)
x8=cos(4*theta)
x9=sin(4*theta)
x10=cos(5*theta)
x11=sin(5*theta)
x12=cos(6*theta)
x13=sin(6*theta)
x0=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13)
p_x=x0%*%(solve(t(x0)%*%x0))%*%(t(x0))  # Least square estimates of the LM
q_data=array(0,dim=c(2,1001,m))
s=1000
radius_dist=array(0,dim=c(s,m))
fda_error=array(0,dim=m)
srvf_error=array(0,dim=m)
para_error=array(0,dim=m)
# Function to compute symmetric difference 
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
  
  sds=c(sds,0)  # to make it sam elength
  return(cbind(sd,sds))
}

actual_area=array(0,dim=m)
star_error=array(0,dim=m)
c1=array(0,dim=m)
theta1=seq(0,2*pi,length.out=(s+1))[-(s+1)]
f1=function(x,y){ifelse(x < 0, atan(y / x) + pi, ifelse(y < 0 , atan(y / x) + 2*pi, atan(y / x)))}
a1=array(0,dim=c(m,s))


for(k in 1:m){ # Main loop
  # SRVF
  radius_dist0=list()
  radius_srvf=list()
  
  S1=array(0,dim=2002)
  for (i in 1:MM1)
  {
    d1=c(K0$qn_out[,,k][1,],K0$qn_out[,,k][2,])  # Change
    m1=c(K0$mu[1,],K0$mu[2,])
    score1=c(t(hv1[,i])%*%(d1-m1))
    S1=S1+score1*hv1[,i]
  }
  m2=S1+m1
  q_data[,,k]=rbind(m2[1:1001],m2[1002:2002])
  beta1=q_to_curve(q_data[,,k])
  T1=1001
  centroid1 = calculatecentroid1(beta1)
  dim(centroid1) = c(length(centroid1),1)
  beta1 = beta1 - repmat1(centroid1,1,T1)
  XY0=rbind(t(beta1)[-1001,],t(beta1)[1,]) # SRVF reconstructed contour
  XY=t(shape_data_out[,,k])                # Actual contour

  theta22=NULL
  r_theta22=NULL
  theta220=NULL
  r_theta220=NULL
  for (i1 in 1:(nrow(XY)-1)){
    # For actual contour 
    xy1=XY[i1,]                       # one vertex of a line segment of a polygon
    xy2=XY[i1+1,]                      # other vertex
    
    t1=f1(xy1[1],xy1[2])              # theta of polar representation of the first vertex
    r1=sqrt(xy1[1]^2+xy1[2]^2)        # length
    
    t2=f1(xy2[1],xy2[2])
    r2=sqrt(xy2[1]^2+xy2[2]^2)
    
    theta11=NULL                      # theta variable for actual contour
    r_theta=NULL                      # corresponding radius
    
    if  (abs(t2-t1)<pi) {theta11=theta1[intersect(which(theta1>=min(t1,t2)),which(theta1<=max(t1,t2)))]} else
    {theta11=theta1[c(which(theta1>=max(t1,t2)),which(theta1<=min(t1,t2)))]}
    
    r_theta=(r1*r2*sin(t1-t2))/(r1*sin(t1-theta11)+r2*sin(theta11-t2))
    
    theta22=c(theta22,theta11)
    r_theta22=c(r_theta22,r_theta)
    
    # Replication for srvf reconstruced contour
    
    xy10=XY0[i1,]          # one vertex of a line segment of a polygon
    xy20=XY0[i1+1,]        # other vertex
    
    t10=f1(xy10[1],xy10[2])         # theta of polar representation of the first vertex
    r10=sqrt(xy10[1]^2+xy10[2]^2)   # length
    
    t20=f1(xy20[1],xy20[2])
    r20=sqrt(xy20[1]^2+xy20[2]^2)
    
    theta110=NULL
    r_theta0=NULL
    
    if  (abs(t20-t10)<pi) {theta110=theta1[intersect(which(theta1>=min(t10,t20)),which(theta1<=max(t10,t20)))]} else
    {theta110=theta1[c(which(theta1>=max(t10,t20)),which(theta1<=min(t10,t20)))]}
    
    r_theta0=(r10*r20*sin(t10-t20))/(r10*sin(t10-theta110)+r20*sin(theta110-t20))
    
    theta220=c(theta220,theta110)
    r_theta220=c(r_theta220,r_theta0)
  }
  
  ###### To detect whether the centoid is inside or outside the contour
  l1=array(0,dim=s)                          # Number of cut points for each directions
  for (i2 in 1:s){
    l1[i2]=length(r_theta22[which(theta22==theta1[i2])])
  }
  if (min(l1)>=1){                          #  Atleast one intersect in all directions
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta22[which(theta22==theta1[i2])])
      radius_dist[i2,k]=max(r_list)            # Star-hull contour
      sn=array(1,dim=length(r_list))
      sn[which(1:length(r_list)%%2==0)]=-1
      radius_dist0[[i2]]=r_list*sn
    }
    if (min(l1)>=2){                          # a typical situation as in example 564 of out of sample
      for (i2 in 1:s){
        r_list=NULL
        r_list=sort(r_theta22[which(theta22==theta1[i2])])
        radius_dist[i2,k]=max(r_list)            # Star-hull contour
        sn=array(1,dim=length(r_list))
        sn[which(1:length(r_list)%%2!=0)]=-1
        radius_dist0[[i2]]=r_list*sn
      }
    }
  }else{                                   # When centroid outside the contour and few directions has no intersection
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta22[which(theta22==theta1[i2])])
      if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
      radius_dist[i2,k]=max(r_list)    # Star-hull contour
      
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

## For SRVF alignment: to detect whethere the centoid is inside or outside the contour
  l1=array(0,dim=s)                # Number of cut points for each directions
  for (i2 in 1:s){
    l1[i2]=length(r_theta220[which(theta220==theta1[i2])])
  }
  if (min(l1)>=1){               #  Atleast one intersect in all directions
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta220[which(theta220==theta1[i2])])
      sn=array(1,dim=length(r_list))
      sn[which(1:length(r_list)%%2==0)]=-1
      radius_srvf[[i2]]=r_list*sn
    }
    if (min(l1)>=2){               # a typical situation as in example 564 of out of sample
      for (i2 in 1:s){
        r_list=NULL
        r_list=sort(r_theta220[which(theta220==theta1[i2])])
        sn=array(1,dim=length(r_list))
        sn[which(1:length(r_list)%%2!=0)]=-1
        radius_srvf[[i2]]=r_list*sn
      }
    }
  }else{  # When centroid outside the contour and few directions has no intersection
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta220[which(theta220==theta1[i2])])
      if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
      
      if (length(r_list)==1) {radius_srvf[[i2]]=r_list}else{
        if (length(r_list)%%2==0){
          sn=array(1,dim=length(r_list))
          sn[which(1:length(r_list)%%2!=0)]=-1
          radius_srvf[[i2]]=r_list*sn
        }else{
          sn=array(1,dim=length(r_list))
          sn[which(1:length(r_list)%%2==0)]=-1
          radius_srvf[[i2]]=r_list*sn
        }
      }
    }
  }
  
  
###### Non-parametric computation
  S3=array(0,dim=c(s,1))
  m2=array(0,dim=c(s,1))
  for (l in 1:MM)
  {
    S3=S3+c((t(c(hv[,l]))%*%(c(data200[,k])-c(mean_x))))*hv[,l]
  }
  m2=S3+mean_x
  m3=pnorm(m2)          # Reverse transform  to Uniform random variables
  retrans_rad=array(0,dim=s)
  for ( i in 1: s){
    retrans_rad[i]=min(data1[r_u>=m3[i]])
  }
  x22=radius_dist[,k]*cos(theta)
  y22=radius_dist[,k]*sin(theta)
  lines(x22,y22,type="l")

###### PARAMETRIC
  x20=(p_x%*%radius_dist[,k])*cos(theta)
  y20=(p_x%*%radius_dist[,k])*sin(theta)
  r_para=sqrt(x20^2+y20^2)

##### Total actual area 
  a0=a1=a2=a3=astar=array(0,dim=length(theta))
  
  for (i in 1:length(theta)){
    c=c(0,unlist(radius_dist0[[i]]))
    cs=sign(unlist(radius_dist0[[i]]))
    a0[i]=aray(c,cs)
    
    d1=symmdist(unlist(radius_dist0[[i]]),retrans_rad[i])  # Distance with non-parametric method
    a1[i]=aray(d1[,1],d1[-nrow(d1),2])
    
    d2=symmdist(unlist(radius_dist0[[i]]),r_para[i])        # With parametric method
    a2[i]=aray(d2[,1],d2[-nrow(d2),2])
    
    d3=symmdist(unlist(radius_dist0[[i]]),unlist(radius_srvf[[i]]))
    a3[i]=aray(d3[,1],d3[-nrow(d3),2])
    
  }
  actual_area[k]=(2*pi/s)*sum(a0)
  fda_error[k]=(2*pi/s)*sum(a1)/actual_area[k]
  para_error[k]=(2*pi/s)*sum(a2)/actual_area[k]
  srvf_error[k]=(2*pi/s)*sum(a3)/actual_area[k]
  print(c(k,para_error[k],srvf_error[k],date()))
  }
c1=cbind(fda_error,para_error,srvf_error)

par(mfcol=c(3,3))
########## With respect to SVF representation
best=which(c1[,3]==min(c1[,3]))
med1=intersect(which(c1[,3]>(median(c1[,3])-0.00005)),which(c1[,3]<(median(c1[,3])+0.00005))) # Median reconstruction
worst=which(c1[,3]==max(c1[,3]))


krs=c(best,med1[1],worst)
for(k in 1:3){ # Main loop
  k=krs[k]
  print(date())
  radius_dist0=list()
   radius_srvf=list()
  
  S1=array(0,dim=2002)
  for (i in 1:MM1)
  {
    d1=c(K0$qn_out[,,k][1,],K0$qn_out[,,k][2,])  
    m1=c(K0$mu[1,],K0$mu[2,])
    score1=c(t(hv1[,i])%*%(d1-m1))
    S1=S1+score1*hv1[,i]
  }
  m2=S1+m1
  q_data[,,k]=rbind(m2[1:1001],m2[1002:2002])
  beta1=q_to_curve(q_data[,,k])
  T1=1001
  centroid1 = calculatecentroid1(beta1)
  dim(centroid1) = c(length(centroid1),1)
  beta1 = beta1 - repmat1(centroid1,1,T1)
  XY0=rbind(t(beta1)[-1001,],t(beta1)[1,]) # SRVF reconstructed contour
  XY=t(shape_data_out[,,k]) # Actual contour
 
  theta22=NULL
  r_theta22=NULL
  theta220=NULL
  r_theta220=NULL
  for (i1 in 1:(nrow(XY)-1)){
    # For actual contour 
    
    xy1=XY[i1,]          # one vertex of a line segment of a polygon
    xy2=XY[i1+1,]        # other vertex
    
    t1=f1(xy1[1],xy1[2])         # theta of polar representation of the first vertex
    r1=sqrt(xy1[1]^2+xy1[2]^2)   # length
    
    t2=f1(xy2[1],xy2[2])
    r2=sqrt(xy2[1]^2+xy2[2]^2)
    
    theta11=NULL                  # theta variable for actual contour
    r_theta=NULL                  # corresponding radius
    
    if  (abs(t2-t1)<pi) {theta11=theta1[intersect(which(theta1>=min(t1,t2)),which(theta1<=max(t1,t2)))]} else
    {theta11=theta1[c(which(theta1>=max(t1,t2)),which(theta1<=min(t1,t2)))]}
    
    r_theta=(r1*r2*sin(t1-t2))/(r1*sin(t1-theta11)+r2*sin(theta11-t2))
    
    theta22=c(theta22,theta11)
    r_theta22=c(r_theta22,r_theta)
    
    # Replication for srvf reconstruced contour
    xy10=XY0[i1,]          # one vertex of a line segment of a polygon
    xy20=XY0[i1+1,]        # other vertex
    
    t10=f1(xy10[1],xy10[2])         # theta of polar representation of the first vertex
    r10=sqrt(xy10[1]^2+xy10[2]^2)   # length
    
    t20=f1(xy20[1],xy20[2])
    r20=sqrt(xy20[1]^2+xy20[2]^2)
    
    theta110=NULL
    r_theta0=NULL
    
    if  (abs(t20-t10)<pi) {theta110=theta1[intersect(which(theta1>=min(t10,t20)),which(theta1<=max(t10,t20)))]} else
    {theta110=theta1[c(which(theta1>=max(t10,t20)),which(theta1<=min(t10,t20)))]}
    
    r_theta0=(r10*r20*sin(t10-t20))/(r10*sin(t10-theta110)+r20*sin(theta110-t20))
    
    theta220=c(theta220,theta110)
    r_theta220=c(r_theta220,r_theta0)
  }
  
  ###### To detect whether the centoid is inside or outside the contour
  l1=array(0,dim=s)  # Number of cut points for each directions
  for (i2 in 1:s){
    l1[i2]=length(r_theta22[which(theta22==theta1[i2])])
  }
  if (min(l1)>=1){             #Atleast one intersect in all directions
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta22[which(theta22==theta1[i2])])
      radius_dist[i2,k]=max(r_list)    # Star-hull contour
      sn=array(1,dim=length(r_list))
      sn[which(1:length(r_list)%%2==0)]=-1
      radius_dist0[[i2]]=r_list*sn
    }
    if (min(l1)>=2){            # a typical situation as in example 564 of out of sample
      for (i2 in 1:s){
        r_list=NULL
        r_list=sort(r_theta22[which(theta22==theta1[i2])])
        radius_dist[i2,k]=max(r_list)    # Star-hull contour
        sn=array(1,dim=length(r_list))
        sn[which(1:length(r_list)%%2!=0)]=-1
        radius_dist0[[i2]]=r_list*sn
      }
    }
  }else{                   # When centroid outside the contour and few directions has no intersection
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta22[which(theta22==theta1[i2])])
      if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
      radius_dist[i2,k]=max(r_list)    # Star-hull contour
      
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

###### For SRVF alignment: to detect whethere the centoid is inside or outside the contour
  l1=array(0,dim=s)  # Number of cut points for each directions
  for (i2 in 1:s){
    l1[i2]=length(r_theta220[which(theta220==theta1[i2])])
  }
  if (min(l1)>=1){ #  Atleast one intersect in all directions
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta220[which(theta220==theta1[i2])])
      sn=array(1,dim=length(r_list))
      sn[which(1:length(r_list)%%2==0)]=-1
      radius_srvf[[i2]]=r_list*sn
    }
    if (min(l1)>=2){ # a typical situation as in example 564 of test-sample 
      for (i2 in 1:s){
        r_list=NULL
        r_list=sort(r_theta220[which(theta220==theta1[i2])])
        sn=array(1,dim=length(r_list))
        sn[which(1:length(r_list)%%2!=0)]=-1
        radius_srvf[[i2]]=r_list*sn
      }
    }
  }else{  # When centroid outside the contour and few directions has no intersection
    for (i2 in 1:s){
      r_list=NULL
      r_list=sort(r_theta220[which(theta220==theta1[i2])])
      if (length(r_list)>=1) {r_list=r_list} else {r_list=0}
      
      if (length(r_list)==1) {radius_srvf[[i2]]=r_list}else{
        if (length(r_list)%%2==0){
          sn=array(1,dim=length(r_list))
          sn[which(1:length(r_list)%%2!=0)]=-1
          radius_srvf[[i2]]=r_list*sn
        }else{
          sn=array(1,dim=length(r_list))
          sn[which(1:length(r_list)%%2==0)]=-1
          radius_srvf[[i2]]=r_list*sn
        }
      }
    }
  }
  
###### Non-parametric computation
  S3=array(0,dim=c(s,1))
  m2=array(0,dim=c(s,1))
  for (l in 1:MM)
  {
    S3=S3+c((t(c(hv[,l]))%*%(c(data200[,k])-c(mean_x))))*hv[,l]
  }
  m2=S3+mean_x
  m3=pnorm(m2)          # Reverse transform  to Uniform random variables
  retrans_rad=array(0,dim=s)
  for ( i in 1: s){
    retrans_rad[i]=min(data1[r_u>=m3[i]])
  }
  xf=retrans_rad*cos(theta)
  yf=retrans_rad*sin(theta)

######## PARAMETRIC ######
  S31=array(0,dim=c(s,1))
  m21=array(0,dim=c(s,1))
  m21=data200[,k]          # Transformed contour
  m31=p_x%*%m21
  m41=pnorm(m31)          # Reverse transform  to Uniform random variables
  retrans_rad1=array(0,dim=s)
  for ( i in 1: s){
    retrans_rad1[i]=min(data1[r_u>=m41[i]])
  }
  xp=retrans_rad1*cos(theta)
  yp=retrans_rad1*sin(theta)  # lines(c(x2,x2[1]),c(y2,y2[1]),col="green",type="l",lty=1)
  r_para=sqrt(xp^2+yp^2)
   
###### Total actual area 
  a0=a1=a2=a3=astar=array(0,dim=length(theta))
  for (i in 1:length(theta)){
    c=c(0,unlist(radius_dist0[[i]]))
    cs=sign(unlist(radius_dist0[[i]]))
    a0[i]=aray(c,cs)
    
    d1=symmdist(unlist(radius_dist0[[i]]),retrans_rad[i]) # Distance with non-parametric method
    a1[i]=aray(d1[,1],d1[-nrow(d1),2])
    
    d2=symmdist(unlist(radius_dist0[[i]]),r_para[i]) # With parametric method
    a2[i]=aray(d2[,1],d2[-nrow(d2),2])
    
    d3=symmdist(unlist(radius_dist0[[i]]),unlist(radius_srvf[[i]]))
    a3[i]=aray(d3[,1],d3[-nrow(d3),2])
    
  }
  actual_area[k]=(2*pi/s)*sum(a0)
  fda_error[k]=(2*pi/s)*sum(a1)/actual_area[k]
  para_error[k]=(2*pi/s)*sum(a2)/actual_area[k]
  srvf_error[k]=(2*pi/s)*sum(a3)/actual_area[k]
  print(c(k,para_error[k],srvf_error[k],date()))
  
  plot(xf,yf, ann=FALSE,xlim=range(XY[,1],xf,XY0[,1],xp),
       ylim=range(XY[,2],yf,XY0[,2],yp),lwd=1.5,lty=2,type="l")
  lines(XY[,1],XY[,2])
  er1=signif(100*fda_error[k],2)
  legend("bottomright",legend=c(paste0=c(er1)))
  
  plot(xp,yp, ann=FALSE,xlim=range(XY[,1],xp,XY0[,1],xf),
       ylim=range(XY[,2],yp,XY0[,2],yf),lwd=1.5,lty=2,type="l")
  lines(XY[,1],XY[,2])
  er1=signif(100*para_error[k],2)
  legend("bottomright",legend=c(paste0=c(er1)))
  
  plot(XY0[,1],XY0[,2], ann=FALSE,xlim=range(XY[,1],xp,XY0[,1],xf),
       ylim=range(XY[,2],yp,XY0[,2],yf),lwd=1.5,lty=2,type="l")
  lines(XY[,1],XY[,2])
  er1=signif(100*srvf_error[k],2)
  legend("bottomright",legend=c(paste0=c(er1)))
  
}
title(xlab = "West-East (in kilometer)",
      ylab = "South-North (in kilometer)",
      outer = TRUE, line = 1,cex.lab=1.5)

print(date())
