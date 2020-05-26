# This code is licensed under the MIT License

# Input: Starhull contours over the years 2002-2012
# Output: Figure 12 of the paper

# execution time: approx 90 minutes(Intel(R) Core(TM) i7-7700HQ CPU @2.81 GHz with RAM 16 GB)

# load libraries
library(sp)
library(geosphere)
library(survival)

# set directory
# setwd("/home")

# input files
input_file <- "data/starhull_complete_contours"
file1 <- list.files(input_file,pattern = "\\.Rda$",full.names = TRUE)
data1=NULL
n1=NULL
for (j in 1:11){
 dt3=NULL
  dt3=readRDS(file=file1[j])           # contours of a year 
  n1[j]=ncol(dt3)
  data1=cbind(data1,dt3)               # all complete larger contour function
}
data_total=data1   
n11=c(0,cumsum(n1))
s=1000
N11=ncol(data1)                         # Number of complete contours
r_u=(rank(data1))/(s*N11)               # ECDF
r_u=replace(r_u,which(r_u==1),0.9999)   
data_tarns=matrix(qnorm(r_u),nrow=s)   # Normal (random variable) transformation

# parametric representation: coefficients
theta=seq(0,2*pi,length=s+1)[-(s+1)] 
x1=seq(1,1,length=s)
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
x_M=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
p_x=(solve(t(x_M)%*%x_M))%*%(t(x_M))  # Least square estimates of the LM  

amp4=NULL
for (m in 1:11)
{
  print(date())
  data2=NULL
  data2=data_tarns[,(n11[m]+1):n11[m+1]]  # Particular year (2011) transformed dat
  amp1=array(0,dim=c(ncol(data2),6))
  amp2=array(0,dim=c(ncol(data2),6))
  for (j in 1:ncol(data2))
  {
    ab1=p_x%*%data2[,j]
    amp1[j,]=c(ab1[1],sqrt(ab1[2]^2+ab1[3]^2),sqrt(ab1[4]^2+ab1[5]^2),sqrt(ab1[6]^2+ab1[7]^2),sqrt(ab1[8]^2+ab1[9]^2),sqrt(ab1[10]^2+ab1[11]^2))
  }
  
amp3=array(0,dim=c(ncol(data2),6))
  for (j in 1:ncol(data2))
  {
    for (i in 1:6){
      amp3[j,i]=min(data_total[r_u>=pnorm(amp1[j,i])])
    }
  }
    amp4=rbind(amp4,amp3)
}

# Histogram of amp2
par(mfrow=c(2,3))
r=(max(amp4[,1])/min(amp4[,1]))^(1/20)
r1=c(min(amp4[,1]),(min(amp4[,1]))*(r^(c(1:19))))
r1=c(r1,150)
s1=c(11.3,13.2,15.4,18,21,24,29,33,39,45,53,62,72,84,98,115,134,156,182,213,250)
for(i in 1:6){
  hist(amp4[,i],col="gray",breaks=r1,xlim=c(0,45),ylim=c(0,.5),xlab="g(Amplitude)",ylab="Density",main="")
  mtext(paste0("(", i, ")"), side = 3, adj = 0.5, line = -1.3)
  }

