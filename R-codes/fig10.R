# This code is licensed under the MIT License

# Input: Starhull contours for the years 2002-2012
# Output: Figure 10 of the paper

# load libraries
library(sp)
library(geosphere)
library(survival)

# set directory
# setwd("/home")

print(date())
# input files
input_file1 <- "data/starhull_complete_contours"
input_file2 <- "data/area_status_combined/area.Rda"
input_file3 <- "data/area_status_combined/status.Rda"

# read input
file1 <- list.files(input_file1,pattern = "\\.Rda$",full.names = TRUE) 
data1 <- NULL
n1 <- NULL
for (j in 1:11){
  dt <- NULL
  dt <- readRDS(file=file1[j])             # complete contours of a year 
  n1[j] <- ncol(dt)
  data1 <- cbind(data1,dt)                 # all complete contour 
}
data_total <- data1   
n11 <- cumsum(n1)
s <- 1000
N11 <- ncol(data1)                        # Number of complete contours
r_u <- (rank(data1))/(s*N11)              # Emperical CDF
r_u <- replace(r_u,which(r_u==1),0.9999)   
data2 <- matrix(qnorm(r_u),nrow=s)   # Normal (random variable) transformation
data_star <- data2

# Computation of Kaplan-Meier weights 
ar <- readRDS(file=input_file2)
st <-  readRDS(file=input_file3)

dt4 <- intersect(which(ar>=200),which(ar<=13500)) # All contours in the range
dt1 <- which(st==1)                                 # censored regions
dt2 <- which(st==2)                                 # complete regions

dt11 <- intersect(dt1,dt4)               # (> 200 sq km && <13500) censored 
dt22 <- intersect(dt2,dt4)               # (> 200 sq km && <13500) complete spell 

N1 <- length(dt22)  # Number of complete observations
N2 <- length(dt11) # Number of censored observations

time1 <- NULL
censor1 <- NULL
time1 <- ar[dt4]    # Region (Area) which lies in the given threshold
censor1 <- st[dt4]
censor1[censor1 < 1.5] <- 0   #Censored is 0
censor1[censor1 > 1.5] <- 1  # Complete observations is 1

time_ind <- order(time1)
area_order <- time1[time_ind]
censor_order <- censor1[time_ind]
area_status <- cbind(area_order,censor_order)

c_max <- max(which(area_status[,2]==1))
area_status1 <- area_status[1:c_max,]
time <- NULL
censor <- NULL
time <- rank(area_status1[,1])
censor <- area_status1[,2]
mini.surv <- survfit(Surv(time, censor)~ 1) # Kaplan-Meier computation
T1 <- mini.surv$time                           # Increasing time period 
F1 <- mini.surv$surv                           # Estimate of survial function
F2 <- c(1,F1[-length(F1)])-F1
F3 <- cbind(T1,F2)                 # (index,prob. mass)
F4 <- c(F3[which(F3[,2]>0),2])     # Weights
order_com <- order(ar[dt22])
x0 <- order_com                   # Order the complete contours
N1 <- length(x0)

mean_x <- c(data_star[,x0]%*%F4)   # Weighted mean
cov1 <- (tcrossprod((data_star[,x0]-mean_x)*(matrix(sqrt(F4), 
                  nrow = s, ncol = N1, byrow = TRUE))))  
h1 <- eigen(cov1)$values  
hv <- eigen(cov1)$vectors  
pc_1 <- ((h1[1:10])/sum(h1))*100;pc_1 # Percentage of the variation explained by the PCS

# plot configuration
par(mfrow=c(2,5),mar = c(0.1, 0.1, 0.1, 0.1), oma = c(4,4,0,0),xpd = NA)
k <- 10139                    # select a random example of contour
S1 <- array(0,dim = c(s,1))
M <- 10                         # Taking out first 10 PCS
mt <- array(0,dim = c(s,M))

#Cumulative sum 
S3 <- array(0,dim=c(s,1))
m2 <- array(0,dim=c(s,1))
for (l in 1:M)
{
  S3 <- S3+(t(hv[,l])%*%(data2[,k]-mean_x))*hv[,l]
  m2 <- S3+mean_x
  m3 <- pnorm(m2)          # Reverse transform  to Uniform random variables
  r2 <- array(0,dim = s)
  for ( i in 1: s){
    r2[i] <- min(data1[r_u>=m3[i]])
  }
   mt[,l] <- r2
}

t2 <- c(1,3,5,7,9)
for(l in 1:length(t2)){
  r2 <-  mt[,t2[l]]
  theta <- seq(0,2*pi,length= s+1)[-(s+1)]
  x22 <- r2*cos(theta)
  y22 <- r2*sin(theta)
  x22 <- c(x22,x22[1])
  y22 <- c(y22,y22[1])
 
  r3 <- data1[,k]
  x33 <- r3*cos(theta)
  y33 <- r3*sin(theta)
  x33 <- c(x33,x33[1])
  y33 <- c(y33,y33[1])
  
  plot(x33,y33, yaxt='n',xaxt='n', ann=FALSE,xlim=range(x22,x33,y22,y33),ylim=range(x22,x33,y22,y33),type="l",cex.main=1,lty=2)
  lines(x22,y22,type="l",lty=1,col="blue")
 
  y=c(-15,0,15)    
  if (l==1) {axis(side = 2,at=y,cex.axis=1,labels = if (l==1) levels(y) else FALSE)}
  if (l==6) {axis(side = 2,at=y,cex.axis=1,labels = if (l==6) levels(y) else FALSE)}
}

# Parametric presentation 
theta <- seq(0,2*pi,length = s+1)[-(s+1)]
x1 <- seq(1,1,length = 1000)
x2 <- cos(theta)
x3 <- sin(theta)
x4 <- cos(2*theta)
x5 <- sin(2*theta)
x6 <- cos(3*theta)
x7 <- sin(3*theta)
x8 <- cos(4*theta)
x9 <- sin(4*theta)
x10 <- cos(5*theta)
x11 <- sin(5*theta)
x_M <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)

#Cumulative sum #######
s2 <- array(0,dim=c(s,1))
x_M1 <- x_M[,1]
p_x <- (solve(t(x_M1)%*%x_M1))%*%(t(x_M1))
s2 <- (x_M1%*%(p_x%*%data2[,k]))
r2 <- NULL
m3 <- pnorm(s2)          # Reverse transform to Uniform random variables
sr2 <- array(0,dim=s)
for ( i in 1: s){
  r2[i] <- min(data1[r_u>=m3[i]])
}
x22 <- r2*cos(theta)
y22 <- r2*sin(theta)
x22 <- c(x22,x22[1])
y22 <- c(y22,y22[1])

r3 <- data1[,k]
x33 <- r3*cos(theta)
y33 <- r3*sin(theta)
x33 <- c(x33,x33[1])
y33 <- c(y33,y33[1])
plot(x33,y33, yaxt='n',xaxt='n', ann=FALSE,xlim=range(x22,x33,y22,y33),
     ylim=range(x22,x33,y22,y33),type="l",cex.main=1,lty=2)
lines(x22,y22,type="l",lty=1,col="blue")
x <- c(-15,0,15)
y <- c(-15,0,15)

axis(side = 1,cex.axis=1,labels = levels(x),at=x)
axis(side = 2,cex.axis=1,labels = levels(y),at=y)

g1=c(2,4,6,8)
for (j in 1:length(g1))
{
  x_M1 <- NULL
  x_M1 <- x_M[,c(g1[j],g1[j]+1)]
  p_x <- (solve(t(x_M1)%*%x_M1))%*%(t(x_M1))
  s2 <- s2+(x_M1%*%(p_x%*%data2[,k]))
  m2 <- s2
  m3 <- pnorm(m2)          # Reverse transform  to Uniform random variables
  r2 <- array(0,dim = s)
  for ( i in 1: s){
    r2[i]=min(data1[r_u>=m3[i]])
  }

  x22 <- r2*cos(theta)
  y22 <- r2*sin(theta)
  x22 <- c(x22,x22[1])
  y22 <- c(y22,y22[1])

  r3 <- data1[,k]
  x33 <- r3*cos(theta)
  y33 <- r3*sin(theta)
  x33 <- c(x33,x33[1])
  y33 <- c(y33,y33[1])
  plot(x33,y33, yaxt='n',xaxt='n', ann=FALSE,xlim=range(x22,x33,y22,y33),ylim=range(x22,x33,y22,y33),type="l",lty=2)
  lines(x22,y22,type="l",lty=1,col="blue")

  x <- seq(min(x22,x33,y22,y33),max(x22,x33,y22,y33),length.out  <- 4)
  y <- seq(min(x22,x33,y22,y33),max(x22,x33,y22,y33),length.out <- 5)

  x <- c(-15,0,15)
  axis(side = 1,cex.axis=0.9, at=x,labels =levels(x))
}
title(xlab = "West-East (in km)",
      ylab = "South-North (in km)",
      outer = TRUE, line = 2,cex.lab=1)

print(date())