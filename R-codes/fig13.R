# This code is licensed under the MIT License

# Input: Starhull contours over the years 2002-2012
# Output: Figure 13 of the paper

# execution time: approx 15 minutes (Intel(R) Core(TM) i7-7700HQ CPU @2.81 GHz with RAM 16 GB)

# setwd("/home/")

# load libraries
library(sp)
library(geosphere)
library(survival)
library(cluster)

# input files
input_file <- "data/starhull_complete_contours"
file1 <- list.files(input_file,pattern = "\\.Rda$",full.names = TRUE) 
data1 <- NULL
n1 <- NULL
for (j in 1:11){
  dt3 <- NULL
  dt3 <- readRDS(file=file1[j])       # contours of a year 
  n1[j] <- ncol(dt3)
  data1 <- cbind(data1,dt3)            # all complete larger contour function
}

data_total <- data1   
n11 <- c(0,cumsum(n1))
s <- 1000                                  # discretization length
N11 <- ncol(data1)                         # number of complete contours
r_u <- (rank(data1))/(s*N11)               # ECDF
r_u <- replace(r_u,which(r_u==1),0.9999)   
data_tarns <- matrix(qnorm(r_u),nrow=s)    # Normal (random variable) transformation

gzero <- min(data_total[r_u>=pnorm(0)])
mean2 <- array(gzero,dim=s)

d1 <- array(0,dim=c(11,4))
d2 <- array(0,dim=c(11,4))

theta <- seq(0,2*pi,length=s+1)[-(s+1)] 
x1 <- seq(1,1,length=1000)
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
p_x <- (solve(t(x_M)%*%x_M))%*%(t(x_M))              # Least square estimates of the LM  

theta <- seq(0,2*pi,length.out = s+1)[-(s+1)]
abhat <- array(0,dim=c(11,2))
gchat <- array(0,dim=11)
chat <- array(0,dim=11)
for (m in 1:11)
{
  print(date())
  data2 <- NULL
  data2 <- data_tarns[,(n11[m]+1):n11[m+1]]  # Particular year (2011) transformed data
  ab <- array(0,dim=c(ncol(data2),2))
  for (j in 1:ncol(data2))
  {
    ab1 <- p_x%*%data2[,j]
    ab[j,] <- c(ab1[4],ab1[5])
  }
  
#First Fourier function 
  x <- ab[,1]
  y <- ab[,2]
  n <- length(x)
  h <- (mean(sd(x),sd(y)))*((1/n)^(1/6))
  x0 <- x
  y0 <- y
  l <- length(x)
  f2 <- array(0,dim=c(l,l))
  f3=function (xx,yy) {(1/n)*(1/h^2)*sum((exp(-(((xx-x)/h)^2)/2))*exp(-(((yy-y)/h)^2)/2))}
  f2 <- outer(x0,y0,Vectorize(f3))
  ij  <-  which(f2==max(f2))
  iy <-  ceiling(ij/l)
  ix <-  ij%%l
  abhat[m,] <- c(x[ix],y[iy])
  chat[m] <- sqrt(x[ix]^2+y[iy]^2)
}
for (i in 1:11){
  gchat[i] <- min(data_total[r_u>=pnorm(chat[i])])
}

plot(range(c(13.5,-13.5)),range(c(13.5,-13.5)),type="n",xlab="West-East (in kilometers)",
     ylab="South-North (in kilometers)",cex.lab=0.8,cex.axis=0.8)
r <- mean2                                 # Mean function 
x1 <- r*cos(theta)
y1 <- r*sin(theta)
lines(x1,y1,lwd=2)

# Plot of one circle with maximum strength of directionality 
m <- 9
r0 <- array(0,dim=1000)
data2 <- NULL
data2 <- data_tarns[,(n11[m]+1):n11[m+1]]  

# Parametric  coefficient 
ab <- array(0,dim=c(ncol(data2),2))
for (j in 1:ncol(data2))
{
  ab1 <- p_x%*%data2[,j]
  ab[j,] <- c(ab1[4],ab1[5])
}

# First Fourier function
x <- ab[,1]
y <- ab[,2]
n <- length(x)
h <- (mean(sd(x),sd(y)))*((1/n)^(1/6))
l <- length(x)
x0 <- x
y0 <- y
f3 <- function (xx,yy) {(1/n)*(1/h^2)*sum((exp(-(((xx-x)/h)^2)/2))*exp(-(((yy-y)/h)^2)/2))}
f2 <- array(0,dim=c(l,l))
f2 <- outer(x0,y0,Vectorize(f3))
ij <-  which(f2==max(f2))
iy <-  ceiling(ij/l)
ix <-  ij%%l
rx <- x0[ix]*cos(2*theta)+y0[iy]*sin(2*theta)
for (i in 1:s){
  r0[i] <- min(data_total[r_u>=pnorm(rx[i])])
}

x_m <- r0*cos(theta)
y_m <- r0*sin(theta)                
lines(x_m,y_m,type="l",col=2)


col_1=c(1,2,4,3,5,6,1,3,2,4,5)
lty_1=c(2,3,4,5,6,8,9,10,1,11,12)
lwd_1=rev(seq(1,1.5,length.out=11))

for (i in 1:11){
  l1=(gchat[i]/chat[i])*(abhat[i,])
  l2=(gchat[i]/chat[i])*(-abhat[i,])
  lines(c(l1[1],l2[1]),c(l1[2],l2[2]),col=col_1[i],lty=lty_1[i],lwd=lwd_1[i])
}

legend(-14.5,-11,c("2002","2003","2004","2005","2006","2007"),col=c(1,2,4,3,5,6),cex=0.6,
         lty=c(2,3,4,5,6,8), box.lwd = 0,box.col = "white",bg = "white",horiz = TRUE,lwd=lwd_1[1:6])
legend(-14.5,-12.4,c("2008","2009","2010","2011","2012","Baseline"),col=c(1,3,2,4,5,1),cex=0.6,
         lty=c(9,10,1,11,12,1), box.lwd = 0,box.col = "white",bg = "white",horiz=TRUE,lwd=c(lwd_1[7:11],1.5))

print(date())
