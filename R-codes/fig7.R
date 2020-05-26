# This code is licensed under the MIT License

# Input: Starhull, area and missing status of regions 
# Output: Figure 7 of the paper.

# Execution time: approx 22 minutes (Intel(R) Core(TM) i7-7700HQ CPU @2.81 GHz with RAM 16 GB)


# load libraries
library(sp)
library(geosphere)
library(survival)

# set directory
# setwd("/home")


# input files
input_file1 <- "data/starhull_complete_contours"    # starhull for 2002-2012
input_file2 <- "data/area_status_combined/area.Rda" # area of the rainfall regions
input_file3 <- "data/area_status_combined/status.Rda" # complete (2) or censored (1)

# read input
file1 <- list.files(input_file1, pattern = "\\.Rda$",full.names = TRUE) 
ar <- readRDS(file=input_file2)
st <- readRDS(file=input_file3)

data1 <- NULL
n1 <- NULL
for (j in 1:11){
  dt <- NULL
  dt <- readRDS(file=file1[j]) # All contours of a year 
  n1[j] <- ncol(dt)
  data1 <- cbind(data1,dt) # all complete larger contour function
}
data_total <- data1   # no: 11140
n11 <- cumsum(n1)
s <- 1000
N11 <- ncol(data1)   # Number of complete contours
r_u <- (rank(data1))/(s*N11) # ECDF
r_u <- replace(r_u,which(r_u==1),0.9999)   ### Check this line again
data2 <- matrix(qnorm(r_u),nrow=s)   # Normal (random variable) transformation
data_star <- data2
dt4 <- intersect(which(ar>=200),which(ar<=13500)) # All contours in the bound
dt1 <- which(st==1)  # Indices of all censored regions
dt2 <- which(st==2)  # Indices of all complete regions

dt11 <- intersect(dt1,dt4)  # Indices of larger (> 200 sq km && <13500) censored spell
dt22 <- intersect(dt2,dt4)  # Indices of larger (> 200 sq km && <13500) complete spell 

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
mini.surv <- survfit(Surv(time, censor)~ 1) # Kaplan Miere computation

T1 <- mini.surv$time  # Increasing time period 
F1 <- mini.surv$surv  # Estimate of survial function

F2 <- c(1,F1[-length(F1)])-F1
F3 <- cbind(T1,F2)               # (index,prob. mass)
F4 <- c(F3[which(F3[,2]>0),2])   # Weights

order_com <- order(ar[dt22])
x0 <- order_com # Order the complete contours 
s=1000
N1 <- length(x0)

mean_x <- c(data_star[,x0]%*%F4)   # Weighted mean
cov1 <- (tcrossprod((data_star[,x0]-mean_x)*(matrix(sqrt(F4), nrow = s, ncol = N1, byrow = TRUE))))
h1 <- eigen(cov1)$values  
hv <- eigen(cov1)$vectors  
h2 <- eigen(cov1)$values[1:20]   
pc_1 <- ((h1[1:10])/sum(h1))*100;pc_1 # Percentage of the variation explained by the 10 PCS

# g(r(\theta))=phi^{-1}(F_n) transformation 
theta <- seq(0,2*pi,length=s+1)[-(s+1)] 
m32 <- pnorm(mean_x)
u41 <- array(0,dim=s)
for ( i in 1: s){
  u41[i] <- min(data_total[r_u>=m32[i]])
}

# Subfigures 
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat  <-  m,heights= c(3,3),widths=c(2,2,2))
par(mar <-  c(0.1, 0.1, 0.1, 0.1), oma = c(4,4,0,0),xpd = NA)       
mm1=6
t1 <- c(1,0,-1)
for (j in 1:mm1){
  plot(0,0, yaxt='n',xaxt='n', ann=FALSE,col="white",type="l",lty=2,cex.main=0.9,
       xlim=range(-28,28),ylim=range(-28,28))  
  x <- c(-20,-10,0,10,20)
  y <- c(-20,-10,0,10,20)
  if (j>=4) {axis(side = 1,cex.axis=0.9,at=x,labels = if (j>= 4) levels(x) else FALSE)}
  if (j==1||j==4){axis(side = 2,cex.axis=0.9,at=y, labels = if (j==1||j==4) levels(y) else FALSE)}
  mtext(paste0("Mode", j), cex=1,side = 1, adj = 0.95, line = -1.3)
  kl <- c(2,1,3)
  cl <- c(4,1,2)
  lwd1 <- c(1.5,1.5,1.5)
  for(k in 1:length(t1))
  {
    s_theta <- NULL
    s_theta <- mean_x+t1[k]*(sqrt(h1[j]))*(hv[,j])
    m31 <- pnorm(s_theta)
    u4 <- array(0,dim=s)
    for ( i in 1:s)
    {
      u4[i] <- min(data_total[r_u>=m31[i]])
    }
    xm3 <- u4*cos(theta)
    ym3 <- u4*sin(theta)
    xm3 <- c(xm3,xm3[1])
    ym3 <- c(ym3,ym3[1])
    lines(xm3,ym3,type="l",col=cl[k],lty=kl[k],lwd=lwd1[k])
    }
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 1,
       legend = c("alpha=-1","alpha=0", "alpha=1"), 
       lwd=c(1.5,1.5,1.5), cex=1, horiz = TRUE,lty=kl,col=cl)

title(xlab = "West-East (in kilometer)",
      ylab = "South-North (in kilometer)",
      outer = TRUE, line = 2,cex=1.5)

